#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "htslib/vcf.h"
#include "htslib/thread_pool.h"
#include "../src/vcf_utils.h"

typedef std::unique_ptr<htsFile, decltype(&hts_close)> vcf_file_t;

vcf_file_t open_vcf(const std::string& filename, const char* mode, htsThreadPool& pool) {
    vcf_file_t file(bcf_open(filename.c_str(), mode), &hts_close);
    if (file == NULL) throw std::runtime_error("Failed to open " + filename + ".");
    if (pool.pool && hts_set_opt(file.get(), HTS_OPT_THREAD_POOL, &pool) != 0) throw std::runtime_error("Failed to attach the shared VCF thread pool to " + filename + ".");
    return file;
}

struct evidence_t {
    float epr;
    bool imprecise;
};

struct epr_buffer_t {
    float* data = NULL;
    int capacity = 0;
    ~epr_buffer_t() { free(data); }
};

void dedup_vcf(const std::string& input_fname, const std::string& output_fname, const std::string& compare_fname, int n_threads) {
    std::ifstream compare(compare_fname);
    if (!compare) throw std::runtime_error("Failed to open " + compare_fname + ".");
    std::vector<std::pair<std::string, std::string> > pairs;
    std::unordered_set<std::string> needed_ids, removed_ids;
    std::string line, id1, id2;
    while (std::getline(compare, line)) {
        std::istringstream fields(line);
        if (!(fields >> id1 >> id2)) throw std::runtime_error("Expected two IDs in " + compare_fname + ".");
        if (id2 == "NONE" || id1 >= id2) continue;
        pairs.push_back({id1, id2});
        needed_ids.insert(id1);
        needed_ids.insert(id2);
    }
    if (compare.bad()) throw std::runtime_error("Failed to read " + compare_fname + ".");

    htsThreadPool pool = {NULL, 0};
    if (n_threads > 1) {
        pool.pool = hts_tpool_init(n_threads - 1);
        if (pool.pool == NULL) throw std::runtime_error("Failed to create the shared VCF thread pool.");
    }
    std::unique_ptr<hts_tpool, decltype(&hts_tpool_destroy)> pool_guard(pool.pool, &hts_tpool_destroy);
    if (!pairs.empty()) {
        vcf_file_t in = open_vcf(input_fname, "r", pool);
        std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr(bcf_hdr_read(in.get()), &bcf_hdr_destroy);
        if (hdr == NULL) throw std::runtime_error("Failed to read the header from " + input_fname + ".");
        std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
        if (record == NULL) throw std::runtime_error("Failed to allocate a VCF record.");
        std::unordered_map<std::string, evidence_t> evidence;
        epr_buffer_t eprs;
        int status;
        while ((status = bcf_read(in.get(), hdr.get(), record.get())) == 0) {
            bcf_unpack(record.get(), BCF_UN_STR);
            std::string id = record->d.id;
            if (!needed_ids.count(id)) continue;
            float epr = 0;
            if (bcf_get_fmt(hdr.get(), record.get(), "EPR") != NULL) {
                int count = bcf_get_format_float(hdr.get(), record.get(), "EPR", &eprs.data, &eprs.capacity);
                if (count <= 0 || bcf_float_is_missing(eprs.data[0]) || bcf_float_is_vector_end(eprs.data[0])) throw std::runtime_error("Missing or invalid FORMAT/EPR for ID " + id + ".");
                epr = eprs.data[0];
            }
            evidence[id] = {epr, bcf_get_info(hdr.get(), record.get(), "IMPRECISE") != NULL};
        }
        if (status != -1) throw std::runtime_error("Failed to read records from " + input_fname + ".");
        if (hts_close(in.release()) != 0) throw std::runtime_error("Failed to finalize reading " + input_fname + ".");
        for (const auto& pair : pairs) {
            for (const auto& id : {pair.first, pair.second}) {
                if (!evidence.count(id)) throw std::runtime_error("Pair ID " + id + " is absent from " + input_fname + ".");
            }
            const evidence_t& first = evidence.at(pair.first);
            const evidence_t& second = evidence.at(pair.second);
            if (first.imprecise != second.imprecise) removed_ids.insert(first.imprecise ? pair.first : pair.second);
            else if (first.epr < second.epr) removed_ids.insert(pair.first);
            else if (second.epr < first.epr) removed_ids.insert(pair.second);
        }
    }

    vcf_file_t in = open_vcf(input_fname, "r", pool);
    std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr(bcf_hdr_read(in.get()), &bcf_hdr_destroy);
    if (hdr == NULL) throw std::runtime_error("Failed to read the header from " + input_fname + ".");
    vcf_file_t out = open_vcf(output_fname, "wz", pool);
    if (bcf_hdr_write(out.get(), hdr.get()) != 0) throw std::runtime_error("Failed to write the VCF header.");
    std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
    if (record == NULL) throw std::runtime_error("Failed to allocate a VCF record.");
    int status;
    while ((status = bcf_read(in.get(), hdr.get(), record.get())) == 0) {
        bcf_unpack(record.get(), BCF_UN_STR);
        if (removed_ids.count(record->d.id)) continue;
        normalize_vcf_end(hdr.get(), record.get());
        if (bcf_write(out.get(), hdr.get(), record.get()) != 0) throw std::runtime_error("Failed to write a VCF record.");
    }
    if (status != -1) throw std::runtime_error("Failed to read records from " + input_fname + ".");
    int in_close_result = hts_close(in.release()), out_close_result = hts_close(out.release());
    if (in_close_result != 0) throw std::runtime_error("Failed to finalize reading " + input_fname + ".");
    if (out_close_result != 0) throw std::runtime_error("Failed to finalize " + output_fname + ".");
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: dedup_vcf <in.vcf.gz> <out.vcf.gz> <compare.txt> [threads]\n";
        return 1;
    }
    try {
        int n_threads = argc == 5 ? std::stoi(argv[4]) : 1;
        if (n_threads < 1) throw std::runtime_error("The number of threads must be positive.");
        dedup_vcf(argv[1], argv[2], argv[3], n_threads);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
