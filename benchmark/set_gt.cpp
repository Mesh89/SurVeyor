#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/thread_pool.h"
#include "htslib/vcf.h"
#include "../libs/cxxopts.h"
#include "../src/vcf_utils.h"

typedef std::unique_ptr<htsFile, decltype(&hts_close)> vcf_file_t;

vcf_file_t open_vcf(const std::string& filename, const char* mode, htsThreadPool& pool) {
    vcf_file_t file(bcf_open(filename.c_str(), mode), &hts_close);
    if (file == NULL) throw std::runtime_error("Failed to open " + filename + ".");
    if (pool.pool && hts_set_opt(file.get(), HTS_OPT_THREAD_POOL, &pool) != 0) throw std::runtime_error("Failed to attach the shared VCF thread pool to " + filename + ".");
    return file;
}

void set_gt(const std::string& input_fname, const std::string& output_fname, int max_copies, int n_threads) {
    // Share BGZF workers between input and output; one thread updates records in order.
    htsThreadPool pool = {NULL, 0};
    if (n_threads > 1) {
        pool.pool = hts_tpool_init(n_threads-1);
        if (pool.pool == NULL) throw std::runtime_error("Failed to create the shared VCF thread pool.");
    }
    std::unique_ptr<hts_tpool, decltype(&hts_tpool_destroy)> pool_guard(pool.pool, &hts_tpool_destroy);
    vcf_file_t in = open_vcf(input_fname, "r", pool);
    std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr(bcf_hdr_read(in.get()), &bcf_hdr_destroy);
    if (hdr == NULL) throw std::runtime_error("Failed to read the header from " + input_fname + ".");
    bool compressed = output_fname.size() >= 3 && output_fname.compare(output_fname.size()-3, 3, ".gz") == 0;
    vcf_file_t out = open_vcf(output_fname, compressed ? "wz" : "w", pool);
    if (bcf_hdr_write(out.get(), hdr.get()) != 0) throw std::runtime_error("Failed to write the VCF header.");

    // Encode each replacement once, then update all samples in one HTSlib call.
    int n_samples = bcf_hdr_nsamples(hdr.get());
    std::vector<int32_t> genotypes(2*n_samples, bcf_gt_missing);
    std::string filters;
    filters.reserve(4*n_samples);
    for (int i = 0; i < n_samples; i++) {
        genotypes[2*i+1] = bcf_gt_unphased(1);
        filters += "PASS";
    }

    int32_t current_rid = -1;
    hts_pos_t current_pos = -1;
    std::unordered_map<std::string, int> allele_counts;
    std::string allele_key;
    std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
    int status;
    while ((status = bcf_read(in.get(), hdr.get(), record.get())) == 0) {
        bcf_unpack(record.get(), BCF_UN_STR);
        if (max_copies > 0) {
            if (record->rid != current_rid || record->pos != current_pos) {
                current_rid = record->rid;
                current_pos = record->pos;
                allele_counts.clear();
            }
            // NUL separators retain the full ordered REF/ALT tuple without collisions.
            allele_key.clear();
            for (int i = 0; i < record->n_allele; i++) {
                allele_key += record->d.allele[i];
                allele_key += '\0';
            }
            int& count = allele_counts[allele_key];
            if (count == max_copies) continue;
            count++;
        }
        if (n_samples > 0) {
            if (record->n_allele < 2) throw std::runtime_error("Cannot set FORMAT/GT to ./1 without an ALT allele.");
            if (bcf_update_genotypes(hdr.get(), record.get(), genotypes.data(), genotypes.size()) != 0) throw std::runtime_error("Failed to update FORMAT/GT.");
            if (bcf_get_fmt(hdr.get(), record.get(), "FT") != NULL && bcf_update_format_char(hdr.get(), record.get(), "FT", filters.data(), filters.size()) != 0) throw std::runtime_error("Failed to update FORMAT/FT.");
        }
        normalize_vcf_end(hdr.get(), record.get());
        if (bcf_write(out.get(), hdr.get(), record.get()) != 0) throw std::runtime_error("Failed to write a VCF record.");
    }
    if (status != -1) throw std::runtime_error("Failed to read records from " + input_fname + ".");
    int in_close_result = hts_close(in.release()), out_close_result = hts_close(out.release());
    if (in_close_result != 0) throw std::runtime_error("Failed to finalize reading " + input_fname + ".");
    if (out_close_result != 0) throw std::runtime_error("Failed to finalize " + output_fname + ".");
}

int main(int argc, char** argv) {
    cxxopts::Options options("set_gt", "Set all sample genotypes to ./1 and existing FORMAT/FT fields to PASS.");
    options.add_options()
        ("input_vcf", "Input VCF or VCF.GZ", cxxopts::value<std::string>())
        ("output_vcf", "Output VCF (BGZF compressed when ending in .gz)", cxxopts::value<std::string>())
        ("deduplicate", "Max records per exact CHROM, POS, REF, ALT in normalized, sorted input (0 disables)", cxxopts::value<int>()->default_value("0"))
        ("t,threads", "Threads including the record-processing thread (1 for serial I/O)", cxxopts::value<int>()->default_value("2"))
        ("h,help", "Print usage");
    options.parse_positional({"input_vcf", "output_vcf"});
    options.positional_help("input_vcf output_vcf");
    try {
        auto args = options.parse(argc, argv);
        if (args.count("help")) { std::cout << options.help(); return 0; }
        if (!args.count("input_vcf") || !args.count("output_vcf")) { std::cerr << options.help(); return 1; }
        int n_threads = args["threads"].as<int>();
        if (n_threads < 1) throw std::runtime_error("--threads must be at least 1");
        int max_copies = args["deduplicate"].as<int>();
        if (max_copies < 0) throw std::runtime_error("--deduplicate must be at least 0");
        set_gt(args["input_vcf"].as<std::string>(), args["output_vcf"].as<std::string>(), max_copies, n_threads);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
