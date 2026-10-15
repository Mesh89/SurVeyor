#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

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

void split_ins_to_dup(const std::string& input_fname, const std::string& ins_to_dup_fname, const std::string& remaining_fname, int n_threads) {
    htsThreadPool pool = {NULL, 0};
    if (n_threads > 1) {
        pool.pool = hts_tpool_init(n_threads - 1);
        if (pool.pool == NULL) throw std::runtime_error("Failed to create the shared VCF thread pool.");
    }
    std::unique_ptr<hts_tpool, decltype(&hts_tpool_destroy)> pool_guard(pool.pool, &hts_tpool_destroy);
    vcf_file_t in = open_vcf(input_fname, "r", pool);
    std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr(bcf_hdr_read(in.get()), &bcf_hdr_destroy);
    if (hdr == NULL) throw std::runtime_error("Failed to read the header from " + input_fname + ".");
    vcf_file_t ins_to_dup = open_vcf(ins_to_dup_fname, "wz", pool);
    vcf_file_t remaining = open_vcf(remaining_fname, "wz", pool);
    if (bcf_hdr_write(ins_to_dup.get(), hdr.get()) != 0 || bcf_hdr_write(remaining.get(), hdr.get()) != 0) throw std::runtime_error("Failed to write a VCF header.");
    std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
    if (record == NULL) throw std::runtime_error("Failed to allocate a VCF record.");
    int status;
    while ((status = bcf_read(in.get(), hdr.get(), record.get())) == 0) {
        bcf_unpack(record.get(), BCF_UN_STR);
        bool converted = bcf_get_info(hdr.get(), record.get(), "INS_TO_DUP") != NULL;
        if (converted) {
            std::string id = record->d.id;
            id.resize(id.size() > 4 ? id.size()-4 : 0);
            if (bcf_update_id(hdr.get(), record.get(), id.c_str()) != 0) throw std::runtime_error("Failed to update a VCF ID.");
        }
        normalize_vcf_end(hdr.get(), record.get());
        if (bcf_write(converted ? ins_to_dup.get() : remaining.get(), hdr.get(), record.get()) != 0) throw std::runtime_error("Failed to write a VCF record.");
    }
    if (status != -1) throw std::runtime_error("Failed to read records from " + input_fname + ".");
    int in_close_result = hts_close(in.release()), ins_to_dup_close_result = hts_close(ins_to_dup.release()), remaining_close_result = hts_close(remaining.release());
    if (in_close_result != 0) throw std::runtime_error("Failed to finalize reading " + input_fname + ".");
    if (ins_to_dup_close_result != 0) throw std::runtime_error("Failed to finalize " + ins_to_dup_fname + ".");
    if (remaining_close_result != 0) throw std::runtime_error("Failed to finalize " + remaining_fname + ".");
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: split_ins_to_dup <in.vcf.gz> <ins_to_dup.vcf.gz> <remaining.vcf.gz> [threads]\n";
        return 1;
    }
    try {
        int n_threads = argc == 5 ? std::stoi(argv[4]) : 1;
        if (n_threads < 1) throw std::runtime_error("The number of threads must be positive.");
        split_ins_to_dup(argv[1], argv[2], argv[3], n_threads);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
