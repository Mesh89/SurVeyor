#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "sam_utils.h"
#include "vcf_utils.h"

void ensure_hpid_header(bcf_hdr_t* hdr) {
    int len = 0;
    bcf_hdr_remove(hdr, BCF_HL_INFO, "HPID");
    const char* hpid_tag = "##INFO=<ID=HPID,Number=1,Type=Integer,Description=\"Identifier of the local haplotype represented by this call.\">";
    bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, hpid_tag, &len));
    if (bcf_hdr_sync(hdr) < 0) {
        throw std::runtime_error("Failed to sync VCF header after adding HPID.");
    }
}

std::shared_ptr<sv_t> make_aux_indel_record(std::shared_ptr<sv_t> parent, std::shared_ptr<sv_t> aux_indel, int idx) {
    aux_indel->id = parent->id + ".AUX_INDEL." + std::to_string(idx);
    aux_indel->source = parent->source;
    aux_indel->sample_info = parent->sample_info;
    return aux_indel;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        throw std::runtime_error("Usage: expand_aux_haplotypes <in.vcf.gz> <out.vcf.gz> <reference.fa>");
    }

    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];
    std::string reference_fname = argv[3];

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in_vcf_file);
    if (hdr == NULL) {
        throw std::runtime_error("Unable to read VCF header from " + in_vcf_fname + ".");
    }
    ensure_hpid_header(hdr);

    std::vector<bcf1_t*> out_records;
    bcf1_t* b = bcf_init();
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, b);
        if (sv == nullptr) continue;
        sv->vcf_entry = bcf_dup(b);

        bcf1_t* main_record = bcf_init();
        sv2bcf(hdr, main_record, sv.get(), chr_seqs.get_seq(sv->chr));
        copy_all_fmt(hdr, b, main_record);
        out_records.push_back(main_record);

        int aux_idx = 0;
        for (const auto& aux_indel : sv->aux_indels) {
            std::shared_ptr<sv_t> aux_record_sv = make_aux_indel_record(sv, aux_indel, aux_idx++);
            bcf1_t* aux_record = bcf_init();
            sv2bcf(hdr, aux_record, aux_record_sv.get(), chr_seqs.get_seq(aux_record_sv->chr));
            copy_all_fmt(hdr, b, aux_record);
            out_records.push_back(aux_record);
        }
    }
    bcf_destroy(b);
    hts_close(in_vcf_file);

    std::sort(out_records.begin(), out_records.end(), [](const bcf1_t* a, const bcf1_t* b) {
        return std::tie(a->rid, a->pos) < std::tie(b->rid, b->pos);
    });

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    if (out_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
    }
    if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
        throw std::runtime_error("Failed to write VCF header to " + out_vcf_fname + ".");
    }
    for (bcf1_t* record : out_records) {
        if (bcf_write(out_vcf_file, hdr, record) != 0) {
            throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname + ".");
        }
    }

    hts_close(out_vcf_file);
    bcf_hdr_destroy(hdr);
    for (bcf1_t* record : out_records) {
        bcf_destroy(record);
    }

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
