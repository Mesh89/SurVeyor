#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "sam_utils.h"
#include "vcf_utils.h"

std::shared_ptr<sv_t> copy_indel_atom(std::shared_ptr<sv_t> sv) {
    std::shared_ptr<sv_t> copy;
    if (sv->svtype() == "DEL") {
        copy = std::make_shared<deletion_t>(sv->chr, sv->start, sv->end, sv->ins_seq, nullptr, nullptr, nullptr, nullptr);
    } else if (sv->svtype() == "INS") {
        copy = std::make_shared<insertion_t>(sv->chr, sv->start, sv->end, sv->ins_seq, nullptr, nullptr, nullptr, nullptr);
    } else {
        throw std::runtime_error("Unsupported AUX_INDEL type in haplotype expansion: " + sv->svtype());
    }
    copy->id = sv->id;
    copy->source = sv->source;
    copy->hpid = sv->hpid;
    return copy;
}

void append_indel_atoms(std::vector<std::shared_ptr<sv_t>>& dest, std::shared_ptr<sv_t> sv) {

    if (sv->svtype() != "DEL" && sv->svtype() != "INS") return;

    if (sv->start != sv->end && !sv->ins_seq.empty()) {
        std::shared_ptr<sv_t> del = std::make_shared<deletion_t>(sv->chr, sv->start, sv->end, "", nullptr, nullptr, nullptr, nullptr);
        del->id = sv->id;
        del->source = sv->source;
        del->hpid = sv->hpid;
        dest.push_back(del);

        std::shared_ptr<sv_t> ins = std::make_shared<insertion_t>(sv->chr, sv->start, sv->start, sv->ins_seq, nullptr, nullptr, nullptr, nullptr);
        ins->id = sv->id;
        ins->source = sv->source;
        ins->hpid = sv->hpid;
        dest.push_back(ins);
        return;
    }

    dest.push_back(copy_indel_atom(sv));
}

std::shared_ptr<sv_t> make_aux_indel_record(std::shared_ptr<sv_t> parent, size_t aux_idx) {
    std::shared_ptr<sv_t> aux_record = copy_indel_atom(parent->aux_indels[aux_idx]);
    aux_record->id = parent->id + ".AUX_INDEL." + std::to_string(aux_idx);
    aux_record->source = parent->source;
    aux_record->sample_info = parent->sample_info;
    aux_record->hpid = parent->hpid;
    aux_record->junction_remap_ref_beg = parent->junction_remap_ref_beg;
    aux_record->junction_remap_ref_end = parent->junction_remap_ref_end;
    aux_record->aux_snps = parent->aux_snps;
    append_indel_atoms(aux_record->aux_indels, parent);
    for (size_t i = 0; i < parent->aux_indels.size(); i++) {
        if (i != aux_idx) {
            aux_record->aux_indels.push_back(copy_indel_atom(parent->aux_indels[i]));
        }
    }
    return aux_record;
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

        for (size_t aux_idx = 0; aux_idx < sv->aux_indels.size(); aux_idx++) {
            std::shared_ptr<sv_t> aux_record_sv = make_aux_indel_record(sv, aux_idx);
            bcf1_t* aux_record = bcf_init();
            sv2bcf(hdr, aux_record, aux_record_sv.get(), chr_seqs.get_seq(aux_record_sv->chr));
            copy_all_fmt(hdr, b, aux_record);
            out_records.push_back(aux_record);
        }
    }
    bcf_destroy(b);
    hts_close(in_vcf_file);

    std::stable_sort(out_records.begin(), out_records.end(), [](const bcf1_t* a, const bcf1_t* b) {
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
