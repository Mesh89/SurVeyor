#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "../libs/cptl_stl.h"
#include "sam_utils.h"
#include "var_utils.h"
#include "vcf_utils.h"

struct aux_haplotype_t {
    std::string sequence;
    std::vector<allele_edit_t> edits;
};

aux_haplotype_t build_aux_haplotype(sv_t& sv, char* chr_seq, hts_pos_t chr_len) {
    size_t n_atoms = sv.aux_snps.size() + sv.aux_indels.size();
    if (n_atoms == 0 || n_atoms > 300) return {};
    std::string svtype = sv.svtype();
    if (svtype != "DEL" && svtype != "INS" && svtype != "RPL") return {};
    if (sv.incomplete_ins_seq() || (svtype != "DEL" && sv.ins_seq.empty())) return {};

    hts_pos_t atom_beg = sv.start, atom_end = sv.end;
    hts_pos_t left_inserted = 0, right_inserted = 0;
    for (const snp_t& snp : sv.aux_snps) {
        atom_beg = std::min(atom_beg, snp.pos);
        atom_end = std::max(atom_end, snp.pos);
    }
    for (const auto& indel : sv.aux_indels) {
        atom_beg = std::min(atom_beg, indel->start);
        atom_end = std::max(atom_end, indel->end);
        if (indel->end < sv.start) left_inserted += indel->ins_seq.length();
        if (indel->start >= sv.end) right_inserted += indel->ins_seq.length();
    }
    if (atom_end-atom_beg > 20000) return {};

    // Flank lengths count haplotype bases, including AUX insertions.
    hts_pos_t left_len = sv.start-atom_beg+201+left_inserted;
    hts_pos_t right_len = atom_end-sv.end+200+right_inserted;
    auto aux_indels = sv.aux_indels;
    auto aux_snps = sv.aux_snps;
    aux_haplotype_t result;
    std::vector<allele_edit_t> right_edits;
    std::unique_ptr<char[]> left, right;
    try {
        left.reset(generate_haplotype_left(chr_seq, sv.start, left_len, aux_indels, aux_snps, &result.edits, nullptr, sv.id, sv.chr));
        right.reset(generate_haplotype_right(chr_seq, chr_len, sv.end+1, right_len, aux_indels, aux_snps, &right_edits, nullptr, sv.id, sv.chr));
    } catch (...) {
        return {};
    }

    result.sequence = std::string(left.get()) + sv.ins_seq;
    int right_offset = result.sequence.length();
    result.sequence += right.get();
    to_uppercase(&result.sequence[0]);
    for (allele_edit_t edit : right_edits) {
        edit.alt_begin += right_offset;
        edit.alt_end += right_offset;
        result.edits.push_back(edit);
    }
    return result;
}

std::vector<std::string> atom_signatures(const std::string& haplotype, int begin, int end, const std::string& other) {
    std::vector<std::pair<int, int>> intervals;
    if (end-begin+20 > 120) {
        intervals.push_back({begin, begin});
        intervals.push_back({end, end});
    } else {
        intervals.push_back({begin, end});
    }

    std::vector<std::string> signatures;
    for (const auto& interval : intervals) {
        for (int flank = 10; flank <= 75; flank += 5) {
            if (interval.first < flank || interval.second+flank > haplotype.size()) break;
            std::string signature = haplotype.substr(interval.first-flank, interval.second-interval.first+2*flank);
            if (other.find(signature) == std::string::npos) {
                signatures.push_back(std::move(signature));
                break;
            }
        }
    }
    return signatures;
}

bool atom_is_contradicted(const aux_haplotype_t& haplotype, const allele_edit_t& edit, const char* chr_seq, const std::vector<std::string>& reads) {
    std::string reference(chr_seq+edit.ref_begin, edit.ref_end-edit.ref_begin);
    std::string reverted = haplotype.sequence;
    reverted.replace(edit.alt_begin, edit.alt_end-edit.alt_begin, reference);
    auto signatures = atom_signatures(haplotype.sequence, edit.alt_begin, edit.alt_end, reverted);
    auto anti_signatures = atom_signatures(reverted, edit.alt_begin, edit.alt_begin+reference.length(), haplotype.sequence);
    if (signatures.empty() || anti_signatures.empty()) return false;

    int signature_reads = 0, anti_signature_reads = 0;
    for (const std::string& read : reads) {
        signature_reads += std::any_of(signatures.begin(), signatures.end(), [&](const std::string& signature) { return read.find(signature) != std::string::npos; });
        anti_signature_reads += std::any_of(anti_signatures.begin(), anti_signatures.end(), [&](const std::string& signature) { return read.find(signature) != std::string::npos; });
    }
    return signature_reads == 0 && anti_signature_reads > 0;
}

void filter_aux_atoms(sv_t& sv, char* chr_seq, hts_pos_t chr_len, open_samFile_t& bam_file, stats_t& stats) {
    aux_haplotype_t haplotype = build_aux_haplotype(sv, chr_seq, chr_len);
    if (haplotype.edits.empty()) return;

    hts_pos_t begin = haplotype.edits.front().ref_begin, end = haplotype.edits.front().ref_end;
    for (const auto& edit : haplotype.edits) {
        begin = std::min(begin, edit.ref_begin);
        end = std::max(end, edit.ref_end);
    }
    hts_pos_t fetch_beg = std::max(hts_pos_t(0), begin-160), fetch_end = std::min(chr_len, end+160);
    int64_t max_reads = int64_t(4) * stats.get_max_depth(sv.chr) * (fetch_end-fetch_beg+stats.read_len) / stats.read_len;
    int tid = sam_hdr_name2tid(bam_file.header, sv.chr.c_str());
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(sam_itr_queryi(bam_file.idx, tid, fetch_beg, fetch_end), &hts_itr_destroy);
    if (!iter) throw std::runtime_error("Unable to fetch reads for " + sv.id + ".");
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), &bam_destroy1);
    std::vector<std::string> reads;
    int status;
    while ((status = sam_itr_next(bam_file.file, iter.get(), read.get())) >= 0) {
        if (!is_primary(read.get()) || read->core.l_qseq == 0) continue;
        std::string sequence(read->core.l_qseq, '\0');
        const uint8_t* encoded = bam_get_seq(read.get());
        for (int i = 0; i < read->core.l_qseq; i++) sequence[i] = seq_nt16_str[bam_seqi(encoded, i)];
        reads.push_back(std::move(sequence));
        if (int64_t(reads.size()) >= max_reads) return;
    }
    if (status < -1) throw std::runtime_error("Unable to read BAM for " + sv.id + ".");

    // Match edits back to the original AUX order, including shared-anchor insertions.
    std::vector<bool> used_edits(haplotype.edits.size(), false);
    auto contradicted = [&](allele_edit_kind_t kind, hts_pos_t ref_begin, hts_pos_t ref_end, const std::string& alt) {
        for (size_t i = 0; i < haplotype.edits.size(); i++) {
            const auto& edit = haplotype.edits[i];
            if (used_edits[i] || edit.kind != kind || edit.ref_begin != ref_begin || edit.ref_end != ref_end) continue;
            if (edit.alt_end-edit.alt_begin != alt.length() || haplotype.sequence.compare(edit.alt_begin, alt.length(), alt) != 0) continue;
            used_edits[i] = true;
            return atom_is_contradicted(haplotype, edit, chr_seq, reads);
        }
        return false;
    };
    std::vector<snp_t> retained_snps;
    for (const auto& snp : sv.aux_snps) {
        if (!contradicted(allele_edit_kind_t::SNP, snp.pos, snp.pos+1, std::string(1, snp.alt_base))) retained_snps.push_back(snp);
    }
    std::vector<std::shared_ptr<sv_t>> retained_indels;
    for (const auto& indel : sv.aux_indels) {
        std::string alt = indel->svtype() == "DEL" ? "" : indel->ins_seq;
        alt.erase(std::remove(alt.begin(), alt.end(), '-'), alt.end());
        if (!contradicted(allele_edit_kind_t::INDEL, indel->start+1, indel->end+1, alt)) retained_indels.push_back(indel);
    }
    sv.aux_snps = std::move(retained_snps);
    sv.aux_indels = std::move(retained_indels);
}

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
    if (parent->svtype() != "INS" || !parent->incomplete_ins_seq()) append_indel_atoms(aux_record->aux_indels, parent);
    for (size_t i = 0; i < parent->aux_indels.size(); i++) {
        if (i != aux_idx) {
            aux_record->aux_indels.push_back(copy_indel_atom(parent->aux_indels[i]));
        }
    }
    return aux_record;
}

int main(int argc, char* argv[]) {
    if (argc != 6 && argc != 7) {
        throw std::runtime_error("Usage: expand_aux_haplotypes <in.vcf.gz> <out.vcf.gz> <reference.fa> <bam> <workdir> [threads]");
    }

    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];
    std::string reference_fname = argv[3];
    std::string bam_fname = argv[4];
    std::string workdir = argv[5];
    int n_threads = argc >= 7 ? std::max(1, std::stoi(argv[6])) : 1;

    config_t config;
    config.parse(workdir + "/config.txt");
    stats_t stats;
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname, true, n_threads);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in_vcf_file);
    if (hdr == NULL) {
        throw std::runtime_error("Unable to read VCF header from " + in_vcf_fname + ".");
    }

    std::vector<std::shared_ptr<sv_t>> svs;
    bcf1_t* b = bcf_init();
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, b);
        if (sv == nullptr) continue;
        sv->vcf_entry = bcf_dup(b);
        svs.push_back(sv);
    }
    bcf_destroy(b);
    hts_close(in_vcf_file);

    std::vector<std::pair<std::shared_ptr<sv_t>, bcf1_t*>> expanded_svs;
    for (const auto& sv : svs) {
        expanded_svs.push_back({sv, sv->vcf_entry});
        for (size_t aux_idx = 0; aux_idx < sv->aux_indels.size(); aux_idx++) {
            expanded_svs.push_back({make_aux_indel_record(sv, aux_idx), sv->vcf_entry});
        }
    }
    std::vector<std::shared_ptr<sv_t>> filter_svs;
    for (const auto& record : expanded_svs) filter_svs.push_back(record.first);
    std::stable_sort(filter_svs.begin(), filter_svs.end(), [](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) { return std::tie(a->chr, a->start) < std::tie(b->chr, b->start); });

    bam_pool_t bam_pool(n_threads, bam_fname, reference_fname);
    ctpl::thread_pool thread_pool(n_threads);
    std::vector<std::future<void>> futures;
    for (size_t block_beg = 0; block_beg < filter_svs.size(); block_beg += 8) {
        size_t block_end = std::min(block_beg+8, filter_svs.size());
        futures.push_back(thread_pool.push([&bam_pool, &chr_seqs, &filter_svs, &stats, block_beg, block_end](int id) {
            open_samFile_t* bam_file = bam_pool.get_bam_reader(id);
            for (size_t i = block_beg; i < block_end; i++) {
                const auto& sv = filter_svs[i];
                filter_aux_atoms(*sv, chr_seqs.get_seq(sv->chr), chr_seqs.get_len(sv->chr), *bam_file, stats);
            }
        }));
    }
    for (auto& future : futures) future.get();

    std::vector<bcf1_t*> out_records;
    for (const auto& record : expanded_svs) {
        bcf1_t* out_record = bcf_init();
        sv2bcf(hdr, out_record, record.first.get(), chr_seqs.get_seq(record.first->chr));
        copy_all_fmt(hdr, record.second, out_record);
        out_records.push_back(out_record);
    }
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
    return 0;
}
