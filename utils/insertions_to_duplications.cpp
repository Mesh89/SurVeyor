#include <string>
#include <unordered_set>

#include "../src/sw_utils.h"
#include "../src/sam_utils.h"
#include "../src/vcf_utils.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "htslib/vcf.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;

StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

int main(int argc, char* argv[]) {

    std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];
    std::string workdir = argv[4];

    chr_seqs.read_fasta_into_map(reference_fname);

    contig_map_t contig_map;
    contig_map.load(workdir);

    config_t config;
    config.parse(workdir + "/config.txt");

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* b = bcf_init();

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln_before, aln_after;

    std::vector<std::shared_ptr<sv_t>> svs;
    std::unordered_set<std::string> new_dup_ids;
    std::unordered_map<std::string, int> unique_key_to_idx;
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, b);
        if (sv == nullptr) continue;

        std::shared_ptr<duplication_t> new_dup = nullptr;

        std::string ins_seq = sv->ins_seq;
        if (sv->svtype() == "INS" && !ins_seq.empty() && ins_seq.length() >= 50 && ins_seq.length() < 16000 
        && ins_seq.find('-') == std::string::npos) {
            std::string contig_name = sv->chr;
            hts_pos_t contig_len = chr_seqs.get_len(contig_name);

            bool skip = false;
            if (sv->start <= ins_seq.length()) skip = true;
            if (sv->start+ins_seq.length() >= contig_len) skip = true;

            if (!skip) {
                aligner.Align(ins_seq.c_str(), chr_seqs.get_seq(contig_name)+sv->start-ins_seq.length()+1, ins_seq.length(), filter, &aln_before, 0);
                aligner.Align(ins_seq.c_str(), chr_seqs.get_seq(contig_name)+sv->start+1, ins_seq.length(), filter, &aln_after, 0);

                std::vector<int> prefix_scores = compute_prefix_scores(aln_after.cigar, ins_seq.length(), 1, -4, -6, -1);
                std::vector<int> suffix_scores = compute_suffix_scores(aln_before.cigar, ins_seq.length(), 1, -4, -6, -1);

                int max_score = 0, best_i = 0;
                for (int i = 0; i <= ins_seq.length(); i++) {
                    if (prefix_scores[i] + suffix_scores[i] > max_score) {
                        max_score = prefix_scores[i] + suffix_scores[i];
                        best_i = i;
                    }
                }

                bool aln_alter_lc = get_left_clip_size(aln_after) >= config.min_clip_len;
                bool aln_before_rc = get_right_clip_size(aln_before) >= config.min_clip_len;
                int aln_len = 0;
                if (!aln_alter_lc) aln_len += aln_after.query_end+aln_after.query_begin+1;
                if (!aln_before_rc) aln_len += aln_before.query_end-aln_before.query_begin+1;
                if (aln_len >= ins_seq.length()) {
                    auto anchor_aln = std::make_shared<sv_t::anchor_aln_t>(sv->start-ins_seq.length()+best_i, sv->start+best_i, ins_seq.length(), 0);
                    new_dup = std::make_shared<duplication_t>(contig_name, anchor_aln->start, anchor_aln->end, "", nullptr, nullptr, anchor_aln, anchor_aln);
                    new_dup->id = sv->id + "_DUP";
                    new_dup->source = sv->source;
                    new_dup->hpid = sv->hpid;
                    new_dup->ins_to_dup_similarity = max_score / double(ins_seq.length() * 1);
                    new_dup_ids.insert(new_dup->id);
                }
            }
        }

        sv->sample_info.gt.clear();
        svs.push_back(sv);
        if (new_dup != NULL) {
            new_dup->sample_info.gt.clear();
            std::string unique_key = new_dup->unique_key();
            if (unique_key_to_idx.count(unique_key)) {
                std::shared_ptr<duplication_t> prev_dup = std::dynamic_pointer_cast<duplication_t>(svs[unique_key_to_idx[unique_key]]);
                if (new_dup->ins_to_dup_similarity > prev_dup->ins_to_dup_similarity) {
                    svs[unique_key_to_idx[unique_key]] = new_dup;
                }
            } else {
                unique_key_to_idx[unique_key] = svs.size();
                svs.push_back(new_dup);
            }
        }
    }

    std::sort(svs.begin(), svs.end(), [&contig_map](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) { 
        size_t chr_a = contig_map.get_id(a->chr);
        size_t chr_b = contig_map.get_id(b->chr);
        return std::tie(chr_a, a->start) < std::tie(chr_b, b->start);
    });
    bcf_hdr_t* out_hdr = generate_vcf_header(chr_seqs, "", config, "");
    bcf_hdr_remove(out_hdr, BCF_HL_INFO, "INS_TO_DUP");
    int len = 0;
    const char* ins_to_dup_tag = "##INFO=<ID=INS_TO_DUP,Number=0,Type=Flag,Description=\"Originally an insertion.\">";
    bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, ins_to_dup_tag, &len));

    if (bcf_hdr_set_samples(out_hdr, NULL, 0) != 0) {
        throw std::runtime_error("Failed to unset samples in VCF header");
    }

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "w");
    if (bcf_hdr_write(out_vcf_file, out_hdr) != 0) {
        throw std::runtime_error("Failed to write VCF header to " +  out_vcf_fname);
    }
    for (const auto& sv : svs) {
        sv2bcf(out_hdr, b, sv.get(), chr_seqs.get_seq(sv->chr), true);
        if (new_dup_ids.count(sv->id) > 0) {
            bcf_update_info_flag(out_hdr, b, "INS_TO_DUP", "", 1);
        }
        if (bcf_write(out_vcf_file, out_hdr, b) != 0) {
            throw std::runtime_error("Failed to write to " + out_vcf_fname);
        }
    }

    bcf_destroy(b);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(out_hdr);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
}
