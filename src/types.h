#ifndef TYPES_H
#define TYPES_H

#include "htslib/hts.h"
#include "utils.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <climits>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <unordered_map>
#include "htslib/vcf.h"

struct consensus_alignment_metrics_t {
    int length = 0, alt_score = 0, ref_score = 0;
    int alt_ref_begin = 0, alt_ref_end = -1;
    std::array<int, 2> split_ref_lengths{{0, 0}};
    std::array<int, 2> split_sizes{{0, 0}};
    std::array<int, 2> split_scores{{0, 0}};
    std::array<int, 2> independent_ref_scores{{0, 0}};
};

struct consensus_t {
    bool left_clipped;
    hts_pos_t start, breakpoint, end;
    hts_pos_t orig_start, orig_end;
    std::string sequence;
    int fwd_reads, rev_reads;
    uint8_t max_mapq;
    hts_pos_t other_bp_lower_boundary = LOWER_BOUNDARY_NON_CALCULATED, other_bp_upper_boundary = UPPER_BOUNDARY_NON_CALCULATED;
    int clip_len, lowq_prefix, lowq_suffix;
    int left_ext_reads = 0, right_ext_reads = 0, hq_left_ext_reads = 0, hq_right_ext_reads = 0;
    bool is_hsr = false;
	bool extended_to_left = false, extended_to_right = false;

    enum : hts_pos_t {
        LOWER_BOUNDARY_NON_CALCULATED = HTS_POS_MIN,
        UPPER_BOUNDARY_NON_CALCULATED = HTS_POS_MAX
    };

    consensus_t(bool left_clipped, hts_pos_t start, hts_pos_t breakpoint, hts_pos_t end,
                const std::string& sequence, int fwd_reads, int rev_reads, int clip_len, uint8_t max_mapq, 
                int lowq_prefix, int lowq_suffix)
                : left_clipped(left_clipped), start(start), breakpoint(breakpoint), end(end),
                orig_start(start), orig_end(end),
                sequence(sequence), fwd_reads(fwd_reads), rev_reads(rev_reads), clip_len(clip_len), max_mapq(max_mapq), 
                lowq_prefix(lowq_prefix), lowq_suffix(lowq_suffix) {}

    consensus_t(std::string& line) {
        std::stringstream ss(line);
        char dir;
        int max_mapq_int;
        ss >> start >> end >> breakpoint >> dir >> sequence >> fwd_reads >> rev_reads
           >> max_mapq_int >> other_bp_lower_boundary >> other_bp_upper_boundary >> lowq_prefix >> lowq_suffix >> is_hsr;
        orig_start = start;
        orig_end = end;
        left_clipped = dir == 'L';
        max_mapq = (uint8_t) max_mapq_int;
        clip_len = left_clipped ? breakpoint - start : end - breakpoint;
    }

    std::string to_string() {
        std::stringstream ss;
        ss << start << " " << end << " " << breakpoint << (left_clipped ? " L " : " R ") << sequence << " ";
        ss << fwd_reads << " " << rev_reads << " " << (int)max_mapq << " " << other_bp_lower_boundary << " " << other_bp_upper_boundary << " " << lowq_prefix << " " << lowq_suffix << " ";
        ss << is_hsr;
        return ss.str();
    }

    int reads() { return fwd_reads + rev_reads; }

    hts_pos_t left_ext_target_start(int max_is, int read_len) {
    	if (!left_clipped) {
    		return start - max_is + read_len;
    	} else {
    		if (other_bp_lower_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint - max_is + read_len;
			} else {
				return other_bp_lower_boundary - max_is + read_len;
			}
    	}
    }
    hts_pos_t left_ext_target_end(int max_is, int read_len) {
		if (!left_clipped) {
			return start;
		} else {
			if (other_bp_upper_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint;
			} else {
				return other_bp_upper_boundary;
			}
		}
	}

    hts_pos_t right_ext_target_start(int max_is, int read_len) {
		if (!left_clipped) {
			if (other_bp_lower_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) {
				return breakpoint;
			} else {
				return other_bp_lower_boundary;
			}
		} else {
			return end;
		}
	}
    hts_pos_t right_ext_target_end(int max_is, int read_len) {
    	if (!left_clipped) {
			if (other_bp_upper_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint + max_is - read_len;
			} else {
				return other_bp_upper_boundary + max_is - read_len;
			}
		} else {
			return end + max_is - read_len;
		}
    }

    hts_pos_t anchor_start() {
        if (left_clipped) {
            return breakpoint;
        } else {
            return start;
        }
    }
    hts_pos_t anchor_end() {
        if (left_clipped) {
            return end;
        } else {
            return breakpoint;
        }
    }

    std::string clip_sequence() {
        if (left_clipped) {
            return sequence.substr(0, clip_len);
        } else {
            return sequence.substr(sequence.length()-clip_len);
        }
    }

    std::string highq_sequence() {
        return sequence.substr(lowq_prefix, sequence.length()-lowq_prefix-lowq_suffix);
    }
};

struct snp_t {
    std::string chr;
    hts_pos_t pos;
    char alt_base;

    snp_t(hts_pos_t pos, char alt_base) : pos(pos), alt_base(alt_base) {}
    snp_t(std::string chr, hts_pos_t pos, char alt_base) : chr(chr), pos(pos), alt_base(alt_base) {}
    snp_t(std::string& snp_str) {
        size_t colon_pos = snp_str.find(':');
        pos = std::stoll(snp_str.substr(0, colon_pos)) - 1;
        alt_base = snp_str[colon_pos+1];
    }

    std::string unique_key() const {
        if (!chr.empty()) return "SNP:" + chr + ":" + std::to_string(pos) + ":" + alt_base;
        return "SNP:" + std::to_string(pos) + ":" + alt_base;
    }
};

struct sv_t {

    struct anchor_aln_t {
        hts_pos_t start, end;
        int seq_len;
        int best_score;

        anchor_aln_t(hts_pos_t start, hts_pos_t end, int seq_len, int best_score) : 
            start(start), end(end), seq_len(seq_len), best_score(best_score) {}

        std::string to_string() {
            std::stringstream ss;
            ss << start << "-" << end;
            return ss.str();
        }
    };

    std::string id;
    std::string chr;
    hts_pos_t start, end;
    std::string ins_seq, inferred_ins_seq;
    int mh_len = 0;
    std::shared_ptr<anchor_aln_t> left_anchor_aln, right_anchor_aln;
    std::shared_ptr<consensus_t> rc_consensus, lc_consensus;
    std::vector<std::shared_ptr<sv_t>> aux_indels;
    std::vector<snp_t> aux_snps;

    std::string source;
    int hpid = 0;
    bool imprecise = false;
    hts_pos_t junction_remap_ref_beg = HTS_POS_MIN, junction_remap_ref_end = HTS_POS_MIN;
    bool hp_genotyped = false;
    hts_pos_t hp_ref_beg = HTS_POS_MIN, hp_ref_end = HTS_POS_MIN;

    static constexpr const double KS_PVAL_NOT_COMPUTED = -1.0;
    static constexpr const int SIZE_NOT_COMPUTED = INT32_MAX;
    static constexpr const double EXPECTED_ALT_READS_FREQ_NOT_COMPUTED = -1.0;

    double ks_pval = KS_PVAL_NOT_COMPUTED;
    int min_conf_size = SIZE_NOT_COMPUTED, max_conf_size = SIZE_NOT_COMPUTED, estimated_size = SIZE_NOT_COMPUTED;

    base_frequencies_t left_anchor_base_freqs, right_anchor_base_freqs;
    base_frequencies_t prefix_ref_base_freqs, suffix_ref_base_freqs;
    base_frequencies_t ins_prefix_base_freqs, ins_suffix_base_freqs;

    struct bp_reads_info_t {
        bool computed = false;

        int reads = 0;
        int consistent_fwd = 0, consistent_rev = 0;
        int consistent_min_mq = INT32_MAX, consistent_max_mq = 0;
        double consistent_avg_mq = 0, consistent_stddev_mq = 0;
        int consistent_high_mq = 0;
        double consistent_avg_score = 0, consistent_stddev_score = 0;
        int consistent_exact_fwd = 0, consistent_exact_rev = 0;
        int consistent_exact_min_mq = INT32_MAX, consistent_exact_max_mq = 0;
        int consistent_exact_high_mq = 0;
        int exact_fwd = 0, exact_rev = 0;
        int exact_min_mq = INT32_MAX, exact_max_mq = 0;
        int exact_high_mq = 0;
        int fwd_mate_cov_bps = 0, rev_mate_cov_bps = 0;
        int fwd_hq_mate_cov_bps = 0, rev_hq_mate_cov_bps = 0;

        int consistent_reads() { return consistent_fwd + consistent_rev; }
        int consistent_exact_reads() { return consistent_exact_fwd + consistent_exact_rev; }
        int exact_reads() { return exact_fwd + exact_rev; }
    };

    struct bp_pairs_info_t {
        bool computed = false;

        int pairs = 0, pos_high_mapq = 0, neg_high_mapq = 0;
        int pos_min_mq = INT32_MAX, pos_max_mq = 0, neg_min_mq = INT32_MAX, neg_max_mq = 0;
        double pos_avg_mq = 0, pos_stddev_mq = 0, neg_avg_mq = 0, neg_stddev_mq = 0;
        
        int lf_span = 0, rf_span = 0;
        double pos_avg_nm = 0, pos_stddev_nm = 0, neg_avg_nm = 0, neg_stddev_nm = 0;
    };

    struct bp_consensus_info_t {

        bp_reads_info_t reads_info;
        bp_pairs_info_t pairs_info;

    };

    struct other_read_info_t {
        bool hq = false;
        bool exact = false;
    };

    struct sample_info_t {
        static const int NOT_COMPUTED = -1;

        std::vector<int> gt;
        double expected_alt1_reads_frac = EXPECTED_ALT_READS_FREQ_NOT_COMPUTED;
        double expected_alt2_reads_frac = EXPECTED_ALT_READS_FREQ_NOT_COMPUTED;
        int max_feasible_alt1_len = NOT_COMPUTED;
        int max_feasible_alt2_len = NOT_COMPUTED;

        bp_consensus_info_t alt_bp1, alt_bp2;
        bp_consensus_info_t ref_bp1, ref_bp2;
        bp_pairs_info_t neutral_bp1_pairs, neutral_bp2_pairs;
        bp_pairs_info_t bp1_stray_pairs, bp2_stray_pairs; // pairs that are discordant and yet do not support the SV

        int oar_bp1_reads = 0, oar_bp2_reads = 0;
        int orr_bp1_reads = 0, orr_bp2_reads = 0;

        // Number of reads assigned away to each destination haplotype, separately by ALT breakpoint.
        std::unordered_map<int, int> oar_bp1_reads_by_hpid;
        std::unordered_map<int, int> oar_bp2_reads_by_hpid;

        // Number of reads attributed to each best-scoring variant within the selected destination haplotype.
        std::unordered_map<std::string, int> oar_bp1_reads_by_vid;
        std::unordered_map<std::string, int> oar_bp2_reads_by_vid;

        // OAR*C and ORR*C sources of truth keyed by suffixed read id; hq and exact derive CHQ and E.
        std::unordered_map<std::string, other_read_info_t> oar_bp1_consistent_reads;
        std::unordered_map<std::string, other_read_info_t> oar_bp2_consistent_reads;
        std::unordered_map<std::string, other_read_info_t> orr_bp1_consistent_reads;
        std::unordered_map<std::string, other_read_info_t> orr_bp2_consistent_reads;

        int alt_ref_equal_reads = 0, alt_ref_equal_reads_highmq = 0;
        int alt_lext_reads = 0, hq_alt_lext_reads = 0, alt_rext_reads = 0, hq_alt_rext_reads = 0;
        consensus_alignment_metrics_t alt_consensus1_metrics, alt_consensus2_metrics;
        consensus_alignment_metrics_t ext_alt_consensus1_metrics, ext_alt_consensus2_metrics;
        int ins_seq_prefix_cov = 0, ins_seq_suffix_cov = 0;
        bool too_deep = false;
        int alt1_hp_len_mode = NOT_COMPUTED, alt1_consistent_hp_len_mode = NOT_COMPUTED;
        double alt1_consistent_hp_len_iqr = NOT_COMPUTED;
        double alt1_hp_5p_mismatch_rate = NOT_COMPUTED, alt1_hp_3p_mismatch_rate = NOT_COMPUTED;
        int alt1_hp_min_mapq = NOT_COMPUTED, alt1_hp_max_mapq = NOT_COMPUTED;
        double alt1_hp_avg_mapq = NOT_COMPUTED, alt1_hp_stddev_mapq = NOT_COMPUTED;
        int ref1_hp_len_mode = NOT_COMPUTED, ref1_consistent_hp_len_mode = NOT_COMPUTED;
        double ref1_consistent_hp_len_iqr = NOT_COMPUTED;
        double ref1_hp_5p_mismatch_rate = NOT_COMPUTED, ref1_hp_3p_mismatch_rate = NOT_COMPUTED;
        int ref1_hp_min_mapq = NOT_COMPUTED, ref1_hp_max_mapq = NOT_COMPUTED;
        double ref1_hp_avg_mapq = NOT_COMPUTED, ref1_hp_stddev_mapq = NOT_COMPUTED;

        int left_flanking_cov = 0, indel_left_cov = 0, indel_right_cov = 0, right_flanking_cov = 0;
        int left_anchor_cov = 0, right_anchor_cov = 0;
        int left_flanking_cov_highmq = 0, indel_left_cov_highmq = 0;
        int indel_right_cov_highmq = 0, right_flanking_cov_highmq = 0;
        int left_anchor_cov_highmq = 0, right_anchor_cov_highmq = 0;

        double alt1_occ_ratio = NOT_COMPUTED, alt2_occ_ratio = NOT_COMPUTED;

        float epr = 0.0;

        std::vector<std::string> filters;

        sample_info_t() {
            gt.push_back(bcf_gt_unphased(1));
        }

        bool is_pass() {
            return filters.empty() || filters[0] == "PASS";
        }
    } sample_info;

    bcf1_t* vcf_entry = NULL;

    sv_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq, 
        std::shared_ptr<consensus_t> rc_consensus, std::shared_ptr<consensus_t> lc_consensus, 
        std::shared_ptr<anchor_aln_t> left_anchor_aln, std::shared_ptr<anchor_aln_t> right_anchor_aln) : 
        chr(chr), start(start), end(end), ins_seq(ins_seq), rc_consensus(rc_consensus), lc_consensus(lc_consensus),
        left_anchor_aln(left_anchor_aln), right_anchor_aln(right_anchor_aln) {
    }

    int rc_fwd_reads() { return rc_consensus ? rc_consensus->fwd_reads : 0; }
    int rc_rev_reads() { return rc_consensus ? rc_consensus->rev_reads : 0; }
    int lc_fwd_reads() { return lc_consensus ? lc_consensus->fwd_reads : 0; }
    int lc_rev_reads() { return lc_consensus ? lc_consensus->rev_reads : 0; }

    bool is_pass() { return sample_info.filters.empty() || sample_info.filters[0] == "PASS"; }

    hts_pos_t start_bp_lower_boundary() {
        if (lc_consensus == NULL) return consensus_t::LOWER_BOUNDARY_NON_CALCULATED;
        return lc_consensus->other_bp_lower_boundary;
    }
    hts_pos_t start_bp_upper_boundary() {
        if (lc_consensus == NULL) return consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
        return lc_consensus->other_bp_upper_boundary;
    }

    hts_pos_t end_bp_lower_boundary() {
        if (rc_consensus == NULL) return consensus_t::LOWER_BOUNDARY_NON_CALCULATED;
        return rc_consensus->other_bp_lower_boundary;
    }
    hts_pos_t end_bp_upper_boundary() {
        if (rc_consensus == NULL) return consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
        return rc_consensus->other_bp_upper_boundary;
    }

    virtual std::string unique_key(bool include_aux = true) {
        std::string key = chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + svtype() + ":" + ins_seq;
        if (!include_aux) return key;
        for (const auto& snp : aux_snps) {
            key += ":" + std::to_string(snp.pos+1) + "," + snp.alt_base;
        }
        for (const auto& sv : aux_indels) {
            key += ":" + sv->unique_key();
        }
        return key;
    }

    virtual std::string svtype() = 0;
    virtual std::string svtype_specific_sort_key() { return ""; }
    virtual hts_pos_t svlen() = 0; // this reflects the SVLEN field in VCF (4.3 and below)
    virtual hts_pos_t svsize() = 0; // this is for filtering, and it is the max number of bases affected either on the ref or in the alt

    std::string left_anchor_aln_string() {
        if (left_anchor_aln == NULL) return "NA";
        return std::to_string(left_anchor_aln->start+1) + "-" + std::to_string(left_anchor_aln->end+1);
    }
    std::string right_anchor_aln_string() {
        if (right_anchor_aln == NULL) return "NA";
        return std::to_string(right_anchor_aln->start+1) + "-" + std::to_string(right_anchor_aln->end+1);
    }

    bool incomplete_ins_seq() { return ins_seq.find("-") != std::string::npos; }

    int known_seq_prefix_len() {
        int d = ins_seq.find("-");
        return d == std::string::npos ? ins_seq.length() : d;
    }
    int known_seq_suffix_len() {
        int d = ins_seq.find("-");
        return d == std::string::npos ? ins_seq.length() : ins_seq.length() - d - 1;
    }

    std::string print_gt() {
        std::stringstream ss;
        for (int i = 0; i < sample_info.gt.size(); i++) {
            if (i > 0) ss << "/";
            ss << (bcf_gt_is_missing(sample_info.gt[i]) ? "." : std::to_string(bcf_gt_allele(sample_info.gt[i])));
        }
        return ss.str();
    }

    int allele_count(int allele) {
        int ac = 0;
        for (int i = 0; i < sample_info.gt.size(); i++) {
            if (!bcf_gt_is_missing(sample_info.gt[i])) ac += (bcf_gt_allele(sample_info.gt[i]) == allele);
        }
        return ac;
    }

    int missing_alleles() {
        int ac = 0;
        for (int i = 0; i < sample_info.gt.size(); i++) {
            ac += bcf_gt_is_missing(sample_info.gt[i]);
        }
        return ac;
    }

    virtual ~sv_t() {
        bcf_destroy1(vcf_entry);
    }
};

inline uint64_t stable_hash64(const std::string& key) {
    uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : key) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    return hash;
}

inline std::string hpid_key(sv_t& sv) {
    std::vector<std::string> atoms;
    atoms.push_back(sv.unique_key(false));
    for (const auto& aux_snp : sv.aux_snps) {
        atoms.push_back(aux_snp.unique_key());
    }
    for (const auto& aux_indel : sv.aux_indels) {
        atoms.push_back(aux_indel->unique_key(false));
    }
    std::sort(atoms.begin(), atoms.end());

    std::stringstream ss;
    for (const auto& atom : atoms) {
        ss << atom << ";";
    }
    return ss.str();
}

inline int make_hpid(sv_t& sv) {
    return 1 + (stable_hash64(hpid_key(sv)) % (uint64_t(INT32_MAX) - 1));
}

inline std::pair<int, std::string> sv_ft_order_key(const std::shared_ptr<sv_t>& sv) {
    if (sv->is_pass()) {
        return std::make_pair(0, std::string());
    }

    std::stringstream ss;
    for (size_t i = 0; i < sv->sample_info.filters.size(); i++) {
        if (i > 0) ss << ";";
        ss << sv->sample_info.filters[i];
    }
    return std::make_pair(1, ss.str());
}

inline bool sv_output_order(const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
    if (a->start != b->start) return a->start < b->start;
    if (a->end != b->end) return a->end < b->end;

    std::string a_svtype = a->svtype();
    std::string b_svtype = b->svtype();
    if (a_svtype != b_svtype) return a_svtype < b_svtype;

    std::string a_svtype_specific_sort_key = a->svtype_specific_sort_key();
    std::string b_svtype_specific_sort_key = b->svtype_specific_sort_key();
    if (a_svtype_specific_sort_key != b_svtype_specific_sort_key) return a_svtype_specific_sort_key < b_svtype_specific_sort_key;

    if (a->ins_seq != b->ins_seq) return a->ins_seq < b->ins_seq;
    if (a->source != b->source) return a->source < b->source;

    std::pair<int, std::string> a_ft_key = sv_ft_order_key(a);
    std::pair<int, std::string> b_ft_key = sv_ft_order_key(b);
    if (a_ft_key != b_ft_key) return a_ft_key < b_ft_key;

    return hpid_key(*a) < hpid_key(*b);
}

struct deletion_t : sv_t {

    bool remapped = false;
    std::string original_range;

    using sv_t::sv_t;

    std::string svtype() { return "DEL"; }
    hts_pos_t svlen() { return start - end + ins_seq.length(); }
    hts_pos_t svsize() { return end - start; }
};

struct duplication_t : sv_t {
    using sv_t::sv_t;
    double ins_to_dup_similarity = 0.0;
    bool cn_unresolved = false;

    std::string svtype() { return "DUP"; }
    hts_pos_t svlen() { return end - start + ins_seq.length(); }
    hts_pos_t svsize() { return end - start + ins_seq.length(); }
};

struct insertion_t : sv_t {
    using sv_t::sv_t;

    std::string svtype() { return "INS"; }
    hts_pos_t svlen() { return ins_seq.length() - (end-start); }
    hts_pos_t svsize() { return ins_seq.length(); }
};

struct replacement_t : sv_t {
    using sv_t::sv_t;

    std::string svtype() { return "RPL"; }
    hts_pos_t svlen() { return 0; }
    hts_pos_t svsize() { return end - start; }
};

struct breakend_t : sv_t {
    bool left_facing;
    std::shared_ptr<consensus_t> leftmost_consensus, rightmost_consensus;
    
    breakend_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq, 
        std::shared_ptr<anchor_aln_t> left_anchor_aln, std::shared_ptr<anchor_aln_t> right_anchor_aln, 
        bool left_facing,
        std::shared_ptr<consensus_t> leftmost_consensus = nullptr, std::shared_ptr<consensus_t> rightmost_consensus = nullptr) :
        sv_t(chr, start, end, ins_seq, nullptr, nullptr, left_anchor_aln, right_anchor_aln),
        left_facing(left_facing), leftmost_consensus(leftmost_consensus), rightmost_consensus(rightmost_consensus) {}

    std::string unique_key(bool include_aux = true) override {
        return sv_t::unique_key(include_aux) + ":" + (left_facing ? "LF" : "RF");
    }

    std::string svtype() override { return "BND"; }
    hts_pos_t svlen() override { return 0; }
    hts_pos_t svsize() override { return end - start; }
};

struct inversion_t : sv_t {

    std::shared_ptr<anchor_aln_t> rbp_left_anchor_aln, rbp_right_anchor_aln;

    hts_pos_t inv_start = 0, inv_end = 0;

    inversion_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq, 
        std::shared_ptr<consensus_t> rc_consensus, std::shared_ptr<consensus_t> lc_consensus,
        std::shared_ptr<anchor_aln_t> lbp_left_anchor_aln, std::shared_ptr<anchor_aln_t> lbp_right_anchor_aln, 
        std::shared_ptr<anchor_aln_t> rbp_left_anchor_aln, std::shared_ptr<anchor_aln_t> rbp_right_anchor_aln) :
    sv_t(chr, start, end, ins_seq, rc_consensus, lc_consensus, lbp_left_anchor_aln, lbp_right_anchor_aln),
    rbp_left_anchor_aln(rbp_left_anchor_aln), rbp_right_anchor_aln(rbp_right_anchor_aln) {
        inv_start = start;
        inv_end = end;
    }

    std::string unique_key(bool include_aux = true) override {
        return sv_t::unique_key(include_aux) + ":INVPOS:" + std::to_string(inv_start) + ":" + std::to_string(inv_end);
    }

    std::string svtype_specific_sort_key() override {
        std::stringstream ss;
        ss << std::setw(20) << std::setfill('0') << inv_start << ":" << std::setw(20) << std::setfill('0') << inv_end;
        return ss.str();
    }

    std::string svtype() override { return "INV"; }
    hts_pos_t svlen() override {
        if (!ins_seq.empty()) {
            return ins_seq.length() - (end-start);
        }
        return (inv_end-inv_start) - (end-start);
    }
    hts_pos_t svsize() override { return end - start; }

    bool is_left_facing() {
        return source[source.length()-2] == 'L';
    }
};

#endif /* TYPES_H */
