#ifndef TYPES_H
#define TYPES_H

#include "htslib/hts.h"
#include "utils.h"
#include <cstdint>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include "htslib/vcf.h"

struct consensus_t {
    bool left_clipped;
    hts_pos_t start, breakpoint, end;
    std::string sequence;
    int fwd_reads, rev_reads;
    uint8_t max_mapq;
    hts_pos_t remap_boundary;
    int clip_len, lowq_clip_portion;
    int left_ext_reads = 0, right_ext_reads = 0, hq_left_ext_reads = 0, hq_right_ext_reads = 0;
    bool is_hsr = false;
	bool extended_to_left = false, extended_to_right = false;

    static const int LOWER_BOUNDARY_NON_CALCULATED = 0, UPPER_BOUNDARY_NON_CALCULATED = INT32_MAX;
    static const int UNKNOWN_CLIP_LEN = INT16_MAX;

    consensus_t(bool left_clipped, hts_pos_t start, hts_pos_t breakpoint, hts_pos_t end,
                const std::string& sequence, int fwd_reads, int rev_reads, int clip_len, uint8_t max_mapq, 
                hts_pos_t remap_boundary, int lowq_clip_portion)
                : left_clipped(left_clipped), start(start), breakpoint(breakpoint), end(end),
                sequence(sequence), fwd_reads(fwd_reads), rev_reads(rev_reads), clip_len(clip_len), max_mapq(max_mapq), 
                remap_boundary(remap_boundary), lowq_clip_portion(lowq_clip_portion) {}
    
    consensus_t(std::string& line, bool is_hsr) : is_hsr(is_hsr) {
        std::stringstream ss(line);
        char dir;
        int max_mapq_int;
        ss >> start >> end >> breakpoint >> dir >> sequence >> fwd_reads >> rev_reads >> max_mapq_int >> remap_boundary >> lowq_clip_portion;
        left_clipped = dir == 'L';
        max_mapq = (uint8_t) max_mapq_int;
        if (is_hsr) {
            clip_len = UNKNOWN_CLIP_LEN;
        } else {
            clip_len = left_clipped ? breakpoint - start : end - breakpoint;
        }
    }

    std::string to_string() {
        std::stringstream ss;
        ss << start << " " << end << " " << breakpoint << (left_clipped ? " L " : " R ") << sequence << " ";
        ss << fwd_reads << " " << rev_reads << " " << (int)max_mapq << " " << remap_boundary << " " << lowq_clip_portion;
        return ss.str();
    }

    int reads() { return fwd_reads + rev_reads; }

    hts_pos_t left_ext_target_start(int max_is, int read_len) {
    	if (!left_clipped) {
    		return start - max_is + read_len;
    	} else {
    		if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint - max_is - 2*sequence.length();
			} else {
				return remap_boundary;
			}
    	}
    }
    hts_pos_t left_ext_target_end(int max_is, int read_len) {
		if (!left_clipped) {
			return start;
		} else {
			if (remap_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
				return breakpoint + max_is + 2*sequence.length();
			} else {
				return remap_boundary + max_is;
			}
		}
	}

    hts_pos_t right_ext_target_start(int max_is, int read_len) {
		if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint - max_is - 2*sequence.length();
			} else {
				return remap_boundary - max_is;
			}
		} else {
			return end;
		}
	}
    hts_pos_t right_ext_target_end(int max_is, int read_len) {
    	if (!left_clipped) {
			if (remap_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
				return breakpoint + max_is + 2*sequence.length();
			} else {
				return remap_boundary;
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

    int anchor_len() { return sequence.length() - clip_len; }

    std::string clip_sequence() {
        if (left_clipped) {
            return sequence.substr(0, clip_len);
        } else {
            return sequence.substr(sequence.length()-clip_len);
        }
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

    std::string source;
    bool imprecise = false;

    static constexpr const double KS_PVAL_NOT_COMPUTED = -1.0;
    static const int SIZE_NOT_COMPUTED = INT32_MAX;
    
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
        int fwd_mate_cov_bps = 0, rev_mate_cov_bps = 0;
        int fwd_hq_mate_cov_bps = 0, rev_hq_mate_cov_bps = 0;

        int consistent_reads() { return consistent_fwd + consistent_rev; }
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

    struct sample_info_t {
        static const int NOT_COMPUTED = -1;

        int* gt;

        bp_consensus_info_t alt_bp1, alt_bp2;
        bp_consensus_info_t ref_bp1, ref_bp2;
        bp_pairs_info_t neutral_bp1_pairs, neutral_bp2_pairs;
        bp_pairs_info_t bp1_stray_pairs, bp2_stray_pairs; // pairs that are discordant and yet do not support the SV

        int assigned_to_other_sv_bp1_reads = 0, assigned_to_other_sv_bp1_consistent = 0, assigned_to_other_sv_bp1_consistent_highmq = 0;
        int assigned_to_other_sv_bp2_reads = 0, assigned_to_other_sv_bp2_consistent = 0, assigned_to_other_sv_bp2_consistent_highmq = 0;
        int alt_ref_equal_reads = 0, alt_ref_equal_reads_highmq = 0;
        int alt_lext_reads = 0, hq_alt_lext_reads = 0, alt_rext_reads = 0, hq_alt_rext_reads = 0;
        int ext_alt_consensus1_length = 0, ext_alt_consensus2_length = 0;
        int ext_alt_consensus1_to_alt_score = 0, ext_alt_consensus1_to_ref_score = 0;
        int ext_alt_consensus2_to_alt_score = 0, ext_alt_consensus2_to_ref_score = 0;
        int alt_consensus1_split_size1 = 0, alt_consensus1_split_size2 = 0;
        int alt_consensus2_split_size1 = 0, alt_consensus2_split_size2 = 0;
        int alt_consensus1_split_score1 = 0, alt_consensus1_split_score2 = 0;
        int alt_consensus1_split_score1_ind_aln = 0, alt_consensus1_split_score2_ind_aln = 0;
        int alt_consensus2_split_score1 = 0, alt_consensus2_split_score2 = 0;
        int alt_consensus2_split_score1_ind_aln = 0, alt_consensus2_split_score2_ind_aln = 0;
        int ins_seq_prefix_cov = 0, ins_seq_suffix_cov = 0;
        bool too_deep = false;

        int left_flanking_cov = 0, indel_left_cov = 0, indel_right_cov = 0, right_flanking_cov = 0;
        int left_anchor_cov = 0, right_anchor_cov = 0;
        int left_flanking_cov_highmq = 0, indel_left_cov_highmq = 0;
        int indel_right_cov_highmq = 0, right_flanking_cov_highmq = 0;
        int left_anchor_cov_highmq = 0, right_anchor_cov_highmq = 0;

        float epr = 0.0;

        std::vector<std::string> filters;

        sample_info_t() {
            gt = (int*)malloc(sizeof(int));
            gt[0] = bcf_gt_unphased(1);
        }

        ~sample_info_t() {
            free(gt);
        }

        bool is_pass() {
            return filters.empty() || filters[0] == "PASS";
        }
    } sample_info;

    int n_gt = 1;
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

    hts_pos_t remap_boundary_upper() {
        if (rc_consensus == NULL) return consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
        return rc_consensus->remap_boundary;
    }

    hts_pos_t remap_boundary_lower() {
        if (lc_consensus == NULL) return consensus_t::LOWER_BOUNDARY_NON_CALCULATED;
        return lc_consensus->remap_boundary;
    }

    std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + svtype() + ":" + ins_seq;
    }

    virtual std::string svtype() = 0;
    virtual hts_pos_t svlen() = 0;

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
        for (int i = 0; i < n_gt; i++) {
            if (i > 0) ss << "/";
            ss << (bcf_gt_is_missing(sample_info.gt[i]) ? "." : std::to_string(bcf_gt_allele(sample_info.gt[i])));
        }
        return ss.str();
    }

    int allele_count(int allele) {
        int ac = 0;
        for (int i = 0; i < n_gt; i++) {
            if (!bcf_gt_is_missing(sample_info.gt[i])) ac += (bcf_gt_allele(sample_info.gt[i]) == allele);
        }
        return ac;
    }

    int missing_alleles() {
        int ac = 0;
        for (int i = 0; i < n_gt; i++) {
            ac += bcf_gt_is_missing(sample_info.gt[i]);
        }
        return ac;
    }

    virtual ~sv_t() {
        bcf_destroy1(vcf_entry);
    }
};

struct deletion_t : sv_t {

    bool remapped = false;
    std::string original_range;

    using sv_t::sv_t;

    std::string svtype() { return "DEL"; }
    hts_pos_t svlen() { return start - end + ins_seq.length(); }
};

struct duplication_t : sv_t {
    using sv_t::sv_t;
    double ins_to_dup_similarity = 0.0;

    std::string svtype() { return "DUP"; }
    hts_pos_t svlen() { return end - start + ins_seq.length(); }
};

struct insertion_t : sv_t {
    using sv_t::sv_t;

    std::string svtype() { return "INS"; }
    hts_pos_t svlen() { return ins_seq.length() - (end-start); }
};

struct breakend_t : sv_t {
    bool left_facing;
    
    breakend_t(std::string chr, hts_pos_t start, hts_pos_t end, std::string ins_seq, 
        std::shared_ptr<consensus_t> rc_consensus, std::shared_ptr<consensus_t> lc_consensus, 
        std::shared_ptr<anchor_aln_t> left_anchor_aln, std::shared_ptr<anchor_aln_t> right_anchor_aln, 
        bool left_facing) : sv_t(chr, start, end, ins_seq, rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln), left_facing(left_facing) {}

    std::string svtype() { return "BND"; }
    hts_pos_t svlen() { return 0; }
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

    std::string svtype() { return "INV"; }
    hts_pos_t svlen() { 
        if (!ins_seq.empty()) {
            return ins_seq.length() - (end-start);
        }
        return (inv_end-inv_start) - (end-start);
    }

    bool is_left_facing() {
        return source[source.length()-2] == 'L';
    }
};

#endif /* TYPES_H */
