#ifndef HP_READ_INFO_H
#define HP_READ_INFO_H

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#include "htslib/sam.h"

#include "../libs/ssw_cpp.h"
#include "genotype.h"
#include "sw_utils.h"

constexpr int UNDEFINED_HP_LEN = -1;

struct hp_read_info_t {
    int hp_len;
    int tail_5p_len, tail_3p_len;
    int aligned_5p_tail_len = 0;
    int tail_5p_mismatches, tail_3p_mismatches;
    bp_support_read_t read;
    bool original_alignment_has_indel_outside_hp = true;
    bool hp_deletion_extends_outside_hp = false;
    bool hp_insertion_has_non_hp_bases = false;
    bool hp_run_extends_into_5p_tail = false;
    bool hp_run_extends_into_3p_tail = false;
    bool ref_hp_has_non_hp_read_base = false;
    bool rescued = false;

    hp_read_info_t(int hp_len = 0, int tail_5p_len = 0, int tail_3p_len = 0,
        int tail_5p_mismatches = 0, int tail_3p_mismatches = 0,
        bp_support_read_t read = bp_support_read_t(), bool original_alignment_has_indel_outside_hp = true,
        bool hp_deletion_extends_outside_hp = false, bool rescued = false, bool hp_insertion_has_non_hp_bases = false) :
        hp_len(hp_len), tail_5p_len(tail_5p_len), tail_3p_len(tail_3p_len),
        tail_5p_mismatches(tail_5p_mismatches), tail_3p_mismatches(tail_3p_mismatches),
        read(read), original_alignment_has_indel_outside_hp(original_alignment_has_indel_outside_hp),
        hp_deletion_extends_outside_hp(hp_deletion_extends_outside_hp),
        hp_insertion_has_non_hp_bases(hp_insertion_has_non_hp_bases), rescued(rescued) {}

    bool is_good_read(int min_tail_len, double max_mismatch_rate) const {
        return tail_5p_len >= min_tail_len && double(tail_5p_mismatches)/tail_5p_len <= max_mismatch_rate &&
               tail_3p_len >= min_tail_len && double(tail_3p_mismatches)/tail_3p_len <= max_mismatch_rate;
    }
};

// TODO: replace with the fast SIMD version in sw_utils.h
// For clipped tails, place the tail at its best ungapped offset in the local reference window.
int best_ungapped_mismatch_count(const std::string& query, char* ref, hts_pos_t ref_len) {
    if (query.empty() || ref_len < (hts_pos_t) query.length()) {
        return query.length();
    }

    int best_mismatches = query.length();
    for (hts_pos_t i = 0; i + query.length() <= ref_len; i++) {
        int mismatches = 0;
        for (int j = 0; j < query.length(); j++) {
            if (query[j] != toupper(ref[i + j])) {
                mismatches++;
            }
        }
        if (mismatches < best_mismatches) {
            best_mismatches = mismatches;
            if (best_mismatches == 0) break;
        }
    }
    return best_mismatches;
}

struct hp_overflow_resolution_t {
    int hp_len, left, right, left_tail_mismatches, right_tail_mismatches;

    hp_overflow_resolution_t(int hp_len = UNDEFINED_HP_LEN, int left = 0, int right = 0,
        int left_tail_mismatches = 0, int right_tail_mismatches = 0) :
        hp_len(hp_len), left(left), right(right),
        left_tail_mismatches(left_tail_mismatches), right_tail_mismatches(right_tail_mismatches) {}
};

hp_overflow_resolution_t resolve_hp_overflow(const std::string& read_seq, int observed_hp_len,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len,
    bool is_rev, double max_5p_mismatch_rate) {

    int ref_hp_len = ref_hp_range.end - ref_hp_range.beg;
    hts_pos_t extend = read_seq.length() - 1;
    hts_pos_t haplotype_beg = std::max(hts_pos_t(0), ref_hp_range.beg - extend);
    hts_pos_t haplotype_end = std::min(contig_len, ref_hp_range.end + extend);
    std::string ref_haplotype(contig_seq + haplotype_beg, haplotype_end - haplotype_beg);
    hts_pos_t hp_beg = ref_hp_range.beg - haplotype_beg;
    hts_pos_t hp_end = ref_hp_range.end - haplotype_beg;

    std::string left_flank = ref_haplotype.substr(0, hp_beg);
    std::string right_flank = ref_haplotype.substr(hp_end);
    hp_overflow_resolution_t best_resolution;
    int best_score = 0;
    bool best_score_is_tied = false;
    for (int hp_len = 0; hp_len <= std::max(ref_hp_len, observed_hp_len); hp_len++) {
        std::string haplotype = left_flank + std::string(hp_len, hp_base) + right_flank;
        ungapped_aln_t aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(), haplotype.c_str(), haplotype.length(), 0);
        if (aln.query_begin != 0 || aln.query_end != (int) read_seq.length()) continue;
        int candidate_hp_beg = left_flank.length();
        int candidate_hp_end = candidate_hp_beg + hp_len;
        if (aln.ref_begin > candidate_hp_beg || aln.ref_end < candidate_hp_end) continue;

        int query_hp_beg = aln.query_begin + candidate_hp_beg - aln.ref_begin;
        int query_hp_end = query_hp_beg + hp_len;
        if (query_hp_beg < 0 || query_hp_end > (int) read_seq.length()) continue;

        if (best_resolution.hp_len == UNDEFINED_HP_LEN || aln.score > best_score) {
            int left_tail_mismatches = 0, right_tail_mismatches = 0;
            for (int qpos = 0; qpos < query_hp_beg; qpos++) {
                int ref_pos = aln.ref_begin + qpos;
                if (toupper(read_seq[qpos]) != toupper(haplotype[ref_pos])) left_tail_mismatches++;
            }
            for (int qpos = query_hp_end; qpos < (int) read_seq.length(); qpos++) {
                int ref_pos = aln.ref_begin + qpos;
                if (toupper(read_seq[qpos]) != toupper(haplotype[ref_pos])) right_tail_mismatches++;
            }
            best_resolution = hp_overflow_resolution_t(hp_len, query_hp_beg, query_hp_end,
                left_tail_mismatches, right_tail_mismatches);
            best_score = aln.score;
            best_score_is_tied = false;
        } else if (aln.score == best_score) {
            best_score_is_tied = true;
        }
    }
    if (best_score_is_tied || best_resolution.hp_len == UNDEFINED_HP_LEN) return hp_overflow_resolution_t();
    int tail_5p_len = is_rev ? read_seq.length() - best_resolution.right : best_resolution.left;
    int tail_5p_mismatches = is_rev ? best_resolution.right_tail_mismatches : best_resolution.left_tail_mismatches;
    return tail_5p_len > 0 && double(tail_5p_mismatches) / tail_5p_len > max_5p_mismatch_rate ?
        hp_overflow_resolution_t() : best_resolution;
}

// If the tail is unclipped, calculate mismatches by simply counting mismatches in its alignment
// If the tail is clipped, force the whole tail to align next to the HP ref region (allowing some extra leeway)
int tail_mismatch_count_simple(bam1_t* read, const std::string& read_seq, const std::vector<hts_pos_t>& qpos_to_rpos,
    char* contig_seq, hts_pos_t contig_len, int q_lo, int q_hi, bool left_side, int leeway, hts_pos_t ref_boundary) {

    if (q_hi <= q_lo || read_seq.empty()) return 0;

    int left_clip = get_left_clip_size(read), right_clip = get_right_clip_size(read);
    if (left_side && left_clip > 0 && q_lo == 0) {
        hts_pos_t ref_hi = std::min(ref_boundary, contig_len);
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, ref_hi - (q_hi - q_lo) - leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }
    if (!left_side && right_clip > 0 && q_hi == (int) read_seq.length()) {
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, std::min(ref_boundary, contig_len));
        hts_pos_t ref_hi = std::min(contig_len, ref_lo + (q_hi - q_lo) + leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }

    int mismatches = 0;
    for (int qpos = q_lo; qpos < q_hi; qpos++) {
        hts_pos_t rpos = qpos_to_rpos[qpos];
        if (rpos == -1 || rpos >= contig_len) continue;
        if (toupper(read_seq[qpos]) != toupper(contig_seq[rpos])) {
            mismatches++;
        }
    }
    return mismatches;
}

// Shared tail-mismatch helper for BAM-backed and SSW-backed HP interpretation.
int tail_mismatch_count_from_mapping(const std::string& read_seq, const std::vector<hts_pos_t>& qpos_to_rpos,
    char* contig_seq, hts_pos_t contig_len, int q_lo, int q_hi, bool left_side,
    bool left_clipped, bool right_clipped, int leeway, hts_pos_t ref_boundary) {

    if (q_hi <= q_lo || read_seq.empty()) return 0;

    if (left_side && left_clipped && q_lo == 0) {
        hts_pos_t ref_hi = std::min(ref_boundary, contig_len);
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, ref_hi - (q_hi - q_lo) - leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }
    if (!left_side && right_clipped && q_hi == (int) read_seq.length()) {
        hts_pos_t ref_lo = std::max((hts_pos_t) 0, std::min(ref_boundary, contig_len));
        hts_pos_t ref_hi = std::min(contig_len, ref_lo + (q_hi - q_lo) + leeway);
        return best_ungapped_mismatch_count(read_seq.substr(q_lo, q_hi - q_lo), contig_seq+ref_lo, ref_hi-ref_lo);
    }

    int mismatches = 0;
    for (int qpos = q_lo; qpos < q_hi; qpos++) {
        hts_pos_t rpos = qpos_to_rpos[qpos];
        if (rpos == -1 || rpos >= contig_len) continue;
        if (toupper(read_seq[qpos]) != toupper(contig_seq[rpos])) {
            mismatches++;
        }
    }
    return mismatches;
}

struct hp_cigar_op_t {
    char op;
    hts_pos_t len;

    hp_cigar_op_t(char op, hts_pos_t len) : op(op), len(len) {}
};

std::vector<hp_cigar_op_t> normalized_cigar(const StripedSmithWaterman::Alignment& aln) {
    std::vector<hp_cigar_op_t> cigar;
    cigar.reserve(aln.cigar.size());
    for (uint32_t cigar_op : aln.cigar) {
        cigar.emplace_back(cigar_int_to_op(cigar_op), cigar_int_to_len(cigar_op));
    }
    return cigar;
}

std::vector<hp_cigar_op_t> normalized_cigar(bam1_t* read) {
    std::vector<hp_cigar_op_t> cigar;
    cigar.reserve(read->core.n_cigar);
    uint32_t* bam_cigar = bam_get_cigar(read);
    for (uint32_t i = 0; i < read->core.n_cigar; i++) {
        cigar.emplace_back(bam_cigar_opchr(bam_cigar[i]), bam_cigar_oplen(bam_cigar[i]));
    }
    return cigar;
}

bool insertion_has_non_hp_bases(const std::string& read_seq, int qpos, int len, char hp_base) {
    if (qpos < 0 || qpos + len > (int) read_seq.length()) return true;
    for (int i = 0; i < len; i++) {
        if (toupper(read_seq[qpos + i]) != hp_base) return true;
    }
    return false;
}

struct hp_alignment_summary_t {
    std::vector<hts_pos_t> qpos_to_rpos;
    std::vector<int> anchors;
    hts_pos_t last_qpos_before_ref_hp = -1;
    hts_pos_t first_qpos_after_ref_hp = -1;
    bool hp_insertion_has_non_hp_bases = false;
};

hp_alignment_summary_t summarize_hp_alignment(const std::vector<hp_cigar_op_t>& cigar,
    hts_pos_t ref_begin, const std::string& read_seq, hts_pair_pos_t ref_hp_range, char hp_base) {

    hp_alignment_summary_t summary;
    summary.qpos_to_rpos.assign(read_seq.length(), -1);
    int qpos = 0;
    hts_pos_t rpos = ref_begin;
    for (const hp_cigar_op_t& cigar_op : cigar) {
        if (cigar_op.op == 'M' || cigar_op.op == '=' || cigar_op.op == 'X') {
            for (int j = 0; j < cigar_op.len; j++) {
                summary.qpos_to_rpos[qpos] = rpos;
                if (rpos < ref_hp_range.beg) {
                    summary.last_qpos_before_ref_hp = qpos;
                } else if (rpos >= ref_hp_range.end && summary.first_qpos_after_ref_hp == -1) {
                    summary.first_qpos_after_ref_hp = qpos;
                } else if (ref_hp_range.beg <= rpos && rpos < ref_hp_range.end) {
                    summary.anchors.push_back(qpos);
                }
                qpos++;
                rpos++;
            }
        } else if (cigar_op.op == 'I') {
            if (ref_hp_range.beg <= rpos && rpos <= ref_hp_range.end &&
                insertion_has_non_hp_bases(read_seq, qpos, cigar_op.len, hp_base)) {
                summary.hp_insertion_has_non_hp_bases = true;
            }
            if (rpos < ref_hp_range.beg) {
                summary.last_qpos_before_ref_hp = qpos + cigar_op.len - 1;
            } else if (rpos >= ref_hp_range.end && summary.first_qpos_after_ref_hp == -1) {
                summary.first_qpos_after_ref_hp =
                    (ref_hp_range.beg == ref_hp_range.end && rpos == ref_hp_range.beg) ?
                    qpos + cigar_op.len : qpos;
            }
            qpos += cigar_op.len;
        } else if (cigar_op.op == 'S') {
            qpos += cigar_op.len;
        } else if (cigar_op.op == 'D' || cigar_op.op == 'N') {
            rpos += cigar_op.len;
        }
    }
    return summary;
}

struct hp_side_indel_info_t {
    int aligned_len = 0;
    int indel_len = 0;
};

struct hp_adjacent_indel_info_t {
    hp_side_indel_info_t left;
    hp_side_indel_info_t right;
};

bool five_p_evidence_permits_iterative_hp_len_estimation(int reads, int indel_reads) {
    return reads > 0 && 10 * indel_reads < reads;
}

hp_adjacent_indel_info_t get_adjacent_indel_info(const std::vector<hp_cigar_op_t>& cigar,
    hts_pos_t ref_begin, hts_pair_pos_t hp_range, int adjacency_window) {

    hp_adjacent_indel_info_t info;
    hts_pos_t rpos = ref_begin;
    for (const hp_cigar_op_t& cigar_op : cigar) {
        if (cigar_op.op == 'I') {
            if (rpos < hp_range.beg && rpos >= hp_range.beg - adjacency_window) info.left.indel_len = cigar_op.len;
            if (rpos > hp_range.end && rpos <= hp_range.end + adjacency_window && info.right.indel_len == 0) info.right.indel_len = cigar_op.len;
        } else if (cigar_op.op == 'D') {
            hts_pos_t deletion_end = rpos + cigar_op.len;
            if (deletion_end <= hp_range.beg && deletion_end >= hp_range.beg - adjacency_window) info.left.indel_len = -cigar_op.len;
            if (rpos >= hp_range.end && rpos <= hp_range.end + adjacency_window && info.right.indel_len == 0) info.right.indel_len = -cigar_op.len;
            rpos = deletion_end;
        } else if (cigar_op.op == 'N') {
            rpos += cigar_op.len;
        } else if (cigar_op.op == 'M' || cigar_op.op == '=' || cigar_op.op == 'X') {
            info.left.aligned_len += std::max<hts_pos_t>(0, std::min(rpos + cigar_op.len, hp_range.beg) - rpos);
            info.right.aligned_len += std::max<hts_pos_t>(0, rpos + cigar_op.len - std::max(rpos, hp_range.end));
            rpos += cigar_op.len;
        }
    }
    if (rpos < hp_range.beg) info.left.aligned_len = 0;
    if (ref_begin > hp_range.end) info.right.aligned_len = 0;
    return info;
}

hp_adjacent_indel_info_t get_adjacent_indel_info(bam1_t* read, hts_pair_pos_t hp_range, int adjacency_window) {
    if (read == NULL) return hp_adjacent_indel_info_t();
    return get_adjacent_indel_info(normalized_cigar(read), read->core.pos, hp_range, adjacency_window);
}

bool cigar_has_indel_outside_hp(const std::vector<hp_cigar_op_t>& cigar, hts_pos_t ref_begin, hts_pair_pos_t hp_range) {
    hts_pos_t rpos = ref_begin;

    for (const hp_cigar_op_t& cigar_op : cigar) {
        if (cigar_op.op == 'I') {
            if (rpos < hp_range.beg || rpos > hp_range.end) {
                return true;
            }
        } else if (cigar_op.op == 'D') {
            hts_pos_t deletion_end = rpos + cigar_op.len;
            bool overlaps_hp = rpos < hp_range.end && deletion_end > hp_range.beg;
            if (!overlaps_hp && (rpos < hp_range.beg || deletion_end > hp_range.end)) {
                return true;
            }
            rpos += cigar_op.len;
        } else if (cigar_op.op == 'N') {
            if (rpos < hp_range.beg || rpos + cigar_op.len > hp_range.end) return true;
            rpos += cigar_op.len;
        } else if (cigar_op.op == 'M' || cigar_op.op == '=' || cigar_op.op == 'X') {
            rpos += cigar_op.len;
        }
    }

    return false;
}

bool cigar_has_hp_deletion_extending_outside_hp(const std::vector<hp_cigar_op_t>& cigar, hts_pos_t ref_begin, hts_pair_pos_t hp_range) {
    hts_pos_t rpos = ref_begin;

    for (const hp_cigar_op_t& cigar_op : cigar) {
        if (cigar_op.op == 'D') {
            hts_pos_t deletion_end = rpos + cigar_op.len;
            bool overlaps_hp = rpos < hp_range.end && deletion_end > hp_range.beg;
            if (overlaps_hp && (rpos < hp_range.beg || deletion_end > hp_range.end)) return true;
            rpos = deletion_end;
        } else if (cigar_op.op == 'N') {
            rpos += cigar_op.len;
        } else if (cigar_op.op == 'M' || cigar_op.op == '=' || cigar_op.op == 'X') {
            rpos += cigar_op.len;
        }
    }
    return false;
}

bool has_unexplained_indel_outside_hp(const hp_read_info_t& hp_read_info, bool allele_has_aux_indels,
    const char* allele_seq, hts_pos_t allele_len, hts_pair_pos_t allele_hp_range,
    StripedSmithWaterman::Aligner& aligner) {

    if (hp_read_info.hp_insertion_has_non_hp_bases) return true;

    // Fast path: without aux indels, the original reference alignment fully
    // determines whether there is an unexplained indel outside the allele HP.
    if (!allele_has_aux_indels) {
        return hp_read_info.original_alignment_has_indel_outside_hp || hp_read_info.hp_deletion_extends_outside_hp;
    }

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter_default;
    aligner.Align(hp_read_info.read.seq.c_str(), allele_seq, allele_len, filter_default, &aln, 0);
    std::vector<hp_cigar_op_t> cigar = normalized_cigar(aln);
    return cigar_has_indel_outside_hp(cigar, aln.ref_begin, allele_hp_range) ||
        cigar_has_hp_deletion_extending_outside_hp(cigar, aln.ref_begin, allele_hp_range);
}

// Shared interpretation core. Callers are responsible for providing the
// query-to-reference mapping summary in the coordinate space of the sequence
// used for HP evaluation.
hp_read_info_t calculate_hp_read_info_core(const std::string& read_seq, const std::vector<hts_pos_t>& qpos_to_rpos,
    const std::vector<int>& anchors, hts_pos_t last_qpos_before_ref_hp, hts_pos_t first_qpos_after_ref_hp,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len,
    bool is_rev, bool left_clipped, bool right_clipped, bp_support_read_t read, int tail_align_leeway,
    int force_resolution_max_hp_len = UNDEFINED_HP_LEN, bool can_use_iterative_hp_len_estimation = false,
    double max_5p_mismatch_rate = 1) {

    if (read_seq.empty()) {
        return hp_read_info_t();
    }

    int aligned_left_tail_len = 0, aligned_right_tail_len = 0;
    for (hts_pos_t rpos : qpos_to_rpos) {
        if (rpos != -1 && rpos < ref_hp_range.beg) aligned_left_tail_len++;
        if (rpos >= ref_hp_range.end) aligned_right_tail_len++;
    }
    int aligned_5p_tail_len = is_rev ? aligned_right_tail_len : aligned_left_tail_len;

    int left, right;
    if (anchors.empty()) {
        if (last_qpos_before_ref_hp == -1 && first_qpos_after_ref_hp == -1) return hp_read_info_t();
        left = last_qpos_before_ref_hp != -1 ? last_qpos_before_ref_hp + 1 : first_qpos_after_ref_hp;
        right = first_qpos_after_ref_hp != -1 ? first_qpos_after_ref_hp : left;
    } else {
        // Start from the aligned bases inside the reference HP and expand through adjacent query HP bases.
        left = anchors.front();
        while (left > 0 && read_seq[left-1] == hp_base) {
            hts_pos_t mapped_rpos = qpos_to_rpos[left-1];
            if (mapped_rpos != -1 && mapped_rpos < ref_hp_range.beg) break;
            left--;
        }

        right = anchors.back() + 1;
        while (right < (int) read_seq.length() && read_seq[right] == hp_base) {
            hts_pos_t mapped_rpos = qpos_to_rpos[right];
            if (mapped_rpos != -1 && mapped_rpos >= ref_hp_range.end) break;
            right++;
        }
    }

    bool hp_run_extends_into_5p_tail = false, hp_run_extends_into_3p_tail = false;
    if (!anchors.empty()) {
        int adjacent_5p_qpos = is_rev ? right : left - 1;
        if (adjacent_5p_qpos >= 0 && adjacent_5p_qpos < (int) read_seq.length() && read_seq[adjacent_5p_qpos] == hp_base) {
            hts_pos_t mapped_rpos = qpos_to_rpos[adjacent_5p_qpos];
            bool maps_outside_ref_hp = mapped_rpos < ref_hp_range.beg || mapped_rpos >= ref_hp_range.end;
            bool maps_to_valid_ref_base = mapped_rpos >= 0 && mapped_rpos < contig_len;
            hp_run_extends_into_5p_tail = maps_outside_ref_hp && maps_to_valid_ref_base && toupper(contig_seq[mapped_rpos]) != hp_base;
        }

        int adjacent_3p_qpos = is_rev ? left - 1 : right;
        if (adjacent_3p_qpos >= 0 && adjacent_3p_qpos < (int) read_seq.length() && read_seq[adjacent_3p_qpos] == hp_base) {
            hts_pos_t mapped_rpos = qpos_to_rpos[adjacent_3p_qpos];
            bool maps_outside_ref_hp = mapped_rpos < ref_hp_range.beg || mapped_rpos >= ref_hp_range.end;
            bool maps_to_valid_ref_base = mapped_rpos >= 0 && mapped_rpos < contig_len;
            hp_run_extends_into_3p_tail = maps_outside_ref_hp && maps_to_valid_ref_base && toupper(contig_seq[mapped_rpos]) != hp_base;
        }
    }

    int resolved_hp_len = right - left;
    if (can_use_iterative_hp_len_estimation && hp_run_extends_into_5p_tail) {
        if (is_rev) {
            while (right < (int) read_seq.length() && read_seq[right] == hp_base) right++;
        } else {
            while (left > 0 && read_seq[left-1] == hp_base) left--;
        }
    }

    if (can_use_iterative_hp_len_estimation && hp_run_extends_into_3p_tail) {
        if (is_rev) {
            while (left > 0 && read_seq[left-1] == hp_base) left--;
        } else {
            while (right < (int) read_seq.length() && read_seq[right] == hp_base) right++;
        }
    }

    hp_overflow_resolution_t resolution;
    if (can_use_iterative_hp_len_estimation && (hp_run_extends_into_5p_tail || hp_run_extends_into_3p_tail || force_resolution_max_hp_len != UNDEFINED_HP_LEN)) {
        int observed_hp_len = std::max(right - left, force_resolution_max_hp_len);
        resolution = resolve_hp_overflow(read_seq, observed_hp_len, ref_hp_range, hp_base, contig_seq, contig_len, is_rev, max_5p_mismatch_rate);
        resolved_hp_len = resolution.hp_len;
        if (resolved_hp_len != UNDEFINED_HP_LEN) {
            left = resolution.left;
            right = resolution.right;
        }
    }

    int left_len = left;
    int right_len = read_seq.length() - right;

    hp_read_info_t hp_read_info;
    hp_read_info.hp_len = resolved_hp_len;
    hp_read_info.hp_run_extends_into_5p_tail = hp_run_extends_into_5p_tail;
    hp_read_info.hp_run_extends_into_3p_tail = hp_run_extends_into_3p_tail;
    bool resolved = resolution.hp_len != UNDEFINED_HP_LEN;
    int left_mismatches = resolved ? resolution.left_tail_mismatches : tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos, contig_seq, contig_len, 0, left, true, left_clipped, right_clipped, tail_align_leeway, ref_hp_range.beg);
    int right_mismatches = resolved ? resolution.right_tail_mismatches : tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos, contig_seq, contig_len, right, read_seq.length(), false, left_clipped, right_clipped, tail_align_leeway, ref_hp_range.end);
    if (is_rev) {
        hp_read_info.tail_5p_len = right_len;
        hp_read_info.tail_3p_len = left_len;
        hp_read_info.tail_5p_mismatches = right_mismatches;
        hp_read_info.tail_3p_mismatches = left_mismatches;
    } else {
        hp_read_info.tail_5p_len = left_len;
        hp_read_info.tail_3p_len = right_len;
        hp_read_info.tail_5p_mismatches = left_mismatches;
        hp_read_info.tail_3p_mismatches = right_mismatches;
    }
    hp_read_info.read = read;
    hp_read_info.aligned_5p_tail_len = resolved ? (is_rev ? right_len : left_len) : aligned_5p_tail_len;

    if (resolved) {
        for (int qpos = left; qpos < right; qpos++) {
            if (toupper(read_seq[qpos]) != hp_base) {
                hp_read_info.ref_hp_has_non_hp_read_base = true;
                break;
            }
        }
    } else {
        for (int anchor : anchors) {
            if (toupper(read_seq[anchor]) != hp_base) {
                hp_read_info.ref_hp_has_non_hp_read_base = true;
                break;
            }
        }
    }

    return hp_read_info;
}

hp_read_info_t calculate_hp_read_info(StripedSmithWaterman::Alignment& aln, const std::string& read_seq,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len, bool is_rev, bp_support_read_t read,
    int tail_align_leeway = 10, bool has_no_left_indel = false, bool has_no_right_indel = false,
    double max_5p_mismatch_rate = 1) {

    if (read_seq.empty() || aln.cigar.empty()) {
        return hp_read_info_t();
    }

    std::vector<hp_cigar_op_t> normalized_aln_cigar = normalized_cigar(aln);
    hp_alignment_summary_t alignment_summary = summarize_hp_alignment(normalized_aln_cigar, aln.ref_begin, read_seq, ref_hp_range, hp_base);

    bool hp_deletion_extends_outside_hp = cigar_has_hp_deletion_extending_outside_hp(normalized_aln_cigar, aln.ref_begin, ref_hp_range);
    int force_resolution_max_hp_len = hp_deletion_extends_outside_hp ? ref_hp_range.end - ref_hp_range.beg : UNDEFINED_HP_LEN;
    bool can_use_iterative_hp_len_estimation = has_no_left_indel && has_no_right_indel;
    if (hp_deletion_extends_outside_hp && !can_use_iterative_hp_len_estimation) return hp_read_info_t(UNDEFINED_HP_LEN);
    if (can_use_iterative_hp_len_estimation) {
        hp_adjacent_indel_info_t adjacent_indel_info = get_adjacent_indel_info(normalized_aln_cigar, aln.ref_begin, ref_hp_range, read_seq.length()/2);
        const hp_side_indel_info_t& three_p_info = is_rev ? adjacent_indel_info.left : adjacent_indel_info.right;
        if (three_p_info.indel_len != 0) {
            int ref_hp_len = ref_hp_range.end - ref_hp_range.beg;
            force_resolution_max_hp_len = std::max(force_resolution_max_hp_len, ref_hp_len + std::max(0, three_p_info.indel_len));
        }
    }
    hp_read_info_t hp_read_info = calculate_hp_read_info_core(read_seq,
        alignment_summary.qpos_to_rpos, alignment_summary.anchors,
        alignment_summary.last_qpos_before_ref_hp, alignment_summary.first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, is_rev,
        get_left_clip_size(aln) > 0, get_right_clip_size(aln) > 0,
        read, tail_align_leeway, force_resolution_max_hp_len, can_use_iterative_hp_len_estimation, max_5p_mismatch_rate);
    hp_read_info.original_alignment_has_indel_outside_hp = cigar_has_indel_outside_hp(normalized_aln_cigar, aln.ref_begin, ref_hp_range);
    hp_read_info.hp_deletion_extends_outside_hp = hp_deletion_extends_outside_hp;
    hp_read_info.hp_insertion_has_non_hp_bases = alignment_summary.hp_insertion_has_non_hp_bases;
    return hp_read_info;
}

hp_read_info_t calculate_hp_read_info(bam1_t* read, hts_pair_pos_t ref_hp_range, char hp_base,
    char* contig_seq, hts_pos_t contig_len, char* ref_allele, hts_pos_t ref_allele_len,
    hts_pair_pos_t ref_allele_hp_range, bool has_no_left_indel = false, bool has_no_right_indel = false,
    int tail_align_leeway = 10, double max_5p_mismatch_rate = 1) {

    if (read == NULL || is_unmapped(read) || !is_primary(read) || read->core.l_qseq <= 0) {
        return hp_read_info_t();
    }

    static const StripedSmithWaterman::Aligner permissive_aligner(2, 2, 6, 1, false);
    std::string read_seq = get_sequence(read);
    std::vector<hp_cigar_op_t> normalized_read_cigar = normalized_cigar(read);
    int force_resolution_max_hp_len = UNDEFINED_HP_LEN;
    if (is_clipped(read, 1)) {
        StripedSmithWaterman::Alignment ref_aln;
        StripedSmithWaterman::Filter filter_default;
        if (permissive_aligner.Align(read_seq.c_str(), ref_allele, ref_allele_len,
            filter_default, &ref_aln, 0) && !ref_aln.cigar.empty() && !is_clipped(ref_aln)) {
            hp_read_info_t hp_read_info = calculate_hp_read_info(ref_aln, read_seq, ref_allele_hp_range,
                hp_base, ref_allele, ref_allele_len, bam_is_rev(read), bp_support_read_t(read), tail_align_leeway,
                has_no_left_indel, has_no_right_indel, max_5p_mismatch_rate);
            hp_read_info.original_alignment_has_indel_outside_hp = cigar_has_indel_outside_hp(normalized_read_cigar, read->core.pos, ref_hp_range);
            return hp_read_info;
        }
        force_resolution_max_hp_len = ref_hp_range.end - ref_hp_range.beg;
    }

    hp_alignment_summary_t alignment_summary = summarize_hp_alignment(normalized_read_cigar, read->core.pos, read_seq, ref_hp_range, hp_base);

    bool hp_deletion_extends_outside_hp = cigar_has_hp_deletion_extending_outside_hp(normalized_read_cigar, read->core.pos, ref_hp_range);
    if (hp_deletion_extends_outside_hp) {
        force_resolution_max_hp_len = std::max<int>(force_resolution_max_hp_len, ref_hp_range.end - ref_hp_range.beg);
    }
    bool is_reverse = bam_is_rev(read);
    bool can_use_iterative_hp_len_estimation = has_no_left_indel && has_no_right_indel;
    if (hp_deletion_extends_outside_hp && !can_use_iterative_hp_len_estimation) return hp_read_info_t(UNDEFINED_HP_LEN);
    if (can_use_iterative_hp_len_estimation) {
        hp_adjacent_indel_info_t adjacent_indel_info = get_adjacent_indel_info(normalized_read_cigar, read->core.pos, ref_hp_range, read->core.l_qseq/2);
        const hp_side_indel_info_t& three_p_info = is_reverse ? adjacent_indel_info.left : adjacent_indel_info.right;
        if (three_p_info.indel_len != 0) {
            int ref_hp_len = ref_hp_range.end - ref_hp_range.beg;
            force_resolution_max_hp_len = std::max(force_resolution_max_hp_len, ref_hp_len + std::max(0, three_p_info.indel_len));
        }
    }

    hp_read_info_t hp_read_info = calculate_hp_read_info_core(read_seq,
        alignment_summary.qpos_to_rpos, alignment_summary.anchors,
        alignment_summary.last_qpos_before_ref_hp, alignment_summary.first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, is_reverse,
        get_left_clip_size(read) > 0, get_right_clip_size(read) > 0,
        bp_support_read_t(read), tail_align_leeway, force_resolution_max_hp_len, can_use_iterative_hp_len_estimation, max_5p_mismatch_rate);
    hp_read_info.original_alignment_has_indel_outside_hp = cigar_has_indel_outside_hp(normalized_read_cigar, read->core.pos, ref_hp_range);
    hp_read_info.hp_deletion_extends_outside_hp = hp_deletion_extends_outside_hp;
    hp_read_info.hp_insertion_has_non_hp_bases = alignment_summary.hp_insertion_has_non_hp_bases;
    return hp_read_info;
}

#endif
