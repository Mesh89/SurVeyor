#ifndef HP_READ_INFO_H
#define HP_READ_INFO_H

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#include "htslib/sam.h"

#include "../libs/ssw_cpp.h"
#include "genotype.h"

struct hp_read_info_t {
    int hp_len;
    int tail_5p_len, tail_3p_len;
    int aligned_5p_tail_len = 0;
    int tail_5p_mismatches, tail_3p_mismatches;
    bp_support_read_t read;
    bool original_alignment_has_indel_outside_hp = true;
    bool hp_deletion_extends_outside_hp = false;
    bool hp_insertion_has_non_hp_bases = false;
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
            if (query[j] != std::toupper(ref[i + j])) {
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
        if (std::toupper(read_seq[qpos]) != std::toupper(contig_seq[rpos])) {
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
        if (std::toupper(read_seq[qpos]) != std::toupper(contig_seq[rpos])) {
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

bool insertion_has_non_hp_bases(const std::string& read_seq, int qpos, int len, char hp_base) {
    if (qpos < 0 || qpos + len > (int) read_seq.length()) return true;
    for (int i = 0; i < len; i++) {
        if (std::toupper(read_seq[qpos + i]) != hp_base) return true;
    }
    return false;
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
    bool is_rev, bool left_clipped, bool right_clipped, bp_support_read_t read, int tail_align_leeway) {

    if (read_seq.empty()) {
        return hp_read_info_t();
    }

    int aligned_left_tail_len = 0, aligned_right_tail_len = 0;
    for (hts_pos_t rpos : qpos_to_rpos) {
        if (rpos != -1 && rpos < ref_hp_range.beg) aligned_left_tail_len++;
        if (rpos >= ref_hp_range.end) aligned_right_tail_len++;
    }
    int aligned_5p_tail_len = is_rev ? aligned_right_tail_len : aligned_left_tail_len;

    // If no query base lands inside the reference HP, derive the run from the
    // query gap between the nearest mapped flanks.
    if (anchors.empty()) {
        int hp_len = 0;
        int left_tail_len, right_tail_len;
        if (last_qpos_before_ref_hp != -1 && first_qpos_after_ref_hp == -1) {
            left_tail_len = last_qpos_before_ref_hp + 1;
            right_tail_len = read_seq.length() - left_tail_len;
        } else if (last_qpos_before_ref_hp == -1 && first_qpos_after_ref_hp != -1) {
            left_tail_len = first_qpos_after_ref_hp;
            right_tail_len = read_seq.length() - first_qpos_after_ref_hp;
        } else if (last_qpos_before_ref_hp != -1 && first_qpos_after_ref_hp != -1) {
            hp_len = std::max<hts_pos_t>(0, first_qpos_after_ref_hp - last_qpos_before_ref_hp - 1);
            left_tail_len = last_qpos_before_ref_hp + 1;
            right_tail_len = read_seq.length() - first_qpos_after_ref_hp;
        } else {
            return hp_read_info_t();
        }

        int left_mismatches = tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos,
            contig_seq, contig_len, 0, left_tail_len, true, left_clipped, right_clipped, tail_align_leeway, ref_hp_range.beg);
        int right_mismatches = tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos,
            contig_seq, contig_len, read_seq.length() - right_tail_len, read_seq.length(), false,
            left_clipped, right_clipped, tail_align_leeway, ref_hp_range.end);

        hp_read_info_t hp_read_info = !is_rev ?
            hp_read_info_t(hp_len, left_tail_len, right_tail_len, left_mismatches, right_mismatches, read) :
            hp_read_info_t(hp_len, right_tail_len, left_tail_len, right_mismatches, left_mismatches, read);
        hp_read_info.aligned_5p_tail_len = aligned_5p_tail_len;
        return hp_read_info;
    }

    // Start from the aligned bases inside the reference HP and expand through
    // adjacent query bases that still look like part of the same HP run.
    int left = anchors.front();
    while (left > 0 && read_seq[left-1] == hp_base) {
        hts_pos_t mapped_rpos = qpos_to_rpos[left-1];
        if (mapped_rpos != -1 && mapped_rpos < ref_hp_range.beg) break;
        left--;
    }

    int right = anchors.back() + 1;
    while (right < (int) read_seq.length() && read_seq[right] == hp_base) {
        hts_pos_t mapped_rpos = qpos_to_rpos[right];
        if (mapped_rpos != -1 && mapped_rpos >= ref_hp_range.end) break;
        right++;
    }

    int left_len = left;
    int right_len = read_seq.length() - right;

    hp_read_info_t hp_read_info;
    hp_read_info.hp_len = right - left;
    int left_mismatches = tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos,
        contig_seq, contig_len, 0, left, true, left_clipped, right_clipped, tail_align_leeway, ref_hp_range.beg);
    int right_mismatches = tail_mismatch_count_from_mapping(read_seq, qpos_to_rpos,
        contig_seq, contig_len, right, read_seq.length(), false, left_clipped, right_clipped, tail_align_leeway, ref_hp_range.end);
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
    hp_read_info.aligned_5p_tail_len = aligned_5p_tail_len;

    for (int anchor : anchors) {
        if (std::toupper(read_seq[anchor]) != hp_base) {
            hp_read_info.ref_hp_has_non_hp_read_base = true;
            break;
        }
    }

    // The inferred run stops at an HP query base aligned just outside the
    // reference HP. If that reference base is non-HP, the query run really
    // extends into the 3' tail and must not count as exact REF support.
    int adjacent_3p_qpos = is_rev ? left - 1 : right;
    if (adjacent_3p_qpos >= 0 && adjacent_3p_qpos < (int) read_seq.length() &&
        read_seq[adjacent_3p_qpos] == hp_base) {
        hts_pos_t mapped_rpos = qpos_to_rpos[adjacent_3p_qpos];
        bool maps_outside_ref_hp = mapped_rpos < ref_hp_range.beg || mapped_rpos >= ref_hp_range.end;
        bool maps_to_valid_ref_base = mapped_rpos >= 0 && mapped_rpos < contig_len;
        hp_read_info.hp_run_extends_into_3p_tail = maps_outside_ref_hp && maps_to_valid_ref_base && toupper(contig_seq[mapped_rpos]) != hp_base;
    }

    return hp_read_info;
}

hp_read_info_t calculate_hp_read_info(bam1_t* read, hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len,
    int tail_align_leeway = 10) {
    if (read == NULL || is_unmapped(read) || !is_primary(read) || read->core.l_qseq <= 0) {
        return hp_read_info_t();
    }

    std::string read_seq = get_sequence(read);
    std::vector<hts_pos_t> qpos_to_rpos(read->core.l_qseq, -1);
    std::vector<int> anchors;

    uint32_t* cigar = bam_get_cigar(read);
    int qpos = 0;
    hts_pos_t rpos = read->core.pos;
    hts_pos_t last_qpos_before_ref_hp = -1, first_qpos_after_ref_hp = -1;
    bool hp_insertion_has_non_hp_bases = false;
    for (uint32_t i = 0; i < read->core.n_cigar; i++) {
        char opchar = bam_cigar_opchr(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);

        if (opchar == 'M' || opchar == '=' || opchar == 'X') {
            for (int j = 0; j < oplen; j++) {
                qpos_to_rpos[qpos] = rpos;
                if (rpos < ref_hp_range.beg) {
                    last_qpos_before_ref_hp = qpos;
                } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                    first_qpos_after_ref_hp = qpos;
                } else if (ref_hp_range.beg <= rpos && rpos < ref_hp_range.end) {
                    anchors.push_back(qpos);
                }
                qpos++;
                rpos++;
            }
        } else if (opchar == 'I') {
            if (ref_hp_range.beg <= rpos && rpos <= ref_hp_range.end && insertion_has_non_hp_bases(read_seq, qpos, oplen, hp_base)) {
                hp_insertion_has_non_hp_bases = true;
            }
            if (rpos < ref_hp_range.beg) {
                last_qpos_before_ref_hp = qpos + oplen - 1;
            } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                first_qpos_after_ref_hp = (ref_hp_range.beg == ref_hp_range.end && rpos == ref_hp_range.beg) ? qpos + oplen : qpos;
            }
            qpos += oplen;
        } else if (opchar == 'S') {
            qpos += oplen;
        } else if (opchar == 'D' || opchar == 'N') {
            rpos += oplen;
        }
    }
    hp_read_info_t hp_read_info = calculate_hp_read_info_core(read_seq, qpos_to_rpos, anchors, last_qpos_before_ref_hp, first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, bam_is_rev(read), get_left_clip_size(read) > 0, get_right_clip_size(read) > 0,
        bp_support_read_t(read), tail_align_leeway);
    std::vector<hp_cigar_op_t> normalized_read_cigar = normalized_cigar(read);
    hp_read_info.original_alignment_has_indel_outside_hp = cigar_has_indel_outside_hp(normalized_read_cigar, read->core.pos, ref_hp_range);
    hp_read_info.hp_deletion_extends_outside_hp = cigar_has_hp_deletion_extending_outside_hp(normalized_read_cigar, read->core.pos, ref_hp_range);
    hp_read_info.hp_insertion_has_non_hp_bases = hp_insertion_has_non_hp_bases;
    return hp_read_info;
}

hp_read_info_t calculate_hp_read_info(StripedSmithWaterman::Alignment& aln, const std::string& read_seq,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len, bool is_rev, bp_support_read_t read,
    int tail_align_leeway = 10) {
    
        if (read_seq.empty() || aln.cigar.empty()) {
        return hp_read_info_t();
    }

    std::vector<hts_pos_t> qpos_to_rpos(read_seq.length(), -1);
    std::vector<int> anchors;

    int qpos = 0;
    hts_pos_t rpos = aln.ref_begin;
    hts_pos_t last_qpos_before_ref_hp = -1, first_qpos_after_ref_hp = -1;
    bool hp_insertion_has_non_hp_bases = false;
    for (uint32_t cigar_op : aln.cigar) {
        char opchar = cigar_int_to_op(cigar_op);
        int oplen = cigar_int_to_len(cigar_op);

        if (opchar == 'M' || opchar == '=' || opchar == 'X') {
            for (int j = 0; j < oplen; j++) {
                qpos_to_rpos[qpos] = rpos;
                if (rpos < ref_hp_range.beg) {
                    last_qpos_before_ref_hp = qpos;
                } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                    first_qpos_after_ref_hp = qpos;
                } else if (ref_hp_range.beg <= rpos && rpos < ref_hp_range.end) {
                    anchors.push_back(qpos);
                }
                qpos++;
                rpos++;
            }
        } else if (opchar == 'I') {
            if (ref_hp_range.beg <= rpos && rpos <= ref_hp_range.end && insertion_has_non_hp_bases(read_seq, qpos, oplen, hp_base)) {
                hp_insertion_has_non_hp_bases = true;
            }
            if (rpos < ref_hp_range.beg) {
                last_qpos_before_ref_hp = qpos + oplen - 1;
            } else if (rpos >= ref_hp_range.end && first_qpos_after_ref_hp == -1) {
                first_qpos_after_ref_hp = qpos;
            }
            qpos += oplen;
        } else if (opchar == 'S') {
            qpos += oplen;
        } else if (opchar == 'D' || opchar == 'N') {
            rpos += oplen;
        }
    }

    hp_read_info_t hp_read_info = calculate_hp_read_info_core(read_seq, qpos_to_rpos, anchors, last_qpos_before_ref_hp, first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, is_rev,
        get_left_clip_size(aln) > 0, get_right_clip_size(aln) > 0, read, tail_align_leeway);
    std::vector<hp_cigar_op_t> normalized_aln_cigar = normalized_cigar(aln);
    hp_read_info.original_alignment_has_indel_outside_hp = cigar_has_indel_outside_hp(normalized_aln_cigar, aln.ref_begin, ref_hp_range);
    hp_read_info.hp_deletion_extends_outside_hp = cigar_has_hp_deletion_extending_outside_hp(normalized_aln_cigar, aln.ref_begin, ref_hp_range);
    hp_read_info.hp_insertion_has_non_hp_bases = hp_insertion_has_non_hp_bases;
    return hp_read_info;
}

#endif
