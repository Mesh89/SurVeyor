#ifndef GENOTYPE_HP_INDELS_H
#define GENOTYPE_HP_INDELS_H

#include <memory>

#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "utils.h"
#include "stat_tests.h"
#include "../libs/ssw_cpp.h"

#include "genotype.h"
#include "var_utils.h"

const double MAX_TAIL_MISMATCH_RATE = 0.2; 

struct hp_read_info_t {
    int hp_len;
    int tail_5p_len, tail_3p_len;
    int tail_5p_mismatches, tail_3p_mismatches;
    bp_support_read_t read;
    bool rescued = false;

    hp_read_info_t(int hp_len = 0, int tail_5p_len = 0, int tail_3p_len = 0,
        int tail_5p_mismatches = 0, int tail_3p_mismatches = 0,
        bp_support_read_t read = bp_support_read_t(), bool rescued = false) :
        hp_len(hp_len), tail_5p_len(tail_5p_len), tail_3p_len(tail_3p_len),
        tail_5p_mismatches(tail_5p_mismatches), tail_3p_mismatches(tail_3p_mismatches),
        read(read), rescued(rescued) {}

    bool is_good_read(int min_tail_len, double max_mismatch_rate) const {
        return tail_5p_len >= min_tail_len && double(tail_5p_mismatches)/tail_5p_len <= max_mismatch_rate &&
               tail_3p_len >= min_tail_len && double(tail_3p_mismatches)/tail_3p_len <= max_mismatch_rate;
    }
};


// Find the mode of HP lengths, optionally restricted to good reads
// Good reads = <20% mismatches on both tails
int find_hp_len_mode(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate, bool good_reads_only) {
    std::unordered_map<int, int> hp_len_counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (!good_reads_only || hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            hp_len_counts[hp_read_info.hp_len]++;
        }
    }

    int mode_hp_len = -1, max_count = 0;
    for (const auto& kv : hp_len_counts) {
        if (kv.second > max_count || (kv.second == max_count && kv.first < mode_hp_len)) {
            mode_hp_len = kv.first;
            max_count = kv.second;
        }
    }
    if (max_count > 0) {
        return mode_hp_len;
    }
    return -1; // no good reads
}

double quantile_linear_interp(std::vector<int>& sorted_values, double q) {
    if (sorted_values.empty()) {
        return -1.0;
    }
    if (sorted_values.size() == 1) {
        return sorted_values[0];
    }

    double pos = q * (sorted_values.size() - 1);
    int lo = std::floor(pos);
    int hi = std::ceil(pos);
    if (lo == hi) {
        return sorted_values[lo];
    }

    double frac = pos - lo;
    return sorted_values[lo] + frac * (sorted_values[hi] - sorted_values[lo]);
}

double find_hp_len_iqr(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate, bool good_reads_only) {
    std::vector<int> hp_lens;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (!good_reads_only || hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            hp_lens.push_back(hp_read_info.hp_len);
        }
    }

    if (hp_lens.empty()) {
        return -1.0;
    }

    std::sort(hp_lens.begin(), hp_lens.end());
    return quantile_linear_interp(hp_lens, 0.75) - quantile_linear_interp(hp_lens, 0.25);
}

void set_hp_tail_mismatch_rates(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate,
    bool good_reads_only, double& avg_5p_mismatch_rate, double& avg_3p_mismatch_rate) {
    
    double sum_5p_mismatch_rate = 0.0;
    double sum_3p_mismatch_rate = 0.0;
    int n_reads = 0;

    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (good_reads_only && !hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) continue;

        sum_5p_mismatch_rate += double(hp_read_info.tail_5p_mismatches) / hp_read_info.tail_5p_len;
        sum_3p_mismatch_rate += double(hp_read_info.tail_3p_mismatches) / hp_read_info.tail_3p_len;
        n_reads++;
    }

    if (n_reads == 0) {
        return;
    }

    avg_5p_mismatch_rate = sum_5p_mismatch_rate / n_reads;
    avg_3p_mismatch_rate = sum_3p_mismatch_rate / n_reads;
}

void set_hp_read_mapq_stats(const std::vector<bp_support_read_t>& reads, int& min_mapq, int& max_mapq, double& avg_mapq, double& stddev_mapq) {
    if (reads.empty()) {
        return;
    }

    std::vector<int> mapqs;
    mapqs.reserve(reads.size());
    for (const bp_support_read_t& read : reads) {
        mapqs.push_back(read.mapq);
    }

    min_mapq = *std::min_element(mapqs.begin(), mapqs.end());
    max_mapq = *std::max_element(mapqs.begin(), mapqs.end());
    avg_mapq = mean(mapqs);
    stddev_mapq = stddev(mapqs);
}

void set_hp_read_mapq_stats(const std::vector<std::shared_ptr<bam1_t>>& reads, int& min_mapq, int& max_mapq, double& avg_mapq, double& stddev_mapq) {
    std::vector<bp_support_read_t> support_reads;
    support_reads.reserve(reads.size());
    for (const std::shared_ptr<bam1_t>& read : reads) {
        support_reads.emplace_back(read.get());
    }
    set_hp_read_mapq_stats(support_reads, min_mapq, max_mapq, avg_mapq, stddev_mapq);
}

std::vector<std::vector<hp_read_info_t>> cluster_reads_by_two_modes(const std::vector<hp_read_info_t>& hp_read_infos,
    int min_tail_len, double max_mismatch_rate) {
    int n = hp_read_infos.size();
    if (n == 0) {
        return {};
    }

    std::unordered_map<int, int> good_counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            good_counts[hp_read_info.hp_len]++;
        }
    }

    if (good_counts.empty()) {
        return {};
    }

    // Pick seed modes from good reads only, breaking count ties toward the smaller length.
    std::vector<std::pair<int, int>> count_items(good_counts.begin(), good_counts.end());
    std::sort(count_items.begin(), count_items.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    if (count_items.size() == 1) {
        return {hp_read_infos};
    }

    int a_mode = count_items[0].first, b_mode = count_items[1].first;
    if (a_mode > b_mode) std::swap(a_mode, b_mode);

    // Absorb every read into the nearest good-read mode; distance ties go to the larger good-read seed.
    std::vector<int> a_idx, b_idx;
    for (int i = 0; i < n; i++) {
        int x = hp_read_infos[i].hp_len;
        int a_dist = abs(x - a_mode);
        int b_dist = abs(x - b_mode);
        if (a_dist < b_dist) {
            a_idx.push_back(i);
        } else if (b_dist < a_dist) {
            b_idx.push_back(i);
        } else if (good_counts[a_mode] >= good_counts[b_mode]) {
            a_idx.push_back(i);
        } else {
            b_idx.push_back(i);
        }
    }

    std::vector<hp_read_info_t> a_reads, b_reads;
    for (int i : a_idx) a_reads.push_back(hp_read_infos[i]);
    for (int i : b_idx) b_reads.push_back(hp_read_infos[i]);
    return {a_reads, b_reads};
}

struct hp_allele_cluster_t {
    int allele_idx;
    int allele_len;
    std::vector<hp_read_info_t> reads;

    hp_allele_cluster_t(int allele_idx, int allele_len) : allele_idx(allele_idx), allele_len(allele_len) {}
};

int nearest_hp_allele_idx(int hp_len, const std::vector<int>& allele_lens, const std::vector<int>* allowed_allele_idxs = NULL) {
    int best_idx = -1;
    int best_dist = INT32_MAX;

    std::vector<int> all_allele_idxs;
    if (allowed_allele_idxs == NULL) {
        for (int i = 0; i < allele_lens.size(); i++) all_allele_idxs.push_back(i);
        allowed_allele_idxs = &all_allele_idxs;
    }

    for (int idx : *allowed_allele_idxs) {
        int dist = abs(hp_len - allele_lens[idx]);
        if (dist < best_dist ||
            (dist == best_dist && (best_idx == -1 || allele_lens[idx] < allele_lens[best_idx])) ||
            (dist == best_dist && allele_lens[idx] == allele_lens[best_idx] && idx < best_idx)) {
            best_idx = idx;
            best_dist = dist;
        }
    }

    return best_idx;
}

std::vector<hp_allele_cluster_t> cluster_reads_by_nearest_allele_len(const std::vector<hp_read_info_t>& hp_read_infos,
    const std::vector<int>& allele_lens) {

    if (hp_read_infos.empty() || allele_lens.empty()) return {};

    // First, assign each read to its nearest allele length
    std::vector<int> allele_counts(allele_lens.size(), 0);
    std::vector<long long> allele_dist_sums(allele_lens.size(), 0);
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        int allele_idx = nearest_hp_allele_idx(hp_read_info.hp_len, allele_lens);
        allele_counts[allele_idx]++;
        allele_dist_sums[allele_idx] += abs(hp_read_info.hp_len - allele_lens[allele_idx]);
    }

    // Then pick the (max 2) alleles with the most assigned reads, breaking count ties toward smaller total distance, then smaller length
    std::vector<int> allele_idxs;
    for (int i = 0; i < allele_lens.size(); i++) {
        if (allele_counts[i] > 0) allele_idxs.push_back(i);
    }
    std::sort(allele_idxs.begin(), allele_idxs.end(), [&](int a, int b) {
        if (allele_counts[a] != allele_counts[b]) return allele_counts[a] > allele_counts[b];
        if (allele_dist_sums[a] != allele_dist_sums[b]) return allele_dist_sums[a] < allele_dist_sums[b];
        if (allele_lens[a] != allele_lens[b]) return allele_lens[a] < allele_lens[b];
        return a < b;
    });
    if (allele_idxs.size() > 2) allele_idxs.resize(2);
    std::sort(allele_idxs.begin(), allele_idxs.end(), [&](int a, int b) {
        if (allele_lens[a] != allele_lens[b]) return allele_lens[a] < allele_lens[b];
        return a < b;
    });

    // Finally, form max 2 clusters of reads by assigning reads to the max 2 retained alleles
    std::vector<hp_allele_cluster_t> clusters;
    for (int allele_idx : allele_idxs) {
        clusters.push_back(hp_allele_cluster_t(allele_idx, allele_lens[allele_idx]));
    }
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        int allele_idx = nearest_hp_allele_idx(hp_read_info.hp_len, allele_lens, &allele_idxs);
        for (hp_allele_cluster_t& cluster : clusters) {
            if (cluster.allele_idx == allele_idx) {
                cluster.reads.push_back(hp_read_info);
                break;
            }
        }
    }

    return clusters;
}


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

// Shared interpretation core. Callers are responsible for providing the
// query-to-reference mapping summary in the coordinate space of the sequence
// used for HP evaluation.
hp_read_info_t calculate_hp_read_info_core(const std::string& read_seq, const std::vector<hts_pos_t>& qpos_to_rpos,
    const std::vector<int>& anchors, hts_pos_t last_qpos_before_ref_hp, hts_pos_t first_qpos_after_ref_hp,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len,
    bool is_rev, bool left_clipped, bool right_clipped, bp_support_read_t read) {

    const int tail_align_leeway = 10;

    if (read_seq.empty()) {
        return hp_read_info_t();
    }

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

        if (!is_rev) {
            return hp_read_info_t(hp_len, left_tail_len, right_tail_len, left_mismatches, right_mismatches, read);
        } else {
            return hp_read_info_t(hp_len, right_tail_len, left_tail_len, right_mismatches, left_mismatches, read);
        }
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

    return hp_read_info;
}

hp_read_info_t calculate_hp_read_info(bam1_t* read, hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len) {
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
    return calculate_hp_read_info_core(read_seq, qpos_to_rpos, anchors, last_qpos_before_ref_hp, first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, bam_is_rev(read),
        get_left_clip_size(read) > 0, get_right_clip_size(read) > 0,
        bp_support_read_t(read));
}

hp_read_info_t calculate_hp_read_info(StripedSmithWaterman::Alignment& aln, const std::string& read_seq,
    hts_pair_pos_t ref_hp_range, char hp_base, char* contig_seq, hts_pos_t contig_len, bool is_rev, bp_support_read_t read) {
    if (read_seq.empty() || aln.cigar.empty()) {
        return hp_read_info_t();
    }

    std::vector<hts_pos_t> qpos_to_rpos(read_seq.length(), -1);
    std::vector<int> anchors;

    int qpos = 0;
    hts_pos_t rpos = aln.ref_begin;
    hts_pos_t last_qpos_before_ref_hp = -1, first_qpos_after_ref_hp = -1;
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

    return calculate_hp_read_info_core(read_seq, qpos_to_rpos, anchors, last_qpos_before_ref_hp, first_qpos_after_ref_hp,
        ref_hp_range, hp_base, contig_seq, contig_len, is_rev,
        get_left_clip_size(aln) > 0, get_right_clip_size(aln) > 0, read);
}

std::vector<int> calculate_aln_scores(std::vector<hp_read_info_t>& hp_read_infos, char* ref_seq, int ref_len, 
    StripedSmithWaterman::Aligner& aligner) {
    std::vector<int> scores;
    scores.reserve(hp_read_infos.size());

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter_score_only(false, false, 0, 32767);
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (hp_read_info.read.seq.empty()) {
            scores.push_back(0);
            continue;
        }

        aligner.Align(hp_read_info.read.seq.c_str(), ref_seq, ref_len, filter_score_only, &aln, 0);
        scores.push_back(aln.sw_score);
    }

    return scores;
}

void genotype_hp_indels_group(std::vector<sv_t*>& hp_indels, hts_pair_pos_t ref_hp_range, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len,
    stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner,
    std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, evidence_logger_t* evidence_logger,
    bool reassign_evidence, evidence_map_t* evidence_map, std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map) {

    if (hp_indels.empty()) return;
    for (sv_t* hp_indel : hp_indels) {
        hp_indel->hp_genotyped = true;
        hp_indel->hp_ref_beg = ref_hp_range.beg;
        hp_indel->hp_ref_end = ref_hp_range.end;
    }

    char hp_base = get_homopolymer_base(hp_indels[0], contig_seq);

    std::vector<hp_read_info_t> hp_read_infos;

    std::stringstream region_ss;
    region_ss << hp_indels[0]->chr << ":" << ref_hp_range.beg-1 << "-" << ref_hp_range.end+1;
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region_ss.str().c_str());

    // Build alt alleles
    hts_pos_t extend = stats.read_len - 1;
    hts_pos_t ref_hp_len = ref_hp_range.end - ref_hp_range.beg;
    hts_pos_t alt_start = std::max(hts_pos_t(0), ref_hp_range.beg - extend);
    hts_pos_t alt_end = std::min(contig_len, ref_hp_range.end + extend);
    hts_pos_t left_flank_len = ref_hp_range.beg - alt_start;
    hts_pos_t right_flank_len = alt_end - ref_hp_range.end;
    hts_pos_t ref_len = left_flank_len + ref_hp_len + right_flank_len;

    std::unique_ptr<char[]> ref_allele(new char[ref_len + 1]);
    strncpy(ref_allele.get(), contig_seq + alt_start, left_flank_len);
    memset(ref_allele.get() + left_flank_len, hp_base, ref_hp_len);
    strncpy(ref_allele.get() + left_flank_len + ref_hp_len, contig_seq + ref_hp_range.end, right_flank_len);
    ref_allele[ref_len] = '\0';
    hts_pair_pos_t ref_allele_hp_range = {left_flank_len, left_flank_len + ref_hp_len};

    std::vector<std::unique_ptr<char[]>> alt_alleles;
    std::vector<int> alt_allele_lens;
    for (sv_t* hp_indel : hp_indels) {
        char* lh_seq = generate_haplotype_left(contig_seq, ref_hp_range.beg - 1, extend,
            hp_indel->aux_indels, hp_indel->aux_snps);
        char* rh_seq = generate_haplotype_right(contig_seq, contig_len, ref_hp_range.end, extend,
            hp_indel->aux_indels, hp_indel->aux_snps);
        hts_pos_t alt_lh_len = strlen(lh_seq);
        hts_pos_t alt_rh_len = strlen(rh_seq);
        hts_pos_t alt_hp_len = ref_hp_len + hp_indel->svlen();
        hts_pos_t alt_len = alt_lh_len + alt_hp_len + alt_rh_len;
        std::unique_ptr<char[]> alt_seq(new char[alt_len + 1]);
        strncpy(alt_seq.get(), lh_seq, alt_lh_len);
        memset(alt_seq.get() + alt_lh_len, hp_base, alt_hp_len);
        strncpy(alt_seq.get() + alt_lh_len + alt_hp_len, rh_seq, alt_rh_len);
        alt_seq[alt_len] = '\0';
        alt_alleles.push_back(std::move(alt_seq));
        alt_allele_lens.push_back(alt_len);
        delete[] lh_seq;
        delete[] rh_seq;
    }

    // For each read overlapping the reference HP, calculate its observed HP length,
    // 5'/3' tail lengths and 5'/3' tail mismatch counts.
    StripedSmithWaterman::Alignment alt_aln;
    StripedSmithWaterman::Filter filter_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Filter filter_default;

    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (!is_proper_pair(read, stats.min_is, stats.max_is)) continue;

        // If this read is assigned to a different, non-HP SV, we don't use it for evidence here
        bool discard_read = reassign_evidence;
        for (int i = 0; i < hp_indels.size(); i++) {
             if (!reassign_evidence || !evidence_map->is_read_assigned_to_different_sv(read, hp_indels[i])) {
                discard_read = false;
                break;
            }
        }
        if (discard_read) {
            for (sv_t* hp_indel : hp_indels) {
                hp_indel->sample_info.assigned_to_other_sv_bp1_reads++;
            }
            continue;
        }

        hp_read_info_t hp_read_info = calculate_hp_read_info(read, ref_hp_range, hp_base, contig_seq, contig_len);
        if (hp_read_info.tail_3p_len < config.min_clip_len || hp_read_info.tail_5p_len < config.min_clip_len) {
            // Even if we are discarding this reads because the tails are too short, we still want to prevent it from being used as evidence for non-HP indels
            // when they are aligned better to one of the HP indels
            std::string seq = get_sequence(read);
            for (int i = 0; i < hp_indels.size(); i++) {
                // if the read is assigned to a different SV, no need to align it, just count and continue

                aligner.Align(seq.c_str(), alt_alleles[i].get(), alt_allele_lens[i], filter_score_only, &alt_aln, 0);
                std::vector<std::shared_ptr<bam1_t>> read_v = { std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1) };
                std::vector<int> alt_aln_scores = { alt_aln.sw_score };
                if (evidence_logger) evidence_logger->log_reads_associations(hp_indels[i]->id, 1, read_v, alt_aln_scores);
            }
            continue;
        }
        hp_read_infos.push_back(hp_read_info);
    }
    hts_itr_destroy(iter);

    // Realign reads that "escaped" to a different locus, but their mate betrays them
    StripedSmithWaterman::Alignment ref_aln;
    std::stringstream possible_mates_ss;
    possible_mates_ss << hp_indels[0]->chr << ":" << std::max(hts_pos_t(0), ref_hp_range.beg - stats.max_is)
        << "-" << std::min(contig_len, ref_hp_range.end + stats.max_is);
    iter = sam_itr_querys(bam_file->idx, bam_file->header, possible_mates_ss.str().c_str());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (!is_dc_pair(read)) continue;

        std::string mate_seq;
        int mate_mapq;
        hts_pos_t endpos = bam_endpos(read);
        if (!bam_is_rev(read) && read->core.pos >= ref_hp_range.beg-stats.max_is && read->core.pos <= ref_hp_range.beg-stats.read_len/2) {
            std::string qname = get_mate_lookup_qname(read);
            if (!mateseqs_w_mapq_chr.count(qname)) continue;
            mate_seq = mateseqs_w_mapq_chr[qname].first;
            mate_mapq = mateseqs_w_mapq_chr[qname].second;
            rc(mate_seq);
        } else if (bam_is_rev(read) && endpos >= ref_hp_range.end+stats.read_len/2 && endpos <= ref_hp_range.end+stats.max_is) {
            std::string qname = get_mate_lookup_qname(read);
            if (!mateseqs_w_mapq_chr.count(qname)) continue;
            mate_seq = mateseqs_w_mapq_chr[qname].first;
            mate_mapq = mateseqs_w_mapq_chr[qname].second;
        } else {
            continue;
        }

        aligner.Align(mate_seq.c_str(), ref_allele.get(), ref_len, filter_default, &ref_aln, 0);

        bool aln_as_rev = !bam_is_rev(read);
        if ((!aln_as_rev && get_left_clip_size(ref_aln) > 0) ||
             (aln_as_rev && get_right_clip_size(ref_aln) > 0)) {
            continue;
        }

        // The rescued read information mostly comes from its mate
        bp_support_read_t rescued_read;
        rescued_read.read_name = bam_get_qname(read);
        rescued_read.mapq = mate_mapq;
        rescued_read.mate_mapq = read->core.qual;
        rescued_read.seq = mate_seq;
        rescued_read.mate_is_reverse = bam_is_rev(read);
        rescued_read.mate_pos = read->core.pos;
        rescued_read.mate_endpos = endpos;
        rescued_read.is_first_in_pair = !is_first_read(read);
        hp_read_info_t hp_read_info = calculate_hp_read_info(ref_aln, mate_seq, ref_allele_hp_range, hp_base, ref_allele.get(), ref_len, aln_as_rev, rescued_read);

        if (hp_read_info.tail_3p_len < config.min_clip_len || hp_read_info.tail_5p_len < config.min_clip_len) {
            continue;
        }

        // We need this to feed into set_bp_consensus_info
        // However, we don't want to feed it to set_hp_read_mapq_stats
        hp_read_info.rescued = true;
        hp_read_infos.push_back(hp_read_info);
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);

    std::vector<int> allele_lens;
    for (sv_t* hp_indel : hp_indels) {
        allele_lens.push_back(ref_hp_range.end - ref_hp_range.beg + hp_indel->svlen());
    }
    allele_lens.push_back(ref_hp_range.end - ref_hp_range.beg); // allele_lens[hp_indels.size()] corresponds to the reference allele

    // Cluster reads into up to two clusters around the candidate allele HP lengths
    std::vector<hp_allele_cluster_t> clusters = cluster_reads_by_nearest_allele_len(hp_read_infos, allele_lens);

    int ref_reads = 0;
    std::vector<bp_support_read_t> ref_good_reads;
    std::vector<bool> ref_is_exact_match;
 
    std::vector<int> alt_reads(hp_indels.size(), 0);
    std::vector<std::vector<hp_read_info_t>> alt_assigned_hp_read_infos(hp_indels.size());
    std::vector<std::vector<bp_support_read_t>> alt_good_reads(hp_indels.size()),alt_good_reads_non_rescued(hp_indels.size());
    std::vector<std::vector<bool>> alt_is_exact_match(hp_indels.size());
    std::unordered_set<std::string> hp_indel_ids;
    for (sv_t* hp_indel : hp_indels) {
        hp_indel_ids.insert(remove_svid_dup_suffix(hp_indel->id));
    }

    // Associate each allele-aware cluster to its selected allele
    for (const hp_allele_cluster_t& allele_cluster : clusters) {
        const std::vector<hp_read_info_t>& cluster = allele_cluster.reads;
        int best_allele_idx = allele_cluster.allele_idx;

        std::vector<hp_read_info_t> good_hp_read_infos;
        std::vector<bp_support_read_t> good_reads, good_reads_non_rescued;
        std::vector<bool> is_exact_match;
        for (const hp_read_info_t& hp_read_info : cluster) {
            if (hp_read_info.is_good_read(config.min_clip_len, MAX_TAIL_MISMATCH_RATE)) {
                good_hp_read_infos.push_back(hp_read_info);
                good_reads.push_back(hp_read_info.read);
                if (!hp_read_info.rescued) {
                    good_reads_non_rescued.push_back(hp_read_info.read);
                }
                is_exact_match.push_back(hp_read_info.hp_len == allele_lens[best_allele_idx]);
            }
        }

        if (best_allele_idx == allele_lens.size() - 1) {
            // Cluster best matches the reference allele
            ref_reads += cluster.size();
            ref_good_reads.insert(ref_good_reads.end(), good_reads.begin(), good_reads.end());
            ref_is_exact_match.insert(ref_is_exact_match.end(), is_exact_match.begin(), is_exact_match.end());
        } else {
            // Cluster best matches an indel allele, so assign its reads to that allele
            alt_reads[best_allele_idx] += cluster.size();
            alt_assigned_hp_read_infos[best_allele_idx].insert(alt_assigned_hp_read_infos[best_allele_idx].end(), cluster.begin(), cluster.end());
            alt_good_reads[best_allele_idx].insert(alt_good_reads[best_allele_idx].end(), good_reads.begin(), good_reads.end());
            alt_good_reads_non_rescued[best_allele_idx].insert(alt_good_reads_non_rescued[best_allele_idx].end(), good_reads_non_rescued.begin(), good_reads_non_rescued.end());
            alt_is_exact_match[best_allele_idx].insert(alt_is_exact_match[best_allele_idx].end(), is_exact_match.begin(), is_exact_match.end());

            std::vector<int> alt_scores = calculate_aln_scores(good_hp_read_infos, alt_alleles[best_allele_idx].get(), alt_allele_lens[best_allele_idx], aligner);
            if (evidence_logger) evidence_logger->log_reads_associations(hp_indels[best_allele_idx]->id, 1, good_reads, alt_scores);
        }
    }

    for (int i = 0; i < hp_indels.size(); i++) {
        // set OR* reads for other indels
        int or1hq = 0, or1e = 0;
        for (int j = 0; j < alt_good_reads[i].size(); j++) {
            if (alt_good_reads[i][j].mate_mapq >= config.high_confidence_mapq) {
                or1hq++;
            }
            if (alt_is_exact_match[i][j]) {
                or1e++;
            }
        }
        for (int j = 0; j < hp_indels.size(); j++) {
            if (i == j) continue;
            hp_indels[j]->sample_info.assigned_to_other_sv_bp1_reads += alt_reads[i];
            hp_indels[j]->sample_info.assigned_to_other_sv_bp1_consistent += alt_good_reads[i].size();
            hp_indels[j]->sample_info.assigned_to_other_sv_bp1_consistent_highmq += or1hq;
            hp_indels[j]->sample_info.assigned_to_other_sv_bp1_consistent_exact += or1e;
        }
        if (reassign_evidence) {
            for (const bp_support_read_t& read : alt_good_reads[i]) {
                for (std::pair<std::string, int>& ov : evidence_map->get_non_chosen_svs_for_read(read)) {
                    if (hp_indel_ids.count(remove_svid_dup_suffix(ov.first))) continue;
                    increase_orc(sv_map, ov.first, ov.second, read.mate_mapq >= config.high_confidence_mapq);
                }
            }
        }

        set_bp_consensus_info(hp_indels[i]->sample_info.alt_bp1.reads_info, alt_reads[i], alt_good_reads[i], alt_is_exact_match[i], 0.0, 0.0);
        hp_indels[i]->sample_info.alt1_hp_len_mode = find_hp_len_mode(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false);
        hp_indels[i]->sample_info.alt1_consistent_hp_len_mode = find_hp_len_mode(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        hp_indels[i]->sample_info.alt1_consistent_hp_len_iqr = find_hp_len_iqr(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        set_hp_tail_mismatch_rates(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false, 
            hp_indels[i]->sample_info.alt1_hp_5p_mismatch_rate, hp_indels[i]->sample_info.alt1_hp_3p_mismatch_rate);

        // We shouldn't use the mapping qualities of rescued reads since they reflect the confidence
        // of the original mapping rather than the remapping to the HP alleles
        set_hp_read_mapq_stats(alt_good_reads_non_rescued[i], hp_indels[i]->sample_info.alt1_hp_min_mapq, hp_indels[i]->sample_info.alt1_hp_max_mapq, 
            hp_indels[i]->sample_info.alt1_hp_avg_mapq, hp_indels[i]->sample_info.alt1_hp_stddev_mapq);
        set_bp_consensus_info(hp_indels[i]->sample_info.ref_bp1.reads_info, ref_reads, ref_good_reads, ref_is_exact_match, 0.0, 0.0);
    }
}

// hp_indels are guaranteed to be on the same chromosome
void genotype_hp_indels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<sv_t*> hp_indels,
    stats_t& stats, config_t& config, contig_map_t& contig_map, bam_pool_t* bam_pool,
    std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr, 
    std::vector<double>* global_crossing_isize_dist,
    evidence_logger_t* evidence_logger, bool reassign_evidence, evidence_map_t* evidence_map,
    std::unordered_map<std::string, std::shared_ptr<sv_t>>* sv_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    
    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    std::unordered_map<std::string, std::vector<sv_t*>> hp_indels_by_ref_hp_range;
    for (sv_t* hp_indel : hp_indels) {
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indel, contig_seq, contig_len);
        std::string ref_hp_range_key = std::to_string(ref_hp_range.beg) + ":" + std::to_string(ref_hp_range.end);
        hp_indels_by_ref_hp_range[ref_hp_range_key].push_back(hp_indel);
    }

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    for (auto& kv : hp_indels_by_ref_hp_range) {
        std::vector<sv_t*>& hp_indels_in_range = kv.second;
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indels_in_range[0], contig_seq, contig_len);
        genotype_hp_indels_group(hp_indels_in_range, ref_hp_range, bam_file, contig_seq, contig_len, stats, config, aligner,
            *mateseqs_w_mapq_chr, evidence_logger, reassign_evidence, evidence_map, *sv_map);
    }
    depth_filter_indel(contig_name, hp_indels, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, hp_indels, bam_file, config, stats);

    release_mates(contig_id);
}

#endif // GENOTYPE_HP_INDELS_H
