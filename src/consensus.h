#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "htslib/hts.h"

struct positional_consensus_t {
    std::string seq, qual;
    std::vector<int> coverage, max_base_freq;
};

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};

inline bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.qual != bs2.qual) return bs1.qual < bs2.qual;
    return bs1.freq < bs2.freq;
}

inline int base_to_index(char base) {
    base = std::toupper((unsigned char) base);
    if (base == 'A') return 0;
    if (base == 'C') return 1;
    if (base == 'G') return 2;
    if (base == 'T') return 3;
    return -1;
}

struct strand_base_scores_t {
    double quals[2][4] = {};
    int coverage[2] = {};

    void add(char base, int qual, bool reverse) {
        int index = base_to_index(base);
        if (index < 0) return;
        coverage[reverse]++;
        quals[reverse][index] += qual;
    }

    char best_base(char current_base) const {
        double scores[4] = {};
        for (int strand = 0; strand < 2; strand++) {
            if (coverage[strand] == 0) continue;
            double denominator = std::sqrt(double(coverage[strand]));
            for (int base = 0; base < 4; base++) scores[base] += quals[strand][base]/denominator;
        }
        int current_index = base_to_index(current_base);
        double best_score = current_index < 0 ? 0 : scores[current_index];
        char best = current_base;
        for (int base = 0; base < 4; base++) {
            if (scores[base] > best_score) {
                best_score = scores[base];
                best = "ACGT"[base];
            }
        }
        return best;
    }
};

inline positional_consensus_t build_positional_consensus(const std::vector<std::string>& seqs, const std::vector<const uint8_t*>& quals,
    const std::vector<hts_pos_t>& read_start_offsets, bool subtract_opposing_qualities = false, const std::vector<bool>& read_is_reverse = {}) {

    positional_consensus_t consensus;
    if (seqs.size() != quals.size() || seqs.size() != read_start_offsets.size()) return consensus;

    hts_pos_t consensus_len = 0;
    for (int i = 0; i < seqs.size(); i++) {
        consensus_len = std::max<hts_pos_t>(consensus_len, read_start_offsets[i] + seqs[i].length());
    }
    consensus.seq.assign(consensus_len, 'N');
    consensus.qual.assign(consensus_len, '!');
    consensus.coverage.assign(consensus_len, 0);
    consensus.max_base_freq.assign(consensus_len, 0);

    for (int i = 0; i < consensus_len; i++) {
        base_score_t base_scores[4] = {base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T')};
        strand_base_scores_t strand_scores;
        for (int j = 0; j < seqs.size(); j++) {
            int qpos = i - read_start_offsets[j];
            if (qpos < 0 || qpos >= seqs[j].length()) continue;

            int base_idx = base_to_index(seqs[j][qpos]);
            if (base_idx < 0) continue;
            consensus.coverage[i]++;
            base_scores[base_idx].freq++;
            int qual = quals[j][qpos];
            // BAM uses 255 for missing qualities; clip consensuses treat these as zero support.
            if (subtract_opposing_qualities && qual == 255) qual = 0;
            base_scores[base_idx].qual += qual;
            if (!read_is_reverse.empty()) strand_scores.add(seqs[j][qpos], qual, read_is_reverse[j]);
        }

        base_score_t best_base_score = std::max(std::max(base_scores[0], base_scores[1]), std::max(base_scores[2], base_scores[3]));
        // Keep the original positional base when its weighted score ties for best.
        if (!read_is_reverse.empty()) best_base_score = base_scores[base_to_index(strand_scores.best_base(best_base_score.base))];
        if (best_base_score.freq > 0) {
            consensus.seq[i] = best_base_score.base;
            int qual = best_base_score.qual;
            if (subtract_opposing_qualities) {
                // Clip consensuses subtract all opposing A/C/G/T support before capping the quality.
                for (const base_score_t& base_score : base_scores) {
                    if (base_score.base != best_base_score.base) qual -= base_score.qual;
                }
            }
            consensus.qual[i] = std::max(0, std::min(qual, 40)) + 33;
        }
        consensus.max_base_freq[i] = best_base_score.freq;
    }
    return consensus;
}

inline std::string build_full_consensus_seq(std::vector<std::string>& seqs, std::vector<uint8_t*>& quals,
    std::vector<hts_pos_t> read_start_offsets, int& lowq_prefix, int& lowq_suffix, std::string& consensus_qual, bool subtract_opposing_qualities = false) {

    std::vector<const uint8_t*> const_quals(quals.begin(), quals.end());
    positional_consensus_t consensus = build_positional_consensus(seqs, const_quals, read_start_offsets, subtract_opposing_qualities);
    consensus_qual = consensus.qual;

    lowq_prefix = 0;
    lowq_suffix = 0;
    bool low_prefix_done = false;
    for (int i = 0; i < consensus.seq.length(); i++) {
        if (consensus.coverage[i] < 3) {
            if (!low_prefix_done && consensus.max_base_freq[i] < 2) {
                lowq_prefix = i + 1;
            } else if (lowq_suffix == 0 && consensus.max_base_freq[i] < 2) {
                lowq_suffix = consensus.seq.length() - i;
            }
        } else {
            low_prefix_done = true;
        }
    }
    return consensus.seq;
}

#endif
