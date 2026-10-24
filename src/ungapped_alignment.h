#ifndef UNGAPPED_ALIGNMENT_H_
#define UNGAPPED_ALIGNMENT_H_

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>

#include "../libs/exact_cache.h"
#include "htslib/hts.h"
#include "simd_macros.h"

struct suffix_prefix_aln_t {
    int overlap, score, mismatches;

    suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}

	double mismatch_rate() { return overlap > 0 ? double(mismatches)/overlap : 1; }
};

suffix_prefix_aln_t aln_suffix_prefix_perfect(const std::string& s1, const std::string& s2, int min_overlap = 1) {
    const char* _s1 = s1.c_str(),* _s2 = s2.c_str();
    int _s1_len = s1.length(), _s2_len = s2.length();

    int max_overlap = std::min(_s1_len, _s2_len);
    for (int i = max_overlap; i >= min_overlap; i--) {
    	if (strncmp(_s1+(_s1_len-i), _s2, i) == 0) {
			return suffix_prefix_aln_t(i, i, 0);
    	}
    }
    return suffix_prefix_aln_t(0, 0, 0);
}

suffix_prefix_aln_t aln_suffix_prefix_qual(const std::string& s1, const std::string& q1, const std::string& s2, const std::string& q2, int min_overlap, double max_seq_error, int serious_q = 40) {
    int s1_len = s1.length(), s2_len = s2.length();
    for (int overlap = std::min(s1_len, s2_len); overlap >= min_overlap; overlap--) {
        int s1_beg = s1_len-overlap;
        int mismatches = 0, max_mismatches = (int) (overlap*max_seq_error);
        bool rejected = false;
        for (int i = 0; i < overlap; i++) {
            if (toupper(s1[s1_beg+i]) == toupper(s2[i])) continue;
            // Missing qualities cannot make a mismatch serious.
            bool serious = s1_beg+i < q1.length() && i < q2.length() && q1[s1_beg+i]-33 >= serious_q && q2[i]-33 >= serious_q;
            if (serious || ++mismatches > max_mismatches) {
                rejected = true;
                break;
            }
        }
        if (!rejected) return suffix_prefix_aln_t(overlap, overlap-mismatches, mismatches);
    }
    return suffix_prefix_aln_t(0, 0, 0);
}

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
int number_of_mismatches_fast(const char* s1, const char* s2, int len, int max_mismatches) {
	// count number of mismatches between the first len characters of s1 and s2
	int n_mismatches = 0;
    const SIMD_INT_16* s1_it = (const SIMD_INT_16*) s1;
    const SIMD_INT_16* s2_it = (const SIMD_INT_16*) s2;
    int scaled_len = len/BYTES_PER_BLOCK_16;
    for (int i = 0; i < scaled_len; i++) {
        SIMD_INT_16 n1 = LOADU_INT_16(s1_it);
        SIMD_INT_16 n2 = LOADU_INT_16(s2_it);
        n_mismatches += BYTES_PER_BLOCK_16 - COUNT_EQUAL_BYTES_16(n1, n2);
        s1_it++;
        s2_it++;
		if (n_mismatches > max_mismatches) return n_mismatches;
    }
    for (int i = scaled_len*BYTES_PER_BLOCK_16; i < len; i++) {
        if (s1[i] != s2[i]) n_mismatches++;
    }
    return n_mismatches;
}

struct ungapped_aln_t {
    int query_begin;
    int query_end;
    int ref_begin;
    int ref_end;
    int mismatches;
    int score;

    ungapped_aln_t(int query_begin, int query_end, int ref_begin, int ref_end, int mismatches, int score)
        : query_begin(query_begin), query_end(query_end), ref_begin(ref_begin), ref_end(ref_end), mismatches(mismatches), score(score) {}

	double mismatch_rate() { return (query_end - query_begin) > 0 ? double(mismatches)/(query_end - query_begin) : 1; }
};

// Score a fixed placement (ref coordinate of query[0]) without gaps or shifting.
// Total query clipping on each side, including reference overhang, cannot exceed its limit.
// Invalid or empty placements return an empty interval with score INT32_MIN.
ungapped_aln_t fixed_ungapped_aln(const char* query, int query_len, const char* ref, int ref_len, hts_pos_t placement,
                                int match_score = 1, int mismatch_score = -4, int max_left_clip = 0, int max_right_clip = 0) {
    ungapped_aln_t best_aln(0, 0, 0, 0, 0, INT32_MIN);
    if (query_len <= 0 || ref_len <= 0 || placement <= -hts_pos_t(query_len) || placement >= ref_len) return best_aln;
    max_left_clip = std::max(0, max_left_clip);
    max_right_clip = std::max(0, max_right_clip);
    int begin = std::max<hts_pos_t>(0, -placement);
    int end = std::min<hts_pos_t>(query_len, ref_len - placement);
    if (begin > max_left_clip || query_len-end > max_right_clip) return best_aln;

    if (begin == max_left_clip && query_len-end == max_right_clip) {
        int mismatches = number_of_mismatches_fast(query + begin, ref + placement + begin, end-begin, INT32_MAX);
        int score = (end-begin-mismatches)*match_score + mismatches*mismatch_score;
        return ungapped_aln_t(begin, end, placement+begin, placement+end, mismatches, score);
    }

    int prefix_score = 0, prefix_mismatches = 0;
    int min_prefix_score = 0, min_prefix_mismatches = 0, best_begin = begin;
    for (int pos = begin; pos < end; pos++) {
        // Eligible starts never expire. On equal prefix scores, keep the earlier (longer) start.
        if (pos <= max_left_clip && prefix_score < min_prefix_score) {
            min_prefix_score = prefix_score;
            min_prefix_mismatches = prefix_mismatches;
            best_begin = pos;
        }
        bool mismatch = query[pos] != ref[placement+pos];
        prefix_score += mismatch ? mismatch_score : match_score;
        prefix_mismatches += mismatch;
        int query_end = pos+1;
        if (query_len-query_end > max_right_clip) continue;
        int score = prefix_score-min_prefix_score;
        if (score > best_aln.score || (score == best_aln.score && query_end-best_begin > best_aln.query_end-best_aln.query_begin)) {
            best_aln = ungapped_aln_t(best_begin, query_end, placement+best_begin, placement+query_end, prefix_mismatches-min_prefix_mismatches, score);
        }
    }
    return best_aln;
}

// Find the best placement, allowing up to max_ref_overflow reference overhang per side.
// With max_query_clip == 0, retain the full overlap and the legacy placement/tie behavior.
// Otherwise each side's TOTAL clipping, including overhang, is limited by max_query_clip.
ungapped_aln_t best_ungapped_aln(const char* query, int query_len, const char* ref, int ref_len,
                                 int max_ref_overflow = 0, int match_score = 1, int mismatch_score = -4, int max_query_clip = 0) {
    if (query_len <= 0 || ref_len <= 0) return ungapped_aln_t(0, 0, 0, 0, 0, 0);

    surveyor_cache::key_builder_t key_builder;
    key_builder.add(query_len);
    key_builder.add_bytes(query, query_len);
    key_builder.add(ref_len);
    key_builder.add_bytes(ref, ref_len);
    key_builder.add(max_ref_overflow);
    key_builder.add(match_score);
    key_builder.add(mismatch_score);
    key_builder.add(max_query_clip);
    std::string cache_key = key_builder.take();
    static thread_local surveyor_cache::exact_cache_t<ungapped_aln_t> cache("SURVEYOR_UNGAPPED_CACHE_LIMIT", 16384);
    const ungapped_aln_t* cached_aln = cache.find(cache_key);
    if (cached_aln) return *cached_aln;

    max_ref_overflow = std::max(0, max_ref_overflow);
    max_query_clip = std::max(0, max_query_clip);

    ungapped_aln_t best_aln(0, 0, 0, 0, 0, INT32_MIN);
    // These are the same placements as the former padded-reference search, including its
    // single placement when the query is longer than the reference plus both overhangs.
    hts_pos_t last_placement = std::max<hts_pos_t>(-max_ref_overflow, hts_pos_t(ref_len)+max_ref_overflow-query_len);
    for (hts_pos_t placement = -max_ref_overflow; placement <= last_placement; placement++) {
        int left_overhang = std::max<hts_pos_t>(0, -placement);
        int right_overhang = std::max<hts_pos_t>(0, placement+query_len-ref_len);
        if (max_query_clip > 0 && (left_overhang > max_ref_overflow || right_overhang > max_ref_overflow)) continue;
        int left_clip = max_query_clip > 0 ? max_query_clip : left_overhang;
        int right_clip = max_query_clip > 0 ? max_query_clip : right_overhang;
        ungapped_aln_t aln = fixed_ungapped_aln(query, query_len, ref, ref_len, placement, match_score, mismatch_score, left_clip, right_clip);
        // Preserve zero-overlap candidates from the legacy search only when clipping is disabled.
        if (max_query_clip == 0 && (placement <= -hts_pos_t(query_len) || placement >= ref_len)) {
            int ref_begin = std::max<hts_pos_t>(0, placement);
            aln = ungapped_aln_t(left_overhang, left_overhang, ref_begin, ref_begin, 0, 0);
        }
        if (aln.score > best_aln.score || (max_query_clip > 0 && aln.score == best_aln.score && aln.query_end-aln.query_begin > best_aln.query_end-best_aln.query_begin)) {
            best_aln = aln;
        }
        if (best_aln.score == query_len*match_score && (max_query_clip == 0 || (match_score >= 0 && mismatch_score <= match_score && best_aln.query_end-best_aln.query_begin == query_len))) break;
    }
    cache.store(std::move(cache_key), best_aln);
    return best_aln;
}

suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {

	if (max_seq_error == 0.0 || max_mismatches == 0) return aln_suffix_prefix_perfect(s1, s2, min_overlap);

    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < (int) s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;	// this means that the suffix of s1 I am considering is larger than the the whole s2,
                                                    // therefore it wouldn't be a suffix-prefix alignment   

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
		int max_acceptable_mm = sp_len * max_seq_error;
		int mismatches = number_of_mismatches_fast(s1_suffix, s2_prefix, sp_len, std::min(max_acceptable_mm, max_mismatches));
		int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;
        if (best_score < score && mismatches <= max_acceptable_mm && mismatches <= max_mismatches) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}


#endif /* UNGAPPED_ALIGNMENT_H_ */
