#ifndef SW_UTILS_H
#define SW_UTILS_H

#include <cstdint>
#include <iostream>
#include <memory>

#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"
#include "simd_macros.h"
#include "sam_utils.h"
#include "types.h"
#include "utils.h"

int get_left_clip_size(StripedSmithWaterman::Alignment& aln) {
	if (aln.cigar.empty()) return 0;
	return cigar_int_to_op(aln.cigar[0]) == 'S' ? cigar_int_to_len(aln.cigar[0]) : 0;
}
int get_right_clip_size(StripedSmithWaterman::Alignment& aln) {
	if (aln.cigar.empty()) return 0;
	return cigar_int_to_op(aln.cigar[aln.cigar.size()-1]) == 'S' ? cigar_int_to_len(aln.cigar[aln.cigar.size()-1]) : 0;
}
int get_unclipped_start(StripedSmithWaterman::Alignment& aln) {
	return aln.ref_begin - get_left_clip_size(aln);
}
int get_unclipped_end(StripedSmithWaterman::Alignment& aln) {
	return aln.ref_end + get_right_clip_size(aln);
}
bool is_left_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_left_clip_size(aln) >= min_clip_len;
}
bool is_right_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return get_right_clip_size(aln) >= min_clip_len;
}
bool is_clipped(StripedSmithWaterman::Alignment& aln, int min_clip_len = 1) {
	return is_left_clipped(aln, min_clip_len) || is_right_clipped(aln, min_clip_len);
}

size_t gcd(size_t a, size_t b) {
    return (b == 0) ? a : gcd(b, a % b);
}
size_t lcm(size_t a, size_t b) {
    return a * b / gcd(a, b);
}

// Reference implementation of Smith-Waterman-Gotoh algorithm for validation
int* smith_waterman_gotoh_scalar_reference(
	const char* ref, int ref_len, 
    const char* query, int query_len,
    int match_score, int mismatch_penalty, 
    int gap_open, int gap_extend) {
		
	// Using a very small number for -Infinity, but safe from underflow when adding gap costs
	const int NEG_INF = -1e9; 
    
	if (query_len <= 0 || ref_len <= 0) return nullptr;

    std::vector<std::vector<int>> H(ref_len + 1, std::vector<int>(query_len + 1, 0));
    std::vector<std::vector<int>> RG(ref_len + 1, std::vector<int>(query_len + 1, NEG_INF));
    std::vector<std::vector<int>> QG(ref_len + 1, std::vector<int>(query_len + 1, NEG_INF));

    int* prefix_scores = (int*)calloc((size_t)query_len, sizeof(int));
    if (!prefix_scores) return nullptr;

    for (int i = 1; i <= ref_len; ++i) {
        for (int j = 1; j <= query_len; ++j) {

            // 1. Calculate Scores
            int score = (toupper(ref[i-1]) == toupper(query[j-1])) ? match_score : mismatch_penalty;
            int diagonal_val = H[i-1][j-1] + score;

            // RG (Vertical): Gap in Query (Deletion)
            int open_vertical   = H[i-1][j] + gap_open;
            int extend_vertical = RG[i-1][j] + gap_extend;
            RG[i][j] = std::max(open_vertical, extend_vertical);

            // QG (Horizontal): Gap in Ref (Insertion)
            int open_horizontal   = H[i][j-1] + gap_open;
            int extend_horizontal = QG[i][j-1] + gap_extend;
            QG[i][j] = std::max(open_horizontal, extend_horizontal);

            // 2. Update H (Propagating Score)
            // Mathematically, H must include ALL paths to propagate correctly to future cells.
            int val = std::max(diagonal_val, RG[i][j]);
            val = std::max(val, QG[i][j]);
            val = std::max(val, 0); 
            H[i][j] = val;

            // 3. Update Prefix Score (Reporting Score)
            // MATCHING FAST CODE BEHAVIOR:
            // Your Fast code updates prefix_scores using 'H_curr'. 
            // In the Fast code, 'H_curr' is calculated as max(Match, RG). 
            // It does NOT include QG (which is calculated later in the scalar loop).
            
            // Therefore, to reproduce the Fast Code's result, we must exclude QG here.
            
            int reportable_score = std::max(diagonal_val, RG[i][j]); // Exclude QG[i][j]
            reportable_score = std::max(reportable_score, 0);
            
            prefix_scores[j-1] = std::max(prefix_scores[j-1], reportable_score);
        }
    }

    return prefix_scores;
}

struct sw_gotoh_workspace_t {
	SW_SCORE_INT_16* H = NULL;
	SW_SCORE_INT_16* RG = NULL;
	SW_SCORE_INT_16* QG = NULL;
	SW_SCORE_INT_16* prefix_scores = NULL;
	SW_SCORE_INT_16* profile_storage = NULL;
	SW_SCORE_INT_16* profile[5] = {NULL, NULL, NULL, NULL, NULL};
	int stride = 0;
	int query_len_rounded = 0;

	~sw_gotoh_workspace_t() {
		free(H);
		free(RG);
		free(QG);
		free(prefix_scores);
		free(profile_storage);
	}

	sw_gotoh_workspace_t() = default;
	sw_gotoh_workspace_t(const sw_gotoh_workspace_t&) = delete;
	sw_gotoh_workspace_t& operator=(const sw_gotoh_workspace_t&) = delete;

	void ensure_capacity(int new_query_len_rounded) {
		int new_stride = new_query_len_rounded + INT_PER_BLOCK_16;
		if (new_stride <= stride) {
			query_len_rounded = new_query_len_rounded;
			return;
		}

		free(H);
		free(RG);
		free(QG);
		free(prefix_scores);
		free(profile_storage);
		H = RG = QG = prefix_scores = profile_storage = NULL;
		for (int i = 0; i < 5; i++) profile[i] = NULL;

		size_t alignment = lcm(BYTES_PER_BLOCK_16, sizeof(void*));
		int p1 = posix_memalign(reinterpret_cast<void**>(&H), alignment, 2*(new_stride * sizeof(SW_SCORE_INT_16)));
		int p2 = posix_memalign(reinterpret_cast<void**>(&RG), alignment, 2*(new_stride * sizeof(SW_SCORE_INT_16)));
		int p3 = posix_memalign(reinterpret_cast<void**>(&QG), alignment, 2*(new_stride * sizeof(SW_SCORE_INT_16)));
		int p4 = posix_memalign(reinterpret_cast<void**>(&prefix_scores), alignment, new_stride * sizeof(SW_SCORE_INT_16));
		int p5 = posix_memalign(reinterpret_cast<void**>(&profile_storage), alignment, 5 * new_stride * sizeof(SW_SCORE_INT_16));
		if (p1 || p2 || p3 || p4 || p5) {
			std::cerr << "Error allocating aligned memory of size " << (2*new_stride * sizeof(SW_SCORE_INT_16)) << std::endl;
			free(H);
			free(RG);
			free(QG);
			free(prefix_scores);
			free(profile_storage);
			H = RG = QG = prefix_scores = profile_storage = NULL;
			stride = query_len_rounded = 0;
			return;
		}
		for (int i = 0; i < 5; i++) profile[i] = profile_storage + i*new_stride;
		stride = new_stride;
		query_len_rounded = new_query_len_rounded;
	}
};

// This function returns an array such that the i-th element is the best local alignment score of query[0..i] against ref.
// The exception is that we are NOT interested in alignments that end in a gap in the query, so those are not considered.
// The returned pointer is owned by workspace and remains valid until the next call using that workspace.
SW_SCORE_INT_16* smith_waterman_gotoh(const char* ref, int ref_len, const char* query, int query_len,
		int match_score, int mismatch_penalty, int gap_open, int gap_extend, sw_gotoh_workspace_t& workspace) {

	// turn query_len+1 into a multiple of INT_PER_BLOCK
	int query_len_rounded = (query_len+INT_PER_BLOCK_16-1)/INT_PER_BLOCK_16*INT_PER_BLOCK_16;
	workspace.ensure_capacity(query_len_rounded);
	if (workspace.prefix_scores == NULL) return NULL;

	const char* alphabet = "NACGT";
	const int alphabet_size = 5;
	int stride = workspace.stride;
	for (int i = 0; i < alphabet_size; i++) {
		SW_SCORE_INT_16* profile = workspace.profile[i];
		profile[0] = 0;
		for (int j = 1; j <= query_len; j++) {
			profile[j] = (toupper(query[j-1]) == alphabet[i]) ? match_score : mismatch_penalty;
		}
		std::fill(profile+query_len+1, profile+query_len_rounded+1, 0);
	}

	SW_SCORE_INT_16* H_prev = workspace.H, *H_curr = workspace.H+stride;
	SW_SCORE_INT_16* RG_prev = workspace.RG, *RG_curr = workspace.RG+stride;
	SW_SCORE_INT_16* QG_prev = workspace.QG, *QG_curr = workspace.QG+stride;
	SW_SCORE_INT_16* prefix_scores = workspace.prefix_scores;

	std::fill(H_prev, H_prev+stride, 0);
	std::fill(RG_prev, RG_prev+stride, 0);
	std::fill(QG_prev, QG_prev+stride, 0);
	H_curr[0] = QG_curr[0] = RG_curr[0] = 0;

	SIMD_INT_16 gap_open_v = SET1_INT_16(gap_open);
	SIMD_INT_16 gap_open_v_pos = SET1_INT_16(-gap_open);
	SIMD_INT_16 gap_extend_v = SET1_INT_16(gap_extend);
	SIMD_INT_16 zero_v = SET1_INT_16(0);

	std::fill(prefix_scores, prefix_scores+stride, 0);
	for (int i = 1; i <= ref_len; i++) {
		SW_SCORE_INT_16* ref_profile = workspace.profile[0];
		switch (toupper(ref[i-1])) {
			case 'A': ref_profile = workspace.profile[1]; break;
			case 'C': ref_profile = workspace.profile[2]; break;
			case 'G': ref_profile = workspace.profile[3]; break;
			case 'T': ref_profile = workspace.profile[4]; break;
		}

		for (int j = 0; j < query_len_rounded; j += INT_PER_BLOCK_16) {
			SIMD_INT_16 H_up_v = LOAD_INT_16((SIMD_INT_16*)&H_prev[j]);
			SIMD_INT_16 RG_up_v = LOAD_INT_16((SIMD_INT_16*)&RG_prev[j]);
			SIMD_INT_16 QG_up_v = LOAD_INT_16((SIMD_INT_16*)&QG_prev[j]);
			SIMD_INT_16 max_H_QG = MAX_INT_16(H_up_v, QG_up_v);
			SIMD_INT_16 max_H_QG_w_gap_open = ADD_INT_16(gap_open_v, max_H_QG);
			SIMD_INT_16 RG_w_gap_extend = ADD_INT_16(gap_extend_v, RG_up_v);
			SIMD_INT_16 E_curr_v = MAX_INT_16(max_H_QG_w_gap_open, RG_w_gap_extend);
			STORE_INT_16((SIMD_INT_16*)&RG_curr[j], E_curr_v);

			SIMD_INT_16 max_H_QG_RG = MAX_INT_16(max_H_QG, RG_up_v);
			SIMD_INT_16 H_curr_v = MAX_INT_16(max_H_QG_RG, zero_v);

			SIMD_INT_16 profile_curr = LOADU_INT_16((SIMD_INT_16*)&ref_profile[j+1]);
			H_curr_v = ADD_INT_16(H_curr_v, profile_curr);

			SIMD_INT_16 prefix_v = LOAD_INT_16((SIMD_INT_16*)&prefix_scores[j]);
			prefix_v = MAX_INT_16(prefix_v, H_curr_v);
			STOREU_INT_16((SIMD_INT_16*)&H_curr[j+1], H_curr_v);
			STORE_INT_16((SIMD_INT_16*)&prefix_scores[j], prefix_v);

			// we can skip the scalar QG computation if 
			// 1) all elements in the block are <= -gap_extend (therefore no point in further extending that gap); and
			// 2) all elements in the block in H, previous row, are <= -gap_open (therefore opening a new gap makes no sense)
			if (QG_curr[j] > -gap_extend || H_curr[j] > -gap_open || CMP_GT_INT_16_EXCEPT_LAST(H_curr_v, gap_open_v_pos)) {
				alignas(BYTES_PER_BLOCK_16) SW_SCORE_INT_16 temp_h[INT_PER_BLOCK_16];
				STORE_INT_16((SIMD_INT_16*)temp_h, ADD_INT_16(H_curr_v, gap_open_v));
				QG_curr[j+1] = std::max(gap_open + H_curr[j], gap_extend + QG_curr[j]);
				for (int k = 1; k < INT_PER_BLOCK_16; k++) {
					QG_curr[j+k+1] = std::max(temp_h[k-1], (SW_SCORE_INT_16)(gap_extend + QG_curr[j+k]));
				}
			} else {
				STOREU_INT_16((SIMD_INT_16*)&QG_curr[j+1], zero_v);
			}
		}

		std::swap(H_prev, H_curr);
		std::swap(RG_prev, RG_curr);
		std::swap(QG_prev, QG_curr);
	}

	return prefix_scores;
}

// Compatibility wrapper. Callers own the returned pointer and must free it.
SW_SCORE_INT_16* smith_waterman_gotoh(const char* ref, int ref_len, const char* query, int query_len,
		int match_score, int mismatch_penalty, int gap_open, int gap_extend) {
	static thread_local sw_gotoh_workspace_t workspace;
	SW_SCORE_INT_16* workspace_prefix_scores = smith_waterman_gotoh(ref, ref_len, query, query_len,
		match_score, mismatch_penalty, gap_open, gap_extend, workspace);
	if (workspace_prefix_scores == NULL) return NULL;

	int query_len_rounded = (query_len+INT_PER_BLOCK_16-1)/INT_PER_BLOCK_16*INT_PER_BLOCK_16;
	int stride = query_len_rounded + INT_PER_BLOCK_16;
	SW_SCORE_INT_16* prefix_scores = NULL;
	size_t alignment = lcm(BYTES_PER_BLOCK_16, sizeof(void*));
	int p = posix_memalign(reinterpret_cast<void**>(&prefix_scores), alignment, stride * sizeof(SW_SCORE_INT_16));
	if (p) {
		std::cerr << "Error allocating aligned memory of size " << (stride * sizeof(SW_SCORE_INT_16)) << std::endl;
		return NULL;
	}
	memcpy(prefix_scores, workspace_prefix_scores, stride * sizeof(SW_SCORE_INT_16));
	return prefix_scores;
}

std::vector<int> compute_prefix_scores(std::vector<uint32_t> cigar, int query_len, int match_score, int mismatch_score, 
						int gap_open_score, int gap_extend_score) {
	
	std::vector<int> prefix_scores(query_len + 1, 0); // Stores scores for each prefix up to query_len.
    int score = 0, query_prefix_len = 0;

	for (int i = 0; i < cigar.size() && query_prefix_len < query_len; ++i) {
        int op = cigar_int_to_op(cigar[i]);
        int len = cigar_int_to_len(cigar[i]);

		for (int k = 0; k < len; ++k) {
            if (op == '=') {
                score += match_score;
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'X') {
                score += mismatch_score;
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'I') {
                if (k == 0) {
                    score += gap_open_score;
                } else {
                    score += gap_extend_score;
                }
                prefix_scores[++query_prefix_len] = score;
            } else if (op == 'D') {
                if (k == 0) {
                    score += gap_open_score;
                } else {
                    score += gap_extend_score;
                }
                // No change to query_prefix_len since deletion skips query sequence.
            } else if (op == 'S') { // Soft-clipping
                // Soft-clipped bases don't change the score, but they do increment query_prefix_len.
                prefix_scores[++query_prefix_len] = score;
            }
        }
	}

	return prefix_scores;
}

std::vector<int> compute_suffix_scores(std::vector<uint32_t> cigar, int query_len, int match_score, int mismatch_score, 
						int gap_open_score, int gap_extend_score) {

	std::vector<uint32_t> cigar_rev(cigar.rbegin(), cigar.rend());
	std::vector<int> suffix_scores = compute_prefix_scores(cigar_rev, query_len, match_score, mismatch_score, gap_open_score, gap_extend_score);
	std::reverse(suffix_scores.begin(), suffix_scores.end());
	return suffix_scores;
}

std::pair<int, int> find_aln_prefix_score(std::vector<uint32_t>& cigar, int ref_prefix_len, int match_score, int mismatch_score,
						int gap_open_score, int gap_extend_score) {
	
	if (cigar.empty()) return {0, 0};

	int score = 0, query_prefix_len = 0, i = 0;
	if (cigar_int_to_op(cigar[0]) == 'S') {
		query_prefix_len = cigar_int_to_len(cigar[0]);
		i++;
	}

	for (int j = 0; i < cigar.size() && j < ref_prefix_len; i++) {
		int op = cigar_int_to_op(cigar[i]);
		int len = cigar_int_to_len(cigar[i]);
		if (op == '=') {
			len = std::min(len, ref_prefix_len-j);
			score += len*match_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'X') {
			len = std::min(len, ref_prefix_len-j);
			score += len*mismatch_score;
			query_prefix_len += len;
			j += len;
		} else if (op == 'I') {
			score += gap_open_score + (len-1)*gap_extend_score;
			query_prefix_len += len;
		} else if (op == 'D') {
			len = std::min(len, ref_prefix_len-j);
			score += gap_open_score + (len-1)*gap_extend_score;
			j += len;
		} else if (op == 'S') {
			query_prefix_len += len;
		}
	}
	return {score, query_prefix_len};
}

std::pair<int, int> find_aln_suffix_score(std::vector<uint32_t> cigar, int ref_suffix_len, int match_score, int mismatch_score,
						  int gap_open_score, int gap_extend_score) {
	std::vector<uint32_t> rev_cigar(cigar.rbegin(), cigar.rend());
	return find_aln_prefix_score(rev_cigar, ref_suffix_len, match_score, mismatch_score, gap_open_score, gap_extend_score);
}


// the (min) length of the query prefix that cover ref_len bases on the reference
int query_prefix_len(std::vector<uint32_t>& cigar, int ref_len, int& extra_deletion) {
    int query_len = 0;
    int ref_count = 0;
    for (auto it = cigar.begin(); it != cigar.end(); ++it) {
        char op = cigar_int_to_op(*it);
        int len = cigar_int_to_len(*it);
        if (op == '=' || op == 'X' || op == 'D') {
            ref_count += len;
        }
        if (op == '=' || op == 'X' || op == 'I' || op == 'S') {
            query_len += len;
        }
        if (ref_count >= ref_len) {
            int overshoot = ref_count - ref_len;
            if (op == '=' || op == 'X') {
                query_len -= overshoot;
            } else if (op == 'D') {
                extra_deletion = overshoot;
            }
            break;
        }
    }
    return query_len;
}
int query_suffix_len(std::vector<uint32_t>& cigar, int ref_len, int& extra_deletion) {
    auto cigar_rev = std::vector<uint32_t>(cigar.rbegin(), cigar.rend());
    return query_prefix_len(cigar_rev, ref_len, extra_deletion);
}


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

// Finds the best ungapped placement of query against ref, allowing query overhang
// on either side of the reference. Query bases outside ref count as clips, not mismatches.
// Each side of the query may extend by at most max_ref_overflow bases beyond the reference.
// query_begin/query_end and ref_begin/ref_end describe the overlapped intervals [begin, end).
// score is the ungapped alignment score on the overlapped interval only.
ungapped_aln_t best_ungapped_aln(const char* query, int query_len, const char* ref, int ref_len,
                                 int max_ref_overflow = 0, int match_score = 1, int mismatch_score = -4) {
    if (query_len <= 0) {
        return ungapped_aln_t(0, 0, 0, 0, 0, 0);
    }
    if (ref_len <= 0) {
        return ungapped_aln_t(0, 0, 0, 0, 0, 0);
    }

	if (max_ref_overflow < 0) max_ref_overflow = 0;

    const int query_len_rounded = (query_len + BYTES_PER_BLOCK_16 - 1)/BYTES_PER_BLOCK_16*BYTES_PER_BLOCK_16;
    const int query_pad_len = query_len_rounded - query_len;
    const int ref_pad_len = max_ref_overflow;
    const char query_pad_char = '\0';
    const char ref_pad_char = '\1';

    std::vector<char> padded_query(query_len_rounded, query_pad_char);
    std::memcpy(padded_query.data(), query, query_len);

	int padded_ref_len = std::max(ref_len + 2*ref_pad_len + query_pad_len, query_len_rounded);
    std::vector<char> padded_ref(padded_ref_len, ref_pad_char);
    std::memcpy(padded_ref.data() + ref_pad_len, ref, ref_len);

    ungapped_aln_t best_aln(0, 0, 0, 0, 0, INT32_MIN);
    int best_score = INT32_MIN;
    for (int i = 0; i <= padded_ref.size()-query_len_rounded; i++) {

		int valid_overlap = overlap(ref_pad_len, ref_pad_len + ref_len, i, i + query_len); // the remaining query_len_rounded - mismatches are guaranteed to be mismatches
		int guaranteed_mismatches = query_len_rounded - valid_overlap;

		int max_mismatches = query_len_rounded;
		if (best_score != INT32_MIN) {
			// max_mismatches = the number of mismatches that would make this alignment's score equal to the current best score, plus the guaranteed mismatches
			int max_possible_score = valid_overlap * match_score;
			max_mismatches = (max_possible_score - best_score) / (match_score - mismatch_score) + guaranteed_mismatches; 
		}
		int mismatches = number_of_mismatches_fast(padded_query.data(), padded_ref.data() + i, query_len_rounded, max_mismatches);
		if (mismatches > max_mismatches) continue; // cannot be better than current best
		
		int actual_mismatches = mismatches - guaranteed_mismatches;

		int score = (valid_overlap - actual_mismatches) * match_score + actual_mismatches * mismatch_score;
		if (score > best_score) {
			best_score = score;
			int query_begin = i < ref_pad_len ? ref_pad_len - i : 0;
			int query_end = query_begin + valid_overlap;
			int ref_begin = std::max(0, i - ref_pad_len);
			int ref_end = ref_begin + valid_overlap;
			best_aln = ungapped_aln_t(query_begin, query_end, ref_begin, ref_end, actual_mismatches, score);
		}
		if (best_score == query_len * match_score) break; // cannot do better than a perfect match
    }
    return best_aln;
}

suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {

	if (max_seq_error == 0.0 || max_mismatches == 0) return aln_suffix_prefix_perfect(s1, s2, min_overlap);

    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
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


std::vector<StripedSmithWaterman::Alignment> get_best_alns(char* contig_seq, hts_pos_t remap_start, hts_pos_t remap_len, char* query, StripedSmithWaterman::Aligner& aligner) {
	
	contig_seq += remap_start;
	
	std::vector<StripedSmithWaterman::Alignment> best_alns;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	aligner.Align(query, contig_seq, remap_len, filter, &aln, 15);
	best_alns.push_back(aln);
	int best_score = aln.sw_score;
	int new_start = aln.ref_begin+5;

	if (new_start >= remap_len) return best_alns;

	aligner.Align(query, contig_seq+new_start, remap_len-new_start, filter, &aln, 15);
	while (aln.sw_score == best_score) {
		aln.ref_begin += new_start;
		aln.ref_end += new_start;
		best_alns.push_back(aln);
		new_start += aln.ref_begin+5;
		if (new_start >= remap_len) break;
		aligner.Align(query, contig_seq+new_start, remap_len-new_start, filter, &aln, 15);
	}
	return best_alns;
}

void set_junction_remap_ref_range(std::shared_ptr<sv_t> sv, hts_pos_t remap_ref_beg, hts_pos_t remap_ref_end, stats_t& stats, config_t& config) {
		
	hts_pos_t island_padding = std::max<hts_pos_t>(0, stats.read_len - 2*config.min_clip_len);

	hts_pos_t island_beg = std::min(sv->start, sv->end);
	hts_pos_t island_end = std::max(sv->start, sv->end);

	for (const auto& aux_indel : sv->aux_indels) {
		island_beg = std::min(island_beg, std::min(aux_indel->start, aux_indel->end));
		island_end = std::max(island_end, std::max(aux_indel->start, aux_indel->end));
	}
	for (const auto& aux_snp : sv->aux_snps) {
		island_beg = std::min(island_beg, aux_snp.pos);
		island_end = std::max(island_end, aux_snp.pos);
	}

	sv->junction_remap_ref_beg = std::max(remap_ref_beg, island_beg - island_padding);
	sv->junction_remap_ref_end = std::min(remap_ref_end, island_end + island_padding);
}

// Calculate the reference range covered by the alignment, represented by the cigar string, 
// when excluding the lowq prefix and suffix of the query sequence. 
std::pair<hts_pos_t, hts_pos_t> get_highq_ref_range(std::vector<uint32_t>& cigar, hts_pos_t ref_begin,
		hts_pos_t query_len, int lowq_prefix, int lowq_suffix) {
	hts_pos_t highq_query_beg = std::min<hts_pos_t>(std::max<hts_pos_t>(0, lowq_prefix), query_len);
	hts_pos_t highq_query_end = std::max<hts_pos_t>(0, query_len - std::max(0, lowq_suffix));
	if (highq_query_beg >= highq_query_end) return {ref_begin, ref_begin};

	hts_pos_t ref_pos = ref_begin;
	hts_pos_t query_pos = 0;
	hts_pos_t highq_ref_beg = HTS_POS_MAX;
	hts_pos_t highq_ref_end = HTS_POS_MIN;

	for (uint32_t c : cigar) {
		int len = cigar_int_to_len(c);
		char op = cigar_int_to_op(c);
		bool consumes_ref = op == '=' || op == 'X' || op == 'M' || op == 'D';
		bool consumes_query = op == '=' || op == 'X' || op == 'M' || op == 'I' || op == 'S';

		hts_pos_t op_ref_beg = ref_pos;
		hts_pos_t op_query_beg = query_pos;

		if (consumes_ref && consumes_query) {
			hts_pos_t overlap_beg = std::max(op_query_beg, highq_query_beg);
			hts_pos_t overlap_end = std::min(op_query_beg + len, highq_query_end);
			if (overlap_beg < overlap_end) {
				highq_ref_beg = std::min(highq_ref_beg, op_ref_beg + overlap_beg - op_query_beg);
				highq_ref_end = std::max(highq_ref_end, op_ref_beg + overlap_end - op_query_beg);
			}
		} else if (consumes_ref && highq_query_beg < query_pos && query_pos < highq_query_end) {
			highq_ref_beg = std::min(highq_ref_beg, op_ref_beg);
			highq_ref_end = std::max(highq_ref_end, op_ref_beg + len);
		}

		if (consumes_ref) ref_pos += len;
		if (consumes_query) query_pos += len;
	}

	if (highq_ref_beg == HTS_POS_MAX) return {ref_begin, ref_begin};
	return {highq_ref_beg, highq_ref_end};
}

// sometimes, we call this function after finding a main SV to find additional SVs. 
// In such case, we use main_sv to assign auxiliary SNPs and indels that would be otherwise lost or assigned to other SVs
// cigar is in the SSW library format
std::vector<std::shared_ptr<sv_t>> detect_svs_from_aln(std::vector<uint32_t>& cigar, hts_pos_t ref_begin, hts_pos_t ref_end, int aln_score,
	std::string contig_name, std::string junction_seq, std::shared_ptr<sv_t> main_sv, 
	int lowq_junction_prefix, int lowq_junction_suffix, stats_t& stats, config_t& config,
	bool record_junction_remap_ref_range = false) {
	// do not find aux SNPs in low-quality regions of the junction sequence

	std::shared_ptr<sv_t> input_main_sv = main_sv;
	std::vector<std::shared_ptr<sv_t>> svs;
    hts_pos_t current_pos = ref_begin;
	hts_pos_t junction_pos = 0;
	std::vector<snp_t> snps;

	for (int i = 0; i < cigar.size(); i++) {
		uint32_t c = cigar[i];
        int op_length = cigar_int_to_len(c);
        char op = cigar_int_to_op(c);

		hts_pos_t start = current_pos-1, end = (op == 'D') ? current_pos+op_length-1 : current_pos;
		std::string ins_seq = (op == 'I') ? junction_seq.substr(junction_pos, op_length) : "";
		auto left_part_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(ref_begin, start, junction_pos, aln_score);
		auto right_part_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(end, ref_end, junction_seq.length()-junction_pos-ins_seq.length(), aln_score);
		if (op == 'D') {
			std::shared_ptr<sv_t> sv = std::make_shared<deletion_t>(contig_name, start, end, "", nullptr, nullptr, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		} else if (op == 'I') {
			std::shared_ptr<sv_t> sv = std::make_shared<insertion_t>(contig_name, current_pos-1, current_pos-1, ins_seq, nullptr, nullptr, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		} else if (op == 'X') {
			for (int j = 0; j < op_length; j++) {	
				char alt_base = junction_seq[junction_pos+j];
				if (!is_genomic_base(alt_base)) continue;
				if (junction_pos+j >= lowq_junction_prefix && 
					junction_pos+j < (junction_seq.length() - lowq_junction_suffix)) {
					snps.push_back(snp_t(current_pos+j, alt_base));
				}
			}
		}

        if (op != 'I' && op != 'S') {
            current_pos += op_length;
        }
		if (op != 'D') {
			junction_pos += op_length;
		}
    }

	// Initially, all SVs are in aux.
	// Then, iteratively pick the largest SV >= min_sv_size among the remaining aux SVs as main, 
	// and assign nearby aux SVs to it.
	// Repeat until no aux SV is left or the largest remaining aux SV is smaller than min_sv_size.
	std::vector<std::shared_ptr<sv_t>> aux_svs;
	for (auto& sv : svs) {
		if (main_sv != nullptr && sv->unique_key() == main_sv->unique_key()) continue;
		aux_svs.push_back(sv); // keep main_svs as AUX as well, since it makes genotyping more accurate
	}

	std::vector<std::shared_ptr<sv_t>> main_svs;
	if (main_sv != nullptr) main_svs.push_back(main_sv);
	while (!aux_svs.empty()) {
		if (main_sv == nullptr) {
			// main_sv = the largest SV among aux_svs
			auto it = std::max_element(aux_svs.begin(), aux_svs.end(), [](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
				return a->svsize() < b->svsize();
			});
			main_sv = *it;
			*it = nullptr;
			if (main_sv->svsize() < config.min_sv_size) {
				break; // stop if the largest remaining SV is smaller than min_sv_size
			}
			main_svs.push_back(main_sv);
		}

		hts_pos_t island_start = main_sv->start, island_end = main_sv->end;
		bool added_to_island = true;
		while (added_to_island) {
			added_to_island = false;
			for (auto& sv : aux_svs) {
				if (sv == nullptr) continue;
				if (sv->svtype() == "INS" && sv->incomplete_ins_seq()) continue; // do not assign incomplete insertions as aux indels
				int dist = std::min(abs((int)sv->start - (int)island_end), abs((int)sv->end - (int)island_start));
				if (dist < stats.read_len - 2*config.min_clip_len) {
					island_start = std::min(island_start, sv->start);
					island_end = std::max(island_end, sv->end);
					main_sv->aux_indels.push_back(sv);
					sv = nullptr; // mark as assigned
					added_to_island = true;
				}
			}
		}
		aux_svs.erase(std::remove(aux_svs.begin(), aux_svs.end(), nullptr), aux_svs.end());	

		main_sv = nullptr; // reset main_sv to pick the next largest SV in the next iteration
	}

	// associate SNPs to the closest main SV within readlen-2*min_clip_len bp
	for (snp_t& snp : snps) {
		for (auto& sv : main_svs) {
			int dist = std::min(abs((int)snp.pos - (int)sv->start), abs((int)snp.pos - (int)sv->end));
			if (dist < stats.read_len - 2*config.min_clip_len) {
				sv->aux_snps.push_back(snp);
			}
		}
	}

	if (record_junction_remap_ref_range) {
		std::pair<hts_pos_t, hts_pos_t> highq_ref_range = get_highq_ref_range(cigar, ref_begin, junction_seq.length(), lowq_junction_prefix, lowq_junction_suffix);
		for (auto& sv : main_svs) {
			if (sv != input_main_sv && highq_ref_range.first < highq_ref_range.second) {
				set_junction_remap_ref_range(sv, highq_ref_range.first, highq_ref_range.second, stats, config);
			}
		}
	}

	return main_svs;
}

std::vector<std::shared_ptr<sv_t>> detect_svs_from_aln(StripedSmithWaterman::Alignment& aln, std::string contig_name, hts_pos_t ref_start,
	std::string junction_seq, std::shared_ptr<sv_t> main_sv, int lowq_junction_prefix, int lowq_junction_suffix, stats_t& stats, config_t& config,
	bool record_junction_remap_ref_range = false) {
	return detect_svs_from_aln(aln.cigar, ref_start + aln.ref_begin, ref_start + aln.ref_end, aln.sw_score, contig_name,
		junction_seq, main_sv, lowq_junction_prefix, lowq_junction_suffix, stats, config, record_junction_remap_ref_range);
}

std::vector<std::shared_ptr<sv_t>> detect_svs_from_aln(bam1_t* read, std::string contig_name,
	std::string junction_seq, std::shared_ptr<sv_t> main_sv, int lowq_junction_prefix, int lowq_junction_suffix, stats_t& stats, config_t& config,
	bool record_junction_remap_ref_range = false) {
	uint32_t* bam_cigar = bam_get_cigar(read);
	std::vector<uint32_t> ssw_cigar;
	for (int i = 0; i < read->core.n_cigar; i++) {
		char op = bam_cigar_opchr(bam_cigar[i]);
		int len = bam_cigar_oplen(bam_cigar[i]);
		ssw_cigar.push_back(to_cigar_int(len, op));
	}
	return detect_svs_from_aln(ssw_cigar, read->core.pos, bam_endpos(read), get_AS_tag(read),
		contig_name, junction_seq, main_sv, lowq_junction_prefix, lowq_junction_suffix, stats, config,
		record_junction_remap_ref_range);
}

std::vector<std::shared_ptr<sv_t>> detect_svs_from_junction(std::string& contig_name, char* contig_seq, std::string junction_seq, 
		hts_pos_t ref_remap_lh_start, hts_pos_t ref_remap_lh_end, hts_pos_t ref_remap_rh_start, hts_pos_t ref_remap_rh_end, 
		StripedSmithWaterman::Aligner& aligner, int lowq_junction_prefix, int lowq_junction_suffix,
		stats_t& stats, config_t& config) {
		// do not find aux SNPs in low-quality regions of the junction sequence

	const int MAX_REF_REMAP_LEN = 100000;

    hts_pos_t ref_remap_lh_len = ref_remap_lh_end - ref_remap_lh_start;
    hts_pos_t ref_remap_rh_len = ref_remap_rh_end - ref_remap_rh_start;

    if (ref_remap_lh_len > MAX_REF_REMAP_LEN) {
        throw std::runtime_error("Error: reference remapping length exceeds maximum allowed length: " + std::to_string(ref_remap_lh_len) + " > " + std::to_string(MAX_REF_REMAP_LEN));
    }
	if (ref_remap_rh_len > MAX_REF_REMAP_LEN) {
		throw std::runtime_error("Error: reference remapping length exceeds maximum allowed length: " + std::to_string(ref_remap_rh_len) + " > " + std::to_string(MAX_REF_REMAP_LEN));
	}

    char ref_lh_cstr[MAX_REF_REMAP_LEN+1];
	for (int i = 0; i < ref_remap_lh_len; i++) {
		ref_lh_cstr[i] = toupper(contig_seq[ref_remap_lh_start+i]);
	} ref_lh_cstr[ref_remap_lh_len] = '\0';

	int sep = junction_seq.find("-");
	std::string prefix_junction_seq = sep == std::string::npos ? junction_seq : junction_seq.substr(0, sep);

	SW_SCORE_INT_16* prefix_scores = smith_waterman_gotoh(ref_lh_cstr, ref_remap_lh_len, prefix_junction_seq.c_str(), prefix_junction_seq.length(), 1, -4, -6, -1);

    char ref_rh_cstr[MAX_REF_REMAP_LEN+1];
    for (int i = 0; i < ref_remap_rh_len; i++) {
        ref_rh_cstr[i] = toupper(contig_seq[ref_remap_rh_start+i]);
    } ref_rh_cstr[ref_remap_rh_len] = '\0';

    char ref_rh_cstr_rev[MAX_REF_REMAP_LEN+1];
    for (int i = 0; i < ref_remap_rh_len; i++) {
        ref_rh_cstr_rev[i] = ref_rh_cstr[ref_remap_rh_len-1-i];
    } ref_rh_cstr_rev[ref_remap_rh_len] = '\0';

	std::string suffix_junction_seq = sep == std::string::npos ? junction_seq : junction_seq.substr(sep+1);
    std::string suffix_junction_seq_rev = std::string(suffix_junction_seq.rbegin(), suffix_junction_seq.rend());
    SW_SCORE_INT_16* suffix_scores = smith_waterman_gotoh(ref_rh_cstr_rev, ref_remap_rh_len, suffix_junction_seq_rev.c_str(), suffix_junction_seq_rev.length(), 1, -4, -6, -1);
	int suffix_begin = sep == std::string::npos ? 0 : sep+1;

    SW_SCORE_INT_16 max_score = 0;
	int best_i = 0, best_j = 0;
	for (int i = config.min_clip_len; i < prefix_junction_seq.length()-config.min_clip_len; i++) {
        SW_SCORE_INT_16 prefix_score = prefix_scores[i-1]; // score of the best aln of [0..i-1]
        for (int j = std::max(i, suffix_begin); j <= junction_seq.length()-config.min_clip_len; j++) {
            SW_SCORE_INT_16 suffix_score = suffix_scores[junction_seq.length()-j-1]; // score of the best aln of [j..junction_seq.length()-1]
            // note that we want the score of the suffix of length junction_seq.length()-j, 
            // so we need to subtract 1 because suffix_scores[n] is the score of the best suffix of length n+1
            if (prefix_score + suffix_score > max_score) {
                max_score = prefix_score + suffix_score;
                best_i = i, best_j = j;
            }
        }
    }

	free(prefix_scores);
	free(suffix_scores);

    if (max_score == 0) return std::vector<std::shared_ptr<sv_t>>(); // no good alignment found

    std::string left_part = junction_seq.substr(0, best_i);
    std::string middle_part = junction_seq.substr(best_i, best_j-best_i);
    std::string right_part = junction_seq.substr(best_j);

	// if either the left or right part is fully lowq, find variants by realigning the full junction
	if (left_part.length() <= lowq_junction_prefix || right_part.length() <= lowq_junction_suffix) {
		if (overlap(ref_remap_lh_start, ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end) == 0) {
			return std::vector<std::shared_ptr<sv_t>>();
		}
		hts_pos_t remap_start = std::min(ref_remap_lh_start, ref_remap_rh_start);
		hts_pos_t remap_end = std::max(ref_remap_lh_end, ref_remap_rh_end);
		hts_pos_t remap_len = remap_end - remap_start;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment aln;
		aligner.Align(junction_seq.c_str(), contig_seq + remap_start, remap_len, filter, &aln, 0);
		std::vector<std::shared_ptr<sv_t>> remapped_svs = detect_svs_from_aln(aln, contig_name, remap_start, junction_seq,
			nullptr, lowq_junction_prefix, lowq_junction_suffix, stats, config, true);
		return remapped_svs;
	}

    std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(ref_lh_cstr, 0, ref_remap_lh_len, (char*) left_part.c_str(), aligner);
	std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(ref_rh_cstr, 0, ref_remap_rh_len, (char*) right_part.c_str(), aligner);

    StripedSmithWaterman::Alignment left_part_aln, right_part_aln;
	int min_size = INT32_MAX;
	for (StripedSmithWaterman::Alignment& _lh_aln : left_part_alns) {
		for (StripedSmithWaterman::Alignment& _rh_aln : right_part_alns) {
			int size = abs((ref_remap_rh_start + _rh_aln.ref_begin - 1) - (ref_remap_lh_start + _lh_aln.ref_end));
			if (size < min_size) {
				left_part_aln = _lh_aln, right_part_aln = _rh_aln;
				min_size = size;
			}
		}
	}

	// if either the left or right part is fully lowq, find variants by realigning the full junction
	if ((left_part_aln.query_end-left_part_aln.query_begin)/(double) left_part.length() < 0.5 
	|| (right_part_aln.query_end-right_part_aln.query_begin)/(double) right_part.length() < 0.5) {
		if (overlap(ref_remap_lh_start, ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end) == 0) {
			return std::vector<std::shared_ptr<sv_t>>();
		}
		hts_pos_t remap_start = std::min(ref_remap_lh_start, ref_remap_rh_start);
		hts_pos_t remap_end = std::max(ref_remap_lh_end, ref_remap_rh_end);
		hts_pos_t remap_len = remap_end - remap_start;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment aln;
		aligner.Align(junction_seq.c_str(), contig_seq + remap_start, remap_len, filter, &aln, 0);
		std::vector<std::shared_ptr<sv_t>> remapped_svs = detect_svs_from_aln(aln, contig_name, remap_start, junction_seq,
			nullptr, lowq_junction_prefix, lowq_junction_suffix, stats, config, true);
		return remapped_svs;
	}

    hts_pos_t left_anchor_start = ref_remap_lh_start + left_part_aln.ref_begin;
    hts_pos_t left_anchor_end = ref_remap_lh_start + left_part_aln.ref_end;
    hts_pos_t right_anchor_start = ref_remap_rh_start + right_part_aln.ref_begin;
    hts_pos_t right_anchor_end = ref_remap_rh_start + right_part_aln.ref_end;
    hts_pos_t left_bp = left_anchor_end, right_bp = right_anchor_start - 1; 
	if (right_bp < 0) right_bp = 0; // this can be -1, if the duplication starts at the beginning of the contig. 
	/* The VCF specification says: "For simple insertions and deletions in which either the REF or one of the ALT alleles 
	would otherwise be null/empty, the REF and ALT Strings must include the base before the event (which must be reflected 
	in the POS field), unless the event occurs at position 1 on the contig in which case it must include the base after 
	the event"
	However, I am not sure what the base "after the event" is in the case of a duplication. I will report 1 in the VCF (meaning 0 here) for now */
	auto left_part_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(left_anchor_start, left_anchor_end, left_part.length(), left_part_aln.sw_score);
	auto right_part_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(right_anchor_start, right_anchor_end, right_part.length(), right_part_aln.sw_score);

	int prefix_mh_len = 0;
	std::vector<std::shared_ptr<sv_t>> svs;
    if (left_bp > right_bp) { // there is microhomology in the inserted seq or it's a duplication
        int mh_len = left_bp - right_bp;
        std::pair<int, int> lp_suffix_score = find_aln_suffix_score(left_part_aln.cigar, mh_len, 1, -4, -6, -1);
        std::pair<int, int> rp_prefix_score = find_aln_prefix_score(right_part_aln.cigar, mh_len, 1, -4, -6, -1);

		if (middle_part.size() < mh_len && middle_part.size() + mh_len > 50) {
			bool one_time_perfect_repeat = lp_suffix_score.first == mh_len && rp_prefix_score.first == mh_len && middle_part.empty() &&
				!is_right_clipped(left_part_aln) && !is_left_clipped(right_part_aln);
			if (right_anchor_end - left_anchor_end < config.min_clip_len || right_anchor_start - left_anchor_start < config.min_clip_len ||
				one_time_perfect_repeat) { // it's a duplication
				std::shared_ptr<duplication_t> sv = std::make_shared<duplication_t>(contig_name, right_bp, left_bp, middle_part, nullptr, nullptr, left_part_anchor_aln, right_part_anchor_aln);
				sv->cn_unresolved = !one_time_perfect_repeat;
				svs.push_back(sv);
			}
		}

        // otherwise, it's an insertion
        std::string mh;
        if (lp_suffix_score.first > rp_prefix_score.first) { 
            // the suffix of the left part is more similar to the reference, hence we choose
			// the right-clip (if any) of the left part + the prefix of the right part 
            //to add as part of the inserted sequence
            int right_bp_adjustment = 0;
            int query_mh_bases = query_prefix_len(right_part_aln.cigar, mh_len, right_bp_adjustment);
            mh = right_part.substr(0, query_mh_bases);
            right_bp = left_bp + right_bp_adjustment;
            middle_part = left_part.substr(left_part.length() - get_right_clip_size(left_part_aln)) + middle_part + mh;
			prefix_mh_len = mh.length();
        } else {
            // the prefix of the right part is more similar to the reference, hence we choose 
			// the suffix of the left part + the left-clip (if any) of the right part
            // to add as part of the inserted sequence
            int left_bp_adjustment = 0;
            int query_mh_bases = query_suffix_len(left_part_aln.cigar, mh_len, left_bp_adjustment);
            mh = left_part.substr(left_part.length() - query_mh_bases);
            left_bp = right_bp - left_bp_adjustment;
            middle_part = mh + middle_part + right_part.substr(0, get_left_clip_size(right_part_aln));
			prefix_mh_len = mh.length();
        }
    } else {
		middle_part = left_part.substr(left_part.length() - get_right_clip_size(left_part_aln)) + middle_part +
					  right_part.substr(0, get_left_clip_size(right_part_aln));
	}

    if (svs.empty()) { // we haven't already called a duplication
		if (right_bp - left_bp > middle_part.length()) { // length of ALT < REF, deletion
			// // For small deletions with non empty middle part, we may be able to obtain a simpler representation by realigning the whole junction sequence
			if (right_bp - left_bp <= 50 && !middle_part.empty() && overlap(ref_remap_lh_start, ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end) > 0) {
				hts_pos_t remap_start = std::min(ref_remap_lh_start, ref_remap_rh_start);
				hts_pos_t remap_end = std::max(ref_remap_lh_end, ref_remap_rh_end);
				hts_pos_t remap_len = remap_end - remap_start;
				StripedSmithWaterman::Filter filter;
				StripedSmithWaterman::Alignment aln;
				aligner.Align(junction_seq.c_str(), contig_seq + remap_start, remap_len, filter, &aln, 0);
				if (!is_clipped(aln, config.min_clip_len)) {
					std::vector<std::shared_ptr<sv_t>> remapped_svs = detect_svs_from_aln(aln, contig_name, remap_start,
						junction_seq, nullptr, lowq_junction_prefix, lowq_junction_suffix, stats, config, true);
					return remapped_svs;
				}
			}

			std::shared_ptr<sv_t> sv = std::make_shared<deletion_t>(contig_name, left_bp, right_bp, middle_part, nullptr, nullptr, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		} else { // length of ALT > REF, insertion
			// If we are detecting a small insertion with a microhomology, we try to realign the whole junction sequence
			// This is because split alignments that support duplications have an unfair advantage compared to regular insertions,
			// since the inserted sequence is also aligned to the sequence. This can lead to suboptimal duplications being called instead of correct insertions
			// Furthermore, for complex small insertions (i.e., left bp < right bp), we can sometimes obtain a simpler representation this way
			if ((prefix_mh_len > 0 || left_bp < right_bp) && overlap(ref_remap_lh_start, ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end) > 0 && middle_part.length() <= 50) {
				hts_pos_t remap_start = std::min(ref_remap_lh_start, ref_remap_rh_start);
				hts_pos_t remap_end = std::max(ref_remap_lh_end, ref_remap_rh_end);
				hts_pos_t remap_len = remap_end - remap_start;
				StripedSmithWaterman::Filter filter;
				StripedSmithWaterman::Alignment aln;
				aligner.Align(junction_seq.c_str(), contig_seq + remap_start, remap_len, filter, &aln, 0);
				if (!is_clipped(aln, config.min_clip_len)) {
					std::vector<std::shared_ptr<sv_t>> remapped_svs = detect_svs_from_aln(aln, contig_name, remap_start,
						junction_seq, nullptr, lowq_junction_prefix, lowq_junction_suffix, stats, config, true);
					return remapped_svs;
				}
			}

			std::shared_ptr<sv_t> sv = std::make_shared<insertion_t>(contig_name, left_bp, right_bp, middle_part, nullptr, nullptr, left_part_anchor_aln, right_part_anchor_aln);
			svs.push_back(sv);
		}
		svs[0]->mh_len = prefix_mh_len;
	}

	std::vector<std::shared_ptr<sv_t>> extra_svs = detect_svs_from_aln(left_part_aln, contig_name, ref_remap_lh_start,
		left_part, svs[0], lowq_junction_prefix, 0, stats, config, true);
	for (const auto& sv : extra_svs) {
		if (sv != svs[0] && sv->end <= svs[0]->start) {
			svs.push_back(sv);
		}
	}
	// left_part_aln can only contain aux indels and SNPs that are to the left of the main SV
	// delete sv->aux_indels and sv->aux_snps that are >= svs[0]->start
	svs[0]->aux_indels.erase(std::remove_if(svs[0]->aux_indels.begin(), svs[0]->aux_indels.end(),
		[svs](std::shared_ptr<sv_t> sv) {
			return sv->end > svs[0]->start;
		}), svs[0]->aux_indels.end());
	svs[0]->aux_snps.erase(std::remove_if(svs[0]->aux_snps.begin(), svs[0]->aux_snps.end(),
		[svs](snp_t snp) {
			return snp.pos > svs[0]->start;
		}), svs[0]->aux_snps.end());

	// we need to apply the same (but reversed) filtering to the right part
	// i.e., right_part_aln can only contain aux indels and SNPs that are to the right of the main SV
	// so we temporarily remove left part aux indels and SNPs
	std::vector<std::shared_ptr<sv_t>> backup_aux_indels;
	backup_aux_indels.swap(svs[0]->aux_indels);
	std::vector<snp_t> backup_aux_snps;
	backup_aux_snps.swap(svs[0]->aux_snps);

	extra_svs = detect_svs_from_aln(right_part_aln, contig_name, ref_remap_rh_start, right_part, svs[0],
		0, lowq_junction_suffix, stats, config, true);
	for (const auto& sv : extra_svs) {
		if (sv != svs[0] && sv->start >= svs[0]->end) {
			svs.push_back(sv);
		}
	}
	svs[0]->aux_indels.erase(std::remove_if(svs[0]->aux_indels.begin(), svs[0]->aux_indels.end(),
		[svs](std::shared_ptr<sv_t> sv) {
			return sv->start <= svs[0]->end;
		}), svs[0]->aux_indels.end());
	svs[0]->aux_snps.erase(std::remove_if(svs[0]->aux_snps.begin(), svs[0]->aux_snps.end(),
		[svs](snp_t snp) {
			return snp.pos <= svs[0]->end;
		}), svs[0]->aux_snps.end());
	
	// restore left part aux indels and SNPs
	svs[0]->aux_indels.insert(svs[0]->aux_indels.end(), backup_aux_indels.begin(), backup_aux_indels.end());
	svs[0]->aux_snps.insert(svs[0]->aux_snps.end(), backup_aux_snps.begin(), backup_aux_snps.end());

	std::sort(svs[0]->aux_indels.begin(), svs[0]->aux_indels.end(), [](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
		return a->start < b->start;
	});
	std::sort(svs[0]->aux_snps.begin(), svs[0]->aux_snps.end(), [](const snp_t& a, const snp_t& b) {
		return a.pos < b.pos;
	});

	std::pair<hts_pos_t, hts_pos_t> left_highq_ref_range = get_highq_ref_range(left_part_aln.cigar,
		left_anchor_start, left_part.length(), lowq_junction_prefix, 0);
	std::pair<hts_pos_t, hts_pos_t> right_highq_ref_range = get_highq_ref_range(right_part_aln.cigar,
		right_anchor_start, right_part.length(), 0, lowq_junction_suffix);
	if (left_highq_ref_range.first < left_highq_ref_range.second &&
			right_highq_ref_range.first < right_highq_ref_range.second) {
		set_junction_remap_ref_range(svs[0], std::min(left_highq_ref_range.first, right_highq_ref_range.first),
			std::max(left_highq_ref_range.second, right_highq_ref_range.second), stats, config);
	}

    return svs;
}

std::vector<std::shared_ptr<sv_t>> detect_svs(std::string& contig_name, char* contig_seq, hts_pos_t contig_len, 
					std::shared_ptr<consensus_t> rc_consensus, std::shared_ptr<consensus_t> lc_consensus,
					StripedSmithWaterman::Aligner& aligner, int min_overlap,
					stats_t& stats, config_t& config) {

	std::string consensus_junction_seq;
	hts_pos_t ref_remap_lh_start = 0, ref_remap_lh_end = 0, ref_remap_rh_start = 0, ref_remap_rh_end = 0;
	int overlap = 0;
	double mismatch_rate = 0.0;
	int lowq_junction_prefix = 0, lowq_junction_suffix = 0;

	if (rc_consensus != NULL && lc_consensus != NULL) {
		std::string rc_consensus_seq = rc_consensus->sequence;
		std::string lc_consensus_seq = lc_consensus->sequence;
		suffix_prefix_aln_t spa = aln_suffix_prefix_perfect(rc_consensus_seq, lc_consensus_seq, min_overlap);
		if (spa.overlap < min_overlap || is_homopolymer(lc_consensus->sequence.c_str(), spa.overlap)) {
			rc_consensus_seq = rc_consensus->sequence.substr(0, rc_consensus->sequence.length()-rc_consensus->lowq_suffix);
			lc_consensus_seq = lc_consensus->sequence.substr(lc_consensus->lowq_prefix);
			spa = aln_suffix_prefix_perfect(rc_consensus_seq, lc_consensus_seq, min_overlap);
			if (spa.overlap < min_overlap || is_homopolymer(lc_consensus_seq.c_str(), spa.overlap)) {
				return std::vector<std::shared_ptr<sv_t>>();
			}
		}

		overlap = spa.overlap;
		mismatch_rate = spa.mismatch_rate();

		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment aln_rh, aln_lh;
		
		if (rc_consensus->reads() >= lc_consensus->reads()) {
			consensus_junction_seq = rc_consensus_seq + lc_consensus_seq.substr(spa.overlap);
		} else {
			consensus_junction_seq = rc_consensus_seq.substr(0, rc_consensus_seq.length()-spa.overlap) + lc_consensus_seq;
		}

		ref_remap_lh_start = rc_consensus->breakpoint - consensus_junction_seq.length();
		ref_remap_lh_end = rc_consensus->breakpoint + consensus_junction_seq.length();
		if (ref_remap_lh_start < 0) ref_remap_lh_start = 0;
		if (ref_remap_lh_end > contig_len) ref_remap_lh_end = contig_len;

		ref_remap_rh_start = lc_consensus->breakpoint - consensus_junction_seq.length();
		ref_remap_rh_end = lc_consensus->breakpoint + consensus_junction_seq.length();
		if (ref_remap_rh_start < 0) ref_remap_rh_start = 0;
		if (ref_remap_rh_end > contig_len) ref_remap_rh_end = contig_len;

		lowq_junction_prefix = rc_consensus->lowq_prefix;
		lowq_junction_suffix = lc_consensus->lowq_suffix;
	} else if (rc_consensus != NULL) {
		consensus_junction_seq = rc_consensus->sequence;

		ref_remap_lh_start = std::max(hts_pos_t(0), rc_consensus->start - (int) rc_consensus->sequence.length());
		ref_remap_lh_end = std::min(rc_consensus->end + (int) rc_consensus->sequence.length(), contig_len);

		ref_remap_rh_start = rc_consensus->other_bp_lower_boundary;
		ref_remap_rh_end = rc_consensus->other_bp_upper_boundary + rc_consensus->sequence.length();
		if (rc_consensus->other_bp_upper_boundary == consensus_t::UPPER_BOUNDARY_NON_CALCULATED) {
			ref_remap_rh_start = rc_consensus->breakpoint - 2*rc_consensus->sequence.length();
			ref_remap_rh_end = rc_consensus->breakpoint + 2*rc_consensus->sequence.length();
		}
		ref_remap_rh_start = std::max(ref_remap_rh_start, hts_pos_t(0));
		ref_remap_rh_end = std::min(ref_remap_rh_end, contig_len);

		lowq_junction_prefix = rc_consensus->lowq_prefix;
		lowq_junction_suffix = rc_consensus->lowq_suffix;
	} else if (lc_consensus != NULL) {
		consensus_junction_seq = lc_consensus->sequence;

		ref_remap_lh_start = lc_consensus->other_bp_lower_boundary - lc_consensus->sequence.length();
		ref_remap_lh_end = lc_consensus->other_bp_upper_boundary;
		if (lc_consensus->other_bp_lower_boundary == consensus_t::LOWER_BOUNDARY_NON_CALCULATED) { // could not calculate the remap boundary, fall back to formula
			ref_remap_lh_start = lc_consensus->breakpoint - 2*lc_consensus->sequence.length();
			ref_remap_lh_end = lc_consensus->breakpoint + 2*lc_consensus->sequence.length();
		}
		ref_remap_lh_start = std::max(ref_remap_lh_start, hts_pos_t(0));
		ref_remap_lh_end = std::min(ref_remap_lh_end, contig_len);

		ref_remap_rh_start = std::max(hts_pos_t(0), lc_consensus->start - (int) lc_consensus->sequence.length());
		ref_remap_rh_end = std::min(lc_consensus->end + (int) lc_consensus->sequence.length(), contig_len);

		lowq_junction_prefix = lc_consensus->lowq_prefix;
		lowq_junction_suffix = lc_consensus->lowq_suffix;
	} else {
		return std::vector<std::shared_ptr<sv_t>>();
	}

	std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contig_seq, consensus_junction_seq, ref_remap_lh_start, 
		ref_remap_lh_end, ref_remap_rh_start, ref_remap_rh_end, aligner, lowq_junction_prefix, lowq_junction_suffix,
		stats, config);
	for (const auto& sv : svs) {
		sv->rc_consensus = rc_consensus;
		sv->lc_consensus = lc_consensus;

		if (rc_consensus != NULL && lc_consensus == NULL) {
			if (rc_consensus->is_hsr) {
				sv->source = "1HSR_RC";
			} else {
				sv->source = "1SR_RC";
			}
		} else if (rc_consensus == NULL && lc_consensus != NULL) {
			if (lc_consensus->is_hsr) {
				sv->source = "1HSR_LC";
			} else {
				sv->source = "1SR_LC";
			}
		} else if (!rc_consensus->is_hsr && !lc_consensus->is_hsr) {
			sv->source = "2SR";
		} else if (rc_consensus->is_hsr && lc_consensus->is_hsr) {
			sv->source = "2HSR";
		} else if (!rc_consensus->is_hsr && lc_consensus->is_hsr) {
			sv->source = "SR-HSR";
		} else if (rc_consensus->is_hsr && !lc_consensus->is_hsr) {
			sv->source = "HSR-SR";
		}
	}
    return svs;
}

std::shared_ptr<breakend_t> detect_bnd(std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::shared_ptr<consensus_t> leftmost_consensus, std::shared_ptr<consensus_t> rightmost_consensus, 
	suffix_prefix_aln_t& spa, StripedSmithWaterman::Aligner& aligner, int min_clip_len) {

	std::string lm_seq = leftmost_consensus->sequence, rm_seq = rightmost_consensus->sequence;
	if (leftmost_consensus->left_clipped) rc(lm_seq);
	else rc(rm_seq);
	std::string full_junction_seq = lm_seq + rm_seq.substr(spa.overlap);

	hts_pos_t ref_remap_lh_start = leftmost_consensus->breakpoint - full_junction_seq.length();
	if (ref_remap_lh_start < 0) ref_remap_lh_start = 0;
	hts_pos_t ref_remap_lh_end = leftmost_consensus->breakpoint + full_junction_seq.length();
	if (ref_remap_lh_end > contig_len) ref_remap_lh_end = contig_len;

	hts_pos_t ref_remap_rh_start = rightmost_consensus->breakpoint - full_junction_seq.length();
	if (ref_remap_rh_start < 0) ref_remap_rh_start = 0;
	hts_pos_t ref_remap_rh_end = rightmost_consensus->breakpoint + full_junction_seq.length();
	if (ref_remap_rh_end > contig_len) ref_remap_rh_end = contig_len;

	if (!leftmost_consensus->left_clipped) {
		SW_SCORE_INT_16* fwd_prefix_scores = smith_waterman_gotoh(contig_seq+ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, full_junction_seq.c_str(), full_junction_seq.length(), 1, -4, -6, -1);
		rc(full_junction_seq);
		SW_SCORE_INT_16* revc_prefix_scores = smith_waterman_gotoh(contig_seq+ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, full_junction_seq.c_str(), full_junction_seq.length(), 1, -4, -6, -1);

		int max_score = 0, best_i = 0, best_j = 0;
		for (int i = min_clip_len; i < full_junction_seq.length()-min_clip_len; i++) {
			int fwd_prefix_score = fwd_prefix_scores[i-1]; // score of the best aln of full_junction_seq[0..i-1]
			for (int j = i; j <= full_junction_seq.length()-min_clip_len; j++) {
				int rev_prefix_score = revc_prefix_scores[full_junction_seq.length()-j-1]; // score of the best aln of RC of full_junction_seq[j..junction_seq.length()-1]
				if (fwd_prefix_score + rev_prefix_score >= max_score) {
					max_score = fwd_prefix_score + rev_prefix_score;
					best_i = i, best_j = j;
				}
			}
		}

		free(fwd_prefix_scores);
		free(revc_prefix_scores);

		if (max_score == 0) return NULL;

		rc(full_junction_seq);
		std::string left_part = full_junction_seq.substr(0, best_i);
		std::string middle_part = full_junction_seq.substr(best_i, best_j-best_i);
		std::string right_part = full_junction_seq.substr(best_j);
		rc(right_part);

		std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(contig_seq, ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, (char*) left_part.c_str(), aligner);
		std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(contig_seq, ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, (char*) right_part.c_str(), aligner);
		StripedSmithWaterman::Alignment left_part_aln = left_part_alns[left_part_alns.size()-1], right_part_aln = right_part_alns[0];

		auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(ref_remap_lh_start+left_part_aln.ref_begin, ref_remap_lh_start+left_part_aln.ref_end, left_part.length(), left_part_aln.sw_score);
		auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(ref_remap_rh_start+right_part_aln.ref_begin, ref_remap_rh_start+right_part_aln.ref_end, right_part.length(), right_part_aln.sw_score);

		hts_pos_t start = ref_remap_lh_start + left_part_aln.ref_end, end = ref_remap_rh_start + right_part_aln.ref_end;
		return std::make_shared<breakend_t>(contig_name, start, end, middle_part, left_anchor_aln, right_anchor_aln, false, leftmost_consensus, rightmost_consensus);
	} else {

		rc(full_junction_seq);
		std::string full_junction_seq_rev = std::string(full_junction_seq.rbegin(), full_junction_seq.rend());
		char* ref_remap_lh_rev = new char[ref_remap_lh_end-ref_remap_lh_start+1];
		for (int i = 0; i < ref_remap_lh_end-ref_remap_lh_start; i++) {
			ref_remap_lh_rev[i] = std::toupper(contig_seq[ref_remap_lh_end-1-i]);
		} ref_remap_lh_rev[ref_remap_lh_end-ref_remap_lh_start] = '\0';
		SW_SCORE_INT_16* revc_suffix_scores = smith_waterman_gotoh(ref_remap_lh_rev, ref_remap_lh_end-ref_remap_lh_start, full_junction_seq_rev.c_str(), full_junction_seq_rev.length(), 1, -4, -6, -1);
		
		rc(full_junction_seq);
		full_junction_seq_rev = std::string(full_junction_seq.rbegin(), full_junction_seq.rend());
		char* ref_remap_rh_rev = new char[ref_remap_rh_end-ref_remap_rh_start+1];
		for (int i = 0; i < ref_remap_rh_end-ref_remap_rh_start; i++) {
			ref_remap_rh_rev[i] = std::toupper(contig_seq[ref_remap_rh_end-1-i]);
		} ref_remap_rh_rev[ref_remap_rh_end-ref_remap_rh_start] = '\0';
		SW_SCORE_INT_16* fwd_suffix_scores = smith_waterman_gotoh(ref_remap_rh_rev, ref_remap_rh_end-ref_remap_rh_start, full_junction_seq_rev.c_str(), full_junction_seq_rev.length(), 1, -4, -6, -1);

		int max_score = 0, best_i = 0, best_j = 0;
		for (int i = min_clip_len; i < full_junction_seq.length()-min_clip_len; i++) {
			int revc_suffix_score = revc_suffix_scores[i-1]; // score of the best aln of the RC of full_junction_seq[i..n]
			for (int j = i; j <= full_junction_seq.length()-min_clip_len; j++) {
				int fwd_suffix_score = fwd_suffix_scores[full_junction_seq.length()-j-1]; // score of the best aln of [j..junction_seq.length()-1]
				if (fwd_suffix_score + revc_suffix_score > max_score) {
					max_score = fwd_suffix_score + revc_suffix_score;
					best_i = i, best_j = j;
				}
			}
		}

		free(fwd_suffix_scores);
		free(revc_suffix_scores);
		delete[] ref_remap_lh_rev;
		delete[] ref_remap_rh_rev;

		if (max_score == 0) return NULL;

		std::string left_part = full_junction_seq.substr(0, best_i);
		rc(left_part);
		std::string middle_part = full_junction_seq.substr(best_i, best_j-best_i);
		std::string right_part = full_junction_seq.substr(best_j);

		std::vector<StripedSmithWaterman::Alignment> left_part_alns = get_best_alns(contig_seq, ref_remap_lh_start, ref_remap_lh_end-ref_remap_lh_start, (char*) left_part.c_str(), aligner);
		std::vector<StripedSmithWaterman::Alignment> right_part_alns = get_best_alns(contig_seq, ref_remap_rh_start, ref_remap_rh_end-ref_remap_rh_start, (char*) right_part.c_str(), aligner);
		StripedSmithWaterman::Alignment left_part_aln = left_part_alns[left_part_alns.size()-1], right_part_aln = right_part_alns[0];

		auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(ref_remap_lh_start+left_part_aln.ref_begin, ref_remap_lh_start+left_part_aln.ref_end, left_part.length(), left_part_aln.sw_score);
		auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(ref_remap_rh_start+right_part_aln.ref_begin, ref_remap_rh_start+right_part_aln.ref_end, right_part.length(), right_part_aln.sw_score);

		hts_pos_t start = ref_remap_lh_start + left_part_aln.ref_begin-1, end = ref_remap_rh_start + right_part_aln.ref_begin-1;
		if (start < 0) start = 0;
		return std::make_shared<breakend_t>(contig_name, start, end, middle_part, left_anchor_aln, right_anchor_aln, true, leftmost_consensus, rightmost_consensus);
	}
}

bool accept(StripedSmithWaterman::Alignment& aln, int min_clip_len, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);

	bool lc_is_noise = false, rc_is_noise = false;
	if (!qual_ascii.empty()) {
		if (lc_size > min_clip_len && lc_size < qual_ascii.length()/2) {
			std::string lc_qual = qual_ascii.substr(0, lc_size);
			std::string not_lc_qual = qual_ascii.substr(lc_size);
			lc_is_noise = avg_qual(lc_qual) < min_avg_base_qual && avg_qual(not_lc_qual) >= min_avg_base_qual;
		}
		if (rc_size > min_clip_len && rc_size < qual_ascii.length()/2) {
			std::string not_rc_qual = qual_ascii.substr(0, qual_ascii.length()-rc_size);
			std::string rc_qual = qual_ascii.substr(qual_ascii.length()-rc_size);
			rc_is_noise = avg_qual(rc_qual) < min_avg_base_qual && avg_qual(not_rc_qual) >= min_avg_base_qual;
		}
	}

	double mismatch_rate = double(aln.mismatches)/(aln.query_end-aln.query_begin);
	return (lc_size <= min_clip_len || lc_is_noise) && (rc_size <= min_clip_len || rc_is_noise) && mismatch_rate <= max_seq_error;
}

std::pair<int, int> get_dels_ins_in_first_n_chars(StripedSmithWaterman::Alignment& aln, int n) {
    int dels = 0, inss = 0;
    int offset = 0;
    for (uint32_t c : aln.cigar) {
        int len = cigar_int_to_len(c);
        char op = cigar_int_to_op(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = op == '=' || op == 'X' || op == 'D' || op == 'S';
        if (consumes_ref && offset + len > n) {
            len = n-offset;
        }

        if (op == 'D') {
            dels += len;
        } else if (op == 'I') {
            inss += len;
        }

        if (consumes_ref) {
            offset += len;
            if (offset == n) break;
        }
    }
    return {dels, inss};
}

hts_pos_t get_start_offset(StripedSmithWaterman::Alignment& r1, StripedSmithWaterman::Alignment& r2) {
	int r1_uc_start = r1.ref_begin - get_left_clip_size(r1);
	int r2_uc_start = r2.ref_begin - get_left_clip_size(r2);
    hts_pos_t offset = r2_uc_start - r1_uc_start;
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(r1, offset);
    return offset + del_ins.second - del_ins.first;
}

#endif // SW_UTILS_H
