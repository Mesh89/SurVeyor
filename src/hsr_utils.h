#ifndef HSR_UTILS_H
#define HSR_UTILS_H

#include <algorithm>
#include <unordered_set>
#include <htslib/sam.h>

#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"
#include "utils.h"
#include "sw_utils.h"
#include "SegTree.h"

int compute_left_half_Ms(bam1_t* r) {
	int border = r->core.l_qseq/2;
	int left_Ms = 0;

	// compute how many Ms (either matches or mismatches are in the first half of the read)
	// in other words, this is readlen/2 - number of insertions in the first half
	int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op_char = bam_cigar_opchr(cigar[i]);
		int op_len = bam_cigar_oplen(cigar[i]);
		if (op_char == 'M' || op_char == '=' || op_char == 'X') {
			left_Ms += std::min(border, qpos+op_len) - qpos;
			qpos += op_len;
		} else if (op_char == 'I' || op_char == 'S') {
			qpos += op_len;
		}

		if (qpos >= border) break;
	}
	return left_Ms;
}

// computes the differences for the left and the right half of the read
std::pair<int, int> compute_left_and_right_differences_indel_as_n_diffs(bam1_t* r) {
    int border = r->core.l_qseq/2;
    int left_Ms = compute_left_half_Ms(r);

    // we computed left_Ms because the MD tag does not include insertions
    // i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
    // but it is in position 125
    // here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
    std::string md_tag = get_md(r);
    int md_len = md_tag.length();
    int m = 0;
    int left_diffs = 0, right_diffs = 0;
    bool del_mode = false;
    for (int i = 0; i < md_len; i++) {
        char c = md_tag[i];
        if (c >= '0' && c <= '9') {
            m = m*10 + c-'0';
            del_mode = false;
        } else if (c >= 'A' && c <= 'T') {
            left_Ms -= m;
            m = 0;
            if (left_Ms > 0) left_diffs++;
            else right_diffs++;
            if (!del_mode) left_Ms--;
        } else if (c == '^') {
			left_Ms -= m;
			m = 0;
            del_mode = true;
        }
    }

    // We also need to count the number of inserted bases, since the MD tag does not report this information
    int qpos = 0;
	uint32_t* cigar = bam_get_cigar(r);
    for (int i = 0; i < r->core.n_cigar; i++) {
        char op_char = bam_cigar_opchr(cigar[i]);
        int op_len = bam_cigar_oplen(cigar[i]);
        if (op_char == 'M' || op_char == '=' || op_char == 'X') {
            qpos += op_len;
        } else if (op_char == 'I') {
            if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
                                                            // half and partially in the right half
                left_diffs += border - qpos;
                right_diffs += qpos + op_len - border;
            } else if (qpos < border) {
                left_diffs += op_len;
            } else if (qpos >= border) {
                right_diffs += op_len;
            }
            qpos += op_len;
        }
    }

    return {left_diffs, right_diffs};
}

enum class consensus_ref_classification_t {
	ALIGNS_WELL,
	LEFT,
	RIGHT,
	AMBIGUOUS
};

// The query is expected to be the high-quality portion of the consensus sequence.
consensus_ref_classification_t classify_consensus_against_ref(
		char* contig_seq, hts_pos_t contig_len, hts_pos_t consensus_start, hts_pos_t consensus_end,
		const std::string& query, const config_t& config) {

	StripedSmithWaterman::Aligner aligner(2, 2, 4, 1, true);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	hts_pos_t ref_start = std::max<hts_pos_t>(0, consensus_start - 10);
	hts_pos_t ref_end = std::min(contig_len, consensus_end + 10);
	aligner.Align(query.c_str(), contig_seq + ref_start, ref_end - ref_start, filter, &aln, 0);

	int dops = 0, dbases = 0, iops = 0, ibases = 0, xbases = 0;
	int left_diffs = 0, right_diffs = 0;
	int qpos = 0, border = query.length()/2;
	auto add_query_diffs = [&](int len) {
		int left_len = std::max(0, std::min(qpos + len, border) - qpos);
		left_diffs += left_len;
		right_diffs += len - left_len;
		qpos += len;
	};

	for (uint32_t op : aln.cigar) {
		char opchar = cigar_int_to_op(op);
		int oplen = cigar_int_to_len(op);
		if (opchar == 'D') {
			dops++;
			dbases += oplen;
			(qpos < border ? left_diffs : right_diffs) += oplen;
		} else if (opchar == 'I') {
			iops++;
			ibases += oplen;
			add_query_diffs(oplen);
		} else if (opchar == 'X') {
			xbases += oplen;
			add_query_diffs(oplen);
		} else if (opchar == '=' || opchar == 'M' || opchar == 'S') {
			qpos += oplen;
		}
	}

	int left_clip_len = get_left_clip_size(aln), right_clip_len = get_right_clip_size(aln);
	bool supports_hsr = dops + iops + xbases >= config.min_diff_hsr || abs(ibases - dbases) >= config.min_sv_size;

	if (left_clip_len == 0 && right_clip_len == 0 && !supports_hsr) {
		return consensus_ref_classification_t::ALIGNS_WELL;
	}
	
	// Infer the direction of the clustr
	// A strong clips (>=min_clip_len) takes the precedence to determine the direction of the cluster
	// If both sides are strong clips, the longer wins
	// If a direction cannot be decided with strong clips alone, we include the number of differences of the left and right halves of the read
	if (left_clip_len > right_clip_len && left_clip_len >= config.min_clip_len) {
		return consensus_ref_classification_t::LEFT;
	} else if (right_clip_len > left_clip_len && right_clip_len >= config.min_clip_len) {
		return consensus_ref_classification_t::RIGHT;
	} else if (left_clip_len + left_diffs > right_clip_len + right_diffs) {
		return consensus_ref_classification_t::LEFT;
	} else if (right_clip_len + right_diffs > left_clip_len + left_diffs) {
		return consensus_ref_classification_t::RIGHT;
	}  else {
		return consensus_ref_classification_t::AMBIGUOUS;
	}
}

// Remove consensues that are completely contained within another consensues and their sequence is a substring of the other consensus
void filter_fully_contained(std::vector<consensus_t*>& consensuses) {
	std::vector<consensus_t*> sorted = consensuses;
	std::sort(sorted.begin(), sorted.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->start < c2->start || (c1->start == c2->start && c1->end > c2->end);
	});

	std::vector<bool> to_delete(sorted.size(), false);
	for (int i = 0; i < sorted.size(); i++) {
		if (to_delete[i]) continue;
		for (int j = i+1; j < sorted.size(); j++) {
			if (to_delete[j]) continue;
			if (sorted[j]->start > sorted[i]->end) break;
			if (sorted[j]->end <= sorted[i]->end) {
				// j is contained within i
				std::string highq_seq = sorted[j]->highq_sequence();
				if (sorted[i]->sequence.find(highq_seq) != std::string::npos) {
					to_delete[j] = true;
				}
			}
		}
	}
	consensuses.clear();
	for (int i = 0; i < sorted.size(); i++) {
		if (!to_delete[i]) consensuses.push_back(sorted[i]);
		else delete sorted[i];
	}
}

// c1 is assumed to be to the left of c2, and neither cluster is not fully contained within the other
void merge_overlapping_pair_of_clusters(consensus_t* c1, consensus_t* c2, consensus_t* target, std::string& c1_seq, std::string& c2_seq, int overlap) {
	if (c1->is_hsr && !c2->is_hsr) {
		target->breakpoint = c2->breakpoint;
	} else if (!c1->is_hsr && c2->is_hsr) {
		target->breakpoint = c1->breakpoint;
	} else {
		target->breakpoint = c1->left_clipped ? c1->breakpoint : c2->breakpoint;
	}
	target->start = c1->start;
	target->end = c2->end;
	target->sequence = c1_seq + c2_seq.substr(overlap);
	target->fwd_reads = c1->fwd_reads + c2->fwd_reads;
	target->rev_reads = c1->rev_reads + c2->rev_reads;
	target->clip_len = 0;
	target->max_mapq = std::max(c1->max_mapq, c2->max_mapq);
	target->other_bp_lower_boundary = std::max(c1->other_bp_lower_boundary, c2->other_bp_lower_boundary);
	target->other_bp_upper_boundary = std::min(c1->other_bp_upper_boundary, c2->other_bp_upper_boundary);
	target->lowq_prefix = c1->lowq_prefix;
	target->lowq_suffix = c2->lowq_suffix;
	target->is_hsr = c1->is_hsr && c2->is_hsr; 
}

bool merge_overlapping_pair_of_clusters(consensus_t* c1, consensus_t* c2, consensus_t* target, int min_overlap) {
	// only merge if the overlap is at least half of the length of the shorter consensus
	hts_pos_t overlap = c1->end - c2->start;
	if (overlap < (int) std::min(c1->sequence.length(), c2->sequence.length())/2) return false;
	
	hts_pos_t c1_hq_end = c1->end - c1->lowq_suffix;
	hts_pos_t c2_hq_start = c2->start + c2->lowq_prefix;
	hts_pos_t hq_overlap = c1_hq_end - c2_hq_start;
	std::string c1_hq_seq = c1->sequence.substr(0, c1->sequence.length()-c1->lowq_suffix);
	std::string c2_hq_seq = c2->sequence.substr(c2->lowq_prefix);
	int min_hq_seq_len = std::min(c1_hq_seq.length(), c2_hq_seq.length());
	if (hq_overlap >= min_hq_seq_len/2 && hq_overlap < min_hq_seq_len) { // first, try to see if the high quality parts are compatible
		suffix_prefix_aln_t spa_hq = aln_suffix_prefix_perfect(c1_hq_seq, c2_hq_seq, min_overlap);
		if (spa_hq.overlap) {
			merge_overlapping_pair_of_clusters(c1, c2, target, c1_hq_seq, c2_hq_seq, spa_hq.overlap);
			return true;
		}
	}
	
	suffix_prefix_aln_t spa = aln_suffix_prefix_perfect(c1->sequence, c2->sequence, min_overlap); // then, try the whole sequences
	if (spa.overlap) {
		merge_overlapping_pair_of_clusters(c1, c2, target, c1->sequence, c2->sequence, spa.overlap);
		return true;
	}
	return false;
}

void merge_overlapping_clusters(std::vector<consensus_t*>& consensuses, int min_overlap) {

	if (consensuses.empty()) return;

	std::vector<consensus_t*> sorted_by_start = consensuses;
	std::sort(sorted_by_start.begin(), sorted_by_start.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->start < c2->start;
	});

	std::vector<consensus_t*> sorted_by_end = consensuses;
	std::sort(sorted_by_end.begin(), sorted_by_end.end(), [](const consensus_t* c1, const consensus_t* c2) {
		return c1->end < c2->end;
	});

	std::unordered_map<consensus_t*, int> start_pos, end_pos;
	for (int i = 0; i < sorted_by_start.size(); i++) {
		start_pos[sorted_by_start[i]] = i;
		end_pos[sorted_by_end[i]] = i;
	}

	std::vector<consensus_t*> sorted_by_reads = consensuses;
	std::sort(sorted_by_reads.begin(), sorted_by_reads.end(), [](consensus_t* c1, consensus_t* c2) {
		return c1->reads() > c2->reads();
	});

	std::unordered_set<consensus_t*> removed;
	for (int i = 0; i < sorted_by_reads.size(); i++) {
		consensus_t* c = sorted_by_reads[i];
		if (removed.count(c) > 0) continue;

		int e_pos = end_pos[c];
		for (int j = e_pos-1; j >= 0; j--) {
			consensus_t* o = sorted_by_end[j];
			if (removed.count(o) > 0) continue;
			if (o->end < c->start) break;
			if (o->start >= c->start || o->end == c->end) continue; // o is fully contained within c

			if (merge_overlapping_pair_of_clusters(o, c, c, min_overlap)) {
				removed.insert(o);
				std::swap(sorted_by_start[start_pos[o]], sorted_by_start[start_pos[c]]);
				std::swap(start_pos[o], start_pos[c]);
				// end positions do not change
			}
		}

		int s_pos = start_pos[c];
		for (int j = s_pos+1; j < sorted_by_start.size(); j++) {
			consensus_t* o = sorted_by_start[j];
			if (removed.count(o) > 0) continue;
			if (o->start > c->end) break;
			if (o->start == c->start || o->end <= c->end) continue; // o is fully contained within c
			
			if (merge_overlapping_pair_of_clusters(c, o, c, min_overlap)) {
				removed.insert(o);
				std::swap(sorted_by_end[end_pos[o]], sorted_by_end[end_pos[c]]);
				std::swap(end_pos[o], end_pos[c]);
				// start positions do not change
			}
		}
	}

	std::vector<consensus_t*> merged;
	for (consensus_t* c : consensuses) {
		if (removed.count(c) == 0) merged.push_back(c);
		else delete c;
	}
	consensuses.swap(merged);
}

void enforce_max_ploidy(std::vector<consensus_t*>& consensuses, int max_ploidy) {

	int contig_len = 0;
	for (consensus_t* c : consensuses) {
		if (c->end > contig_len) contig_len = c->end;
	}

	SegTree segtree(contig_len+1);
	std::vector<consensus_t*> retained;
	std::sort(consensuses.begin(), consensuses.end(), [](consensus_t* c1, consensus_t* c2) {
		return c1->reads() > c2->reads();
	});
	for (consensus_t* c : consensuses) {
		if (!segtree.any_ge(c->start, c->end, max_ploidy)) {
			segtree.add(c->start, c->end, 1);
			retained.push_back(c);
		} else {
			delete c;
		}
	}
	consensuses.swap(retained);
}

#endif
