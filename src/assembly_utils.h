#ifndef ASSEMBLY_UTILS_H_
#define ASSEMBLY_UTILS_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <stack>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "htslib/hts.h"
#include "consensus.h"
#include "ungapped_alignment.h"

struct path_permission_t {
	bool can_start_path, can_end_path;
};
struct seq_w_pp_t {
	std::string seq;
	std::vector<uint8_t> qual;
	path_permission_t clip_pair;

	seq_w_pp_t() : seq(), clip_pair() {}
	seq_w_pp_t(std::string& seq, bool can_start_path, bool can_end_path, std::vector<uint8_t> qual = {}) : seq(seq), qual(std::move(qual)) {
		clip_pair.can_start_path = can_start_path;
		clip_pair.can_end_path = can_end_path;
	}
};

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

std::vector<int> find_rev_topological_order(int n, std::vector<int> out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

bool break_cycle_at_vertex(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj,
		std::vector<std::vector<edge_t> >& l_adj_rev, int i) {

	struct cycle_edge_t {
		int from;
		edge_t edge;
	};

	int n = l_adj.size();
	if (i < 0 || i >= n) return false;

	std::stack<int> s;
	std::vector<bool> visited(n, false);
	std::vector<int> parent(n, -1);
	std::vector<edge_t> parent_edge(n);
	std::vector<cycle_edge_t> cycle;

	parent[i] = i;
	s.push(i);

	while (cycle.empty() && !s.empty()) {
		int curr = s.top();
		s.pop();
		if (visited[curr]) continue;
		visited[curr] = true;

		for (edge_t& e : l_adj[curr]) {
			if (e.next == i) {
				cycle.push_back({curr, e});
				int node = curr;
				while (node != i) {
					cycle.push_back({parent[node], parent_edge[node]});
					node = parent[node];
				}
				break;
			}
			if (!visited[e.next] && parent[e.next] == -1) {
				parent[e.next] = curr;
				parent_edge[e.next] = e;
				s.push(e.next);
			}
		}
	}

	if (cycle.empty()) return false;

	cycle_edge_t weakest = cycle[0];
	for (cycle_edge_t& ce : cycle) {
		if (ce.edge.overlap < weakest.edge.overlap ||
				(ce.edge.overlap == weakest.edge.overlap && ce.edge.score < weakest.edge.score)) {
			weakest = ce;
		}
	}

	std::vector<edge_t>& fwd_edges = l_adj[weakest.from];
	auto fwd_it = std::find_if(fwd_edges.begin(), fwd_edges.end(), [&weakest](edge_t& e) {
		return e.next == weakest.edge.next && e.score == weakest.edge.score && e.overlap == weakest.edge.overlap;
	});
	if (fwd_it != fwd_edges.end()) {
		fwd_edges.erase(fwd_it);
		out_edges[weakest.from]--;
	}

	std::vector<edge_t>& rev_edges = l_adj_rev[weakest.edge.next];
	auto rev_it = std::find_if(rev_edges.begin(), rev_edges.end(), [&weakest](edge_t& e) {
		return e.next == weakest.from && e.score == weakest.edge.score && e.overlap == weakest.edge.overlap;
	});
	if (rev_it != rev_edges.end()) {
		rev_edges.erase(rev_it);
	}

	return true;
}

void break_cycles(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {

	int n = l_adj.size();
	for (int i = 0; i < n; i++) {
		while (break_cycle_at_vertex(out_edges, l_adj, l_adj_rev, i));
	}
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, double max_acceptable_error_rate, int min_clip_len, const std::vector<const uint8_t*>& quals = {}, 
	const std::function<bool(int, const ungapped_aln_t&)>& accept_read = {}, int mismatch_score = -4, const std::vector<bool>& read_is_reverse = {}, const std::string& preferred_bases = {}) {

	std::vector<ungapped_aln_t> ungapped_alns;
	for (std::string& read : reads) {
		ungapped_alns.push_back(best_ungapped_aln(read.c_str(), read.length(), contig.c_str(), contig.length(), min_clip_len - 1, 1, mismatch_score));
	}
	std::vector<int> As(contig.length()), Cs(contig.length()), Gs(contig.length()), Ts(contig.length());
	std::vector<int> A_quals(contig.length()), C_quals(contig.length()), G_quals(contig.length()), T_quals(contig.length());
	std::vector<strand_base_scores_t> strand_scores(read_is_reverse.empty() ? 0 : contig.length());
	for (int i = 0; i < reads.size(); i++) {
		ungapped_aln_t& ungapped_aln = ungapped_alns[i];
		if (ungapped_aln.query_end <= ungapped_aln.query_begin) continue;

		double mismatch_rate = double(ungapped_aln.mismatches)/(ungapped_aln.query_end-ungapped_aln.query_begin);
		if (accept_read ? accept_read(i, ungapped_aln) : mismatch_rate <= max_acceptable_error_rate) {
			for (int j = ungapped_aln.query_begin; j < ungapped_aln.query_end; j++) {
				int ref_pos = j - ungapped_aln.query_begin + ungapped_aln.ref_begin;

				if (ref_pos >= As.size()) {
					// extend As, Cs, Gs, Ts
					int old_size = As.size();
					As.resize(ref_pos+1);
					Cs.resize(ref_pos+1);
					Gs.resize(ref_pos+1);
					Ts.resize(ref_pos+1);
					A_quals.resize(ref_pos+1);
					C_quals.resize(ref_pos+1);
					G_quals.resize(ref_pos+1);
					T_quals.resize(ref_pos+1);
					if (!read_is_reverse.empty()) strand_scores.resize(ref_pos+1);
					for (int k = old_size; k < As.size(); k++) {
						As[k] = Cs[k] = Gs[k] = Ts[k] = A_quals[k] = C_quals[k] = G_quals[k] = T_quals[k] = 0;
					}
				}

				char c = reads[i][j];
				int qual = quals.empty() ? 0 : quals[i][j];
				if (!read_is_reverse.empty()) strand_scores[ref_pos].add(c, qual, read_is_reverse[i]);
				if (c == 'A') { As[ref_pos]++; A_quals[ref_pos] += qual; }
				else if (c == 'C') { Cs[ref_pos]++; C_quals[ref_pos] += qual; }
				else if (c == 'G') { Gs[ref_pos]++; G_quals[ref_pos] += qual; }
				else if (c == 'T') { Ts[ref_pos]++; T_quals[ref_pos] += qual; }
			}
		}
	}

	contig.resize(As.size());

	for (int i = 0; i < contig.length(); i++) {
		if (!read_is_reverse.empty()) {
			contig[i] = strand_scores[i].best_base(contig[i], i < preferred_bases.length() ? preferred_bases[i] : '\0');
			continue;
		}
		if (As[i] + Cs[i] + Gs[i] + Ts[i] == 0) continue;

		auto best_base_score = std::max({
			std::make_tuple(A_quals[i], As[i], 'A'), std::make_tuple(C_quals[i], Cs[i], 'C'),
			std::make_tuple(G_quals[i], Gs[i], 'G'), std::make_tuple(T_quals[i], Ts[i], 'T')
		});
		contig[i] = std::get<2>(best_base_score);
	}
}

bool is_contained(const std::string& s1, const std::string& s2, double max_seq_error) {
	const std::string& shorter = s1.length() <= s2.length() ? s1 : s2;
	const std::string& longer = s1.length() <= s2.length() ? s2 : s1;
	int cap = (int) (shorter.length()*max_seq_error);
	for (int offset = 0; offset <= longer.length()-shorter.length(); offset++) {
		if (number_of_mismatches_fast(shorter.data(), longer.data()+offset, shorter.length(), cap) <= cap) return true;
	}
	return false;
}

// Check an ungapped placement, rejecting mismatches only when both bases have Phred Q >= 40.
bool consensus_placement_compatible(const std::string& seq1, const std::string& qual1, size_t start1, const std::string& seq2, const std::string& qual2, size_t length) {
	for (size_t i = 0; i < length; i++) {
		if (seq1[start1+i] != seq2[i] && qual1[start1+i] >= 40+33 && qual2[i] >= 40+33) return false;
	}
	return true;
}

void make_offsets_nonnegative(std::vector<hts_pos_t>& read_start_offsets) {
    hts_pos_t min_offset = 0;
    for (hts_pos_t offset : read_start_offsets) {
        if (offset < min_offset) min_offset = offset;
    }
    if (min_offset < 0) {
        for (hts_pos_t& offset : read_start_offsets) {
            offset -= min_offset;
        }
    }
}

#endif /* ASSEMBLY_UTILS_H_ */
