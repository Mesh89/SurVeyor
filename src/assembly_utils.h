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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "consensus.h"
#include "ungapped_alignment.h"
#include "utils.h"

bool reject_assembly_overlap(const std::string& left, const std::string& right, int overlap) {
	char hp_base = is_homopolymer(right.data(), overlap);
	if (hp_base == 'C' || hp_base == 'G') return true;
	hp_base = is_homopolymer(left.data()+left.length()-overlap, overlap);
	if (hp_base == 'C' || hp_base == 'G') return true;

	size_t first_non_acgt = right.find_first_not_of("ACGT");
	for (int displacement = 1; displacement <= 6 && displacement < overlap; displacement++) {
		int shorter = overlap-displacement;
		const char* suffix = left.data()+left.length()-shorter;
		if (shorter <= first_non_acgt && number_of_mismatches_fast(suffix, right.data(), shorter, 0) == 0) return true;
	}
	return false;
}

struct path_permission_t {
	bool can_start_path, can_end_path;
};
struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

void add_assembly_edge(int from, const edge_t& edge, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	out_edges[from]++;
	l_adj[from].push_back(edge);
	l_adj_rev[edge.next].push_back({from, edge.score, edge.overlap});
}

struct assembly_path_t {
	std::vector<int> vertices;
	std::vector<int> overlaps; // overlaps[i] joins vertices[i] to vertices[i+1].
};

assembly_path_t select_assembly_path(const std::vector<std::vector<edge_t> >& l_adj_rev, const std::vector<int>& rev_topological_order,
		const std::vector<path_permission_t>& path_permissions, std::vector<bool>& used) {
	int n = l_adj_rev.size();
	// compute longest paths
	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		if (used[i]) continue;
		if (best_scores[i] == 0 && !path_permissions[i].can_end_path) continue;
		for (const edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}

	// Used vertices skip propagation above, but may still be selected as starts.
	int best_score = 0, curr_vertex = 0;
	for (int i = 0; i < best_scores.size(); i++) {
		if (path_permissions[i].can_start_path && best_score < best_scores[i]) {
			best_score = best_scores[i];
			curr_vertex = i;
		}
	}
	assembly_path_t path;
	if (best_score == 0) return path;

	path.vertices.push_back(curr_vertex);
	while (best_edges[curr_vertex].overlap) {
		used[curr_vertex] = true;
		path.overlaps.push_back(best_edges[curr_vertex].overlap);
		curr_vertex = best_edges[curr_vertex].next;
		path.vertices.push_back(curr_vertex);
	}
	used[curr_vertex] = true;
	return path;
}

struct assembly_sequence_t {
	std::string seq;
	std::vector<int> offsets; // One start offset per vertex in path order.
};

struct assembly_consensus_t {
	std::string seq;
	std::vector<base_votes_t> votes;
};

assembly_sequence_t reconstruct_assembly_path(const assembly_path_t& path, const std::vector<const assembly_consensus_t*>& contigs) {
	assembly_sequence_t assembly;
	if (path.vertices.empty()) return assembly;

	assembly.seq = contigs[path.vertices[0]]->seq;
	assembly.offsets.push_back(0);
	for (int i = 0; i < path.overlaps.size(); i++) {
		int overlap = path.overlaps[i];
		assembly.offsets.push_back(assembly.seq.length()-overlap);
		assembly.seq += contigs[path.vertices[i+1]]->seq.substr(overlap);
	}
	return assembly;
}

void add_read_votes(assembly_consensus_t& contig, const uint8_t* quals) {
	for (int i = 0; i < contig.seq.length(); i++) {
		contig.votes[i].add(contig.seq[i], quals ? quals[i] : 0);
	}
}

assembly_consensus_t make_read_contig(const std::string& seq, const uint8_t* quals) {
	assembly_consensus_t contig{seq, std::vector<base_votes_t>(seq.length())};
	add_read_votes(contig, quals);
	return contig;
}

assembly_consensus_t polish_assembly_path(const std::vector<const assembly_consensus_t*>& inputs, assembly_sequence_t assembly) {
	assembly_consensus_t consensus{std::move(assembly.seq), {}};
	consensus.votes.resize(consensus.seq.length());
	for (int v = 0; v < inputs.size(); v++) {
		const std::vector<base_votes_t>& input_votes = inputs[v]->votes;
		for (int i = 0; i < input_votes.size(); i++) {
			consensus.votes[assembly.offsets[v]+i].merge(input_votes[i]);
		}
	}
	for (int i = 0; i < consensus.seq.length(); i++) {
		consensus.seq[i] = consensus.votes[i].best_base(consensus.seq[i]);
	}
	return consensus;
}

std::vector<assembly_consensus_t> assemble_graph_paths(const std::vector<const assembly_consensus_t*>& contigs, const std::vector<std::vector<edge_t> >& l_adj_rev,
		const std::vector<int>& rev_topological_order, const std::vector<path_permission_t>& path_permissions, std::vector<bool>& used) {
	std::vector<assembly_consensus_t> assembled_contigs;
	while (true) {
		assembly_path_t path = select_assembly_path(l_adj_rev, rev_topological_order, path_permissions, used);
		if (path.vertices.empty()) break;

		assembly_sequence_t assembly = reconstruct_assembly_path(path, contigs);
		std::vector<const assembly_consensus_t*> path_inputs;
		for (int vertex : path.vertices) path_inputs.push_back(contigs[vertex]);
		assembled_contigs.push_back(polish_assembly_path(path_inputs, std::move(assembly)));
	}
	return assembled_contigs;
}

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

// Remove every edge from -> e.next for which remove(from, e) is true.
void remove_assembly_edges(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev,
		const std::function<bool(int, const edge_t&)>& remove) {
	for (int i = 0; i < l_adj.size(); i++) {
		l_adj[i].erase(std::remove_if(l_adj[i].begin(), l_adj[i].end(), [&](const edge_t& e) { return remove(i, e); }), l_adj[i].end());
		out_edges[i] = l_adj[i].size();
		l_adj_rev[i].erase(std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(), [&](const edge_t& e) { return remove(e.next, {i, e.score, e.overlap}); }), l_adj_rev[i].end());
	}
}

// Count overlap mismatches where the mean quality of both bases is at least min_qual.
int count_hq_mismatches(const assembly_consensus_t& left, const assembly_consensus_t& right, int overlap, int min_qual) {
	int hq_mismatches = 0, left_start = left.seq.length()-overlap;
	for (int i = 0; i < overlap; i++) {
		char left_base = left.seq[left_start+i], right_base = right.seq[i];
		int left_index = base_to_index(left_base), right_index = base_to_index(right_base);
		if (left_base == right_base || left_index < 0 || right_index < 0) continue;
		const base_votes_t& left_votes = left.votes[left_start+i];
		const base_votes_t& right_votes = right.votes[i];
		int left_qual = left_votes.counts[left_index] ? left_votes.quals[left_index]/left_votes.counts[left_index] : 0;
		int right_qual = right_votes.counts[right_index] ? right_votes.quals[right_index]/right_votes.counts[right_index] : 0;
		if (left_qual >= min_qual && right_qual >= min_qual) hq_mismatches++;
	}
	return hq_mismatches;
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

void correct_contig_with_alignments(std::string& contig, std::vector<std::string>& reads, const std::vector<ungapped_aln_t>& ungapped_alns, const std::vector<const uint8_t*>& quals,
	const std::function<bool(int, const ungapped_aln_t&)>& accept_read, const std::vector<bool>& read_is_reverse, const std::string& preferred_bases, hts_pair_pos_t* supported_bounds = nullptr) {

	// Half-open bounds of positive-quality canonical votes; beg >= end means no support.
	if (supported_bounds != nullptr) *supported_bounds = {HTS_POS_MAX, 0};
	std::vector<base_votes_t> votes(contig.length());
	std::vector<strand_base_scores_t> strand_scores(read_is_reverse.empty() ? 0 : contig.length());
	for (int i = 0; i < reads.size(); i++) {
		const ungapped_aln_t& ungapped_aln = ungapped_alns[i];
		if (ungapped_aln.query_end <= ungapped_aln.query_begin) continue;

		if (accept_read(i, ungapped_aln)) {
			for (int j = ungapped_aln.query_begin; j < ungapped_aln.query_end; j++) {
				int ref_pos = j - ungapped_aln.query_begin + ungapped_aln.ref_begin;

				if (ref_pos >= votes.size()) {
					votes.resize(ref_pos+1);
					if (!read_is_reverse.empty()) strand_scores.resize(ref_pos+1);
				}

				char c = reads[i][j];
				int qual = quals.empty() ? 0 : quals[i][j];
				if (!read_is_reverse.empty()) strand_scores[ref_pos].add(c, qual, read_is_reverse[i]);
				votes[ref_pos].add(c, qual);
				if (supported_bounds != nullptr && qual > 0 && base_to_index(c) >= 0) {
					supported_bounds->beg = std::min<hts_pos_t>(supported_bounds->beg, ref_pos);
					supported_bounds->end = std::max<hts_pos_t>(supported_bounds->end, ref_pos+1);
				}
			}
		}
	}

	contig.resize(votes.size());

	for (int i = 0; i < contig.length(); i++) {
		if (!read_is_reverse.empty()) {
			contig[i] = strand_scores[i].best_base(contig[i], i < preferred_bases.length() ? preferred_bases[i] : '\0');
			continue;
		}
		contig[i] = votes[i].best_base(contig[i]);
	}
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, int min_clip_len, const std::vector<const uint8_t*>& quals,
	const std::function<bool(int, const ungapped_aln_t&)>& accept_read, const std::vector<bool>& read_is_reverse, const std::string& preferred_bases, hts_pair_pos_t* supported_bounds = nullptr) {

	std::vector<ungapped_aln_t> ungapped_alns;
	for (std::string& read : reads) {
		ungapped_alns.push_back(best_ungapped_aln(read.c_str(), read.length(), contig.c_str(), contig.length(), min_clip_len - 1, 1, 0));
	}
	correct_contig_with_alignments(contig, reads, ungapped_alns, quals, accept_read, read_is_reverse, preferred_bases, supported_bounds);
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
