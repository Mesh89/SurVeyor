#ifndef EXTEND_1SR_CONSENSUS_H_
#define EXTEND_1SR_CONSENSUS_H_

#include <set>
#include <unordered_map>
#include <queue>
#include <stack>

#include "../libs/IntervalTree.h"
#include "htslib/sam.h"
#include "utils.h"
#include "sw_utils.h"
#include "sam_utils.h"
#include "types.h"

struct ext_read_t {
	std::string qname;
	hts_pos_t start, end;
	uint8_t mapq;
	bool rev, same_chr, left_clipped;
	int sequence_len = 0;
	char* sequence = NULL;

	ext_read_t(bam1_t* read) {
		init(read);
	}

	void init(bam1_t* read) {
		qname = bam_get_qname(read);
		start = read->core.pos;
		end = bam_endpos(read);
		mapq = read->core.qual;
		rev = bam_is_rev(read);
		same_chr = is_samechr(read);
		left_clipped = is_left_clipped(read, 1);
		if (sequence_len < read->core.l_qseq) {
			if (sequence != NULL) delete[] sequence;
			sequence = new char[read->core.l_qseq+1];
		}
		copy_sequence(read, sequence);
		sequence_len = read->core.l_qseq;
	}

	~ext_read_t() {
		if (sequence != NULL) delete[] sequence;
	}
};

struct ext_read_allocator_t {
    std::queue<ext_read_t*> reads;

    ext_read_allocator_t() {}

    ext_read_t* get(bam1_t* read) {
        ext_read_t* r = nullptr;
		if (!reads.empty()) {
			r = reads.front();
			reads.pop();
		}
        if (!r) return new ext_read_t(read);
        r->init(read);
        return r;
    }

    void release(ext_read_t* read) {
        reads.push(read);
    }

    ~ext_read_allocator_t() {
        while (!reads.empty()) {
            auto* r = reads.front();
            reads.pop();
            delete r;
        }
    }

    ext_read_allocator_t(const ext_read_allocator_t&) = delete;
    ext_read_allocator_t& operator=(const ext_read_allocator_t&) = delete;
};

bool is_vertex_in_cycle(std::vector<std::vector<edge_t> >& l_adj, int i) {
	std::stack<int> s;
	for (edge_t& e : l_adj[i]) {
		s.push(e.next);
	}

	// check if node i can reach itself - i.e., it is part of a cycle
	std::vector<bool> visited(l_adj.size(), false);
	bool i_in_cycle = false;
	while (!s.empty()) {
		int curr = s.top();
		s.pop();
		if (visited[curr]) continue;
		visited[curr] = true;

		if (curr == i) return true;
		for (edge_t& e : l_adj[curr]) {
			s.push(e.next);
		}
	}
	return false;
}

std::vector<bool> find_vertices_in_cycles(std::vector<std::vector<edge_t> >& l_adj) {

	int n = l_adj.size();
	std::vector<bool> in_cycle(n);
	for (int i = 0; i < n; i++) {
		in_cycle[i] = is_vertex_in_cycle(l_adj, i);
	}
	return in_cycle;
}

void build_graph_fwd(std::vector<std::string>& read_seqs, std::vector<hts_pos_t>& read_starts,
		std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap,
		bool strict) {

	l_adj = std::vector<std::vector<edge_t> >(read_seqs.size());
	l_adj_rev = std::vector<std::vector<edge_t> >(read_seqs.size());
	out_edges = std::vector<int>(read_seqs.size(), 0);

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	// allocate memory for homopolymer prefix array
	int max_read_len = 0;
	for (std::string& s : read_seqs) {
		if (s.length() > max_read_len) max_read_len = s.length();
	}
	bool* hp_prefix = new bool[max_read_len];

	std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> kmer_to_idx;
	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		// let this string be s2, and the current string we are trying to extend be s1
		// since we are extending to the right, we are interested in checking whether s1 suffix matches s2 prefix
		// for this reason, we skip kmers corresponding to prefixes of s2 that are shorter than min_overlap
		std::string& seq = read_seqs[i];
		is_homopolymer_prefix(seq.c_str(), seq.length(), hp_prefix);
		for (int j = std::max(0, min_overlap-32); j < seq.length()-1; j++) {
			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= min_overlap && !hp_prefix[j]) {
				kmer_to_idx[kmer].emplace_back(i, j);
			}
		}
	}
	delete[] hp_prefix;

	// sort vectors in kmer_to_idx
	for (auto& kv : kmer_to_idx) {
		std::vector<std::pair<int,int>>& idxs = kv.second;
		// sort by first ascending and second descending
		std::sort(idxs.begin(), idxs.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
			if (a.first != b.first) return a.first < b.first;
			return a.second > b.second;
		});
	}

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int curr_node = bfs.front();
		bfs.pop();

		if (visited[curr_node]) continue;
		visited[curr_node] = true;

		std::string& s1 = read_seqs[curr_node];
		uint64_t kmer = 0;
		for (int j = s1.length()-32; j < s1.length(); j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		int last_accepted_j = -1;
		for (int i = 0; i < kmer_to_idx[kmer].size(); i++) {
			int j = kmer_to_idx[kmer][i].first;
			if (j == last_accepted_j) continue;

			int p = kmer_to_idx[kmer][i].second; // position of kmer in read j
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;
			if (strict && read_starts[curr_node] >= read_starts[j] && read_starts[curr_node] != 0 && read_starts[j] != 0) continue;

			int cmp_len = p+1;
			if (cmp_len < min_overlap || cmp_len > s1.length() || strncmp(s1.c_str()+(s1.length()-cmp_len), s2.c_str(), cmp_len) != 0) continue;

			out_edges[curr_node]++;
			l_adj[curr_node].push_back({j, cmp_len, cmp_len});
			l_adj_rev[j].push_back({curr_node, cmp_len, cmp_len});
			bfs.push(j);
			last_accepted_j = j;
		}
	}
}

void build_graph_rev(std::vector<std::string>& read_seqs, std::vector<hts_pos_t>& read_starts,
		std::vector<int> starting_idxs, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, int min_overlap, 
		bool strict) {

	l_adj = std::vector<std::vector<edge_t> >(read_seqs.size());
	l_adj_rev = std::vector<std::vector<edge_t> >(read_seqs.size());
	out_edges = std::vector<int>(read_seqs.size(), 0);

	uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

	int n = read_seqs.size();
	std::vector<bool> visited(n);

	// allocate memory for homopolymer prefix array
	int max_read_len = 0;
	for (std::string& s : read_seqs) {
		if (s.length() > max_read_len) max_read_len = s.length();
	}
	bool* hp_prefix = new bool[max_read_len];

	std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> kmer_to_idx;
	for (int i = 0; i < n; i++) {
		uint64_t kmer = 0;

		std::string& seq = read_seqs[i];
		is_homopolymer_suffix(seq.c_str(), seq.length(), hp_prefix);
		for (int j = 0; j < seq.length(); j++) {
			if (seq.length()-j+31 < min_overlap) {
				// let this string be s2, and the current string we are trying to extend be s1
				// since we are extending to the left, we are interested in checking whether s2 suffix matches s1 prefix
				// however, kmers beyond this point will provide overlaps < min_overlap, so we skip them
				break;
			}

			uint64_t nv = nucl_bm[seq[j]];
			kmer = ((kmer << 2) | nv);

			if (j >= 31 && !hp_prefix[j-31]) {
				kmer_to_idx[kmer].emplace_back(i, j);
			}
		}
	}
	delete[] hp_prefix;
	// Note: vectors are already sorted by first ascending and second ascending due to the way we built them

	std::queue<int> bfs;
	for (int i : starting_idxs) bfs.push(i);
	while (!bfs.empty()) {
		int curr_node = bfs.front();
		bfs.pop();

		if (visited[curr_node]) continue;
		visited[curr_node] = true;

		std::string& s1 = read_seqs[curr_node];
		if (s1.length() < 32) continue;

		uint64_t kmer = 0;
		for (int j = 0; j < 32; j++) {
			uint64_t nv = nucl_bm[s1[j]];
			kmer = ((kmer << 2) | nv);
		}

		int last_accepted_j = -1;
		for (int i = 0; i < kmer_to_idx[kmer].size(); i++) {
			int j = kmer_to_idx[kmer][i].first;
			if (j == last_accepted_j) continue;

			int p = kmer_to_idx[kmer][i].second; // position of kmer in read j
			std::string& s2 = read_seqs[j];

			if (s1 == s2) continue;
			if (strict && read_starts[curr_node] <= read_starts[j] && read_starts[curr_node] != 0 && read_starts[j] != 0) continue;

			int cmp_len = s2.length()-p+31;
			if (cmp_len < min_overlap || cmp_len > s1.length() || strncmp(s2.c_str()+(s2.length()-cmp_len), s1.c_str(), cmp_len) != 0) {
				continue;
			}

			out_edges[curr_node]++;
			l_adj[curr_node].push_back({j, cmp_len, cmp_len});
			l_adj_rev[j].push_back({curr_node, cmp_len, cmp_len});
			bfs.push(j);
			last_accepted_j = j;
		}
	}
}

void get_extension_read_seqs(IntervalTree<ext_read_t*>& candidate_reads_itree, std::vector<std::string>& read_seqs,
		std::vector<int>& read_mapqs, std::vector<hts_pos_t>& read_starts, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len, int high_confidence_mapq, stats_t& stats, int max_reads = INT32_MAX) {

	hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_start-stats.max_is+stats.read_len);
	hts_pos_t fwd_mates_end = std::max(hts_pos_t(1), target_end-stats.min_is);

	hts_pos_t rev_mates_start = std::min(target_start+stats.min_is, contig_len);
	hts_pos_t rev_mates_end = std::min(target_end+stats.max_is-stats.read_len, contig_len);
	
	std::vector<Interval<ext_read_t*> > target_reads = candidate_reads_itree.findOverlapping(
		std::min(fwd_mates_start, rev_mates_start)-10, std::max(fwd_mates_end, rev_mates_end)+10);
	for (Interval<ext_read_t*> i_read : target_reads) {
		ext_read_t* ext_read = i_read.value;
		if (!ext_read->rev && ext_read->mapq >= high_confidence_mapq &&
			!ext_read->same_chr && ext_read->start >= fwd_mates_start && ext_read->start <= fwd_mates_end) {
			if (mateseqs_w_mapq.count(ext_read->qname) == 0) {
				std::cerr << "ERROR: mateseqs_w_mapq does not contain " << ext_read->qname << std::endl;
				exit(1);
			};
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			rc(mateseq_w_mapq.first);
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
			read_starts.push_back(0); // dummy value
		} else if (ext_read->rev && ext_read->mapq >= high_confidence_mapq &&
				!ext_read->same_chr && ext_read->end >= rev_mates_start && ext_read->end <= rev_mates_end) {
			if (mateseqs_w_mapq.count(ext_read->qname) == 0) {
				std::cerr << "ERROR: mateseqs_w_mapq does not contain " << ext_read->qname << std::endl;
				exit(1);
			};
			std::pair<std::string, int> mateseq_w_mapq = mateseqs_w_mapq[ext_read->qname];
			read_seqs.push_back(mateseq_w_mapq.first);
			read_mapqs.push_back(mateseq_w_mapq.second);
			read_starts.push_back(0); // dummy value
		}

		if (ext_read->end > target_start && ext_read->start < target_end) {
			read_seqs.push_back(ext_read->sequence);
			read_mapqs.push_back(ext_read->mapq);
			read_starts.push_back(ext_read->start);
		}

		if (read_seqs.size() > max_reads) return;
	}
}
std::vector<ext_read_t*> get_extension_reads(std::string contig_name, std::vector<hts_pair_pos_t>& target_ivals, hts_pos_t contig_len,
		config_t& config, stats_t& stats, open_samFile_t* bam_file) {

	std::sort(target_ivals.begin(), target_ivals.end(), [](hts_pair_pos_t& a, hts_pair_pos_t& b) {return a.beg < b.beg;});

	std::vector<hts_pair_pos_t> merged_target_ivals;
	merged_target_ivals.push_back(target_ivals[0]);
	for (int i = 1; i < target_ivals.size(); i++) {
		if (target_ivals[i].beg-merged_target_ivals.back().end < stats.read_len) {
			merged_target_ivals.back().end = std::max(merged_target_ivals.back().end, target_ivals[i].end);
		} else {
			merged_target_ivals.push_back(target_ivals[i]);
		}
	}

	ext_read_allocator_t ext_read_allocator;
	
	std::vector<ext_read_t*> reads;
	bam1_t* read = bam_init1();
	for (hts_pair_pos_t target_ival : merged_target_ivals) {
		hts_pos_t fwd_mates_start = std::max(hts_pos_t(1), target_ival.beg-stats.max_is+stats.read_len);
		hts_pos_t rev_mates_end = std::min(target_ival.end+stats.max_is-stats.read_len, contig_len);

		std::stringstream ss;
		ss << contig_name << ":" << fwd_mates_start << "-" << rev_mates_end;

		std::vector<ext_read_t*> region_reads;
		region_reads.reserve(target_ival.end-target_ival.beg+1);
		ext_read_t* curr_read = nullptr;

		// priority queue to hold left-clipped reads that start at the same position, sorted by mapping quality
		struct ext_read_cmp_mapq {
			bool operator()(ext_read_t* a, ext_read_t* b) {
				return a->mapq > b->mapq;
			}
		};
		std::priority_queue<ext_read_t*, std::vector<ext_read_t*>, ext_read_cmp_mapq> lc_reads_pq;


		hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, ss.str().c_str());
		while (sam_itr_next(bam_file->file, iter, read) >= 0) {
			if (is_unmapped(read) || !is_primary(read)) continue;

			// for same strand pairs that may be due to inversions, add a copy of the unstable end in the the correct orientation
			// same strand is not always due to inversions (transpositions may also cause this), and we do not know which case we are dealing with
			// therefore, when there is no clear stable end, and the read may either be the unstable end of an inversion or the stable end of a transposition,
			// we add a copy of the read in both orientations, and let the extension algorithm decide which one to use
			if (is_samechr(read) && is_samestr(read) && read->core.qual <= get_mq(read)) {
				if ((!bam_is_rev(read) && read->core.mpos < read->core.pos) ||
					(bam_is_mrev(read) && read->core.mpos > read->core.pos)) {
					ext_read_t* ext_read_rc = ext_read_allocator.get(read);
					rc(ext_read_rc->sequence);
					ext_read_rc->rev = !bam_is_rev(read);
					region_reads.push_back(ext_read_rc);
				}
			}

			// we are downsampling by keeping only the highest mapping quality read at each start position 
			// (plus max_depth left-clipped reads because they have a legitimate reason for starting at the same position)
			if (curr_read == nullptr) {
				curr_read = ext_read_allocator.get(read);
			} else if (read->core.pos > curr_read->start) { // new start position, push curr_read and left-clipped reads
				if (!curr_read->left_clipped) {
					region_reads.push_back(curr_read);
				} else {
					ext_read_allocator.release(curr_read);
				}
				while (!lc_reads_pq.empty()) {
					region_reads.push_back(lc_reads_pq.top());
					lc_reads_pq.pop();
				}
				curr_read = ext_read_allocator.get(read);
			} else if (curr_read->mapq < read->core.qual || curr_read->left_clipped && !is_left_clipped(read, 1)) {
				// replace curr_read with the new read if the latter has higher mapping quality or curr_read is left-clipped and the new read is not
				ext_read_allocator.release(curr_read);
				curr_read = ext_read_allocator.get(read);
			} 
			
			if (is_left_clipped(read, 1)) {
				lc_reads_pq.push(ext_read_allocator.get(read));
				if (lc_reads_pq.size() > stats.get_max_depth(contig_name)) { 
					// keep only max_depth left-clipped reads with the highest mapping quality
					ext_read_allocator.release(lc_reads_pq.top());
					lc_reads_pq.pop();
				}
			}
		}
		if (curr_read != nullptr) {
			if (!curr_read->left_clipped) {
				region_reads.push_back(curr_read);
			} else {
				ext_read_allocator.release(curr_read);
			}
		}

		reads.insert(reads.end(), region_reads.begin(), region_reads.end());
		hts_itr_destroy(iter);
	}

	bam_destroy1(read);

	return reads;
}

std::vector<ext_read_t*> get_extension_reads_from_consensuses(std::vector<std::shared_ptr<consensus_t>>& consensuses, std::string contig_name, hts_pos_t contig_len,
		config_t& config, stats_t& stats, open_samFile_t* bam_file) {

	if (consensuses.empty()) return std::vector<ext_read_t*>();

	std::vector<hts_pair_pos_t> target_ivals;
	for (std::shared_ptr<consensus_t> consensus : consensuses) {
		hts_pos_t left_ext_target_start = consensus->left_ext_target_start(stats.max_is, stats.read_len);
		hts_pos_t left_ext_target_end = consensus->left_ext_target_end(stats.max_is, stats.read_len);
		hts_pos_t right_ext_target_start = consensus->right_ext_target_start(stats.max_is, stats.read_len);
		hts_pos_t right_ext_target_end = consensus->right_ext_target_end(stats.max_is, stats.read_len);
		hts_pair_pos_t target_ival;
		target_ival.beg = left_ext_target_start, target_ival.end = left_ext_target_end;
		target_ivals.push_back(target_ival);
		target_ival.beg = right_ext_target_start, target_ival.end = right_ext_target_end;
		target_ivals.push_back(target_ival);
	}
	return get_extension_reads(contig_name, target_ivals, contig_len, config, stats, bam_file);
}

void break_cycles(std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	int n = l_adj.size();
	std::vector<bool> in_cycle = find_vertices_in_cycles(l_adj);
	for (int i = 0; i < n; i++) {
		if (in_cycle[i]) {
			l_adj[i].clear();
			out_edges[i] = 0;
		}
		l_adj_rev[i].erase(
			std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(), [&in_cycle](edge_t& e) {return in_cycle[e.next];}),
			l_adj_rev[i].end()
		);
	}
}

void extend_consensus_to_right(std::shared_ptr<consensus_t> consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len,
		const int high_confidence_mapq, stats_t& stats, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	if (consensus->extended_to_right) return;

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	std::vector<hts_pos_t> read_starts;
	read_seqs.push_back(consensus->sequence);
	read_mapqs.push_back(high_confidence_mapq);
	read_starts.push_back(1); // we want our consensus to be to the left of every other read that extends it to the right (but not 0, since that is reserved for mates, i.e., unknown pos)

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, read_starts, mateseqs_w_mapq, 
		target_start, target_end, contig_len, high_confidence_mapq, stats, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj, l_adj_rev;
	build_graph_fwd(read_seqs, read_starts, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2, false);

	bool cycle_contains_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_contains_0) {
		// remove all edges pointing to 0
		for (int i = 0; i < n; i++) {
			std::vector<edge_t>& l_adj_i = l_adj[i];
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
			out_edges[i] = l_adj_i.size();
		}
		l_adj_rev[0].clear();
	}

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	bool cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) { // try to regenerate the graph enforcing a reference-based order - it might fail due to mates all having pos 0
		build_graph_fwd(read_seqs, read_starts, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2, true);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) { // final fallback - break cycles brutally
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}


	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->sequence;
	while (e.overlap) {
		ext_consensus += read_seqs[e.next].substr(e.overlap);
		e = best_edges[e.next];
		consensus->right_ext_reads++;
		if (read_mapqs[e.next] >= high_confidence_mapq) consensus->hq_right_ext_reads++;
	}

	if (!consensus->left_clipped) {
		if (consensus->clip_len != 0) {
			consensus->clip_len += ext_consensus.length() - consensus->sequence.length();
		}
	} else {
		consensus->end += ext_consensus.length() - consensus->sequence.length();
	}
	consensus->sequence = ext_consensus;
	consensus->extended_to_right = true;
}

void extend_consensus_to_left(std::shared_ptr<consensus_t> consensus, IntervalTree<ext_read_t*>& candidate_reads_itree,
		hts_pos_t target_start, hts_pos_t target_end, hts_pos_t contig_len,
		const int high_confidence_mapq, stats_t& stats, std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq) {

	if (consensus->extended_to_left) return;

	std::vector<std::string> read_seqs;
	std::vector<int> read_mapqs;
	std::vector<hts_pos_t> read_starts;
	read_seqs.push_back(consensus->sequence);
	read_mapqs.push_back(high_confidence_mapq);
	read_starts.push_back(INT32_MAX); // we want our consensus to be to the right of every other read that extends it to the left

	get_extension_read_seqs(candidate_reads_itree, read_seqs, read_mapqs, read_starts, mateseqs_w_mapq, target_start, target_end, contig_len, high_confidence_mapq, stats, 5000);

	if (read_seqs.size() > 5000) return;

	int n = read_seqs.size();
	std::vector<int> out_edges;
	std::vector<std::vector<edge_t> > l_adj, l_adj_rev;
	build_graph_rev(read_seqs, read_starts, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2, false);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	bool cycle_at_0 = is_vertex_in_cycle(l_adj, 0);
	if (cycle_at_0) {
		// remove all edges pointing to 0
		for (int i = 0; i < n; i++) {
			std::vector<edge_t>& l_adj_i = l_adj[i];
			l_adj_i.erase(std::remove_if(l_adj_i.begin(), l_adj_i.end(), [](edge_t& e) {return e.next == 0;}), l_adj_i.end());
			out_edges[i] = l_adj_i.size();
		}
		l_adj_rev[0].clear();
	}

	cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) { // try to regenerate the graph enforcing a reference-based order - it might fail due to mates all having pos 0
		build_graph_rev(read_seqs, read_starts, std::vector<int>(1, 0), out_edges, l_adj, l_adj_rev, stats.read_len/2, true);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	cycle_at_0 = std::find(rev_topological_order.begin(), rev_topological_order.end(), 0) == rev_topological_order.end();
	if (cycle_at_0) { // final fallback - break cycles brutally
		break_cycles(out_edges, l_adj, l_adj_rev);
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
	}

	std::vector<int> best_scores(n);
	std::vector<edge_t> best_edges(n);
	for (int i : rev_topological_order) {
		for (edge_t& e : l_adj_rev[i]) {
			if (best_scores[e.next] < e.score + best_scores[i]) {
				best_scores[e.next] = e.score + best_scores[i];
				best_edges[e.next] = {i, e.score, e.overlap};
			}
		}
	}

	// 0 is the consensus, we start from there
	edge_t e = best_edges[0];
	std::string ext_consensus = consensus->sequence;
	while (e.overlap) {
		ext_consensus = read_seqs[e.next].substr(0, read_seqs[e.next].length()-e.overlap) + ext_consensus;
		e = best_edges[e.next];
		consensus->left_ext_reads++;
		if (read_mapqs[e.next] >= high_confidence_mapq) consensus->hq_left_ext_reads++;
	}

	if (consensus->left_clipped) {
		if (consensus->clip_len != 0) {
			consensus->clip_len += ext_consensus.length() - consensus->sequence.length();
		}
	} else {
		consensus->start -= ext_consensus.length() - consensus->sequence.length();
	}
	consensus->sequence = ext_consensus;
	consensus->extended_to_left = true;
}


#endif /* EXTEND_1SR_CONSENSUS_H_ */
