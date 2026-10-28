#ifndef DE_NOVO_ASSEMBLY_H_
#define DE_NOVO_ASSEMBLY_H_

#include <algorithm>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <mutex>
#include <htslib/sam.h>

#include "../libs/exact_cache.h"
#include "assembly_utils.h"
#include "sw_utils.h"
#include "ungapped_alignment.h"
#include "types.h"
#include "dc_remapper.h"

std::mutex failed_assembly_mtx;

struct build_graph_cache_value_t {
	std::vector<int> out_edges_added;
	std::vector<std::vector<edge_t> > l_adj_added, l_adj_rev_added;
};

void build_graph(std::vector<std::string>& read_seqs, std::vector<int>& order, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev,
		double max_seq_error, int min_overlap) {

	int n = read_seqs.size();
	surveyor_cache::key_builder_t key_builder;
	key_builder.add(max_seq_error);
	key_builder.add(min_overlap);
	uint64_t order_size = order.size();
	key_builder.add(order_size);
	for (int value : order) key_builder.add(value);
	uint64_t reads_size = read_seqs.size();
	key_builder.add(reads_size);
	for (const std::string& read_seq : read_seqs) key_builder.add_string(read_seq);
	std::string cache_key = key_builder.take();
	static thread_local surveyor_cache::exact_cache_t<build_graph_cache_value_t> cache("SURVEYOR_BUILD_GRAPH_CACHE_LIMIT", 256);
	const build_graph_cache_value_t* cached = cache.find(cache_key);
	if (cached) {
		for (int i = 0; i < n; i++) {
			out_edges[i] += cached->out_edges_added[i];
			l_adj[i].insert(l_adj[i].end(), cached->l_adj_added[i].begin(), cached->l_adj_added[i].end());
			l_adj_rev[i].insert(l_adj_rev[i].end(), cached->l_adj_rev_added[i].begin(), cached->l_adj_rev_added[i].end());
		}
		return;
	}

	build_graph_cache_value_t additions;
	additions.out_edges_added = out_edges;
	std::vector<size_t> l_adj_sizes(n), l_adj_rev_sizes(n);
	for (int i = 0; i < n; i++) {
		l_adj_sizes[i] = l_adj[i].size();
		l_adj_rev_sizes[i] = l_adj_rev[i].size();
	}

	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			std::string& s1 = read_seqs[i];
			std::string& s2 = read_seqs[j];
			if (is_contained(s1, s2, max_seq_error)) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(s1, s2, 1, -4, max_seq_error, min_overlap);
			bool spa1_rejected = reject_assembly_overlap(s1, s2, spa1.overlap);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(s2, s1, 1, -4, max_seq_error, min_overlap);
			bool spa2_rejected = reject_assembly_overlap(s2, s1, spa2.overlap);

			bool reverse = spa2.overlap && (!spa1.overlap || spa1.score < spa2.score);
			int from = reverse ? j : i, to = reverse ? i : j;
			const suffix_prefix_aln_t& spa = reverse ? spa2 : spa1;
			bool rejected = reverse ? spa2_rejected : spa1_rejected;
			if (spa.overlap && order[from] <= order[to] && !rejected) {
				add_assembly_edge(from, {to, spa.score, spa.overlap}, out_edges, l_adj, l_adj_rev);
			}
		}
	}

	additions.l_adj_added.resize(n);
	additions.l_adj_rev_added.resize(n);
	for (int i = 0; i < n; i++) {
		additions.out_edges_added[i] = out_edges[i] - additions.out_edges_added[i];
		additions.l_adj_added[i].assign(l_adj[i].begin() + l_adj_sizes[i], l_adj[i].end());
		additions.l_adj_rev_added[i].assign(l_adj_rev[i].begin() + l_adj_rev_sizes[i], l_adj_rev[i].end());
	}
	cache.store(std::move(cache_key), std::move(additions));
}

std::vector<assembly_consensus_t> assemble_contigs(const std::vector<assembly_consensus_t>& read_contigs, const std::vector<path_permission_t>& path_permissions,
		std::vector<int>& order, config_t& config, stats_t& stats) {

	std::vector<std::string> read_seqs;
	for (const assembly_consensus_t& contig : read_contigs) read_seqs.push_back(contig.seq);

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);

	build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, config.max_seq_error, config.min_clip_len);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	if (rev_topological_order.size() < n) {
		// A cycle means a repeat. An overlap with a mismatch where both bases have quality >= 20 likely joins different
		// repeat copies or haplotypes, so remove it; low-quality mismatches are likely sequencing errors, so those overlaps stay.
		remove_assembly_edges(out_edges, l_adj, l_adj_rev, [&read_contigs](int from, const edge_t& e) {
			return e.score != e.overlap && count_hq_mismatches(read_contigs[from], read_contigs[e.next], e.overlap, 20) > 0;
		});
		rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
		auto prune_short_overlaps = [&]() {
			for (int min_overlap = config.min_clip_len+10; rev_topological_order.size() < n && min_overlap <= stats.read_len/2; min_overlap += 10) {
				remove_assembly_edges(out_edges, l_adj, l_adj_rev, [min_overlap](int from, const edge_t& e) { return e.overlap < min_overlap; });
				rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
			}
		};
		prune_short_overlaps();

		// Still cyclic: rebuild with exact overlaps only and prune again.
		if (rev_topological_order.size() < n) {
			std::fill(out_edges.begin(), out_edges.end(), 0);
			for (int i = 0; i < n; i++) {
				l_adj[i].clear();
				l_adj_rev[i].clear();
			}
			build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 0.0, config.min_clip_len);
			rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);
			prune_short_overlaps();
		}
		if (rev_topological_order.size() < n) {
			return {{"HAS_CYCLE", {}}};
		}
	}

	std::vector<const assembly_consensus_t*> graph_contigs;
	for (const assembly_consensus_t& contig : read_contigs) graph_contigs.push_back(&contig);
	std::vector<bool> used(n);
	return assemble_graph_paths(graph_contigs, l_adj_rev, rev_topological_order, path_permissions, used);
}

std::vector<std::string> assemble_sequences(std::string contig_name, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals, config_t& config, stats_t& stats) {

	enum observation_source_t { STABLE_READ, CACHED_MATE, SEMI_MAPPED_READ, CLIP_CONSENSUS };
	struct observation_t {
		bam1_t* read;
		const std::string* clip_qual;
		observation_source_t source;
		bool reverse;
	};
	struct input_t {
		std::string seq;
		path_permission_t permissions;
		observation_t observation;
	};
	std::vector<input_t> left_stable, unstable, right_stable;
	if (r_cluster->clip_consensus) {
		left_stable.push_back({r_cluster->clip_consensus->sequence, {true, false}, {nullptr, &r_cluster->clip_consensus->qual, CLIP_CONSENSUS, false}});
	}
	for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
		std::string seq = get_sequence(read.get());
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		rc(mate_seq);
		left_stable.push_back({seq, {!is_left_clipped(read.get(), 1), true}, {read.get(), nullptr, STABLE_READ, false}});
		unstable.push_back({mate_seq, {true, true}, {read.get(), nullptr, CACHED_MATE, true}});
	}
	for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
		std::string seq = get_sequence(read.get());
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		right_stable.push_back({seq, {true, !is_right_clipped(read.get(), 1)}, {read.get(), nullptr, STABLE_READ, false}});
		unstable.push_back({mate_seq, {true, true}, {read.get(), nullptr, CACHED_MATE, false}});
	}
	if (l_cluster->clip_consensus) {
		right_stable.push_back({l_cluster->clip_consensus->sequence, {false, true}, {nullptr, &l_cluster->clip_consensus->qual, CLIP_CONSENSUS, false}});
	}

	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		unstable.push_back({read_seq, {!is_left_clipped(read, 1), !is_right_clipped(read, 1)}, {read, nullptr, SEMI_MAPPED_READ, true}});
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		unstable.push_back({read_seq, {!is_left_clipped(read, 1), !is_right_clipped(read, 1)}, {read, nullptr, SEMI_MAPPED_READ, false}});
	}

	const int TOO_MANY_READS = stats.get_max_depth(contig_name) * stats.max_is / stats.read_len;
	if (unstable.size() + left_stable.size() + right_stable.size() >= TOO_MANY_READS) return {"TOO_MANY_READS"};

	std::vector<assembly_consensus_t> read_contigs;
	std::vector<path_permission_t> path_permissions;
	std::vector<int> order;
	int group_order = 1;
	for (std::vector<input_t>* group : {&left_stable, &unstable, &right_stable}) {
		for (input_t& input : *group) {
			assembly_consensus_t contig{std::move(input.seq), {}};
			contig.votes.resize(contig.seq.length());
			observation_t& observation = input.observation;
			std::vector<uint8_t> quals;
			if (observation.source == STABLE_READ) quals = get_qual(observation.read);
			else if (observation.source == CACHED_MATE) quals = decode_qualities(get_mate_qual(observation.read, matequals));
			else if (observation.source == SEMI_MAPPED_READ) quals = get_qual(observation.read, true);
			else quals = decode_qualities(*observation.clip_qual);
			if (observation.reverse) std::reverse(quals.begin(), quals.end());
			add_read_votes(contig, quals.empty() ? nullptr : quals.data());
			read_contigs.push_back(std::move(contig));
			path_permissions.push_back(input.permissions);
			order.push_back(group_order);
		}
		group_order++;
	}
	std::vector<assembly_consensus_t> contigs = assemble_contigs(read_contigs, path_permissions, order, config, stats);
	std::vector<std::string> assembled_sequences;
	for (assembly_consensus_t& contig : contigs) assembled_sequences.push_back(std::move(contig.seq));

	return assembled_sequences;
}

std::shared_ptr<sv_t> detect_de_novo_insertion(std::string& contig_name, chr_seqs_map_t& contigs,
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		std::ostream& assembly_failed_no_seq, std::ostream& assembly_failed_cycle_writer, std::ostream& assembly_failed_too_many_reads_writer,
		StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
		config_t& config, stats_t& stats) {

	std::vector<std::string> assembled_sequences = assemble_sequences(contig_name, r_cluster, l_cluster, mateseqs, matequals, config, stats);

	std::string ins_full_id = "NO_ID";
	if (assembled_sequences.empty()) {
		failed_assembly_mtx.lock();
		assembly_failed_no_seq << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_no_seq << contig_name << " " << l_cluster->start << " - INS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "HAS_CYCLE") {
		failed_assembly_mtx.lock();
		assembly_failed_cycle_writer << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_cycle_writer << contig_name << " " << l_cluster->start << " - INS HAS_CYCLE" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "TOO_MANY_READS") {
		failed_assembly_mtx.lock();
		assembly_failed_too_many_reads_writer << ins_full_id << " " << contig_name << " " << r_cluster->end << " + ";
		assembly_failed_too_many_reads_writer << contig_name << " " << l_cluster->start << " - INS TOO_MANY_READS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	}

	// remap assembled sequences to find breakpoints
	char* contig_seq = contigs.get_seq(contig_name);
	int contig_len = contigs.get_len(contig_name);

	std::string assembled_sequence = assembled_sequences[0];

	// divide the assembled sequence into two and remap them
	int extend = 2*assembled_sequence.length();
	int remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
	remap_region_start = std::max(0, remap_region_start);
	int remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
	remap_region_end = std::min(remap_region_end, contig_len-1);

	StripedSmithWaterman::Alignment aln;
	std::string full_assembled_seq;
	std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contig_seq, assembled_sequence, std::string(assembled_sequence.length(), 'I'), remap_region_start, remap_region_end, 
		remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config);
	std::shared_ptr<sv_t> chosen_ins = NULL;
	for (std::shared_ptr<sv_t> sv : svs) {
		if (sv->svtype() == "INS" && sv->svsize() >= config.min_assembled_ins_size) {
			chosen_ins = sv;
			full_assembled_seq = assembled_sequence;
			break;
		}
	}
	if (chosen_ins == NULL) { // cannot identify an SV insertion, check if it is an incomplete assembly
		if (assembled_sequences.size() == 1) return NULL;

		extend = 2*(assembled_sequence.length() + assembled_sequences[1].length());
		remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
		if (remap_region_start < 0) remap_region_start = 0;
		remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
		if (remap_region_end >= contig_len) remap_region_end = contig_len-1;

		std::vector<std::shared_ptr<sv_t>> svs1 = detect_svs_from_junction(contig_name, contig_seq, assembled_sequence + "-" + assembled_sequences[1], std::string(assembled_sequence.length() + 1 + assembled_sequences[1].length(), 'I'),
			remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config);
		std::vector<std::shared_ptr<sv_t>> svs2 = detect_svs_from_junction(contig_name, contig_seq, assembled_sequences[1] + "-" + assembled_sequence, std::string(assembled_sequences[1].length() + 1 + assembled_sequence.length(), 'I'),
			remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config);
		bool sv1_is_ins = !svs1.empty() && svs1[0]->svtype() == "INS" && svs1[0]->svsize() >= config.min_assembled_ins_size;
		bool sv2_is_ins = !svs2.empty() && svs2[0]->svtype() == "INS" && svs2[0]->svsize() >= config.min_assembled_ins_size;
		if (sv1_is_ins && !sv2_is_ins) {
			chosen_ins = svs1[0];
		} else if (!sv1_is_ins && sv2_is_ins) {
			chosen_ins = svs2[0];
		} else if (sv1_is_ins && sv2_is_ins) {
			int sv1_score = svs1[0]->left_anchor_aln->best_score + svs1[0]->right_anchor_aln->best_score;
			int sv2_score = svs2[0]->left_anchor_aln->best_score + svs2[0]->right_anchor_aln->best_score;
			if (sv1_score >= sv2_score) chosen_ins = svs1[0];
			else chosen_ins = svs2[0];
		}

		if (!svs1.empty() && chosen_ins == svs1[0]) full_assembled_seq = assembled_sequence + "-" + assembled_sequences[1];
		else if (!svs2.empty() && chosen_ins == svs2[0]) full_assembled_seq = assembled_sequences[1] + "-" + assembled_sequence;
	}

	if (chosen_ins == NULL) return NULL;

	int sep = chosen_ins->ins_seq.find("-");
	if (sep != std::string::npos && (sep < config.min_clip_len || chosen_ins->ins_seq.length()-sep-1 < config.min_clip_len)) return NULL; // "-" is too close to the edge

	int assembled_reads_count = 0;

	StripedSmithWaterman::Filter filter;
	if (r_cluster->clip_consensus) {
		ungapped_aln_t ungapped_aln = best_ungapped_aln(r_cluster->clip_consensus->sequence.c_str(), r_cluster->clip_consensus->sequence.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			chosen_ins->rc_consensus = r_cluster->clip_consensus;
		}
	}
	for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		std::string mate_qual = get_mate_qual(read.get(), matequals);
		rc(mate_seq);
		mate_qual = std::string(mate_qual.rbegin(), mate_qual.rend());
		StripedSmithWaterman::Alignment aln;
		ungapped_aln_t ungapped_aln = best_ungapped_aln(mate_seq.c_str(), mate_seq.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}
	for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		std::string mate_qual = get_mate_qual(read.get(), matequals);
		StripedSmithWaterman::Alignment aln;
		ungapped_aln_t ungapped_aln = best_ungapped_aln(mate_seq.c_str(), mate_seq.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}
	if (l_cluster->clip_consensus) {
		ungapped_aln_t ungapped_aln = best_ungapped_aln(l_cluster->clip_consensus->sequence.c_str(), l_cluster->clip_consensus->sequence.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			chosen_ins->lc_consensus = l_cluster->clip_consensus;
		}
	}

	// start and end of inserted sequence within the full assembled sequence
	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}

	chosen_ins->source = "DE_NOVO_ASSEMBLY";
	if (assembled_reads_count <= 2) {
		return NULL;
	}
	chosen_ins->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly
	return chosen_ins;
}

#endif /* DE_NOVO_ASSEMBLY_H_ */
