#ifndef DE_NOVO_ASSEMBLY_H_
#define DE_NOVO_ASSEMBLY_H_

#include <algorithm>
#include <cstdint>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <htslib/sam.h>

#include "../libs/exact_cache.h"
#include "assembly_utils.h"
#include "sw_utils.h"
#include "ungapped_alignment.h"
#include "types.h"
#include "dc_remapper.h"
#include "hp_mismatch_rate_thresholds.h"
#include "read_correction.h"

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
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(s2, s1, 1, -4, max_seq_error, min_overlap);

			bool reverse = spa2.overlap && (!spa1.overlap || spa1.score < spa2.score);
			int from = reverse ? j : i, to = reverse ? i : j;
			const suffix_prefix_aln_t& spa = reverse ? spa2 : spa1;
			if (spa.overlap && order[from] <= order[to] && !reject_assembly_overlap(read_seqs[from], read_seqs[to], spa.overlap)) {
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
	// High-quality mismatches likely join different repeat copies or haplotypes; low-quality mismatches may be sequencing errors.
	remove_assembly_edges(out_edges, l_adj, l_adj_rev, [&read_contigs](int from, const edge_t& e) {
		return e.score != e.overlap && count_hq_mismatches(read_contigs[from], read_contigs[e.next], e.overlap, 20) > 0;
	});

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	if (rev_topological_order.size() < n) {
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

// A read as the de novo assembly uses it: after read correction, with the quality of each of its bases.
struct assembly_read_t {
	std::string seq;
	std::vector<uint8_t> quals;
};

std::vector<std::string> assemble_sequences(std::string contig_name, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals, config_t& config, stats_t& stats, const hp_tail_quality_model_t& hp_tail_quality_model, open_samFile_t* flank_bam, std::vector<assembly_read_t>& assembly_reads) {
	assembly_reads.clear();

	const int TOO_MANY_READS = stats.get_max_depth(contig_name) * stats.max_is / stats.read_len;
	if (count_assembly_reads(r_cluster, l_cluster) >= TOO_MANY_READS) return {"TOO_MANY_READS"};

	enum observation_source_t { STABLE_READ, CACHED_MATE, SEMI_MAPPED_READ };
	struct observation_t {
		bam1_t* read;
		observation_source_t source;
		bool reverse;
	};
	struct input_t {
		std::string seq;
		path_permission_t permissions;
		observation_t observation;
	};
	std::vector<input_t> left_stable, unstable, right_stable;
	assembly_read_filter_t read_filter(r_cluster, l_cluster, true);
	if (r_cluster->clip_consensus) {
		for (std::shared_ptr<bam1_t> read : r_cluster->clip_consensus->clipped_reads) {
			if (!read_filter.add_read(read.get(), true)) continue;
			left_stable.push_back({get_sequence(read.get()), {true, false}, {read.get(), STABLE_READ, false}});
		}
	}
	for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
		if (read_filter.add_read(read.get())) {
			left_stable.push_back({get_sequence(read.get()), {!is_left_clipped(read.get(), 1), true}, {read.get(), STABLE_READ, false}});
		}
		if (read_filter.add_mate(read.get())) {
			std::string mate_seq = get_mate_seq(read.get(), mateseqs);
			rc(mate_seq);
			unstable.push_back({mate_seq, {true, true}, {read.get(), CACHED_MATE, true}});
		}
	}
	for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
		if (read_filter.add_read(read.get())) {
			right_stable.push_back({get_sequence(read.get()), {true, !is_right_clipped(read.get(), 1)}, {read.get(), STABLE_READ, false}});
		}
		if (read_filter.add_mate(read.get())) {
			unstable.push_back({get_mate_seq(read.get(), mateseqs), {true, true}, {read.get(), CACHED_MATE, false}});
		}
	}
	if (l_cluster->clip_consensus) {
		for (std::shared_ptr<bam1_t> read : l_cluster->clip_consensus->clipped_reads) {
			if (!read_filter.add_read(read.get(), true)) continue;
			right_stable.push_back({get_sequence(read.get()), {false, true}, {read.get(), STABLE_READ, false}});
		}
	}

	for (bam1_t* read : r_cluster->semi_mapped_mates) {
		if (!read_filter.add_read(read)) continue;
		bool reverse = is_samestr(read);
		std::string read_seq = get_sequence(read);
		if (reverse) rc(read_seq);
		left_stable.push_back({read_seq, {!(reverse ? is_right_clipped(read, 1) : is_left_clipped(read, 1)), true}, {read, STABLE_READ, reverse}});
	}
	for (bam1_t* read : l_cluster->semi_mapped_mates) {
		if (!read_filter.add_read(read)) continue;
		bool reverse = is_samestr(read);
		std::string read_seq = get_sequence(read);
		if (reverse) rc(read_seq);
		right_stable.push_back({read_seq, {true, !(reverse ? is_left_clipped(read, 1) : is_right_clipped(read, 1))}, {read, STABLE_READ, reverse}});
	}

	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		if (!read_filter.add_read(read)) continue;
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		unstable.push_back({read_seq, {!is_left_clipped(read, 1), !is_right_clipped(read, 1)}, {read, SEMI_MAPPED_READ, true}});
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		if (!read_filter.add_read(read)) continue;
		std::string read_seq = get_sequence(read, true);
		unstable.push_back({read_seq, {!is_left_clipped(read, 1), !is_right_clipped(read, 1)}, {read, SEMI_MAPPED_READ, false}});
	}

	std::vector<assembly_consensus_t> read_contigs;
	std::vector<read_correction::correction_input_t> correction_inputs;
	std::vector<path_permission_t> path_permissions;
	std::vector<int> order;
	hp_tail_quality_table_t quality_cache;
	int group_order = 1;
	for (std::vector<input_t>* group : {&left_stable, &unstable, &right_stable}) {
		for (input_t& input : *group) {
			assembly_consensus_t contig{std::move(input.seq), {}};
			contig.votes.resize(contig.seq.length());
			observation_t& observation = input.observation;
			std::vector<uint8_t> quals;
			if (observation.source == STABLE_READ || observation.source == SEMI_MAPPED_READ) {
				quals = recalibrate_clip_read_qualities(observation.read, config, hp_tail_quality_model, quality_cache);
				if ((observation.source == STABLE_READ && observation.reverse) || (observation.source == SEMI_MAPPED_READ && bam_is_rev(observation.read) != observation.reverse)) std::reverse(quals.begin(), quals.end());
			} else {
				quals = decode_qualities(get_mate_qual(observation.read, matequals));
				if (observation.reverse) std::reverse(quals.begin(), quals.end());
				hts_pair_pos_t hp = longest_homopolymer(contig.seq.c_str(), contig.seq.length());
				int hp_len = hp.end-hp.beg;
				if (!quals.empty() && hp_len >= MIN_REF_HP_LEN) {
					std::vector<hp_read_observation_t> observations{{contig.seq, std::move(quals), int(hp.beg), observation.reverse, 0.0, int(contig.seq.length()-hp.end)}};
					recalibrate_hp_3p_tail_qualities(observations, hp_len, contig.seq[hp.beg], hp_tail_quality_model, quality_cache);
					quals = std::move(observations[0].quals);
				}
			}
			add_read_votes(contig, quals.empty() ? nullptr : quals.data());
			correction_inputs.push_back({contig.seq, quals, !quals.empty()});
			read_contigs.push_back(std::move(contig));
			path_permissions.push_back(input.permissions);
			order.push_back(group_order);
		}
		group_order++;
	}
	if (std::any_of(correction_inputs.begin(), correction_inputs.end(), [](const read_correction::correction_input_t& input) { return input.target; })) {
		std::size_t graph_input_count = correction_inputs.size();
		std::unordered_set<std::string> graph_sequences;
		for (const read_correction::correction_input_t& input : correction_inputs) graph_sequences.insert(input.seq);
		hts_pos_t beg = std::max(hts_pos_t(0), std::min(r_cluster->start, l_cluster->start)), end = std::max(r_cluster->end, l_cluster->end);
		hts_pos_t max_flank = stats.get_max_depth(contig_name) * (end-beg+stats.read_len) / stats.read_len, kept_flank = 0;
		int tid = bam_name2id(flank_bam->header, contig_name.c_str());
		hts_itr_t* flank_iter = sam_itr_queryi(flank_bam->idx, tid, beg, end);
		if (flank_iter) {
			bam1_t* flank_read = bam_init1();
			while (sam_itr_next(flank_bam->file, flank_iter, flank_read) >= 0) {
				if (flank_read->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL)) continue;
				std::string seq = get_sequence(flank_read);
				if (graph_sequences.count(seq)) continue;
				std::string reverse_seq = seq;
				rc(reverse_seq);
				if (graph_sequences.count(reverse_seq)) continue;
				correction_inputs.push_back({std::move(seq), get_qual(flank_read), false});
				if (++kept_flank > max_flank) {
					correction_inputs.resize(graph_input_count);
					break;
				}
			}
			bam_destroy1(flank_read);
			sam_itr_destroy(flank_iter);
		}
		std::vector<read_correction::correction_flip_t> flips = read_correction::correct_reads(correction_inputs);
		for (const read_correction::correction_flip_t& flip : flips) {
			const read_correction::correction_input_t& in = correction_inputs[flip.input];
			int qual = in.quals.size() == in.seq.size() ? in.quals[flip.pos] : read_correction::DEFAULT_QUAL;
			read_contigs[flip.input].seq[flip.pos] = flip.new_base;
			read_contigs[flip.input].votes[flip.pos] = base_votes_t();
			read_contigs[flip.input].votes[flip.pos].add(flip.new_base, qual);
		}
	}
	for (const assembly_consensus_t& read : read_contigs) {
		assembly_read_t assembly_read{read.seq, std::vector<uint8_t>(read.seq.length(), 0)};
		for (int i = 0; i < read.seq.length(); i++) {
			int base = base_to_index(read.seq[i]);
			if (base >= 0) assembly_read.quals[i] = std::min(60, std::max(0, read.votes[i].quals[base]));
		}
		assembly_reads.push_back(std::move(assembly_read));
	}
	std::vector<assembly_consensus_t> contigs = assemble_contigs(read_contigs, path_permissions, order, config, stats);
	for (assembly_consensus_t& contig : contigs) {
		if (contig.seq == "HAS_CYCLE") continue;
		std::vector<base_votes_t> votes(contig.seq.length());
		for (const assembly_consensus_t& read : read_contigs) {
			int len = read.seq.length();
			ungapped_aln_t aln = best_ungapped_aln(read.seq.c_str(), len, contig.seq.c_str(), contig.seq.length(), len-std::min(len, 40));
			if (aln.mismatch_rate() > 0.1) continue;
			bool hq_mismatch = false;
			for (int i = aln.query_begin, j = aln.ref_begin; i < aln.query_end; i++, j++) {
				if (read.seq[i] == contig.seq[j]) continue;
				int read_base = base_to_index(read.seq[i]), contig_base = base_to_index(contig.seq[j]);
				if (read_base < 0 || contig_base < 0) continue;
				const base_votes_t& read_votes = read.votes[i];
				const base_votes_t& contig_votes = contig.votes[j];
				int read_qual = read_votes.quals[read_base];
				int contig_qual = contig_votes.diff_qual(contig_base);
				if (read_qual >= 20 && contig_qual >= 20) {
					hq_mismatch = true;
					break;
				}
			}
			if (hq_mismatch) continue;
			for (int i = aln.query_begin, j = aln.ref_begin; i < aln.query_end; i++, j++) votes[j].merge(read.votes[i]);
		}
		for (int i = 0; i < contig.seq.length(); i++) contig.seq[i] = votes[i].best_base(contig.seq[i], true);
		contig.votes = std::move(votes);
	}
	std::vector<std::string> assembled_sequences;
	for (assembly_consensus_t& contig : contigs) assembled_sequences.push_back(std::move(contig.seq));

	return assembled_sequences;
}

std::shared_ptr<sv_t> finalize_de_novo_insertion_from_junction(std::shared_ptr<sv_t> chosen_ins, const std::string& full_assembled_seq,
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals, config_t& config) {

	int sep = chosen_ins->ins_seq.find("-");
	if (sep != std::string::npos && (sep < config.min_clip_len || chosen_ins->ins_seq.length()-sep-1 < config.min_clip_len)) return NULL; // "-" is too close to the edge

	int assembled_reads_count = 0;
	// Reads already counted (through their mate, or as half-mapped reads), so that a clipped read is not counted twice.
	std::set<std::pair<std::string, bool> > counted_reads;
	auto count_clipped_reads = [&](std::shared_ptr<consensus_t> consensus) {
		if (!consensus) return;
		for (std::shared_ptr<bam1_t> read : consensus->clipped_reads) {
			if (!counted_reads.insert({bam_get_qname(read.get()), bool(read->core.flag & BAM_FREAD1)}).second) continue;
			std::string read_seq = get_sequence(read.get());
			ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(),
				full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
			if (ungapped_aln.mismatch_rate() <= config.max_seq_error) assembled_reads_count++;
		}
	};

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
		counted_reads.insert({bam_get_qname(read.get()), !(read->core.flag & BAM_FREAD1)});
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
		counted_reads.insert({bam_get_qname(read.get()), !(read->core.flag & BAM_FREAD1)});
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
		counted_reads.insert({bam_get_qname(read), bool(read->core.flag & BAM_FREAD1)});
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(),
			full_assembled_seq.c_str(), full_assembled_seq.length(), config.min_clip_len - 1);
		counted_reads.insert({bam_get_qname(read), bool(read->core.flag & BAM_FREAD1)});
		if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
			assembled_reads_count++;
		}
	}

	// A complete insertion is also supported by the reads clipped at its breakpoints; for an incomplete one they only show the breakpoints.
	if (sep == std::string::npos) {
		count_clipped_reads(r_cluster->clip_consensus);
		count_clipped_reads(l_cluster->clip_consensus);
	}

	chosen_ins->source = "DE_NOVO_ASSEMBLY";
	if (assembled_reads_count <= 2) {
		return NULL;
	}
	chosen_ins->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly
	return chosen_ins;
}

std::shared_ptr<sv_t> detect_de_novo_insertion_from_junction(std::string& contig_name, chr_seqs_map_t& contigs, const std::string& junction_seq,
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		StripedSmithWaterman::Aligner& aligner_to_base, config_t& config, stats_t& stats) {

	int extend = 2*(junction_seq.length()-(junction_seq.find('-') != std::string::npos));
	int remap_region_start = std::max(0, int(std::min(r_cluster->end, l_cluster->start)-extend));
	int remap_region_end = std::min(int(std::max(r_cluster->end, l_cluster->start)+extend), int(contigs.get_len(contig_name))-1);
	std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contigs.get_seq(contig_name), junction_seq, std::string(junction_seq.length(), 'I'), remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config, -1, 1000, false, nullptr, nullptr, true);
	for (std::shared_ptr<sv_t> sv : svs) {
		if (sv->svtype() == "INS" && sv->svsize() >= config.min_assembled_ins_size) {
			return finalize_de_novo_insertion_from_junction(sv, junction_seq, r_cluster, l_cluster, mateseqs, matequals, config);
		}
	}
	return NULL;
}

std::shared_ptr<sv_t> detect_de_novo_insertion(std::string& contig_name, chr_seqs_map_t& contigs,
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		std::ostream& assembly_failed_no_seq, std::ostream& assembly_failed_cycle_writer, std::ostream& assembly_failed_too_many_reads_writer,
		StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
		config_t& config, stats_t& stats, const hp_tail_quality_model_t& hp_tail_quality_model, open_samFile_t* flank_bam, std::string& junction_seq, std::vector<std::string>& selected_contigs, std::vector<assembly_read_t>& assembly_reads) {

	junction_seq.clear(); // Set only when the selected insertion passes the support checks.
	selected_contigs.clear();
	std::vector<std::string> assembled_sequences = assemble_sequences(contig_name, r_cluster, l_cluster, mateseqs, matequals, config, stats, hp_tail_quality_model, flank_bam, assembly_reads);

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

	selected_contigs = assembled_sequences;

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
		remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config, -1, 1000, false, nullptr, nullptr, true);
	std::shared_ptr<sv_t> chosen_ins = NULL;
	for (std::shared_ptr<sv_t> sv : svs) {
		if (sv->svtype() == "INS" && sv->svsize() >= config.min_assembled_ins_size) {
			chosen_ins = sv;
			full_assembled_seq = assembled_sequence;
			break;
		}
	}
	if (chosen_ins == NULL) {
		if (assembled_sequences.size() == 1) return NULL;

		for (int order = 0; order < 2 && chosen_ins == NULL; order++) {
			const std::string& left = assembled_sequences[order];
			const std::string& right = assembled_sequences[1-order];
			int overlap = std::min(left.length(), right.length());
			for (; overlap >= 30; overlap--) {
				if (left.compare(left.length()-overlap, overlap, right, 0, overlap) == 0) break;
			}
			if (overlap < 30) continue;

			std::string merged = left + right.substr(overlap);
			extend = 2*merged.length();
			remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
			remap_region_start = std::max(0, remap_region_start);
			remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
			remap_region_end = std::min(remap_region_end, contig_len-1);
			std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contig_seq, merged, std::string(merged.length(), 'I'), remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config, -1, 1000, false, nullptr, nullptr, true);
			for (std::shared_ptr<sv_t> sv : svs) {
				if (sv->svtype() == "INS" && sv->svsize() >= config.min_assembled_ins_size) {
					chosen_ins = sv;
					full_assembled_seq = merged;
					break;
				}
			}
		}
	}

	if (chosen_ins == NULL) { // cannot identify an SV insertion, check if it is an incomplete assembly
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

	if (chosen_ins == NULL) {
		extend = 2*assembled_sequences[1].length();
		remap_region_start = std::min(r_cluster->end, l_cluster->start)-extend;
		remap_region_start = std::max(0, remap_region_start);
		remap_region_end = std::max(r_cluster->end, l_cluster->start)+extend;
		remap_region_end = std::min(remap_region_end, contig_len-1);
		std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contig_seq, assembled_sequences[1], std::string(assembled_sequences[1].length(), 'I'), remap_region_start, remap_region_end, remap_region_start, remap_region_end, aligner_to_base, 0, 0, stats, config, -1, 1000, false, nullptr, nullptr, true);
		for (std::shared_ptr<sv_t> sv : svs) {
			if (sv->svtype() == "INS" && sv->svsize() >= config.min_assembled_ins_size) {
				chosen_ins = sv;
				full_assembled_seq = assembled_sequences[1];
				break;
			}
		}
	}

	if (chosen_ins != NULL) chosen_ins = finalize_de_novo_insertion_from_junction(chosen_ins, full_assembled_seq, r_cluster, l_cluster, mateseqs, matequals, config);
	if (chosen_ins == NULL && assembled_sequences.size() > 8) {
		std::shared_ptr<sv_t> third_ins = detect_de_novo_insertion_from_junction(contig_name, contigs, assembled_sequences[2], r_cluster, l_cluster, mateseqs, matequals, aligner_to_base, config, stats);
		if (third_ins != NULL && !third_ins->incomplete_ins_seq()) {
			bool valid_anchors = true;
			for (const auto& anchor : {third_ins->left_anchor_aln, third_ins->right_anchor_aln}) {
				if (!anchor || anchor->best_score < 50 || anchor->seq_len < 50 || anchor->end-anchor->start+1 < 50) valid_anchors = false;
			}
			if (valid_anchors) {
				chosen_ins = third_ins;
				full_assembled_seq = assembled_sequences[2];
			}
		}
	}
	if (chosen_ins) junction_seq = full_assembled_seq;
	return chosen_ins;
}

#endif /* DE_NOVO_ASSEMBLY_H_ */
