#ifndef GUIDED_REFERENCE_ASSEMBLE_H
#define GUIDED_REFERENCE_ASSEMBLE_H

#include <functional>
#include <map>
#include <memory>
#include <queue>
#include <unordered_map>

#include "../libs/ssw_selected.h"
#include "htslib/sam.h"
#include "assembly_utils.h"
#include "sw_utils.h"
#include "ungapped_alignment.h"
#include "utils.h"
#include "sam_utils.h"
#include "dc_remapper.h"
#include "de_novo_assembly.h"
#include "remapping.h"

struct aln_guided_contig_t {
	assembly_consensus_t* contig;
	StripedSmithWaterman::Alignment aln;
};

void build_aln_guided_graph(std::vector<aln_guided_contig_t>& alns, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, config_t& config, bool allow_at_homopolymer_overlaps = false) {
	std::sort(alns.begin(), alns.end(),
			[](aln_guided_contig_t& aln1, aln_guided_contig_t& aln2) {
		return get_unclipped_start(aln1.aln) < get_unclipped_start(aln2.aln);
	});

	for (int i = 0; i < alns.size(); i++) {
		for (int j = i+1; j < alns.size() && get_unclipped_end(alns[i].aln)-get_unclipped_start(alns[j].aln) >= config.min_clip_len; j++) {
			if (is_contained(alns[i].contig->seq, alns[j].contig->seq, config.max_seq_error)) continue;
			int min_overlap = config.min_clip_len;
			// In principle, the suffix-prefix alignment should be the same as the overlap between the two alignments on the reference
			// This may not be the case if there are insertions or deletions in between the reads. However, in that case,
			// we expect the first read to be clipped on the right and/or the second read to be clipped on the left.
			if (!is_right_clipped(alns[i].aln) && !is_left_clipped(alns[j].aln)) {
				int unclipped_i_end = get_unclipped_end(alns[i].aln);
				int unclipped_j_start = get_unclipped_start(alns[j].aln);
				min_overlap = std::max(min_overlap, unclipped_i_end - unclipped_j_start);
			}
			suffix_prefix_aln_t spa = aln_suffix_prefix(alns[i].contig->seq, alns[j].contig->seq, 1, -4, config.max_seq_error, min_overlap);
			bool rejected = reject_assembly_overlap(alns[i].contig->seq, alns[j].contig->seq, spa.overlap);
			// Contig joins may bridge an exact A/T-rich overlap; the completed junction
			// gets its HP length re-estimated from reads after assembly.
			if (rejected && allow_at_homopolymer_overlaps && spa.overlap >= 10 && spa.mismatches == 0) {
				char hp_base = is_homopolymer(alns[j].contig->seq.data(), spa.overlap);
				if (hp_base == 'A' || hp_base == 'T') rejected = false;
			}
			if (spa.overlap && !rejected) {
				add_assembly_edge(i, {j, spa.score, spa.overlap}, out_edges, l_adj, l_adj_rev);
			}
		}
	}

	// two major differences with the regular assembly:
	// 1 - by how the graph is defined, no cycle is possible here
	// 2 - in regular assembly, we only report contigs made of at last 2 reads. Here we report even single reads
}

// remove contigs that are fully contained within another contig
void remove_fully_contained(std::vector<assembly_consensus_t>& assembled_contigs, std::string& reference, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::vector<std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment> > assembled_seqs_w_alns;
	for (assembly_consensus_t& contig : assembled_contigs) {
		aligner.Align(contig.seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		assembled_seqs_w_alns.push_back({std::move(contig), aln});
	}
	std::sort(assembled_seqs_w_alns.begin(), assembled_seqs_w_alns.end(),
			[](const std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment>& p1, const std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment>& p2) {
		if (p1.second.ref_begin != p2.second.ref_begin) {
            return p1.second.ref_begin < p2.second.ref_begin;
        }
        return p1.second.ref_end > p2.second.ref_end;
	});
	int max_end = 0;
	assembled_contigs.clear();
	for (auto& p : assembled_seqs_w_alns) {
		if (max_end < p.second.ref_end) {
			assembled_contigs.push_back(std::move(p.first));
			max_end = p.second.ref_end;
		}
	}
}

std::vector<std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment> > retain_non_overlapping_contigs(
		std::vector<std::tuple<int, assembly_consensus_t, StripedSmithWaterman::Alignment> > contigs_w_score) {
	std::sort(contigs_w_score.begin(), contigs_w_score.end(),
			[](const std::tuple<int, assembly_consensus_t, StripedSmithWaterman::Alignment>& p1, const std::tuple<int, assembly_consensus_t, StripedSmithWaterman::Alignment>& p2) {
		return std::get<0>(p1) > std::get<0>(p2);
	});

	std::vector<bool> remove(contigs_w_score.size(), false);
	for (int i = contigs_w_score.size()-1; i >= 0; i--) {
		for (int j = 0; j < i; j++) {
			StripedSmithWaterman::Alignment& aln1 = std::get<2>(contigs_w_score[i]);
			StripedSmithWaterman::Alignment& aln2 = std::get<2>(contigs_w_score[j]);
			if (overlap(aln1.ref_begin, aln1.ref_end, aln2.ref_begin, aln2.ref_end)) {
				remove[i] = true;
				break;
			}
		}
	}

	std::vector<std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment> > retained;
	for (int i = 0; i < contigs_w_score.size(); i++) {
		if (!remove[i]) {
			retained.push_back({std::move(std::get<1>(contigs_w_score[i])), std::get<2>(contigs_w_score[i])});
		}
	}

	std::sort(retained.begin(), retained.end(),
			[](const std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment>& p1, const std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	return retained;
}

std::vector<assembly_consensus_t> generate_reference_guided_contigs(std::string reference, std::vector<assembly_consensus_t>& read_contigs,
	std::vector<assembly_consensus_t>& discarded_contigs, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, bool allow_at_homopolymer_overlaps = false) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;

	std::vector<aln_guided_contig_t> alns;
	int n = read_contigs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	for (assembly_consensus_t& contig : read_contigs) {
		harsh_aligner.Align(contig.seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		alns.push_back({&contig, aln});
	}
	build_aln_guided_graph(alns, out_edges, l_adj, l_adj_rev, config, allow_at_homopolymer_overlaps);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	std::vector<const assembly_consensus_t*> graph_contigs;
	for (const aln_guided_contig_t& aln_contig : alns) graph_contigs.push_back(aln_contig.contig);
	std::vector<bool> used(n);
	std::vector<path_permission_t> path_permissions(n, {true, true});
	std::vector<assembly_consensus_t> assembled_contigs = assemble_graph_paths(graph_contigs, l_adj_rev, rev_topological_order, path_permissions, used);
	for (int i = 0; i < n; i++) {
		if (!used[i]) discarded_contigs.push_back(*alns[i].contig);
	}
	return assembled_contigs;
}

// Keep joining compatible paths while preserving contigs that were not used by a path.
void join_reference_guided_contigs(std::string& reference, std::vector<assembly_consensus_t>& contigs, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, bool allow_at_homopolymer_overlaps = false) {
	while (true) {
		remove_fully_contained(contigs, reference, aligner);
		std::vector<assembly_consensus_t> temp;
		int prev_ass_seqs = contigs.size();
		contigs = generate_reference_guided_contigs(reference, contigs, temp, harsh_aligner, config, allow_at_homopolymer_overlaps);
		contigs.insert(contigs.end(), temp.begin(), temp.end());
		if (prev_ass_seqs <= contigs.size()) break;
	}
}

// Soft-clipped contig ends may contain HP-tail errors rather than insertion sequence.
// Remove their weak outer portions using quality margins, retaining supported overhangs
// and leaving the reference-aligned interior intact.
void trim_low_quality_clipped_ends(std::vector<assembly_consensus_t>& contigs, const std::string& reference, StripedSmithWaterman::Aligner& harsh_aligner, int min_clip_len, int min_qual) {
	StripedSmithWaterman::Filter filter;
	for (assembly_consensus_t& contig : contigs) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(contig.seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		const std::vector<base_votes_t>& votes = contig.trim_votes.empty() ? contig.votes : contig.trim_votes;
		auto margin = [&](int i) {
			int base = base_to_index(contig.seq[i]);
			return base < 0 ? 0 : std::max(0, votes[i].diff_qual(base));
		};
		int begin = 0, end = contig.seq.length();
		int left_clip = get_left_clip_size(aln), right_clip = get_right_clip_size(aln);
		// Keep the maximum low-quality deficit before reliable support outweighs it.
		if (left_clip >= min_clip_len) {
			int sum = 0, max_sum = 0;
			for (int i = 0; i < left_clip; i++) {
				sum += min_qual-margin(i);
				if (sum < 0) break;
				if (sum > max_sum) { max_sum = sum; begin = i+1; }
			}
		}
		if (right_clip >= min_clip_len) {
			int sum = 0, max_sum = 0;
			for (int i = contig.seq.length()-1; i >= int(contig.seq.length())-right_clip; i--) {
				sum += min_qual-margin(i);
				if (sum < 0) break;
				if (sum > max_sum) { max_sum = sum; end = i; }
			}
		}
		if ((begin == 0 && end == contig.seq.length()) || end-begin < min_clip_len) continue;
		contig.seq = contig.seq.substr(begin, end-begin);
		contig.votes = std::vector<base_votes_t>(contig.votes.begin()+begin, contig.votes.begin()+end);
		if (!contig.trim_votes.empty()) contig.trim_votes = std::vector<base_votes_t>(contig.trim_votes.begin()+begin, contig.trim_votes.begin()+end);
	}
}

std::vector<std::string> generate_reference_guided_consensus(std::string reference, 
		std::vector<std::string>& seqs_lf, std::vector<std::string>& seqs_is, std::vector<std::string>& seqs_rf,
		const std::vector<const uint8_t*>& quals_lf, const std::vector<const uint8_t*>& quals_is, const std::vector<const uint8_t*>& quals_rf,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, config_t& config, stats_t& stats,
		bool scaffold = true, bool retain_initial_contigs = false, const std::vector<const uint8_t*>* trim_quals = nullptr) {

	std::vector<std::string> read_seqs;
	for (std::string& seq : seqs_lf) {
		read_seqs.push_back(seq);
	}
	for (std::string& seq : seqs_is) {
		read_seqs.push_back(seq);
	}
	for (std::string& seq : seqs_rf) {
		read_seqs.push_back(seq);
	}
	std::vector<const uint8_t*> read_quals = quals_lf;
	read_quals.insert(read_quals.end(), quals_is.begin(), quals_is.end());
	read_quals.insert(read_quals.end(), quals_rf.begin(), quals_rf.end());

	std::vector<assembly_consensus_t> read_contigs;
	for (int i = 0; i < read_seqs.size(); i++) {
		read_contigs.push_back(make_read_contig(read_seqs[i], read_quals.empty() ? nullptr : read_quals[i], trim_quals == nullptr ? nullptr : (*trim_quals)[i]));
	}
	std::vector<assembly_consensus_t> discarded_contigs;
	std::vector<assembly_consensus_t> assembled_contigs = generate_reference_guided_contigs(reference, read_contigs, discarded_contigs, harsh_aligner, config);
	std::vector<std::string> initial_sequences;
	if (retain_initial_contigs) {
		for (const assembly_consensus_t& contig : assembled_contigs) initial_sequences.push_back(contig.seq);
	}

	join_reference_guided_contigs(reference, assembled_contigs, aligner, harsh_aligner, config, trim_quals != nullptr);
	// Join first so overlapping reads can protect an overhang, then trim before scoring.
	// Callers that omit trim qualities, including genotype, retain their existing behavior.
	if (trim_quals != nullptr) {
		trim_low_quality_clipped_ends(assembled_contigs, reference, harsh_aligner, config.min_clip_len, 35);
		// Removing unreliable ends can expose exact overlaps between previously separate contigs.
		join_reference_guided_contigs(reference, assembled_contigs, aligner, harsh_aligner, config, true);
	}
	// score the retained contigs
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::vector<int> scores(assembled_contigs.size());
	for (std::string& read_seq : read_seqs) {
		int best_score = 0;
		std::vector<int> read_scores(assembled_contigs.size());
		for (int i = 0; i < assembled_contigs.size(); i++) {
			std::string& contig = assembled_contigs[i].seq;
			harsh_aligner.Align(read_seq.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
			if (accept(aln, config.min_clip_len) && aln.sw_score > best_score) {
				best_score = aln.sw_score;
				read_scores[i] = aln.sw_score;
			}
		}
		for (int i = 0; i < assembled_contigs.size(); i++) {
			if (read_scores[i] == best_score) {
				scores[i] += best_score;
			}
		}
	}

	// sort by descending score
	std::vector<std::tuple<int, assembly_consensus_t, StripedSmithWaterman::Alignment> > contigs_w_score;
	for (int i = 0; i < assembled_contigs.size(); i++) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(assembled_contigs[i].seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		contigs_w_score.push_back({scores[i], std::move(assembled_contigs[i]), aln});
	}
	std::vector<std::pair<assembly_consensus_t, StripedSmithWaterman::Alignment>> retained = retain_non_overlapping_contigs(std::move(contigs_w_score));
	std::vector<assembly_consensus_t> retained_contigs;
	consensus_contigs_alns.clear();
	for (auto& p : retained) {
		retained_contigs.push_back(std::move(p.first));
		consensus_contigs_alns.push_back(p.second);
	}
	auto append_initial_sequences = [&](std::vector<std::string>& sequences) {
		for (std::string& seq : initial_sequences) {
			if (std::find(sequences.begin(), sequences.end(), seq) != sequences.end()) continue;
			harsh_aligner.Align(seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
			sequences.push_back(std::move(seq));
			consensus_contigs_alns.push_back(aln);
		}
	};

	if (!scaffold || retained_contigs.empty()) {
		std::vector<std::string> sequences;
		for (assembly_consensus_t& contig : retained_contigs) sequences.push_back(std::move(contig.seq));
		append_initial_sequences(sequences);
		return sequences;
	}

	// we assembled as much as possible using guidance from the reference
	// however, reads may be misaligned for in the case of an incomplete or rearranged reference
	// we can try to "scaffold" the contigs using the reads that were rejected during the assembly

	std::vector<path_permission_t> bridge_permissions(discarded_contigs.size(), {true, true});
	std::vector<int> bridge_order(discarded_contigs.size(), 2);
	std::vector<assembly_consensus_t> scaffolding_contigs = assemble_contigs(discarded_contigs, bridge_permissions, bridge_order, config, stats);
	
	std::vector<int> linked(retained_contigs.size()-1, -1);
	std::vector<std::pair<int, int> > link_overlaps(retained_contigs.size()-1);
	for (int i = 0; i < scaffolding_contigs.size(); i++) {
		std::string& scaffold = scaffolding_contigs[i].seq;
		int best_link = -1, best_link_w = 0;
		std::pair<int, int> best_link_overlap;
		for (int j = 0; j < retained_contigs.size()-1; j++) {
			if (linked[j] >= 0) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(retained_contigs[j].seq, scaffold, 1, -4, config.max_seq_error, config.min_clip_len);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(scaffold, retained_contigs[j+1].seq, 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa1.overlap && spa2.overlap && best_link_w < spa1.score+spa2.score) {
				if (reject_assembly_overlap(retained_contigs[j].seq, scaffold, spa1.overlap) || reject_assembly_overlap(scaffold, retained_contigs[j+1].seq, spa2.overlap)) continue;
				best_link = j, best_link_w = spa1.score+spa2.score;
				best_link_overlap = {spa1.overlap, spa2.overlap};
			}
		}

		if (best_link >= 0) {
			linked[best_link] = i;
			link_overlaps[best_link] = best_link_overlap;
		}
	}

	std::vector<assembly_consensus_t> scaffolded_contigs;
	assembly_sequence_t assembly{retained_contigs[0].seq, {0}};
	std::vector<const assembly_consensus_t*> path_inputs{&retained_contigs[0]};
	for (int i = 1; i < retained_contigs.size(); i++) {
		if (linked[i-1] == -1) {
			scaffolded_contigs.push_back(polish_assembly_path(path_inputs, std::move(assembly)));
			assembly = {retained_contigs[i].seq, {0}};
			path_inputs = {&retained_contigs[i]};
		} else {
			const assembly_consensus_t& bridge = scaffolding_contigs[linked[i-1]];
			auto& lo = link_overlaps[i-1];
			int bridge_offset = int(assembly.seq.length())-lo.first;
			int right_offset = bridge_offset+int(bridge.seq.length())-lo.second;
			assembly.seq += bridge.seq.substr(lo.first);
			// Retained bases supply the fallback; the right contig wins intersecting spans.
			assembly.seq.resize(right_offset);
			assembly.seq += retained_contigs[i].seq;
			assembly.offsets.push_back(bridge_offset);
			assembly.offsets.push_back(right_offset);
			path_inputs.push_back(&bridge);
			path_inputs.push_back(&retained_contigs[i]);
		}
	}
	scaffolded_contigs.push_back(polish_assembly_path(path_inputs, std::move(assembly)));

	std::vector<std::tuple<int, assembly_consensus_t, StripedSmithWaterman::Alignment> > scaffolded_w_score;
	for (assembly_consensus_t& contig : scaffolded_contigs) {
		harsh_aligner.Align(contig.seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		scaffolded_w_score.push_back({aln.sw_score, std::move(contig), aln});
	}

	retained = retain_non_overlapping_contigs(std::move(scaffolded_w_score));
	std::vector<std::string> sequences;
	consensus_contigs_alns.clear();
	for (auto& p : retained) {
		sequences.push_back(std::move(p.first.seq));
		consensus_contigs_alns.push_back(p.second);
	}

	append_initial_sequences(sequences);
	return sequences;
}

std::vector<std::string> generate_reference_guided_consensus(std::string reference, 
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, config_t& config, stats_t& stats, const hp_tail_quality_model_t& hp_tail_quality_model) {

	// Build raw base-voting and recalibrated trimming evidence for the same inputs.
	// Each quality stream follows its sequence's assembly orientation, so HP 3' tails
	// stay on the correct end.
	std::vector<std::string> seqs_lf, seqs_is, seqs_rf;
	std::vector<const uint8_t*> quals_lf, quals_is, quals_rf;
	std::vector<const uint8_t*> trim_quals_lf, trim_quals_is, trim_quals_rf;
	std::vector<std::vector<uint8_t>> quals_storage;
	std::vector<std::vector<uint8_t>> trim_quals_storage;
	size_t n_inputs = count_assembly_reads(r_cluster, l_cluster);
	quals_storage.reserve(n_inputs);
	trim_quals_storage.reserve(n_inputs);
	hp_tail_quality_table_t quality_cache;
	assembly_read_filter_t read_filter(r_cluster, l_cluster);
	for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
		if (read_filter.add_read(read.get())) {
			std::string read_seq = get_sequence(read.get());
			seqs_lf.push_back(read_seq);
			quals_storage.push_back(get_qual(read.get()));
			quals_lf.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(recalibrate_clip_read_qualities(read.get(), config, hp_tail_quality_model, quality_cache));
			trim_quals_lf.push_back(trim_quals_storage.back().data());
		}
		if (read_filter.add_mate(read.get())) {
			std::string mate_seq = get_mate_seq(read.get(), mateseqs);
			quals_storage.push_back(decode_qualities(get_mate_qual(read.get(), matequals)));
			rc(mate_seq);
			std::reverse(quals_storage.back().begin(), quals_storage.back().end());
			seqs_is.push_back(mate_seq);
			quals_is.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(quals_storage.back());
			hts_pair_pos_t hp = longest_homopolymer(mate_seq.c_str(), mate_seq.length());
			int hp_len = hp.end-hp.beg;
			if (!trim_quals_storage.back().empty() && hp_len >= MIN_REF_HP_LEN) {
				std::vector<hp_read_observation_t> observations{{mate_seq, std::move(trim_quals_storage.back()), int(hp.beg), true, 0.0, int(mate_seq.length()-hp.end)}};
				recalibrate_hp_3p_tail_qualities(observations, hp_len, mate_seq[hp.beg], hp_tail_quality_model, quality_cache);
				trim_quals_storage.back() = std::move(observations[0].quals);
			}
			trim_quals_is.push_back(trim_quals_storage.back().data());
		}
	}
	for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
		if (read_filter.add_read(read.get())) {
			std::string read_seq = get_sequence(read.get());
			seqs_rf.push_back(read_seq);
			quals_storage.push_back(get_qual(read.get()));
			quals_rf.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(recalibrate_clip_read_qualities(read.get(), config, hp_tail_quality_model, quality_cache));
			trim_quals_rf.push_back(trim_quals_storage.back().data());
		}
		if (read_filter.add_mate(read.get())) {
			std::string mate_seq = get_mate_seq(read.get(), mateseqs);
			quals_storage.push_back(decode_qualities(get_mate_qual(read.get(), matequals)));
			seqs_is.push_back(mate_seq);
			quals_is.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(quals_storage.back());
			hts_pair_pos_t hp = longest_homopolymer(mate_seq.c_str(), mate_seq.length());
			int hp_len = hp.end-hp.beg;
			if (!trim_quals_storage.back().empty() && hp_len >= MIN_REF_HP_LEN) {
				std::vector<hp_read_observation_t> observations{{mate_seq, std::move(trim_quals_storage.back()), int(hp.beg), false, 0.0, int(mate_seq.length()-hp.end)}};
				recalibrate_hp_3p_tail_qualities(observations, hp_len, mate_seq[hp.beg], hp_tail_quality_model, quality_cache);
				trim_quals_storage.back() = std::move(observations[0].quals);
			}
			trim_quals_is.push_back(trim_quals_storage.back().data());
		}
	}
	if (r_cluster->clip_consensus) {
		for (std::shared_ptr<bam1_t> read : r_cluster->clip_consensus->clipped_reads) {
			if (!read_filter.add_read(read.get(), true)) continue;
			seqs_lf.push_back(get_sequence(read.get()));
			quals_storage.push_back(get_qual(read.get()));
			quals_lf.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(recalibrate_clip_read_qualities(read.get(), config, hp_tail_quality_model, quality_cache));
			trim_quals_lf.push_back(trim_quals_storage.back().data());
		}
	}
	if (l_cluster->clip_consensus) {
		for (std::shared_ptr<bam1_t> read : l_cluster->clip_consensus->clipped_reads) {
			if (!read_filter.add_read(read.get(), true)) continue;
			seqs_rf.push_back(get_sequence(read.get()));
			quals_storage.push_back(get_qual(read.get()));
			quals_rf.push_back(quals_storage.back().data());
			trim_quals_storage.push_back(recalibrate_clip_read_qualities(read.get(), config, hp_tail_quality_model, quality_cache));
			trim_quals_rf.push_back(trim_quals_storage.back().data());
		}
	}
	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		if (!read_filter.add_read(read)) continue;
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		quals_storage.push_back(get_qual(read, true));
		std::reverse(quals_storage.back().begin(), quals_storage.back().end());
		seqs_is.push_back(read_seq);
		quals_is.push_back(quals_storage.back().data());
		trim_quals_storage.push_back(recalibrate_clip_read_qualities(read, config, hp_tail_quality_model, quality_cache));
		if (!bam_is_rev(read)) std::reverse(trim_quals_storage.back().begin(), trim_quals_storage.back().end());
		trim_quals_is.push_back(trim_quals_storage.back().data());
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		if (!read_filter.add_read(read)) continue;
		std::string read_seq = get_sequence(read, true);
		quals_storage.push_back(get_qual(read, true));
		seqs_is.push_back(read_seq);
		quals_is.push_back(quals_storage.back().data());
		trim_quals_storage.push_back(recalibrate_clip_read_qualities(read, config, hp_tail_quality_model, quality_cache));
		if (bam_is_rev(read)) std::reverse(trim_quals_storage.back().begin(), trim_quals_storage.back().end());
		trim_quals_is.push_back(trim_quals_storage.back().data());
	}

	std::vector<const uint8_t*> trim_quals = trim_quals_lf;
	trim_quals.insert(trim_quals.end(), trim_quals_is.begin(), trim_quals_is.end());
	trim_quals.insert(trim_quals.end(), trim_quals_rf.begin(), trim_quals_rf.end());
	return generate_reference_guided_consensus(reference, seqs_lf, seqs_is, seqs_rf, quals_lf, quals_is, quals_rf, aligner, harsh_aligner, consensus_contigs_alns, config, stats, true, false, &trim_quals);
}

// Re-estimate junction HP lengths from reads anchored on both sides. The sequencing-3'
// flank tolerates HP-tail errors, while the 5' flank must identify the run reliably.
void reestimate_junction_homopolymers(std::string& junction, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::vector<hts_pair_pos_t>& assembled_contig_intervals) {
	std::vector<std::pair<std::string, bool>> voters; // Sequence and whether its 3' end is on the left.
	assembly_read_filter_t read_filter(r_cluster, l_cluster);
	for (auto cluster : {r_cluster, l_cluster}) {
		bool right_cluster = cluster == r_cluster;
		if (cluster->clip_consensus) {
			for (auto read : cluster->clip_consensus->clipped_reads) {
				if (read_filter.add_read(read.get(), true)) voters.push_back({get_sequence(read.get()), bam_is_rev(read.get()) != 0});
			}
		}
		for (const auto& read : cluster->cluster->reads) {
			if (read_filter.add_read(read.get())) voters.push_back({get_sequence(read.get()), bam_is_rev(read.get()) != 0});
			if (read_filter.add_mate(read.get())) {
				std::string mate = get_mate_seq(read.get(), mateseqs);
				if (right_cluster) rc(mate);
				voters.push_back({std::move(mate), right_cluster});
			}
		}
	}
	for (auto cluster : {r_cluster, l_cluster}) {
		bool right_cluster = cluster == r_cluster;
		for (bam1_t* read : cluster->semi_mapped_reads) {
			if (!read_filter.add_read(read)) continue;
			std::string seq = get_sequence(read, true);
			if (right_cluster) rc(seq);
			voters.push_back({std::move(seq), right_cluster});
		}
	}

	std::string corrected;
	std::vector<hts_pos_t> boundary_map(junction.length()+1);
	for (int begin = 0; begin < junction.length();) {
		int end = begin+1;
		char base = junction[begin];
		while (end < junction.length() && junction[end] == base) end++;
		int length = end-begin, new_length = length;
		if (base_to_index(base) >= 0 && length >= 10 && begin >= 10 && junction.length()-end >= 10) {
			std::map<int, int> votes;
			int total_votes = 0;
			for (const auto& voter : voters) {
				const std::string& seq = voter.first;
				for (int read_begin = 0; read_begin < seq.length();) {
					int read_end = read_begin+1;
					while (read_end < seq.length() && seq[read_end] == seq[read_begin]) read_end++;
					int read_length = read_end-read_begin;
					if (seq[read_begin] == base && read_length >= std::max(1, length/2)) {
						int left_len = std::min({20, begin, read_begin});
						int right_len = std::min({20, int(junction.length())-end, int(seq.length())-read_end});
						if (left_len >= 10 && right_len >= 10) {
							int left_limit = voter.second ? 2*left_len/5 : 1, right_limit = voter.second ? 1 : 2*right_len/5;
							int left_mm = number_of_mismatches_fast(seq.data()+read_begin-left_len, junction.data()+begin-left_len, left_len, left_limit);
							int right_mm = number_of_mismatches_fast(seq.data()+read_end, junction.data()+end, right_len, right_limit);
							if (left_mm <= left_limit && right_mm <= right_limit) {
								votes[read_length]++;
								total_votes++;
								break; // At most one vote per read for this junction run.
							}
						}
					}
					read_begin = read_end;
				}
			}
			for (const auto& vote : votes) {
				if (vote.second > total_votes/2) new_length = vote.first;
			}
		}
		// Preserve boundary order and exact run endpoints; scale boundaries inside a resized HP.
		for (int i = 0; i <= length; i++) boundary_map[begin+i] = corrected.length()+hts_pos_t(i)*new_length/length;
		corrected.append(new_length, base);
		begin = end;
	}
	for (auto& interval : assembled_contig_intervals) {
		interval.beg = boundary_map[interval.beg];
		interval.end = boundary_map[interval.end];
	}
	junction = std::move(corrected);
}

void remap_reference_guided_junction(const std::string& junction_seq, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		remap_side_t& rc_remap_side, remap_side_t& lc_remap_side, std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats) {
	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	rc_remap_side.read_remap_infos.clear();
	rc_remap_side.stable_read_spans.clear();
	rc_remap_side.clip_remap_info.reset();
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		std::shared_ptr<bam1_t> read = r_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read.get());
		std::string read_qual = get_qual_ascii(read.get());
		std::string mateseq = get_mate_seq(read.get(), mateseqs);
		std::string matequal = get_mate_qual(read.get(), matequals);
		rc(mateseq);
		matequal = std::string(matequal.rbegin(), matequal.rend());
		
		aligner.Align(read_seq.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		std::pair<int, int> stable_span{aln.ref_begin, aln.ref_end};
		aligner.Align(mateseq.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		rc_remap_side.read_remap_infos.push_back(remap_info_t(aln, accepted, config.min_clip_len));
		if (accepted) rc_remap_side.stable_read_spans.push_back(stable_span);
	}
	if (r_cluster->clip_consensus) {
		harsh_aligner.Align(r_cluster->clip_consensus->sequence.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		rc_remap_side.clip_remap_info.reset(new remap_info_t(aln, accepted, config.min_clip_len));
	}

	lc_remap_side.read_remap_infos.clear();
	lc_remap_side.stable_read_spans.clear();
	lc_remap_side.clip_remap_info.reset();
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		std::shared_ptr<bam1_t> read = l_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read.get());
		std::string read_qual = get_qual_ascii(read.get());
		std::string mateseq = get_mate_seq(read.get(), mateseqs);
		std::string matequal = get_mate_qual(read.get(), matequals);

		aligner.Align(read_seq.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		std::pair<int, int> stable_span{aln.ref_begin, aln.ref_end};
		aligner.Align(mateseq.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		lc_remap_side.read_remap_infos.push_back(remap_info_t(aln, accepted, config.min_clip_len));
		if (accepted) lc_remap_side.stable_read_spans.push_back(stable_span);
	}
	if (l_cluster->clip_consensus) {
		harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), junction_seq.c_str(), junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		lc_remap_side.clip_remap_info.reset(new remap_info_t(aln, accepted, config.min_clip_len));
	}

}

std::string generate_consensus_sequences(std::string contig_name, chr_seqs_map_t& contigs, contig_map_t& contig_map, 
			std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, 
	        remap_side_t& rc_remap_side, remap_side_t& lc_remap_side, region_t& best_region, bool is_rc, 
			std::vector<hts_pair_pos_t>& assembled_contig_intervals, std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
			StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats, const hp_tail_quality_model_t& hp_tail_quality_model) {

	assembled_contig_intervals.clear();

	if (best_region.score.remap_start >= best_region.score.remap_end) return "";

	insertion_cluster_t refined_r_cluster(std::make_shared<cluster_t>());
	insertion_cluster_t refined_l_cluster(std::make_shared<cluster_t>());
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		if (rc_remap_side.read_remap_infos[i].accepted) {
			refined_r_cluster.add_stable_read(r_cluster->cluster->reads[i]);
		}
	}
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		if (lc_remap_side.read_remap_infos[i].accepted) {
			refined_l_cluster.add_stable_read(l_cluster->cluster->reads[i]);
		}
	}
	if (r_cluster->clip_consensus && rc_remap_side.clip_remap_info && rc_remap_side.clip_remap_info->accepted) {
		refined_r_cluster.add_clip_cluster(r_cluster->clip_consensus);
	}
	if (l_cluster->clip_consensus && lc_remap_side.clip_remap_info && lc_remap_side.clip_remap_info->accepted) {
		refined_l_cluster.add_clip_cluster(l_cluster->clip_consensus);
	}

	if (refined_r_cluster.empty() || refined_l_cluster.empty()) return "";
	if (refined_r_cluster.end-refined_r_cluster.start <= 0 || refined_l_cluster.end-refined_l_cluster.start <= 0) return "";

    // build reference-guide as left flanking + predicted transposed region + right flanking
	char* left_flanking = new char[refined_r_cluster.end-refined_r_cluster.start+2];
	strncpy(left_flanking, contigs.get_seq(contig_name)+refined_r_cluster.start, refined_r_cluster.end-refined_r_cluster.start+1);
	left_flanking[refined_r_cluster.end-refined_r_cluster.start+1] = '\0';
	to_uppercase(left_flanking);

	int ins_seq_start = best_region.start + best_region.score.remap_start,
		ins_seq_len = best_region.score.remap_end - best_region.score.remap_start + 1;
	char* pred_ins_seq = new char[ins_seq_len+1];
	strncpy(pred_ins_seq, contigs.get_seq(contig_map.get_name(best_region.contig_id))+ins_seq_start, ins_seq_len);
	pred_ins_seq[ins_seq_len] = '\0';
	if (is_rc) rc(pred_ins_seq);
	to_uppercase(pred_ins_seq);

	char* right_flanking = new char[refined_l_cluster.end-refined_l_cluster.start+2];
	strncpy(right_flanking, contigs.get_seq(contig_name)+refined_l_cluster.start, refined_l_cluster.end-refined_l_cluster.start+1);
	right_flanking[refined_l_cluster.end-refined_l_cluster.start+1] = '\0';
	to_uppercase(right_flanking);

	// NOTE: prefixes and suffixes may be duplicated - i.e., suffix of left flanking also appearing as prefix of predicted ins seq
	// This is a problem because it predicts an insertion where there is none. One such case is fragment LINE insertions.
	// There are small regions in LINE/L1P* that accumulate a lot of mutations compared to reference, and it may be predicted as a block-swap
	// (which is not an insertion as the length of alt allele is not changed)
	// However, because of the duplicated suffix/prefix, a small insertion is predicted
	std::string left_flanking_str = left_flanking, pred_ins_seq_str = pred_ins_seq, right_flanking_str = right_flanking;

	suffix_prefix_aln_t spa12 = aln_suffix_prefix(left_flanking_str, pred_ins_seq_str, 1, -4, 1.0, config.min_clip_len);
	int overlap12 = spa12.overlap;
	if (!overlap12) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(pred_ins_seq_str.c_str(), left_flanking_str.c_str(), left_flanking_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == left_flanking_str.length()-1) overlap12 = aln.query_end;
	}

	suffix_prefix_aln_t spa23 = aln_suffix_prefix(pred_ins_seq_str, right_flanking_str, 1, -4, 1.0, config.min_clip_len);
	int overlap23 = spa23.overlap;
	if (!overlap23) {
		StripedSmithWaterman::Alignment aln;
		StripedSmithWaterman::Filter filter;
		aligner.Align(right_flanking_str.c_str(), pred_ins_seq_str.c_str(), pred_ins_seq_str.length(), filter, &aln, 0);
		if (aln.query_begin == 0 && aln.ref_end == pred_ins_seq_str.length()-1) overlap23 = aln.query_end;
	}

	std::string full_junction_sequence = left_flanking_str;
	full_junction_sequence += pred_ins_seq_str.substr(overlap12);
	full_junction_sequence += right_flanking_str.substr(overlap23);

	delete[] left_flanking;
	delete[] pred_ins_seq;
	delete[] right_flanking;


	// assembled contigs
	std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;
	std::vector<std::string> consensus_contigs = generate_reference_guided_consensus(full_junction_sequence, r_cluster, l_cluster, mateseqs, matequals, aligner, harsh_aligner, consensus_contigs_alns, config, stats, hp_tail_quality_model);

	std::string corrected_junction_seq;

	if (consensus_contigs.empty()) {
		corrected_junction_seq = full_junction_sequence;
	} else {
		std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
		StripedSmithWaterman::Filter filter;
		for (int i = 0; i < consensus_contigs.size(); i++) {
			StripedSmithWaterman::Alignment aln;
			aligner.Align(consensus_contigs[i].c_str(), full_junction_sequence.c_str(), full_junction_sequence.length(), filter, &aln, 0);
			contigs_sorted_by_pos.push_back({consensus_contigs[i], aln});
		}
		std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
				[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
			return p1.second.ref_begin-p1.second.query_begin < p2.second.ref_begin-p2.second.query_begin;
		});
	
		// complement holes in contigs with reference
		int cur = 0;
		for (auto& p : contigs_sorted_by_pos) {
			int span_begin = std::max(0, p.second.ref_begin-p.second.query_begin);
			int span_end = std::min(int(full_junction_sequence.length())-1, p.second.ref_end+int(p.first.length())-1-p.second.query_end);
			if (span_begin > cur) corrected_junction_seq += full_junction_sequence.substr(cur, span_begin-cur);
			hts_pos_t contig_begin = corrected_junction_seq.length();
			corrected_junction_seq += p.first;
			// Half-open coordinates in the corrected junction, excluding template-filled gaps.
			assembled_contig_intervals.push_back({contig_begin, hts_pos_t(corrected_junction_seq.length())});
			cur = std::max(cur, span_end+1);
		}
		corrected_junction_seq += full_junction_sequence.substr(cur);
	}

	// Resolve HP-length errors in the completed junction before recomputing read placements.
	reestimate_junction_homopolymers(corrected_junction_seq, r_cluster, l_cluster, mateseqs, assembled_contig_intervals);

	remap_reference_guided_junction(corrected_junction_seq, r_cluster, l_cluster, rc_remap_side, lc_remap_side, mateseqs, matequals, aligner, harsh_aligner, config, stats);

	return corrected_junction_seq;
}


void update_read(bam1_t* read, region_t& chosen_region, remap_info_t& remap_info, bool region_is_rc) {
    read->core.mtid = chosen_region.original_bam_id;
    read->core.mpos = chosen_region.start + remap_info.start;
    if (region_is_rc == bam_is_rev(read)) {
        read->core.flag |= BAM_FMREVERSE; //sets flag to true
    } else {
        read->core.flag &= ~BAM_FMREVERSE; //sets flag to false
    }
    bam_aux_update_str(read, "MC", remap_info.cigar.length() + 1, remap_info.cigar.c_str());
    char ok = remap_info.accepted ? 'T' : 'F';
    bam_aux_append(read, "OK", 'A', 1, (const uint8_t*) &ok);
    bam_aux_append(read, "MS", 'i', 4, (const uint8_t*) &remap_info.score);
}

std::shared_ptr<insertion_t> detect_reference_guided_insertion_from_junction(std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::string& junction_seq,
			std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, 
			remap_side_t& ro_remap_side, remap_side_t& lo_remap_side,
			const std::vector<hts_pair_pos_t>& assembled_contig_intervals,
			StripedSmithWaterman::Aligner& aligner, stats_t& stats, config_t& config) {

	int remap_start = std::max(hts_pos_t(0), r_cluster->start-50);
	int remap_end = std::min(l_cluster->end+50, contig_len-1);
	if (remap_start >= remap_end) return NULL;

	 std::vector<std::shared_ptr<sv_t>> insertions = detect_svs_from_junction(contig_name, contig_seq, junction_seq, std::string(junction_seq.length(), 'I'),
                remap_start, remap_end, remap_start, remap_end, aligner, 0, 0, stats, config, -1, 1000, false, nullptr, &assembled_contig_intervals);

	if (insertions.empty() || insertions[0]->svtype() != "INS" || insertions[0]->svsize() < config.min_assembled_ins_size) return NULL;

	std::shared_ptr<insertion_t> insertion = std::dynamic_pointer_cast<insertion_t>(insertions[0]);

	std::shared_ptr<consensus_t> rc_consensus = r_cluster->clip_consensus && ro_remap_side.clip_remap_info && ro_remap_side.clip_remap_info->accepted ? r_cluster->clip_consensus : nullptr;
	std::shared_ptr<consensus_t> lc_consensus = l_cluster->clip_consensus && lo_remap_side.clip_remap_info && lo_remap_side.clip_remap_info->accepted ? l_cluster->clip_consensus : nullptr;
	if (!rc_consensus && std::none_of(ro_remap_side.read_remap_infos.begin(), ro_remap_side.read_remap_infos.end(), [](const remap_info_t& ri) { return ri.accepted; })) return NULL;
	if (!lc_consensus && std::none_of(lo_remap_side.read_remap_infos.begin(), lo_remap_side.read_remap_infos.end(), [](const remap_info_t& ri) { return ri.accepted; })) return NULL;

	// Get transposed sequence coverage
	std::vector<std::pair<int,int>> covered_segments;
	for (remap_info_t ri : ro_remap_side.read_remap_infos) {
		if (ri.accepted) covered_segments.push_back({ri.start, ri.end});
	}
	for (remap_info_t ri : lo_remap_side.read_remap_infos) {
		if (ri.accepted) covered_segments.push_back({ri.start, ri.end});
	}
	if (ro_remap_side.clip_remap_info && ro_remap_side.clip_remap_info->accepted) {
		covered_segments.push_back({ro_remap_side.clip_remap_info->start, ro_remap_side.clip_remap_info->end});
	}
	if (lo_remap_side.clip_remap_info && lo_remap_side.clip_remap_info->accepted) {
		covered_segments.push_back({lo_remap_side.clip_remap_info->start, lo_remap_side.clip_remap_info->end});
	}

	size_t found = junction_seq.find(insertion->ins_seq);
	int ins_seq_start = (found != std::string::npos) ? (int) found : insertion->left_anchor_aln->seq_len - insertion->mh_len;
	int ins_seq_end = ins_seq_start + insertion->ins_seq.length();

	int i = 0;
	std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.first < p2.first;});
	while (i < covered_segments.size() && covered_segments[i].second <= ins_seq_start) i++; // find left-most segment that either overlaps or is fully contained in the inserted sequence
	int prefix_cov_start = 0, prefix_cov_end = 0;
	if (i < covered_segments.size()) {
		prefix_cov_start = std::max(ins_seq_start, covered_segments[i].first);
		prefix_cov_end = covered_segments[i].second;
		while (i < covered_segments.size() && covered_segments[i].first <= prefix_cov_end) {
			prefix_cov_end = std::max(prefix_cov_end, covered_segments[i].second);
			i++;
		}
		prefix_cov_start -= ins_seq_start;
		prefix_cov_end -= ins_seq_start;
	}
	if (prefix_cov_end >= insertion->ins_seq.length()) {
		prefix_cov_end = insertion->ins_seq.length()-1;
	}
	if (prefix_cov_end == 0 || prefix_cov_start >= insertion->ins_seq.length()) {
		return NULL;
	}

	i = covered_segments.size()-1;
	std::sort(covered_segments.begin(), covered_segments.end(), [] (std::pair<int,int>& p1, std::pair<int,int>& p2) {return p1.second < p2.second;});
	while (i >= 0 && covered_segments[i].first >= ins_seq_end) i--; // find right-most segment that either overlaps or is fully contained in the inserted sequence
	int suffix_cov_start = 0, suffix_cov_end = 0;
	if (i >= 0) {
		suffix_cov_end = std::min(ins_seq_end-1, covered_segments[i].second);
		suffix_cov_start = covered_segments[i].first;
		while (i >= 0 && covered_segments[i].second >= suffix_cov_start) {
			suffix_cov_start = std::min(suffix_cov_start, covered_segments[i].first);
			i--;
		}
		suffix_cov_start -= ins_seq_start;
		if (suffix_cov_start < 0) suffix_cov_start = 0;
		suffix_cov_end -= ins_seq_start;
	}

	// Extend already-covered outer ends through adjacent stable reads, without changing
	// the mate/clip support that determines whether the insertion has an internal gap.
	std::vector<std::pair<int, int>> stable_read_spans = ro_remap_side.stable_read_spans;
	stable_read_spans.insert(stable_read_spans.end(), lo_remap_side.stable_read_spans.begin(), lo_remap_side.stable_read_spans.end());
	int covered_start = ins_seq_start+prefix_cov_start, covered_end = ins_seq_start+suffix_cov_end;
	int previous_start, previous_end;
	do {
		previous_start = covered_start;
		previous_end = covered_end;
		for (const auto& span : stable_read_spans) {
			if (span.first < covered_start && span.second >= covered_start-1) covered_start = std::max(ins_seq_start, span.first);
			if (span.second > covered_end && span.first <= covered_end+1) covered_end = std::min(ins_seq_end-1, span.second);
		}
	} while (covered_start != previous_start || covered_end != previous_end);
	prefix_cov_start = covered_start-ins_seq_start;
	suffix_cov_end = covered_end-ins_seq_start;
	prefix_cov_end = std::max(prefix_cov_end, prefix_cov_start);
	suffix_cov_start = std::min(suffix_cov_start, suffix_cov_end);

	if (prefix_cov_end < suffix_cov_start) {
		insertion->inferred_ins_seq = insertion->ins_seq;
		std::string ins_seq_prefix = insertion->ins_seq.substr(prefix_cov_start, prefix_cov_end-prefix_cov_start+1);
		std::string ins_seq_suffix = insertion->ins_seq.substr(suffix_cov_start, suffix_cov_end-suffix_cov_start+1);
		insertion->ins_seq = ins_seq_prefix + "-" + ins_seq_suffix;
	} else {
		insertion->ins_seq = insertion->ins_seq.substr(prefix_cov_start, suffix_cov_end-prefix_cov_start+1);
	}

	if (rc_consensus) insertion->rc_consensus = rc_consensus;
	if (lc_consensus) insertion->lc_consensus = lc_consensus;
	insertion->source = "REFERENCE_GUIDED_ASSEMBLY";
	insertion->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly

	return insertion;
}

std::shared_ptr<insertion_t> detect_reference_guided_assembly_insertion(std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::string& junction_seq,
			std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
			remap_side_t& ro_remap_side, remap_side_t& lo_remap_side,
			region_t& best_region, bool is_rc, const std::vector<hts_pair_pos_t>& assembled_contig_intervals,
			StripedSmithWaterman::Aligner& aligner, stats_t& stats, config_t& config) {
	return detect_reference_guided_insertion_from_junction(contig_name, contig_seq, contig_len, junction_seq, r_cluster, l_cluster, ro_remap_side, lo_remap_side, assembled_contig_intervals, aligner, stats, config);
}

#endif // GUIDED_REFERENCE_ASSEMBLE_H
