#ifndef GUIDED_REFERENCE_ASSEMBLE_H
#define GUIDED_REFERENCE_ASSEMBLE_H

#include <functional>
#include <memory>
#include <queue>

#include "../libs/ssw_cpp.h"
#include "htslib/sam.h"
#include "sw_utils.h"
#include "utils.h"
#include "sam_utils.h"
#include "dc_remapper.h"
#include "assemble.h"
#include "remapping.h"

void build_aln_guided_graph(std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& alns, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev, config_t& config) {
	std::sort(alns.begin(), alns.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& aln1, std::pair<std::string, StripedSmithWaterman::Alignment>& aln2) {
		return get_unclipped_start(aln1.second) < get_unclipped_start(aln2.second);
	});

	for (int i = 0; i < alns.size(); i++) {
		for (int j = i+1; j < alns.size() && get_unclipped_end(alns[i].second)-get_unclipped_start(alns[j].second) >= config.min_clip_len; j++) {
			int min_overlap = config.min_clip_len;
			// In principle, the suffix-prefix alignment should be the same as the overlap between the two alignments on the reference
			// This may not be the case if there are insertions or deletions in between the reads. However, in that case,
			// we expect the first read to be clipped on the right and/or the second read to be clipped on the left.
			if (!is_right_clipped(alns[i].second) && !is_left_clipped(alns[j].second)) {
				int unclipped_i_end = get_unclipped_end(alns[i].second);
				int unclipped_j_start = get_unclipped_start(alns[j].second);
				min_overlap = std::max(min_overlap, unclipped_i_end - unclipped_j_start);
			}
			suffix_prefix_aln_t spa = aln_suffix_prefix(alns[i].first, alns[j].first, 1, -4, config.max_seq_error, min_overlap);
			bool spa_i_homopolymer = is_homopolymer(alns[i].first.c_str()+alns[i].first.length()-spa.overlap, spa.overlap);
			bool spa_j_homopolymer = is_homopolymer(alns[j].first.c_str(), spa.overlap);
			if (spa.overlap && !spa_i_homopolymer && !spa_j_homopolymer) {
				out_edges[i]++;
				l_adj[i].push_back({j, spa.score, spa.overlap});
				l_adj_rev[j].push_back({i, spa.score, spa.overlap});
			}
		}
	}

	// two major differences with the regular assembly:
	// 1 - by how the graph is defined, no cycle is possible here
	// 2 - in regular assembly, we only report contigs made of at last 2 reads. Here we report even single reads
}

int calculate_score(std::vector<std::string> contigs, std::vector<std::string> read_seqs) {
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	int score = 0;
	for (std::string& read : read_seqs) {
		int best_score = 0;
		for (std::string& contig : contigs) {
			aligner.Align(read.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
			if (accept(aln, 0) && aln.sw_score > best_score) {
				best_score = aln.sw_score;
			}
		}
		score += best_score;
	}
	return score;
}



struct seq_pair_overlap_t {
	int i1, i2;
	int overlap, score;

	seq_pair_overlap_t(int i1, int i2, int overlap, int score) : i1(i1), i2(i2), overlap(overlap), score(score) {}
};

// remove contigs that are fully contained within another contig
void remove_fully_contained(std::vector<std::string>& assembled_sequences, std::string& reference, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > assembled_seqs_w_alns;
	for (std::string& assembled_sequence : assembled_sequences) {
		aligner.Align(assembled_sequence.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		assembled_seqs_w_alns.push_back({assembled_sequence, aln});
	}
	std::sort(assembled_seqs_w_alns.begin(), assembled_seqs_w_alns.end(),
			[](const std::pair<std::string, StripedSmithWaterman::Alignment>& p1, const std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		if (p1.second.ref_begin != p2.second.ref_begin) {
            return p1.second.ref_begin < p2.second.ref_begin;
        }
        return p1.second.ref_end > p2.second.ref_end;
	});
	int max_end = 0;
	assembled_sequences.clear();
	for (auto& p : assembled_seqs_w_alns) {
		if (max_end < p.second.ref_end) {
			assembled_sequences.push_back(p.first);
			max_end = p.second.ref_end;
		}
	}
}

std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > retain_non_overlapping_contigs(
		std::vector<std::tuple<int, std::string, StripedSmithWaterman::Alignment> > contigs_w_score) {
	std::sort(contigs_w_score.begin(), contigs_w_score.end(),
			[](const std::tuple<int, std::string, StripedSmithWaterman::Alignment>& p1, const std::tuple<int, std::string, StripedSmithWaterman::Alignment>& p2) {
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

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > retained;
	for (int i = 0; i < contigs_w_score.size(); i++) {
		if (!remove[i]) {
			retained.push_back({std::get<1>(contigs_w_score[i]), std::get<2>(contigs_w_score[i])});
		}
	}

	std::sort(retained.begin(), retained.end(),
			[](const std::pair<std::string, StripedSmithWaterman::Alignment>& p1, const std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	return retained;
}

std::vector<std::string> generate_reference_guided_contigs(std::string reference, std::vector<std::string>& read_seqs, 
	std::vector<std::string>& discarded_seqs, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > alns;
	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	for (std::string& seq : read_seqs) {
		harsh_aligner.Align(seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		alns.push_back({seq, aln});
	}
	build_aln_guided_graph(alns, out_edges, l_adj, l_adj_rev, config);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = alns[curr_vertex].first;
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		std::vector<int> path;
		used_reads.push_back(alns[curr_vertex].first);
		path.push_back(curr_vertex);
		while (best_edges[curr_vertex].overlap) {
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += alns[curr_vertex].first.substr(overlap);
			used_reads.push_back(alns[curr_vertex].first);
			path.push_back(curr_vertex);
		}
		used[curr_vertex] = true;

		std::string corrected_assembled_sequence = assembled_sequence;
		correct_contig(corrected_assembled_sequence, used_reads, config.max_seq_error, config.min_clip_len);
		assembled_sequences.push_back(corrected_assembled_sequence);
	}
	for (int i = 0; i < n; i++) {
		if (!used[i]) discarded_seqs.push_back(alns[i].first);
	}
	return assembled_sequences;
}

std::vector<std::string> generate_reference_guided_consensus(std::string reference, 
		std::vector<std::string>& seqs_lf, std::vector<std::string>& seqs_is, std::vector<std::string>& seqs_rf,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, config_t& config, stats_t& stats,
		bool scaffold = true) {

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

	std::vector<std::string> discarded_reads;
	std::vector<std::string> assembled_sequences = generate_reference_guided_contigs(reference, read_seqs, discarded_reads, harsh_aligner, config);

	// try and join contigs
	while (true) {
		remove_fully_contained(assembled_sequences, reference, aligner);
		std::vector<std::string> temp;
		int prev_ass_seqs = assembled_sequences.size();
		assembled_sequences = generate_reference_guided_contigs(reference, assembled_sequences, temp, harsh_aligner, config);
		assembled_sequences.insert(assembled_sequences.end(), temp.begin(), temp.end());
		if (prev_ass_seqs <= assembled_sequences.size()) break;
	}

	// score the retained contigs
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::vector<int> scores(assembled_sequences.size());
	for (std::string& read_seq : read_seqs) {
		int best_score = 0;
		std::vector<int> read_scores(assembled_sequences.size());
		for (int i = 0; i < assembled_sequences.size(); i++) {
			std::string& contig = assembled_sequences[i];
			harsh_aligner.Align(read_seq.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
			if (accept(aln, config.min_clip_len) && aln.sw_score > best_score) {
				best_score = aln.sw_score;
				read_scores[i] = aln.sw_score;
			}
		}
		for (int i = 0; i < assembled_sequences.size(); i++) {
			if (read_scores[i] == best_score) {
				scores[i] += best_score;
			}
		}
	}

	// sort by descending score
	std::vector<std::tuple<int, std::string, StripedSmithWaterman::Alignment> > contigs_w_score;
	for (int i = 0; i < assembled_sequences.size(); i++) {
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(assembled_sequences[i].c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		contigs_w_score.push_back({scores[i], assembled_sequences[i], aln});
	}
	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment>> retained = retain_non_overlapping_contigs(contigs_w_score);
	std::vector<std::string> retained_assembled_seqs;
	consensus_contigs_alns.clear();
	for (auto& p : retained) {
		retained_assembled_seqs.push_back(p.first);
		consensus_contigs_alns.push_back(p.second);
	}

	// we assembled as much as possible using guidance from the reference
	// however, reads may be misaligned for in the case of an incomplete or rearranged reference
	// we can try to "scaffold" the contigs using the reads that were rejected during the assembly

	std::vector<seq_w_pp_t> seqs_w_pp, temp1, temp2;
	for (std::string& seq : discarded_reads) {
		seqs_w_pp.push_back({seq, true, true});
	}
	std::vector<std::string> scaffolding_sequences = assemble_reads(temp1, seqs_w_pp, temp2, config, stats);

	if (!scaffold || retained_assembled_seqs.empty()) return retained_assembled_seqs;
	
	std::vector<int> linked(retained_assembled_seqs.size()-1, -1);
	std::vector<std::pair<int, int> > link_overlaps(retained_assembled_seqs.size()-1);
	for (int i = 0; i < scaffolding_sequences.size(); i++) {
		std::string& scaffold = scaffolding_sequences[i];
		int best_link = -1, best_link_w = 0;
		std::pair<int, int> best_link_overlap;
		for (int j = 0; j < retained_assembled_seqs.size()-1; j++) {
			if (linked[j] >= 0) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(retained_assembled_seqs[j], scaffold, 1, -4, config.max_seq_error, config.min_clip_len);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(scaffold, retained_assembled_seqs[j+1], 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa1.overlap && spa2.overlap && best_link_w < spa1.score+spa2.score) {
				best_link = j, best_link_w = spa1.score+spa2.score;
				best_link_overlap = {spa1.overlap, spa2.overlap};
			}
		}

		if (best_link >= 0) {
			linked[best_link] = i;
			link_overlaps[best_link] = best_link_overlap;
		}
	}

	std::vector<std::string> scaffolded_sequences;
	std::string curr_seq = retained_assembled_seqs[0];
	for (int i = 1; i < retained_assembled_seqs.size(); i++) {
		if (linked[i-1] == -1) {
			scaffolded_sequences.push_back(curr_seq);
			curr_seq = retained_assembled_seqs[i];
		} else {
			std::string link = scaffolding_sequences[linked[i-1]];
			auto& lo = link_overlaps[i-1];
			link = link.substr(lo.first, link.length()-lo.first-lo.second);
			curr_seq += link + retained_assembled_seqs[i];
		}
	}
	scaffolded_sequences.push_back(curr_seq);

	std::vector<std::tuple<int, std::string, StripedSmithWaterman::Alignment> > scaffolded_w_score;
	for (std::string& seq : scaffolded_sequences) {
		harsh_aligner.Align(seq.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		scaffolded_w_score.push_back({aln.sw_score, seq, aln});
	}

	retained = retain_non_overlapping_contigs(scaffolded_w_score);
	retained_assembled_seqs.clear();
	consensus_contigs_alns.clear();
	for (auto& p : retained) {
		retained_assembled_seqs.push_back(p.first);
		consensus_contigs_alns.push_back(p.second);
	}

	return retained_assembled_seqs;
}

std::vector<std::string> generate_reference_guided_consensus(std::string reference, 
		std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, config_t& config, stats_t& stats) {

	std::vector<std::string> seqs_lf, seqs_is, seqs_rf;
	for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
		std::string read_seq = get_sequence(read.get());
		seqs_lf.push_back(read_seq);
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		rc(mate_seq);
		seqs_is.push_back(mate_seq);
	}
	for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
		std::string read_seq = get_sequence(read.get());
		seqs_rf.push_back(read_seq);
		std::string mate_seq = get_mate_seq(read.get(), mateseqs);
		seqs_is.push_back(mate_seq);
	}
	if (r_cluster->clip_consensus) {
		seqs_lf.push_back(r_cluster->clip_consensus->sequence);
	}
	if (l_cluster->clip_consensus) {
		seqs_rf.push_back(l_cluster->clip_consensus->sequence);
	}

	return generate_reference_guided_consensus(reference, seqs_lf, seqs_is, seqs_rf, aligner, harsh_aligner, consensus_contigs_alns, config, stats);
}

std::string generate_consensus_sequences(std::string contig_name, chr_seqs_map_t& contigs, contig_map_t& contig_map, 
			std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, 
	        remap_side_t& rc_remap_side, remap_side_t& lc_remap_side, region_t& best_region, bool is_rc, 
			bool& left_bp_precise, bool& right_bp_precise, std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
			StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner, config_t& config, stats_t& stats) {

	left_bp_precise = false, right_bp_precise = false;

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
	std::vector<std::string> consensus_contigs = generate_reference_guided_consensus(full_junction_sequence, r_cluster, l_cluster, mateseqs,
			aligner, harsh_aligner, consensus_contigs_alns, config, stats);

	std::string corrected_junction_seq;

	if (consensus_contigs.empty()) {
		corrected_junction_seq = full_junction_sequence;
	} else {
		int l_bp_in_seq = refined_r_cluster.end-refined_r_cluster.start, r_bp_in_seq = full_junction_sequence.length()-(refined_l_cluster.end-refined_l_cluster.start);

		const int OFFSET = 50;
		for (StripedSmithWaterman::Alignment& aln : consensus_contigs_alns) {
			left_bp_precise |= (aln.ref_begin < l_bp_in_seq-OFFSET && aln.ref_end >= l_bp_in_seq+OFFSET);
			right_bp_precise |= (aln.ref_begin < r_bp_in_seq-OFFSET && aln.ref_end >= r_bp_in_seq+OFFSET);
		}
	
		std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
		for (int i = 0; i < consensus_contigs.size(); i++) {
			contigs_sorted_by_pos.push_back({consensus_contigs[i], consensus_contigs_alns[i]});
		}
		std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
				[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
			return p1.second.ref_begin < p2.second.ref_begin;
		});
	
		// complement holes in contigs with reference
		corrected_junction_seq = full_junction_sequence.substr(0, contigs_sorted_by_pos[0].second.ref_begin);
		corrected_junction_seq += contigs_sorted_by_pos[0].first;
		for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
			corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[i-1].second.ref_end+1, contigs_sorted_by_pos[i].second.ref_begin-contigs_sorted_by_pos[i-1].second.ref_end-1);
			corrected_junction_seq += contigs_sorted_by_pos[i].first;
		}
		corrected_junction_seq += full_junction_sequence.substr(contigs_sorted_by_pos[contigs_sorted_by_pos.size()-1].second.ref_end+1);
	}

	StripedSmithWaterman::Alignment aln;
	StripedSmithWaterman::Filter filter;
	rc_remap_side.read_remap_infos.clear();
	rc_remap_side.clip_remap_info.reset();
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
		std::shared_ptr<bam1_t> read = r_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read.get());
		std::string read_qual = get_qual_ascii(read.get());
		std::string mateseq = get_mate_seq(read.get(), mateseqs);
		std::string matequal = get_mate_qual(read.get(), matequals);
		rc(mateseq);
		matequal = std::string(matequal.rbegin(), matequal.rend());
		
		aligner.Align(read_seq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		aligner.Align(mateseq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		rc_remap_side.read_remap_infos.push_back(remap_info_t(aln, accepted, config.min_clip_len));
	}
	if (r_cluster->clip_consensus) {
		harsh_aligner.Align(r_cluster->clip_consensus->sequence.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		rc_remap_side.clip_remap_info.reset(new remap_info_t(aln, accepted, config.min_clip_len));
	}

	lc_remap_side.read_remap_infos.clear();
	lc_remap_side.clip_remap_info.reset();
	for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
		std::shared_ptr<bam1_t> read = l_cluster->cluster->reads[i];
		std::string read_seq = get_sequence(read.get());
		std::string read_qual = get_qual_ascii(read.get());
		std::string mateseq = get_mate_seq(read.get(), mateseqs);
		std::string matequal = get_mate_qual(read.get(), matequals);

		aligner.Align(read_seq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error, read_qual, stats.min_avg_base_qual);
		aligner.Align(mateseq.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		accepted &= accept(aln, config.min_clip_len, config.max_seq_error, matequal, stats.min_avg_base_qual);
		lc_remap_side.read_remap_infos.push_back(remap_info_t(aln, accepted, config.min_clip_len));
	}
	if (l_cluster->clip_consensus) {
		harsh_aligner.Align(l_cluster->clip_consensus->sequence.c_str(), corrected_junction_seq.c_str(), corrected_junction_seq.length(), filter, &aln, 0);
		bool accepted = accept(aln, config.min_clip_len, config.max_seq_error);
		lc_remap_side.clip_remap_info.reset(new remap_info_t(aln, accepted, config.min_clip_len));
	}

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

std::shared_ptr<insertion_t> detect_reference_guided_assembly_insertion(std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::string& junction_seq,
			std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, 
			remap_side_t& ro_remap_side, remap_side_t& lo_remap_side,
			region_t& best_region, bool is_rc, bool left_bp_precise, bool right_bp_precise, 
			StripedSmithWaterman::Aligner& aligner, stats_t& stats, config_t& config) {

	int remap_start = std::max(hts_pos_t(0), r_cluster->start-50);
	int remap_end = std::min(l_cluster->end+50, contig_len-1);
	if (remap_start >= remap_end) return NULL;

	 std::vector<std::shared_ptr<sv_t>> insertions = detect_svs_from_junction(contig_name, contig_seq, junction_seq, std::string(junction_seq.length(), 'I'),
                remap_start, remap_end, remap_start, remap_end, aligner, 0, 0, stats, config);

	if (insertions.empty() || insertions[0]->svtype() != "INS" || insertions[0]->svsize() < config.min_assembled_ins_size) return NULL;

	std::shared_ptr<insertion_t> insertion = std::dynamic_pointer_cast<insertion_t>(insertions[0]);

	std::shared_ptr<insertion_cluster_t> refined_r_cluster = std::make_shared<insertion_cluster_t>();
    for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
        std::shared_ptr<bam1_t> r = r_cluster->cluster->reads[i];
        if (ro_remap_side.read_remap_infos[i].accepted) {
            refined_r_cluster->add_stable_read(r);
        }
    }
	if (r_cluster->clip_consensus && ro_remap_side.clip_remap_info && ro_remap_side.clip_remap_info->accepted) {
        refined_r_cluster->add_clip_cluster(r_cluster->clip_consensus);
    }

	std::shared_ptr<insertion_cluster_t> refined_l_cluster = std::make_shared<insertion_cluster_t>();
    for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
        std::shared_ptr<bam1_t> r = l_cluster->cluster->reads[i];
        if (lo_remap_side.read_remap_infos[i].accepted) {
            refined_l_cluster->add_stable_read(r);
        }
    }
    if (l_cluster->clip_consensus && lo_remap_side.clip_remap_info && lo_remap_side.clip_remap_info->accepted) {
        refined_l_cluster->add_clip_cluster(l_cluster->clip_consensus);
    }

	if (refined_r_cluster->empty() || refined_l_cluster->empty()) {
		return NULL;
	}

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

	int ins_seq_start = insertion->left_anchor_aln->seq_len - insertion->mh_len;
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

	if (prefix_cov_end < suffix_cov_start) {
		insertion->inferred_ins_seq = insertion->ins_seq;
		std::string ins_seq_prefix = insertion->ins_seq.substr(prefix_cov_start, prefix_cov_end-prefix_cov_start+1);
		std::string ins_seq_suffix = insertion->ins_seq.substr(suffix_cov_start, suffix_cov_end-suffix_cov_start+1);
		insertion->ins_seq = ins_seq_prefix + "-" + ins_seq_suffix;
	} else {
		insertion->ins_seq = insertion->ins_seq.substr(prefix_cov_start, suffix_cov_end-prefix_cov_start+1);
	}

	if (refined_r_cluster->clip_consensus) insertion->rc_consensus = refined_r_cluster->clip_consensus;
	if (refined_l_cluster->clip_consensus) insertion->lc_consensus = refined_l_cluster->clip_consensus;
	insertion->imprecise = !left_bp_precise || !right_bp_precise;
	insertion->source = "REFERENCE_GUIDED_ASSEMBLY";
	insertion->mh_len = 0; // current value is just a temporary approximation, reset it. genotype will calculate it correctly

	return insertion;
}

#endif // GUIDED_REFERENCE_ASSEMBLE_H
