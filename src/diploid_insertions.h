#ifndef DIPLOID_INSERTIONS_H_
#define DIPLOID_INSERTIONS_H_

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "dc_remapper.h"
#include "../libs/ssw_selected.h"

struct insertion_read_evidence_t {
	int read_end; // BAM read-end flags from get_read_key.
	std::string seq;
	std::vector<uint8_t> quals;
	bool reverse; // Stored sequence is reverse-complemented relative to the sequenced read.
};

struct insertion_fragment_evidence_t {
	std::string qname;
	std::vector<insertion_read_evidence_t> reads;
};

std::vector<insertion_fragment_evidence_t> collect_insertion_read_evidence(std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals) {

	std::vector<insertion_fragment_evidence_t> fragments;
	std::unordered_map<std::string, int> fragment_indices;
	assembly_read_filter_t read_filter(r_cluster, l_cluster);
	auto add_read = [&](const std::pair<std::string, int>& key, std::string seq, std::vector<uint8_t> quals, bool reverse) {
		auto inserted = fragment_indices.emplace(key.first, fragments.size());
		if (inserted.second) fragments.push_back({key.first, {}});
		fragments[inserted.first->second].reads.push_back({key.second, std::move(seq), std::move(quals), reverse});
	};
	for (auto cluster : {r_cluster, l_cluster}) {
		bool right_cluster = cluster == r_cluster;
		if (cluster->clip_consensus) {
			for (auto read : cluster->clip_consensus->clipped_reads) {
				if (read_filter.add_read(read.get(), true)) add_read(get_read_key(read.get()), get_sequence(read.get()), get_qual(read.get()), bam_is_rev(read.get()) != 0);
			}
		}
		for (auto read : cluster->cluster->reads) {
			if (read_filter.add_read(read.get())) add_read(get_read_key(read.get()), get_sequence(read.get()), get_qual(read.get()), bam_is_rev(read.get()) != 0);
			if (read_filter.add_mate(read.get())) {
				std::string seq = get_mate_seq(read.get(), mateseqs);
				std::vector<uint8_t> quals = decode_qualities(get_mate_qual(read.get(), matequals));
				if (right_cluster) {
					rc(seq);
					std::reverse(quals.begin(), quals.end());
				}
				add_read(get_read_key(read.get(), true), std::move(seq), std::move(quals), right_cluster);
			}
		}
	}
	for (auto cluster : {r_cluster, l_cluster}) {
		bool right_cluster = cluster == r_cluster;
		for (bam1_t* read : cluster->semi_mapped_reads) {
			if (!read_filter.add_read(read)) continue;
			std::string seq = get_sequence(read, true);
			std::vector<uint8_t> quals = get_qual(read, true);
			if (right_cluster) {
				rc(seq);
				std::reverse(quals.begin(), quals.end());
			}
			add_read(get_read_key(read), std::move(seq), std::move(quals), right_cluster);
		}
	}
	return fragments;
}

struct junction_edit_t {
	int begin, end; // Zero-based, half-open; equal for an insertion.
	std::string alt;
	bool correction; // Application priority; not part of the edit identity.

	bool operator<(const junction_edit_t& other) const { return std::tie(begin, end, alt) < std::tie(other.begin, other.end, other.alt); }
	bool operator==(const junction_edit_t& other) const { return begin == other.begin && end == other.end && alt == other.alt; }
};

struct junction_read_alignment_t {
	const insertion_read_evidence_t* read;
	int ref_begin, ref_end; // Zero-based, half-open.
	std::vector<int> query_positions; // Relative to ref_begin; -1 for a deleted reference base.
	std::vector<junction_edit_t> differences;
	int score, query_begin, query_end;
	bool clean;
};

inline bool junction_canonical_base(char base) {
	return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

void left_normalize_junction_indel(const std::string& junction, junction_edit_t& edit, int limit) {
	if (edit.begin == edit.end) {
		while (edit.begin > limit && junction_canonical_base(junction[edit.begin-1]) && junction[edit.begin-1] == edit.alt.back()) {
			std::rotate(edit.alt.begin(), edit.alt.end()-1, edit.alt.end());
			edit.begin--;
			edit.end--;
		}
	} else {
		while (edit.begin > limit && junction_canonical_base(junction[edit.begin-1]) && junction[edit.begin-1] == junction[edit.end-1]) {
			edit.begin--;
			edit.end--;
		}
	}
}

bool align_insertion_junction_read(const std::string& junction, const insertion_read_evidence_t& read, StripedSmithWaterman::Aligner& aligner, junction_read_alignment_t& result) {
	if (read.seq.length() < 50 || read.quals.size() != read.seq.length()) return false;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	// Exclude shorter prefixes of the best placement from the runner-up score.
	if (!aligner.Align(read.seq.c_str(), junction.c_str(), junction.length(), filter, &alignment, std::max(15, int(read.seq.length()/2)))) return false;
	int aligned_bases = 0;
	for (uint32_t encoded_op : alignment.cigar) {
		if (bam_cigar_type(bam_cigar_op(encoded_op)) == 3) aligned_bases += bam_cigar_oplen(encoded_op);
	}
	if (aligned_bases < 50) return false;
	if (alignment.sw_score_next_best*10 > alignment.sw_score*9) return false;

	StripedSmithWaterman::Filter score_filter(false, false, 0, 0);
	StripedSmithWaterman::Alignment alternative;
	int suffix_begin = alignment.ref_begin+1;
	if (aligner.Align(read.seq.c_str(), junction.c_str()+suffix_begin, junction.length()-suffix_begin, score_filter, &alternative, 0) && alternative.sw_score >= alignment.sw_score) return false;

	result = {&read, alignment.ref_begin, alignment.ref_end+1, std::vector<int>(alignment.ref_end-alignment.ref_begin+1, -1), {}, alignment.sw_score, alignment.query_begin, alignment.query_end+1, true};
	if (alignment.ref_begin > 0) {
		for (int pos = 0; pos < result.query_begin; pos++) if (read.quals[pos] >= 20) result.clean = false;
	}
	if (result.ref_end < junction.length()) {
		for (int pos = result.query_end; pos < read.seq.length(); pos++) if (read.quals[pos] >= 20) result.clean = false;
	}
	int ref_pos = alignment.ref_begin, query_pos = alignment.query_begin, previous_difference_end = 0;
	auto high_quality_base = [&](int pos) { return read.quals[pos] >= 20 && junction_canonical_base(read.seq[pos]); };
	for (int i = 0; i < alignment.cigar.size(); i++) {
		uint32_t encoded_op = alignment.cigar[i];
		int op = bam_cigar_op(encoded_op), len = bam_cigar_oplen(encoded_op);
		if (op == BAM_CSOFT_CLIP) continue;
		if (bam_cigar_type(op) == 3) {
			for (int j = 0; j < len; j++) {
				result.query_positions[ref_pos+j-result.ref_begin] = query_pos+j;
				if (read.quals[query_pos+j] >= 20 && (read.seq[query_pos+j] != junction[ref_pos+j] || !junction_canonical_base(read.seq[query_pos+j]) || !junction_canonical_base(junction[ref_pos+j]))) result.clean = false;
				if (read.seq[query_pos+j] == junction[ref_pos+j]) continue;
				if (high_quality_base(query_pos+j) && junction_canonical_base(junction[ref_pos+j])) result.differences.push_back({ref_pos+j, ref_pos+j+1, read.seq.substr(query_pos+j, 1)});
				previous_difference_end = ref_pos+j+1;
			}
		} else if (op == BAM_CINS || op == BAM_CDEL) {
			result.clean = false;
			bool anchored = i > 0 && i+1 < alignment.cigar.size() && bam_cigar_type(bam_cigar_op(alignment.cigar[i-1])) == 3 && bam_cigar_type(bam_cigar_op(alignment.cigar[i+1])) == 3;
			int query_end = query_pos+(op == BAM_CINS ? len : 0);
			bool high_quality = anchored && high_quality_base(query_pos-1) && high_quality_base(query_end);
			if (anchored) high_quality &= junction_canonical_base(junction[ref_pos-1]) && junction_canonical_base(junction[ref_pos+(op == BAM_CDEL ? len : 0)]);
			for (int j = 0; j < len && high_quality; j++) {
				if (op == BAM_CINS) high_quality = high_quality_base(query_pos+j);
				else high_quality = junction_canonical_base(junction[ref_pos+j]);
			}
			if (high_quality) {
				junction_edit_t edit{ref_pos, ref_pos+(op == BAM_CDEL ? len : 0), op == BAM_CINS ? read.seq.substr(query_pos, len) : ""};
				left_normalize_junction_indel(junction, edit, previous_difference_end);
				result.differences.push_back(std::move(edit));
			}
			previous_difference_end = ref_pos+(op == BAM_CDEL ? len : 0);
		}
		if (bam_cigar_type(op)&1) query_pos += len;
		if (bam_cigar_type(op)&2) ref_pos += len;
	}
	return true;
}

struct insertion_junction_difference_t {
	junction_edit_t edit;
	std::vector<int> supporting_fragments, informative_fragments;
	bool correction;
};

struct insertion_junction_differences_t {
	std::vector<insertion_junction_difference_t> differences;
	bool too_many_differences = false;
};

// -1: no high-quality observation; 0: another allele; 1: this difference.
int junction_difference_read_state(const std::string& junction, const junction_edit_t& edit, const junction_read_alignment_t& alignment) {
	if (std::find(alignment.differences.begin(), alignment.differences.end(), edit) != alignment.differences.end()) return 1;
	for (const auto& other : alignment.differences) {
		if ((edit.begin < other.end && other.begin < edit.end) || (edit.begin == edit.end && other.begin == other.end && edit.begin == other.begin)) return 0;
	}
	bool snp = edit.end-edit.begin == 1 && edit.alt.length() == 1;
	int left = snp ? edit.begin : edit.begin-1, right = snp ? edit.begin : edit.end;
	// An indel can move through a repeat; absence is informative only across all equivalent placements.
	if (edit.begin == edit.end) {
		while (right < junction.length() && junction_canonical_base(junction[right]) && junction[right] == edit.alt[(right-edit.begin)%edit.alt.length()]) right++;
	} else if (edit.alt.empty()) {
		int begin = edit.begin;
		while (right < junction.length() && junction_canonical_base(junction[right]) && junction[begin] == junction[right]) {
			begin++;
			right++;
		}
	}
	if (left < alignment.ref_begin || right >= alignment.ref_end) return -1;
	if (!junction_canonical_base(junction[left]) || !junction_canonical_base(junction[right])) return -1;
	int query_begin = alignment.query_positions[left-alignment.ref_begin], query_end = alignment.query_positions[right-alignment.ref_begin];
	if (query_begin < 0 || query_end < 0) return -1;
	for (int pos = query_begin; pos <= query_end; pos++) {
		if (alignment.read->quals[pos] < 20 || !junction_canonical_base(alignment.read->seq[pos])) return -1;
	}
	return 0;
}

insertion_junction_differences_t identify_insertion_junction_differences(const std::string& junction, const std::vector<insertion_fragment_evidence_t>& fragments, StripedSmithWaterman::Aligner& aligner, int min_support_percent = 20) {
	std::vector<std::vector<junction_read_alignment_t>> alignments(fragments.size());
	std::set<junction_edit_t> candidates;
	for (int i = 0; i < fragments.size(); i++) {
		for (const auto& read : fragments[i].reads) {
			junction_read_alignment_t alignment;
			if (!align_insertion_junction_read(junction, read, aligner, alignment)) continue;
			candidates.insert(alignment.differences.begin(), alignment.differences.end());
			alignments[i].push_back(std::move(alignment));
		}
	}

	insertion_junction_differences_t result;
	for (const auto& edit : candidates) {
		insertion_junction_difference_t difference;
		difference.edit = edit;
		for (int i = 0; i < alignments.size(); i++) {
			int state = -1;
			bool conflicting = false;
			for (const auto& alignment : alignments[i]) {
				int read_state = junction_difference_read_state(junction, edit, alignment);
				if (read_state < 0) continue;
				if (state >= 0 && state != read_state) conflicting = true;
				state = read_state;
			}
			if (conflicting || state < 0) continue;
			difference.informative_fragments.push_back(i);
			if (state == 1) difference.supporting_fragments.push_back(i);
		}
		int support = difference.supporting_fragments.size(), informative = difference.informative_fragments.size();
		if (support < 2 || support*100 < min_support_percent*informative) continue;
		difference.correction = support*10 >= informative*9;
		result.differences.push_back(std::move(difference));
		if (result.differences.size() > 12) {
			result.differences.clear();
			result.too_many_differences = true;
			return result;
		}
	}
	return result;
}

struct insertion_phase_link_t {
	int first, second;
	bool opposite;
	int support, opposition;
};

std::vector<insertion_phase_link_t> phase_insertion_differences(const insertion_junction_differences_t& differences) {
	std::vector<insertion_phase_link_t> links;
	for (int i = 0; i < differences.differences.size(); i++) {
		const auto& first = differences.differences[i];
		if (first.correction) continue;
		for (int j = i+1; j < differences.differences.size(); j++) {
			const auto& second = differences.differences[j];
			if (second.correction) continue;
			std::vector<int> shared;
			std::set_intersection(first.informative_fragments.begin(), first.informative_fragments.end(), second.informative_fragments.begin(), second.informative_fragments.end(), std::back_inserter(shared));
			int same = 0, opposite = 0;
			for (int fragment : shared) {
				bool first_alt = std::binary_search(first.supporting_fragments.begin(), first.supporting_fragments.end(), fragment);
				bool second_alt = std::binary_search(second.supporting_fragments.begin(), second.supporting_fragments.end(), fragment);
				if (first_alt == second_alt) same++;
				else opposite++;
			}
			if (same > 0 && same >= 4*opposite) links.push_back({i, j, false, same, opposite});
			else if (opposite > 0 && opposite >= 4*same) links.push_back({i, j, true, opposite, same});
		}
	}
	std::sort(links.begin(), links.end(), [](const insertion_phase_link_t& first, const insertion_phase_link_t& second) {
		int first_margin = first.support-first.opposition, second_margin = second.support-second.opposition;
		if (first_margin != second_margin) return first_margin > second_margin;
		if (first.support != second.support) return first.support > second.support;
		return std::tie(first.first, first.second) < std::tie(second.first, second.second);
	});
	std::vector<int> parent(differences.differences.size()), size(parent.size(), 1), parity(parent.size());
	for (int i = 0; i < parent.size(); i++) parent[i] = i;
	auto find_root = [&](int site) {
		int phase = 0;
		while (parent[site] != site) {
			phase ^= parity[site];
			site = parent[site];
		}
		return std::make_pair(site, phase);
	};
	std::vector<insertion_phase_link_t> forest;
	for (const auto& link : links) {
		auto first = find_root(link.first), second = find_root(link.second);
		if (first.first == second.first) continue; // Existing stronger links already fix this relationship.
		int root = first.first, child = second.first;
		if (size[root] < size[child]) std::swap(root, child);
		parent[child] = root;
		parity[child] = first.second^second.second^link.opposite;
		size[root] += size[child];
		forest.push_back(link);
	}
	return forest;
}

struct reconstructed_insertion_junction_t {
	std::string seq;
	std::vector<hts_pos_t> before, after; // Original boundaries before/after inserted bases.
	std::vector<hts_pair_pos_t> assembled_contig_intervals;
	bool haplotype;
};

bool build_reconstructed_insertion_junction(const std::string& primary, std::vector<junction_edit_t> edits, const std::vector<hts_pair_pos_t>& intervals, reconstructed_insertion_junction_t& result) {
	std::sort(edits.begin(), edits.end(), [](const junction_edit_t& first, const junction_edit_t& second) {
		if (first.begin != second.begin) return first.begin < second.begin;
		if (first.correction != second.correction) return first.correction;
		return first < second;
	});
	result.seq.clear();
	result.before.assign(primary.length()+1, 0);
	result.after.assign(primary.length()+1, 0);
	result.assembled_contig_intervals.clear();
	result.haplotype = false;
	int pos = 0, previous_insertion = -1;
	for (const auto& edit : edits) {
		if (edit.begin < pos || (edit.begin == edit.end && edit.begin == previous_insertion)) continue;
		while (pos < edit.begin) {
			result.seq += primary[pos++];
			result.before[pos] = result.after[pos] = result.seq.length();
		}
		hts_pos_t output_begin = result.seq.length();
		result.seq += edit.alt;
		if (edit.begin == edit.end) {
			result.after[pos] = result.seq.length();
			previous_insertion = pos;
		} else {
			for (int boundary = edit.begin+1; boundary <= edit.end; boundary++) result.before[boundary] = result.after[boundary] = output_begin+hts_pos_t(boundary-edit.begin)*edit.alt.length()/(edit.end-edit.begin);
			pos = edit.end;
		}
	}
	while (pos < primary.length()) {
		result.seq += primary[pos++];
		result.before[pos] = result.after[pos] = result.seq.length();
	}
	// Insertions on outer contig boundaries do not extend assembled coverage.
	for (const auto& interval : intervals) {
		hts_pos_t begin = result.after[interval.beg], end = result.before[interval.end];
		if (begin < end) result.assembled_contig_intervals.push_back({begin, end});
	}
	return true;
}

std::vector<reconstructed_insertion_junction_t> reconstruct_insertion_junctions(const std::string& primary, const insertion_junction_differences_t& differences, const std::vector<hts_pair_pos_t>& intervals = {}) {
	if (differences.too_many_differences) return {};
	std::vector<insertion_phase_link_t> links = phase_insertion_differences(differences);
	int corrections = 0;
	std::vector<int> variables;
	for (int i = 0; i < differences.differences.size(); i++) {
		if (differences.differences[i].correction) corrections |= 1<<i;
		else variables.push_back(i);
	}
	std::vector<int> masks{corrections};
	if (!variables.empty()) {
		std::vector<int> phase(differences.differences.size(), -1);
		for (int seed : variables) {
			if (phase[seed] >= 0) continue;
			phase[seed] = 0; // Retain the original junction's phase in each disconnected group.
			for (int pass = 0; pass < variables.size(); pass++) {
				for (const auto& link : links) {
					if (phase[link.first] >= 0) phase[link.second] = phase[link.first]^link.opposite;
					else if (phase[link.second] >= 0) phase[link.first] = phase[link.second]^link.opposite;
				}
			}
		}
		int original_phase = corrections, complement = corrections;
		for (int site : variables) {
			if (phase[site]) original_phase |= 1<<site;
			else complement |= 1<<site;
		}
		masks.push_back(original_phase);
		masks.push_back(complement);
	}
	std::vector<reconstructed_insertion_junction_t> junctions;
	std::set<std::string> sequences;
	for (int mask : masks) {
		std::vector<junction_edit_t> edits;
		for (int i = 0; i < differences.differences.size(); i++) {
			if (mask&(1<<i)) {
				junction_edit_t edit = differences.differences[i].edit;
				edit.correction = differences.differences[i].correction;
				edits.push_back(std::move(edit));
			}
		}
		reconstructed_insertion_junction_t junction;
		if (!build_reconstructed_insertion_junction(primary, std::move(edits), intervals, junction)) continue;
		junction.haplotype = (mask&~corrections) != 0;
		if (junction.seq != primary && sequences.insert(junction.seq).second) junctions.push_back(std::move(junction));
	}
	return junctions;
}

bool polish_insertion_junction(const std::string& primary, const std::vector<insertion_fragment_evidence_t>& fragments, reconstructed_insertion_junction_t& alternative, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Filter score_filter(false, false, 0, 0);
	for (int round = 0; round < 2; round++) {
		std::vector<insertion_fragment_evidence_t> selected;
		for (const auto& fragment : fragments) {
			insertion_fragment_evidence_t supporting{fragment.qname, {}};
			for (const auto& read : fragment.reads) {
				junction_read_alignment_t alignment;
				if (!align_insertion_junction_read(alternative.seq, read, aligner, alignment)) continue;
				StripedSmithWaterman::Alignment primary_alignment;
				if (aligner.Align(read.seq.c_str(), primary.c_str(), primary.length(), score_filter, &primary_alignment, 0) && alignment.score > primary_alignment.sw_score) supporting.reads.push_back(read);
			}
			if (!supporting.reads.empty()) selected.push_back(std::move(supporting));
		}
		insertion_junction_differences_t residual = identify_insertion_junction_differences(alternative.seq, selected, aligner, 75);
		if (residual.too_many_differences) return false;
		std::vector<junction_edit_t> edits;
		for (const auto& difference : residual.differences) {
			junction_edit_t edit = difference.edit;
			edit.correction = difference.correction;
			edits.push_back(std::move(edit));
		}
		if (edits.empty()) break;
		reconstructed_insertion_junction_t polished;
		if (!build_reconstructed_insertion_junction(alternative.seq, std::move(edits), alternative.assembled_contig_intervals, polished)) return false;
		for (hts_pos_t& boundary : alternative.before) boundary = polished.before[boundary];
		for (hts_pos_t& boundary : alternative.after) boundary = polished.after[boundary];
		alternative.seq = std::move(polished.seq);
		alternative.assembled_contig_intervals = std::move(polished.assembled_contig_intervals);
	}

	int support = 0, clean = 0;
	for (const auto& fragment : fragments) {
		for (const auto& read : fragment.reads) {
			junction_read_alignment_t alignment;
			if (!align_insertion_junction_read(alternative.seq, read, aligner, alignment)) continue;
			StripedSmithWaterman::Alignment primary_alignment;
			if (!aligner.Align(read.seq.c_str(), primary.c_str(), primary.length(), score_filter, &primary_alignment, 0) || alignment.score <= primary_alignment.sw_score) continue;
			support++;
			clean += alignment.clean;
		}
	}
	return clean >= 3 && (!alternative.haplotype || clean*5 >= support*2);
}

std::vector<reconstructed_insertion_junction_t> reconstruct_polished_insertion_junctions(const std::string& primary, const std::vector<insertion_fragment_evidence_t>& fragments, StripedSmithWaterman::Aligner& aligner, const std::vector<hts_pair_pos_t>& intervals = {}) {
	insertion_junction_differences_t differences = identify_insertion_junction_differences(primary, fragments, aligner);
	std::vector<reconstructed_insertion_junction_t> alternatives = reconstruct_insertion_junctions(primary, differences, intervals), polished;
	std::set<std::string> sequences;
	for (auto& alternative : alternatives) {
		if (polish_insertion_junction(primary, fragments, alternative, aligner) && alternative.seq != primary && sequences.insert(alternative.seq).second) polished.push_back(std::move(alternative));
	}
	return polished;
}

bool equivalent_insertion_alleles(const sv_t& first, const sv_t& second, const char* reference) {
	if (first.chr != second.chr) return false;
	hts_pos_t begin = std::min(first.start, second.start)+1, end = std::max(first.end, second.end)+1;
	auto allele = [&](const sv_t& insertion) {
		return std::string(reference+begin, insertion.start+1-begin)+insertion.ins_seq+std::string(reference+insertion.end+1, end-insertion.end-1);
	};
	return allele(first) == allele(second);
}

void append_distinct_insertion_candidate(std::vector<std::shared_ptr<sv_t>>& candidates, std::shared_ptr<sv_t> candidate, const char* reference) {
	for (const auto& existing : candidates) {
		if (equivalent_insertion_alleles(*existing, *candidate, reference)) return;
	}
	candidates.push_back(std::move(candidate));
}

#endif /* DIPLOID_INSERTIONS_H_ */
