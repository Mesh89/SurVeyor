#ifndef READ_CORRECTION_H_
#define READ_CORRECTION_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace read_correction {

const int K = 15;
const int MIN_ANCHORS = 2;
const int MISMATCH_HQ = 4;
const int MISMATCH_LQ = 2;
const double OVERLAP_HQ_ERR = 0.01; // Strict error < threshold: q >= 21.
const int MIN_OVERLAP = 30;
const int REPEAT_MARGIN = 10;
const int DEFAULT_QUAL = 16;
const double SAME_LOCUS_PRIOR = 0.9;
const double LOCUS_HQ_ERR = 0.005; // Strict error < threshold: q >= 24.
const int LAMBDA_HQ = 8;
const double POLY_PRIOR = 0.001;
const int POLY_GRID = 20;
const double FLIP_THRESHOLD = 0.99;

inline double base_error(uint8_t qual) {
	static const std::array<double, 256> errors = []() {
		std::array<double, 256> table;
		for (int q = 0; q < 256; q++) table[q] = std::min(0.5, std::max(1e-4, std::pow(10.0, -q/10.0)));
		return table;
	}();
	return errors[qual];
}

struct correction_input_t {
	std::string seq;            // assembly orientation; any characters, only A/C/G/T are informative
	std::vector<uint8_t> quals; // same length as seq; empty (or another length) = DEFAULT_QUAL for every base
	bool target;               // true = may be corrected
};
struct correction_flip_t { int input, pos; char new_base; double posterior; };

namespace detail {

inline int base_index(char base) {
	switch (base) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return -1;
	}
}

inline double input_error(const correction_input_t& input, int pos) {
	return base_error(input.quals.size() == input.seq.size() ? input.quals[pos] : DEFAULT_QUAL);
}

struct kmer_entry_t {
	uint64_t kmer;
	int input, pos;
};

struct kmer_range_t { int begin, end; };

const uint64_t EMPTY_KMER = ~0ULL;

struct target_kmer_index_t {
	std::vector<kmer_entry_t> entries;
	std::vector<uint64_t> table;
	std::vector<kmer_range_t> ranges;

	explicit target_kmer_index_t(const std::vector<correction_input_t>& inputs) {
		const uint64_t mask = (uint64_t(1) << (2*K))-1;
		for (int i = 0; i < int(inputs.size()); i++) {
			if (!inputs[i].target) continue;
			const std::string& seq = inputs[i].seq;
			uint64_t kmer = 0;
			int valid = 0;
			for (int p = 0; p < int(seq.length()); p++) {
				int base = base_index(seq[p]);
				if (base < 0) {
					kmer = 0;
					valid = 0;
					continue;
				}
				kmer = ((kmer << 2) | base) & mask;
				if (++valid >= K) entries.push_back({kmer, i, p+1-K});
			}
		}
		std::sort(entries.begin(), entries.end(), [](const kmer_entry_t& a, const kmer_entry_t& b) {
			if (a.kmer != b.kmer) return a.kmer < b.kmer;
			if (a.input != b.input) return a.input < b.input;
			return a.pos < b.pos;
		});
		size_t capacity = 1;
		while (capacity < 2*entries.size()) capacity *= 2;
		table.assign(capacity, EMPTY_KMER);
		ranges.resize(capacity);
		for (int begin = 0; begin < int(entries.size());) {
			int end = begin+1;
			while (end < int(entries.size()) && entries[end].kmer == entries[begin].kmer) end++;
			uint64_t kmer = entries[begin].kmer;
			size_t slot = ((kmer*0x9E3779B97F4A7C15ULL) >> 20) & (capacity-1);
			while (table[slot] != EMPTY_KMER) slot = (slot+1) & (capacity-1);
			table[slot] = kmer;
			ranges[slot] = {begin, end};
			begin = end;
		}
	}

	kmer_range_t find(uint64_t kmer) const {
		size_t slot = ((kmer*0x9E3779B97F4A7C15ULL) >> 20) & (table.size()-1);
		while (table[slot] != EMPTY_KMER) {
			if (table[slot] == kmer) return ranges[slot];
			slot = (slot+1) & (table.size()-1);
		}
		return {0, 0};
	}
};

struct anchored_offset_t { int input, offset; };

inline std::vector<std::vector<anchored_offset_t> > find_target_offsets_sorted(const std::vector<correction_input_t>& inputs, const target_kmer_index_t& target_index) {
	std::vector<std::vector<anchored_offset_t> > offsets(inputs.size());
	std::vector<std::pair<int, int> > hits;
	const uint64_t mask = (uint64_t(1) << (2*K))-1;
	for (int j = 0; j < int(inputs.size()); j++) {
		hits.clear();
		const std::string& seq = inputs[j].seq;
		uint64_t kmer = 0;
		int valid = 0;
		for (int p = 0; p < int(seq.length()); p++) {
			int base = base_index(seq[p]);
			if (base < 0) {
				kmer = 0;
				valid = 0;
				continue;
			}
			kmer = ((kmer << 2) | base) & mask;
			if (++valid < K) continue;
			int pj = p+1-K;
			kmer_range_t range = target_index.find(kmer);
			for (int hit = range.begin; hit < range.end; hit++) {
				const kmer_entry_t& entry = target_index.entries[hit];
				if (entry.input != j) hits.push_back({entry.input, entry.pos-pj});
			}
		}
		std::sort(hits.begin(), hits.end());
		for (int begin = 0; begin < int(hits.size());) {
			int end = begin+1;
			while (end < int(hits.size()) && hits[end] == hits[begin]) end++;
			if (end-begin >= MIN_ANCHORS) offsets[hits[begin].first].push_back({j, hits[begin].second});
			begin = end;
		}
	}
	return offsets;
}

inline std::vector<std::vector<anchored_offset_t> > find_target_offsets(const std::vector<correction_input_t>& inputs, const target_kmer_index_t& target_index) {
	std::vector<int> targets, ranks(inputs.size());
	int max_len = 0;
	for (int i = 0; i < int(inputs.size()); i++) {
		max_len = std::max(max_len, int(inputs[i].seq.length()));
		if (inputs[i].target) {
			ranks[i] = targets.size();
			targets.push_back(i);
		}
	}
	size_t width = 2*size_t(max_len)+1;
	if (targets.size()*width > (size_t(1) << 26)) return find_target_offsets_sorted(inputs, target_index);
	std::vector<uint32_t> counts(targets.size()*width, 0), touched, collected;
	std::vector<std::vector<anchored_offset_t> > offsets(inputs.size());
	const uint64_t mask = (uint64_t(1) << (2*K))-1;
	for (int j = 0; j < int(inputs.size()); j++) {
		const std::string& seq = inputs[j].seq;
		uint64_t kmer = 0;
		int valid = 0;
		for (int p = 0; p < int(seq.length()); p++) {
			int base = base_index(seq[p]);
			if (base < 0) {
				kmer = 0;
				valid = 0;
				continue;
			}
			kmer = ((kmer << 2) | base) & mask;
			if (++valid < K) continue;
			int pj = p+1-K;
			kmer_range_t range = target_index.find(kmer);
			for (int hit = range.begin; hit < range.end; hit++) {
				const kmer_entry_t& entry = target_index.entries[hit];
				if (entry.input == j) continue;
				uint32_t cell = ranks[entry.input]*width+(entry.pos-pj+max_len);
				if (counts[cell]++ == 0) touched.push_back(cell);
			}
		}
		collected.clear();
		for (uint32_t cell : touched) {
			if (counts[cell] >= MIN_ANCHORS) collected.push_back(cell);
			counts[cell] = 0;
		}
		touched.clear();
		std::sort(collected.begin(), collected.end());
		for (uint32_t cell : collected) offsets[targets[cell/width]].push_back({j, int(cell%width)-max_len});
	}
	return offsets;
}

struct overlap_t {
	int input, offset, begin, end, locus_term;
};

struct target_overlaps_t {
	std::vector<overlap_t> overlaps;
	std::vector<bool> disagreements;
};

inline target_overlaps_t find_target_overlaps(int i, const std::vector<correction_input_t>& inputs, const std::vector<anchored_offset_t>& offsets, const std::vector<std::vector<uint8_t> >& quality_classes) {
	target_overlaps_t result;
	result.disagreements.resize(inputs[i].seq.length(), false);
	const correction_input_t& target = inputs[i];
	for (int begin = 0; begin < int(offsets.size());) {
		int j = offsets[begin].input, end = begin+1;
		while (end < int(offsets.size()) && offsets[end].input == j) end++;
		const correction_input_t& other = inputs[j];
		overlap_t best{j, 0, 0, 0, 0};
		int best_score = 0, second_score = 0, candidates = 0;
		for (int k = begin; k < end; k++) {
			int d = offsets[k].offset;
			int a = std::max(0, d), b = std::min(int(target.seq.length()), int(other.seq.length())+d);
			int cur = 0, score = 0, cs = a, seg_begin = a, seg_end = a;
			for (int p = a; p < b; p++) {
				if (target.seq[p] == other.seq[p-d]) cur++;
				else cur -= (quality_classes[i][p] & 1) && (quality_classes[j][p-d] & 1) ? MISMATCH_HQ : MISMATCH_LQ;
				if (cur <= 0) {
					cur = 0;
					cs = p+1;
				} else if (cur > score) {
					score = cur;
					seg_begin = cs;
					seg_end = p+1;
				}
			}
			if (seg_end-seg_begin < MIN_OVERLAP) continue;
			candidates++;
			if (score > best_score) {
				second_score = best_score;
				best_score = score;
				best = {j, d, seg_begin, seg_end, 0};
			} else if (score > second_score) second_score = score;
		}
		if (candidates > 0 && (candidates == 1 || best_score-second_score >= REPEAT_MARGIN)) {
			for (int p = best.begin; p < best.end; p++) {
				int other_pos = p-best.offset;
				if (target.seq[p] == other.seq[other_pos]) continue;
				if ((quality_classes[i][p] & 2) && (quality_classes[j][other_pos] & 2)) best.locus_term -= LAMBDA_HQ;
				if (base_index(other.seq[other_pos]) >= 0) result.disagreements[p] = true;
			}
			result.overlaps.push_back(best);
		}
		begin = end;
	}
	return result;
}

struct pile_observation_t {
	int base, memo_index;
};

inline double base_likelihood(int observed, int allele, double error) {
	return observed == allele ? 1.0-error : error/3.0;
}

struct likelihood_terms_t {
	double match, mismatch;
	// Grid cases: observed base matches u, matches v, or matches neither.
	std::array<std::array<double, POLY_GRID>, 3> grid;
};

struct likelihood_memo_t {
	struct weight_t {
		double weight;
		uint64_t bits;
		std::array<int, 256> indices;
	};
	std::vector<weight_t> weights;
	std::vector<likelihood_terms_t> terms;

	int weight_index(double weight) {
		uint64_t bits;
		std::memcpy(&bits, &weight, sizeof(bits));
		for (int i = 0; i < int(weights.size()); i++) if (weights[i].bits == bits) return i;
		weight_t entry;
		entry.weight = weight;
		entry.bits = bits;
		entry.indices.fill(-1);
		weights.push_back(entry);
		return weights.size()-1;
	}

	int find(int weight_index, uint8_t qual) {
		int& index = weights[weight_index].indices[qual];
		if (index >= 0) return index;
		double weight = weights[weight_index].weight, error = base_error(qual);
		likelihood_terms_t entry;
		entry.match = std::log(weight*base_likelihood(0, 0, error)+(1.0-weight)/4.0);
		entry.mismatch = std::log(weight*base_likelihood(0, 1, error)+(1.0-weight)/4.0);
		for (int observed = 0; observed < 3; observed++) {
			double pu = base_likelihood(observed, 0, error), pv = base_likelihood(observed, 1, error);
			for (int g = 0; g < POLY_GRID; g++) {
				double f = (g+0.5)/POLY_GRID;
				entry.grid[observed][g] = std::log(weight*(f*pu+(1.0-f)*pv)+(1.0-weight)/4.0);
			}
		}
		index = terms.size();
		terms.push_back(entry);
		return index;
	}
};

inline double flip_posterior(int x, int y, double error, const std::vector<pile_observation_t>& pile, const std::array<double, 4>& ll, const likelihood_memo_t& memo) {
	static const double mono_prior = std::log((1.0-3.0*POLY_PRIOR)/4.0);
	static const double poly_prior = std::log(POLY_PRIOR);
	std::array<double, 7> hypotheses;
	for (int b = 0; b < 4; b++) hypotheses[b] = ll[b]+mono_prior;
	std::array<int, 4> ranked{{0, 1, 2, 3}};
	std::stable_sort(ranked.begin(), ranked.end(), [&ll](int a, int b) { return ll[a] > ll[b]; });
	const int pairs[3][2] = {{ranked[0], ranked[1]}, {ranked[0], ranked[2]}, {ranked[1], ranked[2]}};
	std::array<double, 3> shares;
	for (int pair = 0; pair < 3; pair++) {
		int u = pairs[pair][0], v = pairs[pair][1];
		std::array<double, POLY_GRID> grid_ll;
		grid_ll.fill(0.0);
		for (const pile_observation_t& observation : pile) {
			int observed = observation.base == u ? 0 : observation.base == v ? 1 : 2;
			const std::array<double, POLY_GRID>& terms = memo.terms[observation.memo_index].grid[observed];
			for (int g = 0; g < POLY_GRID; g++) grid_ll[g] += terms[g];
		}
		double max_ll = *std::max_element(grid_ll.begin(), grid_ll.end());
		double pu = base_likelihood(x, u, error), pv = base_likelihood(x, v, error);
		double s = 0.0, su = 0.0;
		for (int g = 0; g < POLY_GRID; g++) {
			double f = (g+0.5)/POLY_GRID, scale = std::exp(grid_ll[g]-max_ll);
			s += scale*(f*pu+(1.0-f)*pv);
			su += scale*f*pu;
		}
		s /= POLY_GRID;
		su /= POLY_GRID;
		hypotheses[4+pair] = max_ll+std::log(s)+poly_prior;
		shares[pair] = su/s;
	}
	double max_h = *std::max_element(hypotheses.begin(), hypotheses.end()), total = 0.0;
	std::array<double, 7> post;
	for (int h = 0; h < 7; h++) {
		post[h] = std::exp(hypotheses[h]-max_h);
		total += post[h];
	}
	for (double& probability : post) probability /= total;
	double probability = post[y];
	for (int pair = 0; pair < 3; pair++) {
		if (pairs[pair][0] == y) probability += post[4+pair]*shares[pair];
		else if (pairs[pair][1] == y) probability += post[4+pair]*(1.0-shares[pair]);
	}
	return probability;
}

inline void append_target_flips(int i, const std::vector<correction_input_t>& inputs, const target_overlaps_t& target, const std::vector<std::vector<uint8_t> >& quality_classes, likelihood_memo_t& memo, std::vector<correction_flip_t>& flips) {
	static const double locus_prior = std::log(SAME_LOCUS_PRIOR/(1.0-SAME_LOCUS_PRIOR));
	const correction_input_t& input = inputs[i];
	std::vector<std::array<int, 2> > weight_indices(target.overlaps.size());
	for (int k = 0; k < int(target.overlaps.size()); k++) {
		const overlap_t& overlap = target.overlaps[k];
		for (int high_quality = 0; high_quality < 2; high_quality++) {
			double term = high_quality ? -LAMBDA_HQ : 0.0;
			double llr = locus_prior+overlap.locus_term-term;
			double weight = llr <= -50.0 ? 0.0 : 1.0/(1.0+std::exp(-llr));
			weight_indices[k][high_quality] = memo.weight_index(weight);
		}
	}
	std::vector<pile_observation_t> pile;
	pile.reserve(target.overlaps.size());
	for (int p = 0; p < int(input.seq.length()); p++) {
		if (!target.disagreements[p]) continue;
		int x = base_index(input.seq[p]);
		if (x < 0) continue;
		double error = input_error(input, p);
		std::array<double, 4> ll;
		double match = std::log(1.0-error), mismatch = std::log(error/3.0);
		for (int b = 0; b < 4; b++) ll[b] = b == x ? match : mismatch;
		pile.clear();
		for (int k = 0; k < int(target.overlaps.size()); k++) {
			const overlap_t& overlap = target.overlaps[k];
			if (p < overlap.begin || p >= overlap.end) continue;
			const correction_input_t& other = inputs[overlap.input];
			int other_pos = p-overlap.offset, observed = base_index(other.seq[other_pos]);
			if (observed < 0) continue;
			int high_quality = observed != x && (quality_classes[i][p] & 2) && (quality_classes[overlap.input][other_pos] & 2);
			uint8_t qual = other.quals.size() == other.seq.size() ? other.quals[other_pos] : DEFAULT_QUAL;
			int memo_index = memo.find(weight_indices[k][high_quality], qual);
			pile.push_back({observed, memo_index});
			const likelihood_terms_t& terms = memo.terms[memo_index];
			for (int b = 0; b < 4; b++) ll[b] += b == observed ? terms.match : terms.mismatch;
		}
		int y = std::max_element(ll.begin(), ll.end())-ll.begin();
		if (y == x) continue;
		double posterior = flip_posterior(x, y, error, pile, ll, memo);
		if (posterior >= FLIP_THRESHOLD) flips.push_back({i, p, "ACGT"[y], posterior});
	}
}

} // namespace detail

inline std::vector<correction_flip_t> correct_reads(const std::vector<correction_input_t>& inputs) {
	if (std::none_of(inputs.begin(), inputs.end(), [](const correction_input_t& input) { return input.target; })) return {};
	std::vector<std::vector<uint8_t> > quality_classes(inputs.size());
	for (int i = 0; i < int(inputs.size()); i++) {
		quality_classes[i].resize(inputs[i].seq.size());
		for (int p = 0; p < int(inputs[i].seq.size()); p++) {
			double error = detail::input_error(inputs[i], p);
			quality_classes[i][p] = (error < OVERLAP_HQ_ERR ? 1 : 0) | (error < LOCUS_HQ_ERR ? 2 : 0);
		}
	}
	detail::target_kmer_index_t target_index(inputs);
	std::vector<std::vector<detail::anchored_offset_t> > offsets = detail::find_target_offsets(inputs, target_index);
	std::vector<correction_flip_t> flips;
	detail::likelihood_memo_t memo;
	for (int i = 0; i < int(inputs.size()); i++) {
		if (!inputs[i].target) continue;
		detail::target_overlaps_t target = detail::find_target_overlaps(i, inputs, offsets[i], quality_classes);
		detail::append_target_flips(i, inputs, target, quality_classes, memo, flips);
	}
	return flips;
}

} // namespace read_correction

#endif // READ_CORRECTION_H_
