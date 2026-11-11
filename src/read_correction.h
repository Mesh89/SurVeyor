#ifndef READ_CORRECTION_H_
#define READ_CORRECTION_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
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

struct target_kmer_index_t {
	std::vector<kmer_entry_t> entries;
	std::unordered_map<uint64_t, kmer_range_t> table;

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
		for (int begin = 0; begin < int(entries.size());) {
			int end = begin+1;
			while (end < int(entries.size()) && entries[end].kmer == entries[begin].kmer) end++;
			table[entries[begin].kmer] = {begin, end};
			begin = end;
		}
	}

	kmer_range_t find(uint64_t kmer) const {
		auto it = table.find(kmer);
		return it == table.end() ? kmer_range_t{0, 0} : it->second;
	}
};

struct anchored_offset_t { int input, offset; };

inline std::vector<std::vector<anchored_offset_t> > find_target_offsets(const std::vector<correction_input_t>& inputs, const target_kmer_index_t& target_index) {
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

struct overlap_t {
	int input, offset, begin, end, locus_term;
};

struct target_overlaps_t {
	std::vector<overlap_t> overlaps;
	std::vector<bool> disagreements;
};

inline target_overlaps_t find_target_overlaps(int i, const std::vector<correction_input_t>& inputs, const std::vector<anchored_offset_t>& offsets) {
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
				else cur -= input_error(target, p) < OVERLAP_HQ_ERR && input_error(other, p-d) < OVERLAP_HQ_ERR ? MISMATCH_HQ : MISMATCH_LQ;
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
				if (input_error(target, p) < LOCUS_HQ_ERR && input_error(other, other_pos) < LOCUS_HQ_ERR) best.locus_term -= LAMBDA_HQ;
				if (base_index(other.seq[other_pos]) >= 0) result.disagreements[p] = true;
			}
			result.overlaps.push_back(best);
		}
		begin = end;
	}
	return result;
}

struct pile_observation_t {
	int base;
	double error, weight;
};

inline double base_likelihood(int observed, int allele, double error) {
	return observed == allele ? 1.0-error : error/3.0;
}

inline double flip_posterior(int x, int y, double error, const std::vector<pile_observation_t>& pile, const std::array<double, 4>& ll) {
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
		for (int g = 0; g < POLY_GRID; g++) {
			double f = (g+0.5)/POLY_GRID, sum = 0.0;
			for (const pile_observation_t& observation : pile) {
				double pu = base_likelihood(observation.base, u, observation.error), pv = base_likelihood(observation.base, v, observation.error);
				sum += std::log(observation.weight*(f*pu+(1.0-f)*pv)+(1.0-observation.weight)/4.0);
			}
			grid_ll[g] = sum;
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

inline void append_target_flips(int i, const std::vector<correction_input_t>& inputs, const target_overlaps_t& target, std::vector<correction_flip_t>& flips) {
	static const double locus_prior = std::log(SAME_LOCUS_PRIOR/(1.0-SAME_LOCUS_PRIOR));
	const correction_input_t& input = inputs[i];
	std::vector<pile_observation_t> pile;
	pile.reserve(target.overlaps.size());
	for (int p = 0; p < int(input.seq.length()); p++) {
		if (!target.disagreements[p]) continue;
		int x = base_index(input.seq[p]);
		if (x < 0) continue;
		double error = input_error(input, p);
		std::array<double, 4> ll;
		for (int b = 0; b < 4; b++) ll[b] = std::log(base_likelihood(x, b, error));
		pile.clear();
		for (const overlap_t& overlap : target.overlaps) {
			if (p < overlap.begin || p >= overlap.end) continue;
			const correction_input_t& other = inputs[overlap.input];
			int other_pos = p-overlap.offset, observed = base_index(other.seq[other_pos]);
			if (observed < 0) continue;
			double other_error = input_error(other, other_pos);
			double term = observed != x && error < LOCUS_HQ_ERR && other_error < LOCUS_HQ_ERR ? -LAMBDA_HQ : 0.0;
			double llr = locus_prior+overlap.locus_term-term;
			double weight = llr <= -50.0 ? 0.0 : 1.0/(1.0+std::exp(-llr));
			pile.push_back({observed, other_error, weight});
			for (int b = 0; b < 4; b++) ll[b] += std::log(weight*base_likelihood(observed, b, other_error)+(1.0-weight)/4.0);
		}
		int y = std::max_element(ll.begin(), ll.end())-ll.begin();
		if (y == x) continue;
		double posterior = flip_posterior(x, y, error, pile, ll);
		if (posterior >= FLIP_THRESHOLD) flips.push_back({i, p, "ACGT"[y], posterior});
	}
}

} // namespace detail

inline std::vector<correction_flip_t> correct_reads(const std::vector<correction_input_t>& inputs) {
	if (std::none_of(inputs.begin(), inputs.end(), [](const correction_input_t& input) { return input.target; })) return {};
	detail::target_kmer_index_t target_index(inputs);
	std::vector<std::vector<detail::anchored_offset_t> > offsets = detail::find_target_offsets(inputs, target_index);
	std::vector<correction_flip_t> flips;
	for (int i = 0; i < int(inputs.size()); i++) {
		if (!inputs[i].target) continue;
		detail::target_overlaps_t target = detail::find_target_overlaps(i, inputs, offsets[i]);
		detail::append_target_flips(i, inputs, target, flips);
	}
	return flips;
}

} // namespace read_correction

#endif // READ_CORRECTION_H_
