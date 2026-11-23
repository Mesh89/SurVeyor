/* Included by each backend to reuse its existing traceback and CIGAR validation. */

static int consensus_exact_start(const int8_t* read, int nq, const int8_t* ref, int nr) {
	if (nr < nq) return -1;
	if (!memcmp(ref, read, nq)) return 0;
	if (nr == nq) return -1;

	/* A prefix table keeps the search linear even for highly repetitive targets. */
	int* prefix = (int*)malloc((size_t)nq * sizeof(int));
	if (!prefix) return -1;
	prefix[0] = 0;
	for (int i = 1, matched = 0; i < nq; ++i) {
		while (matched && read[i] != read[matched]) matched = prefix[matched - 1];
		if (read[i] == read[matched]) ++matched;
		prefix[i] = matched;
	}
	int start = -1;
	for (int i = 0, matched = 0; i < nr; ++i) {
		while (matched && ref[i] != read[matched]) matched = prefix[matched - 1];
		if (ref[i] == read[matched]) ++matched;
		if (matched == nq) {
			start = i - nq + 1;
			break;
		}
	}
	free(prefix);
	return start;
}

enum { CONSENSUS_MAX_LOSS = 96, CONSENSUS_SEED_SLOTS = 64 };

typedef struct {
	uint32_t prefix;
	int first;
} consensus_seed_slot_t;

static consensus_seed_slot_t* consensus_seed_slot(consensus_seed_slot_t* slots, uint32_t prefix) {
	unsigned i = (prefix * UINT32_C(0x9e3779b1)) >> 26;
	while (slots[i].first && slots[i].prefix != prefix) i = (i + 1) % CONSENSUS_SEED_SLOTS;
	return slots + i;
}

/* The rolling prefix only locates candidates. Verify the full seed. */
static int consensus_seed_diagonals(const int8_t* read, int nq, const int8_t* ref, int nr,
		int budget, int* low, int* high) {
	const int nseeds = budget / 5 + 1;
	const int separation = 4;
	const int len = (nq - 2 * budget - (nseeds - 1) * separation) / nseeds;
	if (len < 8 || len > nr) return 0;
	const int stride = len + separation;
	const int prefix_len = len < 16 ? len : 16;
	const uint32_t prefix_mask = UINT32_MAX >> (32 - 2 * prefix_len);
	consensus_seed_slot_t slots[CONSENSUS_SEED_SLOTS] = {{0, 0}};
	int next[CONSENSUS_MAX_LOSS / 5 + 1];
	for (int seed = 0; seed < nseeds; ++seed) {
		const int8_t* bases = read + budget + seed * stride;
		uint32_t prefix = 0;
		for (int j = 0; j < prefix_len; ++j) prefix = (prefix << 2) | bases[j];
		consensus_seed_slot_t* slot = consensus_seed_slot(slots, prefix);
		next[seed] = slot->first;
		slot->prefix = prefix;
		slot->first = seed + 1;
	}

	*low = INT_MAX;
	*high = INT_MIN;
	uint32_t prefix = 0;
	for (int j = 0; j < prefix_len - 1; ++j) prefix = (prefix << 2) | ref[j];
	for (int i = 0; i <= nr - len; ++i) {
		prefix = ((prefix << 2) | ref[i + prefix_len - 1]) & prefix_mask;
		for (int id = consensus_seed_slot(slots, prefix)->first; id; id = next[id - 1]) {
			int q = budget + (id - 1) * stride;
			if (memcmp(read + q, ref + i, len)) continue;
			int diagonal = i - q;
			if (diagonal < *low) *low = diagonal;
			if (diagonal > *high) *high = diagonal;
			if ((int64_t)*high - *low > 2 * budget) return 0;
		}
	}
	return *low != INT_MAX;
}

typedef struct {
	int16_t match, ins, del;
} consensus_wave_cell_t;

/* Advance the furthest query offset on each diagonal in increasing loss.
 * Relative to a +1 match, a mismatch costs 5, an insertion costs 7 then 2
 * per extra base, and a deletion costs 6 then 1. The initial offset s
 * accounts for left clipping. Reaching query offset q gives score q - s. */
static alignment_end consensus_wavefront_end(const int8_t* read, int nq, const int8_t* ref, int nr,
		int low, int high, int budget, int reverse) {
	enum { WAVE_ROWS = 8, WAVE_COLUMNS = 4 * CONSENSUS_MAX_LOSS + 3 };
	consensus_wave_cell_t waves[WAVE_ROWS][WAVE_COLUMNS];
	memset(waves, -1, sizeof(waves));
	alignment_end best = {0, 0, 0};
	int limit = budget;
	for (int s = 0; s <= limit; ++s) {
		consensus_wave_cell_t* row = waves[s % WAVE_ROWS];
		const consensus_wave_cell_t* mismatch = waves[(s + WAVE_ROWS - 5) % WAVE_ROWS];
		const consensus_wave_cell_t* open_ins = waves[(s + WAVE_ROWS - 7) % WAVE_ROWS];
		const consensus_wave_cell_t* extend_ins = waves[(s + WAVE_ROWS - 2) % WAVE_ROWS];
		const consensus_wave_cell_t* open_del = waves[(s + WAVE_ROWS - 6) % WAVE_ROWS];
		const consensus_wave_cell_t* extend_del = waves[(s + WAVE_ROWS - 1) % WAVE_ROWS];
		for (int d = low, i = 1; d <= high; ++d, ++i) {
			int end = nr - d < nq ? nr - d : nq;
			int q = s <= end && s + d >= 0 ? s : -1;
			int pos = mismatch[i].match;
			if (pos >= 0 && pos < end && pos + 1 > q) q = pos + 1;

			pos = open_ins[i + 1].match;
			int ins = pos >= 0 && pos < end ? pos + 1 : -1;
			pos = extend_ins[i + 1].ins;
			if (pos >= 0 && pos < end && pos + 1 > ins) ins = pos + 1;
			row[i].ins = ins;
			if (ins > q) q = ins;

			pos = open_del[i - 1].match;
			int del = pos <= end ? pos : -1;
			pos = extend_del[i - 1].del;
			if (pos <= end && pos > del) del = pos;
			row[i].del = del;
			if (del > q) q = del;

			if (q >= 0) {
				if (reverse) {
					while (q + 8 <= end && !memcmp(read + (nq - q - 8), ref + (nr - q - d - 8), 8)) q += 8;
					while (q < end && read[nq - q - 1] == ref[nr - q - d - 1]) ++q;
				} else {
					while (q + 8 <= end && !memcmp(read + q, ref + (q + d), 8)) q += 8;
					while (q < end && read[q] == ref[q + d]) ++q;
				}
				int score = q - s, r = q + d - 1;
				/* Include the last loss row: a full-query end can tie an earlier clipped end. */
				if (score >= nq - limit && (score > best.score || (r < best.ref || (r == best.ref && q - 1 < best.read)))) {
					best.score = score;
					best.ref = r;
					best.read = q - 1;
					limit = nq - score;
				}
			}
			row[i].match = q;
		}
	}
	return best;
}

enum ssw_consensus_path ssw_try_align_consensus(const int8_t* read, int32_t nq,
		const int8_t* ref, int32_t nr, const int8_t* mat, int32_t n,
		uint8_t go, uint8_t ge, s_align* result) {
	/* Stay below the SIMD word limit and use the standard genotype scoring. */
	if (n < 4 || n > 128 || go != 6 || ge != 1 || nq <= 0 || nq >= 32767 || nr <= 0) return SSW_CONSENSUS_FULL;
	for (int r = 0; r < 4; ++r) for (int q = 0; q < 4; ++q)
		if (mat[r * n + q] != (r == q ? 1 : -4)) return SSW_CONSENSUS_FULL;
	for (int i = 0; i < nq; ++i) if ((unsigned)read[i] >= 4) return SSW_CONSENSUS_FULL;
	for (int i = 0; i < nr; ++i) if ((unsigned)ref[i] >= 4) return SSW_CONSENSUS_FULL;

	int start = consensus_exact_start(read, nq, ref, nr);
	s_align a = {0};
	a.ref_end2 = -1;
	if (start >= 0) {
		/* Only a full ungapped match can score nq. The first occurrence ends first. */
		a.cigar = (uint32_t*)malloc(sizeof(uint32_t));
		if (!a.cigar) return SSW_CONSENSUS_FULL;
		a.cigar[0] = to_cigar_int(nq, 'M');
		a.cigarLen = 1;
		a.score1 = nq;
		a.ref_begin1 = start;
		a.ref_end1 = start + nq - 1;
		a.read_end1 = nq - 1;
		*result = a;
		return SSW_CONSENSUS_EXACT;
	}

	const int budget = nq / 5 < CONSENSUS_MAX_LOSS ? nq / 5 : CONSENSUS_MAX_LOSS;
	if (nr < nq - budget) return SSW_CONSENSUS_FULL;
	int seed_low, seed_high;
	if (!consensus_seed_diagonals(read, nq, ref, nr, budget, &seed_low, &seed_high)) return SSW_CONSENSUS_FULL;
	/* Clipping within the budget cannot touch a seed. A mismatch costs 5,
	 * a deletion breaks at most one seed, and an insertion covering k
	 * seeds must span at least 5*k-4 bases. Thus each broken seed costs at
	 * least 5, and every qualifying path contains a verified seed. Its
	 * diagonal can move by at most the loss, including on tied paths. */
	int low = seed_low - budget, high = seed_high + budget;
	int offset = low > 0 ? low : 0;
	int64_t last = (int64_t)nq + high;
	int span = (last < nr ? (int)last : nr) - offset;
	alignment_end end = consensus_wavefront_end(read, nq, ref + offset, span, low - offset, high - offset, budget, 0);
	if (end.score < nq - budget) return SSW_CONSENSUS_FULL;

	/* Reflect the band around the winning end to search the reversed prefixes. */
	int loss = nq - end.score;
	low = seed_low - loss - offset;
	high = seed_high + loss - offset;
	int end_diagonal = end.ref - end.read;
	alignment_end begin = consensus_wavefront_end(read, end.read + 1, ref + offset, end.ref + 1, end_diagonal - high, end_diagonal - low, end.read + 1 - end.score, 1);
	if (begin.score != end.score) return SSW_CONSENSUS_FULL;

	a.score1 = end.score;
	a.ref_end1 = offset + end.ref;
	a.read_end1 = end.read;
	a.ref_begin1 = a.ref_end1 - begin.ref;
	a.read_begin1 = end.read - begin.read;
	int ref_span = begin.ref + 1, read_span = begin.read + 1;
	cigar* path = banded_sw(ref + a.ref_begin1, read + a.read_begin1, ref_span, read_span, a.score1, go, ge, abs(ref_span - read_span) + 1, mat, n);
	if (!path) return SSW_CONSENSUS_FULL;
	a.cigar = path->seq;
	a.cigarLen = path->length;
	free(path);
	s_profile profile = {0};
	profile.read = read;
	profile.mat = mat;
	profile.n = n;
	if (!ssw_valid_cigar(&profile, ref, &a, go, ge)) {
		free(a.cigar);
		return SSW_CONSENSUS_FULL;
	}
	*result = a;
	return SSW_CONSENSUS_BANDED;
}
