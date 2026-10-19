/* Instantiated for byte/word lanes and SSE2/AVX2 by ssw_clip_impl.h.
 * H stores raw score plus clip for paths starting at query base zero.
 * Local restarts stay at zero. Subtract 2*clip at interior endpoints and
 * clip at the last query base to obtain the penalized local score. */
static int CP_KERNEL(const s_profile* p, const int8_t* ref, int nr, int reverse,
        int go, int ge, int clip, int stop, s_align* a, int* columns) {
    const int lanes = CP_BYTES/sizeof(CP_INT), nq = p->readLen;
    const int seg = (nq+lanes-1)/lanes, offset = (nq-1)%seg, lane = (nq-1)/seg;
    if (clip+CP_BIAS >= CP_LIMIT) return 0;
    CP_VEC* mem = (CP_VEC*)CP_ALLOC((4+(a->prefix_scores != NULL))*(size_t)seg+(columns ? 2 : 0));
    if (!mem) return -1;
    CP_VEC *h = mem, *old = h+seg, *e = old+seg, *best_h = e+seg;
    CP_VEC* prefix = a->prefix_scores ? best_h+seg : NULL;
    CP_VEC* valid = columns ? best_h+seg+(prefix ? seg : 0) : NULL;
    const CP_VEC zero = CP_ZERO(), gap_o = CP_SET(go), gap_e = CP_SET(ge);
    if (valid) for (int k = 0; k < lanes; ++k) {
        ((CP_INT*)valid)[k] = k < nq/seg ? (CP_INT)-1 : 0;
        ((CP_INT*)(valid+1))[k] = k == nq/seg ? (CP_INT)-1 : 0;
    }
    const CP_VEC boundary = CP_LOW(clip);
#if CP_LIMIT == 255
    const CP_VEC bias = CP_SET(CP_BIAS);
#endif
    const CP_VEC boundary_f = CP_LOW(clip > go ? clip-go : 0);
    const CP_VEC* profile = CP_PROFILE;
    CP_VEC* temporary = NULL;
    if (!profile) profile = temporary = CP_MAKE_PROFILE;
    if (!profile) { free(mem); return -1; }
#ifdef CP_FULL_BYTE
    CP_VEC* unbiased = (CP_VEC*)CP_ALLOC((size_t)p->n*seg);
    if (!unbiased) { free(temporary); free(mem); return -1; }
    const CP_VEC sign = CP_SET(-128), bias = CP_SET(p->bias);
    for (int j = 0; j < p->n*seg; ++j) unbiased[j] = CP_UNBIAS(profile[j], bias);
    free(temporary); profile = temporary = unbiased;
#endif
    for (int j = 0; j < seg; ++j) {
        for (int k = 0; k < lanes; ++k) {
            int64_t value = clip-go-(int64_t)(j+k*seg)*ge;
            ((CP_INT*)(h+j))[k] = value > 0 ? value : 0;
        }
        e[j] = CP_DEC(h[j], gap_o);
        if (prefix) prefix[j] = zero;
    }
    int best = 0, end = -1, raw_max = 0, exact = 1;
    CP_VEC maximum = zero;
    for (int i = 0; i < nr; ++i) {
        CP_VEC f = boundary_f, v = CP_OR(CP_SHIFT(h[seg-1]), boundary), column = zero;
        CP_VEC* swap = old; old = h; h = swap;
        const CP_VEC* scores = profile+ref[reverse ? nr-1-i : i]*seg;
        for (int j = 0; j < seg; ++j) {
            v = CP_DIAG(v, scores[j]);
            if (prefix) prefix[j] = CP_MAX(prefix[j], v);
            v = CP_MAX(v, CP_MAX(e[j], f));
            column = CP_MAX(column, v);
            h[j] = v;
            CP_VEC opened = CP_DEC(v, gap_o);
            e[j] = CP_MAX(CP_DEC(e[j], gap_e), opened);
            f = CP_MAX(CP_DEC(f, gap_e), opened);
            v = old[j];
        }
        for (int sweep = 0; sweep < 2; ++sweep) {
            if (sweep) f = CP_PREFIX(f, (int64_t)seg*ge);
            f = CP_OR(CP_SHIFT(f), boundary_f);
            for (int j = 0; j < seg; ++j) {
                v = CP_MAX(h[j], f);
                h[j] = v;
                column = CP_MAX(column, v);
                CP_VEC opened = CP_DEC(v, gap_o);
                e[j] = CP_MAX(e[j], opened);
                f = CP_DEC(f, gap_e);
                if (go > ge && CP_EQUAL(CP_DEC(f, opened), zero)) goto repaired;
            }
        }
repaired:
        {
            CP_VEC next = CP_MAX(maximum, column);
            if (!CP_EQUAL(maximum, next)) {
                raw_max = CP_REDUCE(next);
                maximum = CP_SET(raw_max);
                if (raw_max+CP_BIAS >= CP_LIMIT) { exact = 0; break; }
            }
            int last = ((CP_INT*)(h+offset))[lane]-clip;
            int score = raw_max-2*clip;
            if (last > score) score = last;
            if (score > best) {
                best = score; end = i;
                memcpy(best_h, h, (size_t)seg*sizeof(CP_VEC));
            }
            if (columns) {
                /* Padding can retain scores after the real query has ended.
                 * It cannot raise the global maximum, but must be excluded
                 * when comparing individual reference columns. */
                column = CP_AND(column, valid[0]);
                if (nq%seg) {
                    CP_VEC partial = zero;
                    for (int j = 0; j < nq%seg; ++j) partial = CP_MAX(partial, h[j]);
                    column = CP_MAX(column, CP_AND(partial, valid[1]));
                }
                int value = CP_REDUCE(column)-2*clip;
                if (last > value) value = last;
                columns[i] = value > 0 ? value : 0;
            }
            if (best == stop) break;
        }
    }
    a->score1 = best; a->ref_end1 = end; a->read_end1 = -1;
    if (best) for (int j = 0; j < nq; ++j) {
        int value = ((CP_INT*)(best_h+j%seg))[j/seg]-clip-(j == nq-1 ? 0 : clip);
        if (value == best) { a->read_end1 = j; break; }
    }
    if (prefix) for (int j = 0; j < nq; ++j) {
        int value = ((CP_INT*)(prefix+j%seg))[j/seg]-clip-(j == nq-1 ? 0 : clip);
        a->prefix_scores[j] = value > 0 ? value : 0;
    }
    free(temporary); free(mem);
    return exact;
}
