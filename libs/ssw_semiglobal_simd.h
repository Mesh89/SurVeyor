/* Instantiated for byte and word lanes by ssw_semiglobal_impl.h.
 * Store B + H(i,j) - i*match_max. Every transition is a nonnegative
 * subtraction, so a saturated zero can never recover and corrupt a positive
 * result. Only a zero at a requested endpoint requires a wider retry. */
static int SG_KERNEL(const s_profile* p, const int8_t* ref, int nr, int reverse,
        int go, int ge, int match_max, int stop, int mask, s_align_semiglobal* a, int32_t* columns) {
    const int lanes = SG_BYTES / sizeof(SG_INT), nq = p->readLen;
    const int seg = (nq + lanes - 1) / lanes;
    const int offset = (nq-1) % seg, lane = (nq-1) / seg;
    SG_VEC* mem = (SG_VEC*)SG_ALLOC((p->n+3+(a->prefix_scores != NULL))*(size_t)seg);
    if (!mem) return -1;
    SG_VEC *profile = mem, *h = profile+p->n*seg, *old = h+seg, *e = old+seg;
    SG_VEC* prefix = a->prefix_scores ? e+seg : NULL;
    const SG_VEC zero = SG_ZERO();
    const SG_VEC gap_o = SG_SET(go), gap_e = SG_SET(ge);
    const SG_VEC ins_o = SG_SET(go+match_max > SG_LIMIT ? SG_LIMIT : go+match_max);
    const SG_VEC ins_e = SG_SET(ge+match_max > SG_LIMIT ? SG_LIMIT : ge+match_max);
    const SG_VEC boundary_h = SG_LOW(SG_LIMIT);
    const SG_VEC boundary_f = SG_LOW(go+match_max < SG_LIMIT ? SG_LIMIT-go-match_max : 0);
    const SG_VEC* qp = SG_PROFILE;
    SG_VEC* temporary = NULL;
    if (!qp) qp = temporary = SG_MAKE_PROFILE;
    const SG_VEC bias = SG_SET(SG_PROFILE_BIAS);
    for (int j = 0; j < p->n*seg; ++j) profile[j] = SG_SUB(bias, qp[j]);
    free(temporary);
    for (int j = 0; j < seg; ++j) {
        for (int k = 0; k < lanes; ++k) {
            int64_t i = j+k*seg+1;
            int64_t value = SG_LIMIT-go-(i-1)*ge-i*match_max;
            ((SG_INT*)(h+j))[k] = value > 0 ? value : 0;
        }
        e[j] = SG_DEC(h[j], gap_o);
        if (prefix) prefix[j] = zero;
    }
    int best = ((SG_INT*)(h+offset))[lane], end = -1;
    for (int i = 0; i < nr; ++i) {
        SG_VEC f = boundary_f, v = SG_OR(SG_SHIFT(h[seg-1]), boundary_h);
        SG_VEC* swap = old; old = h; h = swap;
        const SG_VEC* scores = profile + ref[reverse ? nr-1-i : i]*seg;
        for (int j = 0; j < seg; ++j) {
            v = SG_DEC(v, scores[j]);
            if (prefix) prefix[j] = SG_MAX(prefix[j], v);
            v = SG_MAX(v, SG_MAX(e[j], f));
            h[j] = v;
            e[j] = SG_MAX(SG_DEC(e[j], gap_e), SG_DEC(v, gap_o));
            f = SG_MAX(SG_DEC(f, ins_e), SG_DEC(v, ins_o));
            v = old[j];
        }
        /* Close lane crossings before the repair sweep. The discount includes
         * the query bias; reopening a gap cannot improve its extension. */
        {
            f = SG_PREFIX(f, (int64_t)seg*(ge+match_max));
            f = SG_OR(SG_SHIFT(f), boundary_f);
            for (int j = 0; j < seg; ++j) {
                v = SG_MAX(h[j], f);
                h[j] = v;
                e[j] = SG_MAX(e[j], SG_DEC(v, gap_o));
                f = SG_DEC(f, ins_e);
            }
        }
        {
            int value = ((SG_INT*)(h+offset))[lane];
            if (value > best) { best = value; end = i; }
            if (columns) columns[i] = value ? value-SG_LIMIT+nq*match_max : INT32_MIN;
            if (best && (best-SG_LIMIT+nq*match_max == stop || (best == SG_LIMIT && !columns && !prefix))) break;
        }
    }
    int exact = best != 0;
    a->score1 = best-SG_LIMIT+nq*match_max;
    a->ref_end1 = end;
    if (prefix) for (int i = 0; i < nq; ++i) {
        int value = ((SG_INT*)(prefix+i%seg))[i/seg];
        a->prefix_scores[i] = value-SG_LIMIT+(i+1)*match_max;
        if (!value) exact = 0;
    }
    /* A masked second-best endpoint may need more precision than the best. */
    if (columns) {
        int next = INT32_MIN, available = 0;
        for (int i = 0; i < nr; ++i) if (i < end-(int64_t)mask || i >= end+(int64_t)mask) {
            available = 1;
            if (columns[i] > next) next = columns[i];
        }
        if (available && next == INT32_MIN) exact = 0;
    }
    free(mem);
    return exact;
}
