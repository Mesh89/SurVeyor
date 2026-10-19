/* Fixed penalty per clipped query end. Shared by SSE2 and AVX2. */
#ifdef CP_AVX2
#define CP_VEC __m256i
#define CP_BYTES 32
#define CP_ALLOC(n) ssw_avx2_alloc(n, 0)
#define CP_ZERO _mm256_setzero_si256
#define CP_OR _mm256_or_si256
#define CP_AND _mm256_and_si256
#define CP_LOW(x) _mm256_inserti128_si256(_mm256_setzero_si256(), _mm_cvtsi32_si128(x), 0)
#define CP_SHIFT ssw_avx2_shift1
#define CP_SET(x) _mm256_set1_epi8((char)(x))
#define CP_DEC _mm256_subs_epu8
#define CP_MAX _mm256_max_epu8
#define CP_DIAG(v, s) CP_DEC(_mm256_adds_epu8(v, s), bias)
#define CP_EQUAL(x, y) (_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)) == -1)
#define CP_PREFIX ssw_avx2_prefix_byte
#else
#define CP_VEC __m128i
#define CP_BYTES 16
#define CP_ALLOC(n) malloc((n)*sizeof(__m128i))
#define CP_ZERO _mm_setzero_si128
#define CP_OR _mm_or_si128
#define CP_AND _mm_and_si128
#define CP_LOW _mm_cvtsi32_si128
#define CP_SHIFT(x) _mm_slli_si128(x, 1)
#define CP_SET(x) _mm_set1_epi8((char)(x))
#define CP_DEC _mm_subs_epu8
#define CP_MAX _mm_max_epu8
#define CP_DIAG(v, s) CP_DEC(_mm_adds_epu8(v, s), bias)
#define CP_EQUAL(x, y) (_mm_movemask_epi8(_mm_cmpeq_epi8(x, y)) == 0xffff)
#define CP_PREFIX ssw_prefix_byte
#endif

static inline int cp_max_byte(CP_VEC v) {
#ifdef CP_AVX2
    __m128i x = _mm_max_epu8(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
#else
    __m128i x = v;
#endif
    x = _mm_max_epu8(x, _mm_srli_si128(x, 8));
    x = _mm_max_epu8(x, _mm_srli_si128(x, 4));
    x = _mm_max_epu8(x, _mm_srli_si128(x, 2));
    x = _mm_max_epu8(x, _mm_srli_si128(x, 1));
    return _mm_cvtsi128_si32(x)&255;
}
#define CP_INT uint8_t
#define CP_LIMIT 255
#define CP_BIAS p->bias
#define CP_REDUCE cp_max_byte
#define CP_KERNEL cp_byte
#define CP_PROFILE p->profile_byte
#define CP_MAKE_PROFILE qP_byte(p->read, p->mat, nq, p->n, p->bias)
#include "ssw_clip_simd.h"
#undef CP_INT
#undef CP_LIMIT
#undef CP_BIAS
#undef CP_KERNEL
#undef CP_DIAG

/* Near the byte ceiling, a biased add can overflow even when H itself fits.
 * A signed add to H-128 retains all 256 byte values. Use this
 * specialization only when query_length*max_match + clip <= 255, so even a
 * score of 255 is known to be exact. Normal byte/word kernels stay unchanged. */
#ifdef CP_AVX2
#define CP_DIAG(v, s) _mm256_xor_si256(_mm256_adds_epi8(_mm256_xor_si256(v, sign), s), sign)
#define CP_UNBIAS _mm256_sub_epi8
#else
#define CP_DIAG(v, s) _mm_xor_si128(_mm_adds_epi8(_mm_xor_si128(v, sign), s), sign)
#define CP_UNBIAS _mm_sub_epi8
#endif
#define CP_FULL_BYTE
#define CP_INT uint8_t
#define CP_LIMIT 256
#define CP_BIAS 0
#define CP_KERNEL cp_byte_full
#include "ssw_clip_simd.h"
#undef CP_FULL_BYTE
#undef CP_UNBIAS
#undef CP_INT
#undef CP_LIMIT
#undef CP_BIAS
#undef CP_REDUCE
#undef CP_KERNEL
#undef CP_PROFILE
#undef CP_MAKE_PROFILE
#undef CP_SHIFT
#undef CP_SET
#undef CP_DEC
#undef CP_MAX
#undef CP_DIAG
#undef CP_EQUAL
#undef CP_PREFIX

#ifdef CP_AVX2
#define CP_SHIFT ssw_avx2_shift2
#define CP_SET(x) _mm256_set1_epi16(x)
#define CP_DEC _mm256_subs_epu16
#define CP_MAX _mm256_max_epi16
#define CP_DIAG _mm256_adds_epi16
#define CP_EQUAL(x, y) (_mm256_movemask_epi8(_mm256_cmpeq_epi16(x, y)) == -1)
#define CP_PREFIX ssw_avx2_prefix_word
#else
#define CP_SHIFT(x) _mm_slli_si128(x, 2)
#define CP_SET(x) _mm_set1_epi16(x)
#define CP_DEC _mm_subs_epu16
#define CP_MAX _mm_max_epi16
#define CP_DIAG _mm_adds_epi16
#define CP_EQUAL(x, y) (_mm_movemask_epi8(_mm_cmpeq_epi16(x, y)) == 0xffff)
#define CP_PREFIX ssw_prefix_word
#endif

static inline int cp_max_word(CP_VEC v) {
#ifdef CP_AVX2
    __m128i x = _mm_max_epi16(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
#else
    __m128i x = v;
#endif
    x = _mm_max_epi16(x, _mm_srli_si128(x, 8));
    x = _mm_max_epi16(x, _mm_srli_si128(x, 4));
    x = _mm_max_epi16(x, _mm_srli_si128(x, 2));
    return _mm_extract_epi16(x, 0);
}
#define CP_INT uint16_t
#define CP_LIMIT 32767
#define CP_BIAS 0
#define CP_REDUCE cp_max_word
#define CP_KERNEL cp_word
#define CP_PROFILE p->profile_word
#define CP_MAKE_PROFILE qP_word(p->read, p->mat, nq, p->n)
#include "ssw_clip_simd.h"
#undef CP_INT
#undef CP_LIMIT
#undef CP_BIAS
#undef CP_REDUCE
#undef CP_KERNEL
#undef CP_PROFILE
#undef CP_MAKE_PROFILE
#undef CP_SHIFT
#undef CP_SET
#undef CP_DEC
#undef CP_MAX
#undef CP_DIAG
#undef CP_EQUAL
#undef CP_PREFIX
#undef CP_VEC
#undef CP_BYTES
#undef CP_ALLOC
#undef CP_ZERO
#undef CP_OR
#undef CP_AND
#undef CP_LOW

/* Same exceptional cases as local SSW: word overflow or open < extend. */
static int cp_scalar(const s_profile* p, const int8_t* ref, int nr, int go, int ge, int clip,
        uint8_t flag, int filters, int filterd, s_align* a, int* columns) {
    int nq = p->readLen;
    int* h = (int*)malloc((nq+1)*sizeof(int));
    int* e = (int*)calloc(nq+1, sizeof(int));
    int* m = go < ge ? (int*)calloc(nq+1, sizeof(int)) : NULL;
    int* ins = go < ge ? (int*)malloc((nq+1)*sizeof(int)) : NULL;
    uint8_t* dir = flag & 15 ? (uint8_t*)malloc((size_t)nq*nr) : NULL;
    if (!h || !e || (go < ge && (!m || !ins)) || ((flag&15) && !dir)) { free(h); free(e); free(m); free(ins); free(dir); return 0; }
    h[0] = clip;
    for (int j = 1; j <= nq; ++j) {
        int64_t value = clip-go-(int64_t)(j-1)*ge;
        h[j] = value > 0 ? value : 0;
        if (ins) ins[j] = h[j];
        if (a->prefix_scores) a->prefix_scores[j-1] = 0;
    }
    a->score1 = 0; a->ref_end1 = a->read_end1 = -1;
    int ok = 1;
    for (int i = 0; i < nr; ++i) {
        int diagonal = clip, f = 0, ml = clip, el = 0, column = 0;
        for (int j = 1; j <= nq; ++j) {
            int old = h[j], eo = old-go, ee = e[j]-ge, fo = h[j-1]-go, fe = f-ge;
            uint8_t ep = 0, fp = 0;
            if (m) {
                eo = (m[j] >= ins[j] ? m[j] : ins[j])-go; ep = m[j] >= ins[j] ? 1 : 3;
                fo = (ml >= el ? ml : el)-go; fp = ml >= el ? 1 : 2;
            }
            if (eo <= ee) ep = 2;
            if (fo <= fe) fp = 3;
            e[j] = eo > ee ? eo : ee; f = fo > fe ? fo : fe;
            int v = diagonal+p->mat[ref[i]*p->n+p->read[j-1]], cost = clip+(j == nq ? 0 : clip);
            uint8_t d = v > 0 ? 1 : 0;
            if (v < 0) v = 0;
            if (a->prefix_scores && v-cost > a->prefix_scores[j-1]) {
                if (v-cost > 65535) { errno = ERANGE; ok = 0; goto done; }
                a->prefix_scores[j-1] = v-cost;
            }
            if (m) { m[j] = ml = v; ins[j] = f; el = e[j]; }
            if (e[j] > v) { v = e[j]; d = 2; }
            if (f > v) { v = f; d = 3; }
            if (dir) dir[(size_t)i*nq+j-1] = d|(ep<<2)|(fp<<4);
            h[j] = v; diagonal = old;
            int score = v-cost;
            if (score > 65535) { errno = ERANGE; ok = 0; goto done; }
            if (score > column) column = score;
            if (score > a->score1) { a->score1 = score; a->ref_end1 = i; a->read_end1 = j-1; }
        }
        if (columns) columns[i] = column;
    }
    if (a->score1 && dir && ((flag&15) != 2 || a->score1 >= filters)) {
        int i = a->ref_end1, j = a->read_end1, state = 0, count = 0;
        uint32_t* ops = (uint32_t*)malloc(((size_t)i+j+2)*sizeof(uint32_t));
        if (!ops) { ok = 0; goto done; }
        while (i >= 0 && j >= 0) {
            uint8_t d = dir[(size_t)i*nq+j];
            int move = state ? state : d&3;
            if (!move) break;
            if (move == 1) { sg_push(ops, &count, 1, 0); --i; --j; state = 0; }
            else if (move == 2) { sg_push(ops, &count, 1, 2); --i; state = (d>>2)&3; }
            else { sg_push(ops, &count, 1, 1); --j; state = (d>>4)&3; }
        }
        if (i < 0 && j >= 0 && clip-go-(int64_t)j*ge > 0) { sg_push(ops, &count, j+1, 1); j = -1; }
        a->read_begin1 = j+1; a->ref_begin1 = i+1;
        if ((flag&7) && (!(flag&2) || a->score1 >= filters) && (!(flag&4) || (a->ref_end1-a->ref_begin1 <= filterd && a->read_end1-a->read_begin1 <= filterd))) {
            for (int k = 0; k < count/2; ++k) { uint32_t t = ops[k]; ops[k] = ops[count-1-k]; ops[count-1-k] = t; }
            a->cigar = ops; a->cigarLen = count;
        } else free(ops);
    }
done:
    free(h); free(e); free(m); free(ins); free(dir);
    return ok;
}

s_align* ssw_align_clip(const s_profile* p, const int8_t* ref, int32_t nr,
        uint8_t go, uint8_t ge, uint8_t clip, uint8_t flag, uint16_t filters, int32_t filterd, int32_t mask) {
    if (!clip) return ssw_align(p, ref, nr, go, ge, flag, filters, filterd, mask);
    if (!p || !ref || nr <= 0) { errno = EINVAL; return NULL; }
    s_align* a = (s_align*)calloc(1, sizeof(*a));
    if (!a) return NULL;
    a->ref_begin1 = a->read_begin1 = a->ref_end1 = a->read_end1 = -1;
    a->ref_end2 = mask >= 15 ? 0 : -1;
    int* columns = mask >= 15 ? (int*)malloc((size_t)nr*sizeof(int)) : NULL;
    if (flag & SSW_REPORT_PREFIX_SCORES) a->prefix_scores = (uint16_t*)calloc(p->readLen, sizeof(uint16_t));
    if ((mask >= 15 && !columns) || ((flag & SSW_REPORT_PREFIX_SCORES) && !a->prefix_scores)) goto error;
    int exact = 0, word = !p->profile_byte, stop = -1, match_max = 0;
    for (int i = 0; i < p->n*p->n; ++i) if (p->mat[i] > match_max) match_max = p->mat[i];
    int64_t ceiling = (int64_t)p->readLen*match_max+clip;
    int full_byte = ceiling <= 255 && ceiling+p->bias >= 255;
    if (!columns && !a->prefix_scores) {
        /* No later endpoint can improve a perfect full-query alignment. */
        if (ceiling-clip <= 65535) stop = p->readLen*match_max;
    }
    if (go >= ge) {
        if (!word) exact = full_byte ? cp_byte_full(p, ref, nr, 0, go, ge, clip, stop, a, columns) : cp_byte(p, ref, nr, 0, go, ge, clip, stop, a, columns);
        if (!exact) { word = 1; exact = cp_word(p, ref, nr, 0, go, ge, clip, stop, a, columns); }
        if (exact < 0) goto error;
    }
    if (!exact) {
scalar:
        if (!cp_scalar(p, ref, nr, go, ge, clip, flag, filters, filterd, a, columns)) goto error;
    } else if (a->score1 && (flag&15) && ((flag&15) != 2 || a->score1 >= filters)) {
        int tail_cost = a->read_end1 == p->readLen-1 ? 0 : clip;
        int8_t* reversed = seq_reverse(p->read, a->read_end1);
        s_profile* reverse = ssw_init(reversed, a->read_end1+1, p->mat, p->n, word ? 1 : 0);
        s_align back = {0};
        int ok = 0;
        if (!word) ok = full_byte ? cp_byte_full(reverse, ref, a->ref_end1+1, 1, go, ge, clip, a->score1+tail_cost, &back, NULL) : cp_byte(reverse, ref, a->ref_end1+1, 1, go, ge, clip, a->score1+tail_cost, &back, NULL);
        if (!ok) ok = cp_word(reverse, ref, a->ref_end1+1, 1, go, ge, clip, a->score1+tail_cost, &back, NULL);
        init_destroy(reverse); free(reversed);
        if (ok < 0) goto error;
        if (!ok) goto scalar;
        /* Earliest forward endpoints ensure this reversed optimum is anchored
         * to the chosen end. A clipped forward tail contributes no reverse cost. */
        a->ref_begin1 = a->ref_end1-back.ref_end1;
        a->read_begin1 = a->read_end1-back.read_end1;
        if ((flag&7) && (!(flag&2) || a->score1 >= filters) && (!(flag&4) || (a->ref_end1-a->ref_begin1 <= filterd && a->read_end1-a->read_begin1 <= filterd))) {
            s_profile span = *p;
            span.read += a->read_begin1; span.readLen = a->read_end1-a->read_begin1+1;
            s_align_semiglobal trace = {0};
            trace.score1 = a->score1+tail_cost+(a->read_begin1 ? clip : 0);
            if (!sg_trace(&span, ref+a->ref_begin1, a->ref_end1-a->ref_begin1+1, go, ge, &trace)) goto error;
            a->cigar = trace.cigar; a->cigarLen = trace.cigarLen;
        }
    }
    if (columns) for (int i = 0; i < nr; ++i) {
        if ((i < a->ref_end1-(int64_t)mask || i > a->ref_end1+(int64_t)mask) && columns[i] > a->score2) { a->score2 = columns[i]; a->ref_end2 = i; }
    }
    free(columns);
    return a;
error:
    free(columns); align_destroy(a);
    return NULL;
}
