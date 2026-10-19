/* Shared by the SSE2 and AVX2 translation units; local SW is unchanged. */
#include <limits.h>
#include <errno.h>

#ifdef SG_AVX2
#define SG_VEC __m256i
#define SG_BYTES 32
#define SG_ALLOC(n) ssw_avx2_alloc(n, 0)
#define SG_ZERO _mm256_setzero_si256
#define SG_OR _mm256_or_si256
#define SG_LOW(x) _mm256_inserti128_si256(_mm256_setzero_si256(), _mm_cvtsi32_si128(x), 0)
#define SG_SHIFT ssw_avx2_shift1
#define SG_SET(x) _mm256_set1_epi8((char)(x))
#define SG_SUB _mm256_sub_epi8
#define SG_DEC _mm256_subs_epu8
#define SG_MAX _mm256_max_epu8
#define SG_PREFIX ssw_avx2_prefix_byte
#else
#define SG_VEC __m128i
#define SG_BYTES 16
#define SG_ALLOC(n) malloc((n)*sizeof(__m128i))
#define SG_ZERO _mm_setzero_si128
#define SG_OR _mm_or_si128
#define SG_LOW _mm_cvtsi32_si128
#define SG_SHIFT(x) _mm_slli_si128(x, 1)
#define SG_SET(x) _mm_set1_epi8((char)(x))
#define SG_SUB _mm_sub_epi8
#define SG_DEC _mm_subs_epu8
#define SG_MAX _mm_max_epu8
#define SG_PREFIX ssw_prefix_byte
#endif
#define SG_INT uint8_t
#define SG_LIMIT 255
#define SG_KERNEL sg_byte
#define SG_PROFILE p->profile_byte
#define SG_MAKE_PROFILE qP_byte(p->read, p->mat, nq, p->n, p->bias)
#define SG_PROFILE_BIAS (match_max+p->bias)
#include "ssw_semiglobal_simd.h"
#undef SG_INT
#undef SG_LIMIT
#undef SG_KERNEL
#undef SG_PROFILE
#undef SG_MAKE_PROFILE
#undef SG_PROFILE_BIAS
#undef SG_SHIFT
#undef SG_SET
#undef SG_SUB
#undef SG_DEC
#undef SG_MAX
#undef SG_PREFIX

#ifdef SG_AVX2
#define SG_SHIFT ssw_avx2_shift2
#define SG_SET(x) _mm256_set1_epi16(x)
#define SG_SUB _mm256_sub_epi16
#define SG_DEC _mm256_subs_epu16
#define SG_MAX _mm256_max_epi16
#define SG_PREFIX ssw_avx2_prefix_word
#else
#define SG_SHIFT(x) _mm_slli_si128(x, 2)
#define SG_SET(x) _mm_set1_epi16(x)
#define SG_SUB _mm_sub_epi16
#define SG_DEC _mm_subs_epu16
#define SG_MAX _mm_max_epi16
#define SG_PREFIX ssw_prefix_word
#endif
#define SG_INT uint16_t
#define SG_LIMIT 32767
#define SG_KERNEL sg_word
#define SG_PROFILE p->profile_word
#define SG_MAKE_PROFILE qP_word(p->read, p->mat, nq, p->n)
#define SG_PROFILE_BIAS match_max
#include "ssw_semiglobal_simd.h"
#undef SG_INT
#undef SG_LIMIT
#undef SG_KERNEL
#undef SG_PROFILE
#undef SG_MAKE_PROFILE
#undef SG_PROFILE_BIAS
#undef SG_SHIFT
#undef SG_SET
#undef SG_SUB
#undef SG_DEC
#undef SG_MAX
#undef SG_PREFIX
#undef SG_VEC
#undef SG_BYTES
#undef SG_ALLOC
#undef SG_ZERO
#undef SG_OR
#undef SG_LOW

void align_destroy_semiglobal(s_align_semiglobal* a) {
    if (!a) return;
    free(a->prefix_scores);
    free(a->cigar);
    free(a);
}

static void sg_push(uint32_t* ops, int* count, int len, int op) {
    if (!len) return;
    if (*count && (ops[*count-1]&15) == (unsigned)op) ops[*count-1] += (uint32_t)len<<4;
    else ops[(*count)++] = ((uint32_t)len<<4)|op;
}

/* Exact fallback for unusual gap penalties or exhausted word precision.
 * H/M/I are separate when opening a gap is cheaper than extending it. */
static int sg_scalar(const s_profile* p, const int8_t* ref, int nr, int go, int ge,
        uint8_t flag, int filters, int filterd, s_align_semiglobal* a, int32_t* columns) {
    int nq = p->readLen, neg = INT32_MIN/2;
    int *h = (int*)malloc((nq+1)*sizeof(int)), *e = (int*)malloc((nq+1)*sizeof(int));
    int *m = go < ge ? (int*)malloc((nq+1)*sizeof(int)) : NULL;
    int *ins = go < ge ? (int*)malloc((nq+1)*sizeof(int)) : NULL;
    uint8_t* dir = flag & 15 ? (uint8_t*)malloc((size_t)nq*nr) : NULL;
    if (!h || !e || (go < ge && (!m || !ins)) || ((flag&15) && !dir)) { free(h); free(e); free(m); free(ins); free(dir); return 0; }
    h[0] = 0;
    for (int j = 1; j <= nq; ++j) {
        h[j] = -go-(j-1)*ge; e[j] = neg;
        if (m) { m[j] = neg; ins[j] = h[j]; }
        if (a->prefix_scores) a->prefix_scores[j-1] = INT32_MIN;
    }
    a->score1 = h[nq]; a->ref_end1 = -1;
    for (int i = 0; i < nr; ++i) {
        int diagonal = 0, f = neg, ml = 0, el = neg;
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
            int v = diagonal+p->mat[ref[i]*p->n+p->read[j-1]];
            if (a->prefix_scores && v > a->prefix_scores[j-1]) a->prefix_scores[j-1] = v;
            if (m) { m[j] = ml = v; ins[j] = f; el = e[j]; }
            uint8_t d = 1;
            if (e[j] > v) { v = e[j]; d = 2; }
            if (f > v) { v = f; d = 3; }
            if (dir) dir[(size_t)i*nq+j-1] = d|(ep<<2)|(fp<<4);
            h[j] = v; diagonal = old;
        }
        if (columns) columns[i] = h[nq];
        if (h[nq] > a->score1) { a->score1 = h[nq]; a->ref_end1 = i; }
    }
    if (dir && ((flag&15) != 2 || a->score1 >= filters)) {
        int i = a->ref_end1, j = nq-1, state = 0, count = 0;
        uint32_t* ops = (uint32_t*)malloc(((size_t)nq+i+2)*sizeof(uint32_t));
        if (!ops) { free(h); free(e); free(m); free(ins); free(dir); return 0; }
        while (i >= 0 && j >= 0) {
            uint8_t d = dir[(size_t)i*nq+j];
            int move = state ? state : d&3;
            if (move == 1) { sg_push(ops, &count, 1, 0); --i; --j; state = 0; }
            else if (move == 2) { sg_push(ops, &count, 1, 2); --i; state = (d>>2)&3; }
            else { sg_push(ops, &count, 1, 1); --j; state = (d>>4)&3; }
        }
        sg_push(ops, &count, j+1, 1);
        a->read_begin1 = 0; a->ref_begin1 = i+1;
        if ((flag&7) && (!(flag&2) || a->score1 >= filters) && (!(flag&4) || (a->ref_end1-a->ref_begin1 <= filterd && nq-1 <= filterd))) {
            for (int k = 0; k < count/2; ++k) { uint32_t t = ops[k]; ops[k] = ops[count-1-k]; ops[count-1-k] = t; }
            a->cigar = ops; a->cigarLen = count;
        } else free(ops);
    }
    free(h); free(e); free(m); free(ins); free(dir);
    return 1;
}

/* Global traceback between the endpoints established by the SIMD passes.
 * Widen until the endpoint reproduces the known optimum, including negatives. */
static int sg_trace(const s_profile* p, const int8_t* ref, int nr, int go, int ge, s_align_semiglobal* a) {
    int nq = p->readLen, band = abs(nq-nr)+1, limit = nq > nr ? nq : nr, neg = INT32_MIN/2;
    int* h = (int*)malloc((nr+1)*sizeof(int));
    int* e = (int*)malloc((nr+1)*sizeof(int));
    uint8_t* dir = NULL;
    if (!h || !e) { free(h); free(e); return 0; }
    int width;
    for (;;) {
        if (band > limit) band = limit;
        width = 2*band+1 < nr+1 ? 2*band+1 : nr+1;
        dir = (uint8_t*)malloc((size_t)(nq+1)*width);
        if (!dir) { free(h); free(e); return 0; }
        for (int j = 0; j <= nr; ++j) { h[j] = j == 0 ? 0 : (j <= band ? -go-(j-1)*ge : neg); e[j] = neg; }
        for (int i = 1; i <= nq; ++i) {
            int beg = i-band > 1 ? i-band : 1, end = i+band < nr ? i+band : nr;
            int diagonal = h[beg-1], left = beg == 1 ? -go-(i-1)*ge : neg, f = neg;
            h[beg-1] = left;
            for (int j = beg; j <= end; ++j) {
                int old = h[j], eo = old-go, ee = e[j]-ge, fo = left-go, fe = f-ge;
                e[j] = eo > ee ? eo : ee; f = fo > fe ? fo : fe;
                int v = diagonal+p->mat[ref[j-1]*p->n+p->read[i-1]];
                uint8_t d = 1;
                if (f > v) { v = f; d = 2; }
                if (e[j] > v) { v = e[j]; d = 3; }
                dir[(size_t)i*width+j-(i > band ? i-band : 0)] = d|(eo <= ee ? 4 : 0)|(fo <= fe ? 8 : 0);
                h[j] = left = v; diagonal = old;
            }
            if (end < nr) h[end+1] = e[end+1] = neg;
        }
        if (h[nr] == a->score1) break;
        free(dir); dir = NULL;
        if (band == limit) { free(h); free(e); return 0; }
        band = band > limit/2 ? limit : band*2;
    }
    uint32_t* ops = (uint32_t*)malloc(((size_t)nq+nr+1)*sizeof(uint32_t));
    if (!ops) { free(h); free(e); free(dir); return 0; }
    int i = nq, j = nr, state = 0, count = 0;
    while (i && j) {
        uint8_t d = dir[(size_t)i*width+j-(i > band ? i-band : 0)];
        int move = state ? state : d&3;
        if (move == 1) { sg_push(ops, &count, 1, 0); --i; --j; state = 0; }
        else if (move == 2) { sg_push(ops, &count, 1, 2); --j; state = d&8 ? 2 : 0; }
        else { sg_push(ops, &count, 1, 1); --i; state = d&4 ? 3 : 0; }
    }
    sg_push(ops, &count, i, 1); sg_push(ops, &count, j, 2);
    for (int k = 0; k < count/2; ++k) { uint32_t t = ops[k]; ops[k] = ops[count-1-k]; ops[count-1-k] = t; }
    a->cigar = ops; a->cigarLen = count;
    free(h); free(e); free(dir);
    return 1;
}

s_align_semiglobal* ssw_align_semiglobal(const s_profile* p, const int8_t* ref, int32_t nr,
        uint8_t go, uint8_t ge, uint8_t flag, int32_t filters, int32_t filterd, int32_t mask) {
    if (!p || !ref || nr <= 0) { errno = EINVAL; return NULL; }
    int nq = p->readLen, match_max = 0;
    for (int i = 0; i < p->n*p->n; ++i) if (p->mat[i] > match_max) match_max = p->mat[i];
    if ((int64_t)nq*(match_max+ge)+go > INT32_MAX/2) { errno = ERANGE; return NULL; }
    s_align_semiglobal* a = (s_align_semiglobal*)calloc(1, sizeof(*a));
    if (!a) return NULL;
    a->ref_begin1 = a->read_begin1 = -1; a->read_end1 = nq-1;
    a->score2 = INT32_MIN; a->ref_end2 = -1;
    int32_t* columns = mask >= 15 ? (int32_t*)malloc((size_t)nr*sizeof(int32_t)) : NULL;
    if (flag & SSW_REPORT_PREFIX_SCORES) a->prefix_scores = (int32_t*)malloc((size_t)nq*sizeof(int32_t));
    if ((mask >= 15 && !columns) || ((flag & SSW_REPORT_PREFIX_SCORES) && !a->prefix_scores)) goto error;
    int exact = 0, word = !p->profile_byte;
    if (go >= ge) {
        if (!word) exact = sg_byte(p, ref, nr, 0, go, ge, match_max, INT32_MIN, mask, a, columns);
        if (!exact) { word = 1; exact = sg_word(p, ref, nr, 0, go, ge, match_max, INT32_MIN, mask, a, columns); }
        if (exact < 0) goto error;
    }
    if (!exact) {
        if (!sg_scalar(p, ref, nr, go, ge, flag, filters, filterd, a, columns)) goto error;
    } else if ((flag&15) && ((flag&15) != 2 || a->score1 >= filters)) {
        a->read_begin1 = 0;
        if (a->ref_end1 < 0) a->ref_begin1 = 0;
        else {
            int8_t* reversed = seq_reverse(p->read, nq-1);
            s_profile* reverse = ssw_init(reversed, nq, p->mat, p->n, word ? 1 : 0);
            s_align_semiglobal back = {0};
            int ok = word ? sg_word(reverse, ref, a->ref_end1+1, 1, go, ge, match_max, a->score1, 0, &back, NULL) : sg_byte(reverse, ref, a->ref_end1+1, 1, go, ge, match_max, a->score1, 0, &back, NULL);
            init_destroy(reverse); free(reversed);
            if (ok <= 0) goto error;
            a->ref_begin1 = a->ref_end1-back.ref_end1;
        }
        if ((flag&7) && (!(flag&2) || a->score1 >= filters) && (!(flag&4) || (a->ref_end1-a->ref_begin1 <= filterd && nq-1 <= filterd))) {
            if (a->ref_end1 < 0) {
                a->cigar = (uint32_t*)malloc(sizeof(uint32_t));
                if (!a->cigar) goto error;
                a->cigar[0] = to_cigar_int(nq, 'I'); a->cigarLen = 1;
            } else if (!sg_trace(p, ref+a->ref_begin1, a->ref_end1-a->ref_begin1+1, go, ge, a)) goto error;
        }
    }
    a->ref_end2 = -1;
    if (columns) for (int i = 0; i < nr; ++i) {
        if ((i < a->ref_end1-(int64_t)mask || i >= a->ref_end1+(int64_t)mask) && columns[i] > a->score2) { a->score2 = columns[i]; a->ref_end2 = i; }
    }
    free(columns);
    return a;
error:
    free(columns); align_destroy_semiglobal(a);
    return NULL;
}
