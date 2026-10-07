#ifndef SSW_AVX2_SIMD_H
#define SSW_AVX2_SIMD_H
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>
#ifndef __AVX2__
#error "Compile ssw_avx2.c with AVX2 enabled (-mavx2)."
#endif

/* AVX2 byte shifts are lane-local: explicitly carry the low half into high. */
static inline __m256i ssw_avx2_shift1(__m256i x) {
    return _mm256_alignr_epi8(x, _mm256_permute2x128_si256(x,x,0x08),15);
}
static inline __m256i ssw_avx2_shift2(__m256i x) {
    return _mm256_alignr_epi8(x, _mm256_permute2x128_si256(x,x,0x08),14);
}
static inline __m256i ssw_avx2_shift4(__m256i x) {
    return _mm256_alignr_epi8(x, _mm256_permute2x128_si256(x,x,0x08),12);
}
static inline __m256i ssw_avx2_shift8(__m256i x) {
    return _mm256_alignr_epi8(x, _mm256_permute2x128_si256(x,x,0x08),8);
}
static inline __m256i ssw_avx2_shift16(__m256i x) {
    return _mm256_permute2x128_si256(x,x,0x08);
}
/* Exact closure of F[lane] -> F[lane+1] - segLen*gap_extend.
 * Logarithmic prefix maximum replaces up to 31/15 serial lazy sweeps. */
static inline __m256i ssw_avx2_prefix_byte(__m256i f,int cost) {
#define STEP(S,D) f=_mm256_max_epu8(f,_mm256_subs_epu8(ssw_avx2_shift##S(f),_mm256_set1_epi8(cost*(D)>255?255:cost*(D))))
    STEP(1,1);STEP(2,2);STEP(4,4);STEP(8,8);STEP(16,16);
#undef STEP
    return f;
}
static inline __m256i ssw_avx2_prefix_word(__m256i f,int cost) {
#define STEP(S,D) f=_mm256_max_epi16(f,_mm256_subs_epu16(ssw_avx2_shift##S(f),_mm256_set1_epi16(cost*(D)>32767?32767:cost*(D))))
    STEP(2,1);STEP(4,2);STEP(8,4);STEP(16,8);
#undef STEP
    return f;
}
static inline uint8_t ssw_avx2_hmax_byte(__m256i x) {
    __m128i v=_mm_max_epu8(_mm256_castsi256_si128(x),_mm256_extracti128_si256(x,1));
    v=_mm_max_epu8(v,_mm_srli_si128(v,8));
    v=_mm_max_epu8(v,_mm_srli_si128(v,4));
    v=_mm_max_epu8(v,_mm_srli_si128(v,2));
    v=_mm_max_epu8(v,_mm_srli_si128(v,1));
    return (uint8_t)_mm_extract_epi16(v,0);
}
static inline uint16_t ssw_avx2_hmax_word(__m256i x) {
    __m128i v=_mm_max_epi16(_mm256_castsi256_si128(x),_mm256_extracti128_si256(x,1));
    v=_mm_max_epi16(v,_mm_srli_si128(v,8));
    v=_mm_max_epi16(v,_mm_srli_si128(v,4));
    v=_mm_max_epi16(v,_mm_srli_si128(v,2));
    return (uint16_t)_mm_extract_epi16(v,0);
}
static inline void* ssw_avx2_alloc(size_t count, int zero) {
    void* p=NULL;
    if (posix_memalign(&p,32,count*sizeof(__m256i))!=0) return NULL;
    if (zero) memset(p,0,count*sizeof(__m256i));
    return p;
}
#endif
