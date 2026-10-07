#ifndef SIMD_MACROS_H_
#define SIMD_MACROS_H_

#include <cstdint>
#include <cstring>

#if !defined(USE_SCALAR) && !defined(USE_SSE) && !defined(USE_AVX2) && !defined(USE_NEON) //&& !defined(USE_AVX512)
  // #if defined(__AVX512F__) && defined(__AVX512BW__)  CURRENTLY I HAVE NO WAY TO TEST THIS AS I HAVE NO ACCESS TO AVX512 HARDWARE
    // AVX512BW is required for byte/word (16-bit) operations
    // #define USE_AVX512
  #if defined(__AVX2__)
    #define USE_AVX2
  #elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define USE_NEON
  #elif defined(__SSE4_2__)
    #define USE_SSE
  #else
    #define USE_SCALAR
  #endif
#endif

#ifdef USE_SCALAR

typedef int16_t SIMD_INT_16;

static inline SIMD_INT_16 load_int_unaligned_16(const void* p) {
    SIMD_INT_16 v;
    std::memcpy(&v, p, sizeof v); 
    return v;
}

#define LOADU_INT_16(p) load_int_unaligned_16(p)

// Count bytes macro
#define COUNT_EQUAL_BYTES_16(a, b) (\
    ((a & 0xFF00) == (b & 0xFF00)) + \
    ((a & 0x00FF) == (b & 0x00FF)) )

#define BYTES_PER_BLOCK_16 (sizeof(SIMD_INT_16))

#elif defined(USE_NEON)
#include <arm_neon.h>

typedef int16x8_t SIMD_INT_16;

static inline SIMD_INT_16 load_int_16_neon(const void* p) {
    return vreinterpretq_s16_u8(vld1q_u8(reinterpret_cast<const uint8_t*>(p)));
}

static inline int count_equal_bytes_16_neon(SIMD_INT_16 a, SIMD_INT_16 b) {
    uint8x16_t eq = vceqq_u8(vreinterpretq_u8_s16(a), vreinterpretq_u8_s16(b));
    uint8x16_t ones = vshrq_n_u8(eq, 7);
    uint16x8_t sum16 = vpaddlq_u8(ones);
    uint32x4_t sum32 = vpaddlq_u16(sum16);
    uint64x2_t sum64 = vpaddlq_u32(sum32);
    return (int)(vgetq_lane_u64(sum64, 0) + vgetq_lane_u64(sum64, 1));
}

#define LOADU_INT_16(p) load_int_16_neon(p)

#define COUNT_EQUAL_BYTES_16(a, b) count_equal_bytes_16_neon((a), (b))

#define BYTES_PER_BLOCK_16 16

#elif defined(USE_SSE)
#include <immintrin.h>

typedef __m128i SIMD_INT_16;

#define LOADU_INT_16 _mm_loadu_si128

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u32(_mm_movemask_epi8(_mm_cmpeq_epi8(a, b)))

#define BYTES_PER_BLOCK_16 16

#elif defined(USE_AVX2)
#include <immintrin.h>

typedef __m256i SIMD_INT_16;

#define LOADU_INT_16 _mm256_loadu_si256

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)))

#define BYTES_PER_BLOCK_16 32

#elif defined(USE_AVX512)
#include <immintrin.h>

typedef __m512i SIMD_INT_16;

#define LOADU_INT_16 _mm512_loadu_si512

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(a, b))

#define BYTES_PER_BLOCK_16 64

#else
#error "Either USE_SCALAR, USE_NEON, USE_SSE, USE_AVX2, or USE_AVX512 must be defined"
#endif

#endif /* SIMD_MACROS_H_ */
