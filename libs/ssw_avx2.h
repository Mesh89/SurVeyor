#ifndef SSW_AVX2_H
#define SSW_AVX2_H

/* Shares the result ABI and CIGAR conventions with the bundled SSW API.
 * Select this header instead of ssw.h; compile/link ssw_avx2.c with -mavx2.
 * Define SSW_AVX2_NO_REMAP to access both C backends in the same translation unit.
 * Positive sequence lengths and encoded residues in [0,n) are required.
 * Unusual gap_open < gap_extend uses exact scalar affine DP.
 * Scores above 65535 return NULL with errno=ERANGE
 * (the shared result ABI cannot represent them); invalid arguments: EINVAL.
 * Query/matrix storage must outlive the profile, as with SSW.
 */
#include "ssw.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct _ssw_avx2_profile ssw_avx2_profile;
extern const uint8_t ssw_avx2_encoded_ops[];
ssw_avx2_profile* ssw_avx2_init(const int8_t*, int32_t, const int8_t*, int32_t, int8_t);
void ssw_avx2_init_destroy(ssw_avx2_profile*);
s_align* ssw_avx2_align(const ssw_avx2_profile*, const int8_t*, int32_t,
                      uint8_t, uint8_t, uint8_t, uint16_t, int32_t, int32_t);
void ssw_avx2_align_destroy(s_align*);
int32_t ssw_avx2_mark_mismatch(int32_t, int32_t, int32_t, const int8_t*,
                             const int8_t*, int32_t, uint32_t**, int32_t*);
static inline uint32_t ssw_avx2_to_cigar_int(uint32_t len, unsigned char op) {
    return (len << BAM_CIGAR_SHIFT) | ssw_avx2_encoded_ops[op];
}
#ifdef __cplusplus
}
#endif
#ifndef SSW_AVX2_NO_REMAP
#define s_profile ssw_avx2_profile
#define ssw_init ssw_avx2_init
#define ssw_align ssw_avx2_align
#define init_destroy ssw_avx2_init_destroy
#define align_destroy ssw_avx2_align_destroy
#define mark_mismatch ssw_avx2_mark_mismatch
#define to_cigar_int ssw_avx2_to_cigar_int
#define encoded_ops ssw_avx2_encoded_ops
#endif
#endif
