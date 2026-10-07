// Compile the implementation matching ssw_selected.h, with the same target flags.
#ifdef __AVX2__
#include "ssw_avx2.c"
#else
#include "ssw.c"
#endif
