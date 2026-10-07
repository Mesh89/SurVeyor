#ifndef SURVEYOR_SSW_SELECTED_H
#define SURVEYOR_SSW_SELECTED_H

#ifdef __AVX2__
#include "ssw_avx2.h"
#ifdef __cplusplus
#include "ssw_avx2_cpp.h"
#endif
#else
#include "ssw.h"
#ifdef __cplusplus
#include "ssw_cpp.h"
#endif
#endif

#endif
