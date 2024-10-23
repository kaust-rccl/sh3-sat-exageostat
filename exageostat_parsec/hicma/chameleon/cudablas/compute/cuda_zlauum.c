/**
 *
 * @file cuda_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlauum GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

#if defined(CHAMELEON_USE_MAGMA)
int CUDA_zlauum(
        char uplo, magma_int_t n,
        cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info)
{
    int ret;
    ret = magma_zlauum_gpu( uplo, n, dA, ldda, info);
    if (ret != MAGMA_SUCCESS) {
        fprintf(stderr, "Error in MAGMA: %d\n", ret);
        exit(-1);
    }
    return MORSE_SUCCESS;
}
#endif
