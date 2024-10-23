/**
 *
 * @file cuda_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztrmm GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_ztrmm(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUBLAS_STREAM_PARAM)
{

#if defined(CHAMELEON_USE_CUBLAS_V2)

    cublasZtrmm(
        CUBLAS_HANDLE
        morse_cublas_const(side), morse_cublas_const(uplo),
        morse_cublas_const(transa), morse_cublas_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
        B, ldb,
        B, ldb);

#else

    cublasZtrmm(
        CUBLAS_HANDLE
        morse_cublas_const(side), morse_cublas_const(uplo),
        morse_cublas_const(transa), morse_cublas_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
                             B, ldb);
#endif

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}

