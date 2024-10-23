/**
 *
 * @file cuda_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zhemm GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c
 *
 */
#include "cudablas.h"

int CUDA_zhemm(MORSE_enum side, MORSE_enum uplo,
               int m, int n,
               cuDoubleComplex *alpha,
               const cuDoubleComplex *A, int lda,
               const cuDoubleComplex *B, int ldb,
               cuDoubleComplex *beta,
               cuDoubleComplex *C, int ldc,
               CUBLAS_STREAM_PARAM)
{
    cublasZhemm(CUBLAS_HANDLE
                morse_cublas_const(side), morse_cublas_const(uplo),
                m, n,
                CUBLAS_VALUE(alpha), A, lda,
                                     B, ldb,
                CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );
    return MORSE_SUCCESS;
}
