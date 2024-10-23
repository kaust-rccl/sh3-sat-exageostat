/**
 *
 * @file cuda_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherk GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c
 *
 */
#include "cudablas.h"

int CUDA_zherk( MORSE_enum uplo, MORSE_enum trans,
                int n, int k,
                double *alpha,
                const cuDoubleComplex *A, int lda,
                double *beta,
                cuDoubleComplex *B, int ldb,
                CUBLAS_STREAM_PARAM)
{
    cublasZherk( CUBLAS_HANDLE
                 morse_cublas_const(uplo), morse_cublas_const(trans),
                 n, k,
                 CUBLAS_VALUE(alpha), A, lda,
                 CUBLAS_VALUE(beta),  B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
