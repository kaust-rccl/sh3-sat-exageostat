/**
 *
 * @file cuda_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgeadd GPU kernel
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2015-09-17
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

#if !defined(CHAMELEON_USE_CUBLAS_V2)
#error "This file requires cublas api v2 support"
#endif

/**
 ******************************************************************************
 *
 * @ingroup CUDA_MORSE_Complex64_t
 *
 *  CUDA_zgeadd adds to matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = MorseNoTrans:   op(A) = A
 *          = MorseTrans:     op(A) = A'
 *          = MorseConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = MorseNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = MorseNoTrans, and k=N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
int CUDA_zgeadd(MORSE_enum trans,
                int m, int n,
                const cuDoubleComplex *alpha,
                const cuDoubleComplex *A, int lda,
                const cuDoubleComplex *beta,
                cuDoubleComplex *B, int ldb,
                CUBLAS_STREAM_PARAM)
{
    cublasZgeam(CUBLAS_HANDLE
                morse_cublas_const(trans), morse_cublas_const(MorseNoTrans),
                m, n,
                CUBLAS_VALUE(alpha), A, lda,
                CUBLAS_VALUE(beta),  B, ldb,
                B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
