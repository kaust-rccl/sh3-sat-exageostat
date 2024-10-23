/**
 *
 * @file core_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgeadd CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @date 2015-11-03
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 ******************************************************************************
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgeadd adds to matrices together as in PBLAS pzgeadd.
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
#if defined(MORSE_HAVE_WEAK)
#pragma weak CORE_zgeadd = PCORE_zgeadd
#define CORE_zgeadd PCORE_zgeadd
#endif
int CORE_zgeadd(MORSE_enum trans, int M, int N,
                      MORSE_Complex64_t  alpha,
                const MORSE_Complex64_t *A, int LDA,
                      MORSE_Complex64_t  beta,
                      MORSE_Complex64_t *B, int LDB)
{
    int i, j;

    if ((trans < MorseNoTrans) || (trans > MorseConjTrans))
    {
        coreblas_error(1, "illegal value of trans");
        return -1;
    }

    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if ( ((trans == MorseNoTrans) && (LDA < chameleon_max(1,M)) && (M > 0)) ||
         ((trans != MorseNoTrans) && (LDA < chameleon_max(1,N)) && (N > 0)) )
    {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case MorseConjTrans:
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B = beta * (*B) + alpha * conj(A[LDA*i]);
            }
            B += LDB-M;
        }
        break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    case MorseTrans:
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B = beta * (*B) + alpha * A[LDA*i];
            }
            B += LDB-M;
        }
        break;

    case MorseNoTrans:
    default:
        for (j=0; j<N; j++) {
            for(i=0; i<M; i++, B++, A++) {
                *B = beta * (*B) + alpha * (*A);
            }
            A += LDA-M;
            B += LDB-M;
        }
    }
    return MORSE_SUCCESS;
}
