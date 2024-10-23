/**
 *
 * @file core_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztradd CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
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
 *  CORE_ztradd adds to matrices together as in PBLAS pztradd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          = MorseUpperLower: A and B are general matrices.
 *          = MorseUpper: op(A) and B are upper trapezoidal matrices.
 *          = MorseLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = MorseNoTrans:   op(A) = A
 *          = MorseTrans:     op(A) = A'
 *          = MorseConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
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
#pragma weak CORE_ztradd = PCORE_ztradd
#define CORE_ztradd PCORE_ztradd
#define CORE_zgeadd PCORE_zgeadd
int
CORE_zgeadd(MORSE_enum trans, int M, int N,
                  MORSE_Complex64_t alpha,
            const MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t beta,
                  MORSE_Complex64_t *B, int LDB);
#endif
int CORE_ztradd(MORSE_enum uplo, MORSE_enum trans, int M, int N,
                      MORSE_Complex64_t  alpha,
                const MORSE_Complex64_t *A, int LDA,
                      MORSE_Complex64_t  beta,
                      MORSE_Complex64_t *B, int LDB)
{
    int i, j;

    if (uplo == MorseUpperLower){
        int rc = CORE_zgeadd( trans, M, N, alpha, A, LDA, beta, B, LDB );
        if (rc != MORSE_SUCCESS)
            return rc-1;
        else
            return rc;
    }

    if ((uplo != MorseUpper) &&
        (uplo != MorseLower))
    {
        coreblas_error(1, "illegal value of uplo");
        return -1;
    }

    if ((trans < MorseNoTrans) || (trans > MorseConjTrans))
    {
        coreblas_error(2, "illegal value of trans");
        return -2;
    }

    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ( ((trans == MorseNoTrans) && (LDA < chameleon_max(1,M)) && (M > 0)) ||
         ((trans != MorseNoTrans) && (LDA < chameleon_max(1,N)) && (N > 0)) )
    {
        coreblas_error(7, "Illegal value of LDA");
        return -7;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(9, "Illegal value of LDB");
        return -9;
    }

    /**
     * MorseLower
     */
    if (uplo == MorseLower) {
        switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        case MorseConjTrans:
            for (j=0; j<N; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-M+j+1;
            }
            break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

        case MorseTrans:
            for (j=0; j<N; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * A[LDA*i];
                }
                B += LDB-M+j+1;
            }
            break;

        case MorseNoTrans:
        default:
            for (j=0; j<N; j++) {
                for(i=j; i<M; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                B += LDB-M+j+1;
                A += LDA-M+j+1;
            }
        }
    }
    /**
     * MorseUpper
     */
    else {
        switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        case MorseConjTrans:
            for (j=0; j<N; j++, A++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-mm;
            }
            break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

        case MorseTrans:
            for (j=0; j<N; j++, A++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B = beta * (*B) + alpha * (A[LDA*i]);
                }
                B += LDB-mm;
            }
            break;

        case MorseNoTrans:
        default:
            for (j=0; j<N; j++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                B += LDB-mm;
                A += LDA-mm;
            }
        }
    }
    return MORSE_SUCCESS;
}
