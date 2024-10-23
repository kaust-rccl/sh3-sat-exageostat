/**
 *
 * @file core_ztsmlq_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztsmlq_hetra1 CPU kernel
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <coreblas/lapacke.h>
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 * CORE_ztsmlq_hetra1: see CORE_ztsmlq
 *
 * This kernel applies a Right transformation on | A1' A2 |
 * and does not handle the transpose of A1.
 * Needs therefore to make the explicit transpose of A1 before
 * and after the application of the block of reflectors
 * Can be further optimized by changing accordingly the underneath
 * kernel ztsrfb!
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg MorseLeft  : apply Q or Q**H from the Left;
 *         @arg MorseRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg MorseNoTrans   :  No transpose, apply Q;
 *         @arg MorseConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] m1
 *         The number of rows of the tile A1. m1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. n1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. m2 >= 0.
 *         m2 = m1 if side == MorseRight.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. n2 >= 0.
 *         n2 = n1 if side == MorseLeft.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A1
 *         On entry, the m1-by-n1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,m1).
 *
 * @param[in,out] A2
 *         On entry, the m2-by-n2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTSLQT in the first k rows of its array argument V.
 *!
 * @param[in] ldv
 *         The leading dimension of the array V. ldv >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-n1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-m1 if side == MorseLeft
 *             LDWORK-by-IB if side == MorseRight
 *
 * @param[in] ldwork
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == MorseLeft
 *             LDWORK >= max(1,n1) if side == MorseRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
#if defined(MORSE_HAVE_WEAK)
#pragma weak CORE_ztsmlq_hetra1 = PCORE_ztsmlq_hetra1
#define CORE_ztsmlq_hetra1 PCORE_ztsmlq_hetra1
#define CORE_ztsmlq PCORE_ztsmlq
int  CORE_ztsmlq(MORSE_enum side, MORSE_enum trans,
                 int m1, int n1, int m2, int n2, int K, int IB,
                 MORSE_Complex64_t *A1, int lda1,
                 MORSE_Complex64_t *A2, int lda2,
                 const MORSE_Complex64_t *V, int ldv,
                 const MORSE_Complex64_t *T, int ldt,
                 MORSE_Complex64_t *WORK, int LDWORK);
#endif
int CORE_ztsmlq_hetra1( MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        MORSE_Complex64_t *A1, int lda1,
                        MORSE_Complex64_t *A2, int lda2,
                        const MORSE_Complex64_t *V, int ldv,
                        const MORSE_Complex64_t *T, int ldt,
                        MORSE_Complex64_t *WORK, int ldwork)
{
    int i, j;

    if ( (m1 != n1) ) {
        coreblas_error(3, "Illegal value of M1, N1");
        return -3;
    }

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
#if defined(PRECISION_z) || defined(PRECISION_c)
        A1[j + j*lda1] = conj(A1[j + j*lda1]);
#endif
        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = conj(*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = conj(*WORK);
        }
    }

    CORE_ztsmlq(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2,
                V,  ldv,  T,  ldt,
                WORK, ldwork);

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
#if defined(PRECISION_z) || defined(PRECISION_c)
        A1[j + j*lda1] = conj(A1[j + j*lda1]);
#endif
        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = conj(*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = conj(*WORK);
        }
    }

    return MORSE_SUCCESS;
}
