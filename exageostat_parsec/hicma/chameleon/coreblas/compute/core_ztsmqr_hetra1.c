/**
 *
 * @file core_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztsmqr_hetra1 CPU kernel
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
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
 *  CORE_ztsmqr_hetra1: see CORE_ztsmqr
 *
 *  This kernel applies a left transformation on | A1'|
 *                                               | A2 |
 *
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
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. M2 >= 0.
 *         M2 = M1 if side == MorseRight.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. N2 >= 0.
 *         N2 = N1 if side == MorseLeft.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTSQRT in the first k columns of its array argument V.
 *
 * @param[in] ldv
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-N1 if side == MorseLeft
 *             LDWORK-by-IB if side == MorseRight
 *
 * @param[in] ldwork
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == MorseLeft
 *             LDWORK >= max(1,M1) if side == MorseRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
#if defined(MORSE_HAVE_WEAK)
#pragma weak CORE_ztsmqr_hetra1 = PCORE_ztsmqr_hetra1
#define CORE_ztsmqr_hetra1 PCORE_ztsmqr_hetra1
#define CORE_ztsmqr PCORE_ztsmqr
int  CORE_ztsmqr(MORSE_enum side, MORSE_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK, int LDWORK);
#endif
int CORE_ztsmqr_hetra1( MORSE_enum side, MORSE_enum trans,
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

    CORE_ztsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

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
