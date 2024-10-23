/**
 *
 * @file core_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztpmlqt CPU kernel
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 * CORE_ztpmlqt applies a complex orthogonal matrix Q obtained from a
 * "triangular-pentagonal" complex block reflector H to a general complex matrix
 * C, which consists of two blocks A and B.
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
 * @param[in] M
 *         The number of rows of the tile B. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile B. N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of V.
 *          K >= L >= 0.  See Further Details.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTPQRT in the first k rows of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] A
 *         A is COMPLEX*16 array, dimension (LDA,N) if side = MorseLeft
 *         or (LDA,K) if SIDE = MorseRight
 *         On entry, the K-by-N or M-by-K matrix A.
 *         On exit, A is overwritten by the corresponding block of
 *         Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
 *
 * @param[in] LDA
 *         The leading dimension of the array A. LDA >= max(1,M).
 *         If side = MorseLeft,  LDA >= max(1,K);
 *         If side = Morseright, LDA >= max(1,M).
 *
 * @param[in,out] B
 *         On entry, the M-by-N tile B.
 *         On exit, B is overwritten by the corresponding block of
 *         Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
 *
 * @param[in] LDB
 *         The leading dimension of the tile B. LDB >= max(1,M).
 *
 * @param[out] WORK
 *         Workspace array of size LDWORK-by-NB.
 *         LDWORK = N if side = MorseLeft, or  M if side = MorseRight.
 *
 *******************************************************************************
 *
 * @par Further Details:
 * =====================
 *
 *  The columns of the pentagonal matrix V contain the elementary reflectors
 *  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
 *  trapezoidal block V2:
 *
 *        V = [V1] [V2].
 *
 *  The size of the trapezoidal block V2 is determined by the parameter L,
 *  where 0 <= L <= K; V2 is lower trapezoidal, consisting of the first L
 *  rows of a K-by-K upper triangular matrix.  If L=K, V2 is lower triangular;
 *  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
 *
 *  If side = MorseLeft:  C = [A]  where A is K-by-N,  B is M-by-N and V is K-by-M.
 *                            [B]
 *
 *  If side = MorseRight: C = [A B]  where A is M-by-K, B is M-by-N and V is K-by-N.
 *
 *  The complex orthogonal matrix Q is formed from V and T.
 *
 *  If trans='N' and side='L', C is on exit replaced with Q * C.
 *
 *  If trans='C' and side='L', C is on exit replaced with Q**H * C.
 *
 *  If trans='N' and side='R', C is on exit replaced with C * Q.
 *
 *  If trans='C' and side='R', C is on exit replaced with C * Q**H.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_ztpmlqt( MORSE_enum side, MORSE_enum trans,
                  int M, int N, int K, int L, int IB,
                  const MORSE_Complex64_t *V, int LDV,
                  const MORSE_Complex64_t *T, int LDT,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB,
                  MORSE_Complex64_t *WORK )
{
    int m1, n1, ldwork;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        coreblas_error(1, "Illegal value of side");
        return -1;
    }

    if ( side == MorseLeft ) {
        m1 = K;
        n1 = N;
        ldwork = IB;
    }
    else {
        m1 = M;
        n1 = K;
        ldwork = chameleon_max( n1, N );
    }

    /* TS case */
    if (L == 0) {
        CORE_ztsmlq( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, ldwork );
    }
    /* TT case */
    else if( L == N ) {
        CORE_zttmlq( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, ldwork );
    }
    else {
        //LAPACKE_ztpmlqt_work( LAPACK_COL_MAJOR, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK );
        coreblas_error( 6, "Illegal value of L (only 0 or M handled for now)");
        return -6;
    }

    return MORSE_SUCCESS;
}
