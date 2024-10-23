/**
 *
 * @file cuda_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztpmqrt GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int
CUDA_ztpmqrt( MORSE_enum side, MORSE_enum trans,
              int M, int N, int K, int L, int IB,
              const cuDoubleComplex *V, int LDV,
              const cuDoubleComplex *T, int LDT,
                    cuDoubleComplex *A, int LDA,
                    cuDoubleComplex *B, int LDB,
                    cuDoubleComplex *WORK,
              CUBLAS_STREAM_PARAM )
{
    int m1, n1, ldwork, ldworkc, ws;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        cudablas_error(1, "Illegal value of side");
        return -1;
    }

    if ( side == MorseLeft ) {
        m1 = K;
        n1 = N;
        ldwork  = IB;
        ldworkc = M;
        ws = ldwork * n1;
    }
    else {
        m1 = M;
        n1 = K;
        ldwork  = m1;
        ldworkc = IB;
        ws = ldwork * IB;
    }

    /* TS case */
    if (L == 0) {
        CUDA_ztsmqr( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, ldwork, WORK + ws, ldworkc,
                     CUBLAS_STREAM_VALUE );
    }
    /* TT case */
    else  if( L == M ) {
        CUDA_zttmqr( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, ldwork, WORK + ws, ldworkc,
                     CUBLAS_STREAM_VALUE );
    }
    else {
        cudablas_error(-6, "TPMQRT not available on GPU for general cases yet\n" );
        return -6;
    }

    return MORSE_SUCCESS;
}
