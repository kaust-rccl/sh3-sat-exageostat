/**
 *
 * @file cuda_zparfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zparfb GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

/**
 *****************************************************************************
 *
 * @ingroup CUDA_MORSE_Complex64_t
 *
 *  CUDA_zparfb applies a complex upper triangular block reflector H
 *  or its transpose H' to a complex rectangular matrix formed by
 *  coupling two tiles A1 and A2. Matrix V is:
 *
 *          COLUMNWISE                    ROWWISE
 *
 *         |     K     |                 |      N2-L     |   L  |
 *      __ _____________ __           __ _________________        __
 *         |    |      |                 |               | \
 *         |    |      |                 |               |   \    L
 *    M2-L |    |      |              K  |_______________|_____\  __
 *         |    |      | M2              |                      |
 *      __ |____|      |                 |                      | K-L
 *         \    |      |              __ |______________________| __
 *       L   \  |      |
 *      __     \|______| __              |          N2          |
 *
 *         | L |  K-L  |
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg MorseLeft  : apply Q or Q**H from the Left;
 *         @arg MorseRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg MorseNoTrans   : No transpose, apply Q;
 *         @arg MorseConjTrans : ConjTranspose, apply Q**H.
 *
 * @param[in] direct
 *         Indicates how H is formed from a product of elementary
 *         reflectors
 *         @arg MorseForward  : H = H(1) H(2) . . . H(k) (Forward)
 *         @arg MorseBackward : H = H(k) . . . H(2) H(1) (Backward)
 *
 * @param[in] storev
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored:
 *         @arg MorseColumnwise
 *         @arg MorseRowwise
 *
 * @param[in] M1
 *         The number of columns of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of rows of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of columns of the tile A2. M2 >= 0.
 *
 * @param[in] N2
 *         The number of rows of the tile A2. N2 >= 0.
 *
 * @param[in] K
 *         The order of the matrix T (= the number of elementary
 *         reflectors whose product defines the block reflector).
 *
 * @param[in] L
 *         The size of the triangular part of V
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,N2).
 *
 * @param[in] V
 *         (LDV,K) if STOREV = 'C'
 *         (LDV,M2) if STOREV = 'R' and SIDE = 'L'
 *         (LDV,N2) if STOREV = 'R' and SIDE = 'R'
 *         Matrix V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *         If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M2);
 *         if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N2);
 *         if STOREV = 'R', LDV >= K.
 *
 * @param[out] T
 *         The triangular K-by-K matrix T in the representation of the
 *         block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= K.
 *
 * @param[in,out] WORK
 *         Workspace of dimension LDWORK-by-N1 if side == MorseLeft, LDWORK-by-K
 *         otherwise.
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK: LDWORK >= K, if side ==
 *         MorseLeft, LDWORK >= M1 otehrwise.
 *
 * @param[in,out] WORKC
 *         Optionnal additional workspace to replace the TRMM operation by a GEMM kernel.
 *         This workspace is of dimension LDWORK-by-K if side == MorseLeft, LDWORK-by-N2
 *         otherwise.
 *
 * @param[in] LDWORKC
 *         The leading dimension of the array WORKC: LDWORKC >= M2, if side ==
 *         MorseLeft, LDWORK >= K otehrwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CUDA_zparfb(MORSE_enum side, MORSE_enum trans,
            MORSE_enum direct, MORSE_enum storev,
            int M1, int N1, int M2, int N2, int K, int L,
                  cuDoubleComplex *A1, int LDA1,
                  cuDoubleComplex *A2, int LDA2,
            const cuDoubleComplex *V, int LDV,
            const cuDoubleComplex *T, int LDT,
                  cuDoubleComplex *WORK, int LDWORK,
                  cuDoubleComplex *WORKC, int LDWORKC,
            CUBLAS_STREAM_PARAM )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zzero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
#else
    double zzero = 0.0;
    double zone  = 1.0;
    double mzone = -1.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    int j;
    MORSE_enum transW;
    MORSE_enum transA2;

    CUBLAS_GET_STREAM;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        return -1;
    }
    if ((trans != MorseNoTrans) && (trans != MorseConjTrans)) {
        return -2;
    }
    if ((direct != MorseForward) && (direct != MorseBackward)) {
        return -3;
    }
    if ((storev != MorseColumnwise) && (storev != MorseRowwise)) {
        return -4;
    }
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == MorseRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == MorseLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }
    if ( ((LDWORK < K ) && (side == MorseLeft )) ||
         ((LDWORK < M1) && (side == MorseRight)) ) {
        return -20;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return MORSE_SUCCESS;

    if (direct == MorseForward) {

        if (side == MorseLeft) {

            /*
             * Column or Rowwise / Forward / Left
             * ----------------------------------
             *
             * Form  H * A  or  H' * A  where  A = ( A1 )
             *                                     ( A2 )
             */

            /*
             * W = A1 + V' * A2:
             *      W = A1
             *      W = W + V' * A2
             *
             */
            cudaMemcpy2DAsync( WORK, LDWORK * sizeof(cuDoubleComplex),
                               A1,   LDA1   * sizeof(cuDoubleComplex),
                               K * sizeof(cuDoubleComplex), N1,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == MorseColumnwise ? MorseConjTrans : MorseNoTrans;
            transA2 = storev == MorseColumnwise ? MorseNoTrans : MorseConjTrans;

            cublasZgemm(CUBLAS_HANDLE
                        morse_cublas_const(transW), morse_cublas_const(MorseNoTrans),
                        K, N1, M2,
                        CUBLAS_SADDR(zone),
                        V     /* K*M2  */, LDV,
                        A2    /* M2*N1 */, LDA2,
                        CUBLAS_SADDR(zone),
                        WORK  /* K*N1  */, LDWORK);

            if (WORKC == NULL) {
                /* W = op(T) * W */
                CUDA_ztrmm( MorseLeft, MorseUpper, trans, MorseNonUnit,
                            K, N2,
                            CUBLAS_SADDR(zone), T,    LDT,
                                                WORK, LDWORK,
                            CUBLAS_STREAM_VALUE );

                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    cublasZaxpy(CUBLAS_HANDLE
                                K, CUBLAS_SADDR(mzone),
                                (WORK + LDWORK*j), 1,
                                (A1 + LDA1*j),     1);
                }

                /* A2 = A2 - op(V) * W  */
                cublasZgemm(CUBLAS_HANDLE
                            morse_cublas_const(transA2), morse_cublas_const(MorseNoTrans),
                            M2, N2, K,
                            CUBLAS_SADDR(mzone), V    /* M2*K  */, LDV,
                                                 WORK /* K*N2  */, LDWORK,
                            CUBLAS_SADDR(zone),  A2   /* m2*N2 */, LDA2);

            } else {
                /* Wc = V * op(T) */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(transA2), morse_cublas_const(trans),
                             M2, K, K,
                             CUBLAS_SADDR(zone),  V, LDV,
                                                  T, LDT,
                             CUBLAS_SADDR(zzero), WORKC, LDWORKC );

                /* A1 = A1 - opt(T) * W */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(trans), morse_cublas_const(MorseNoTrans),
                             K, N1, K,
                             CUBLAS_SADDR(mzone), T,    LDT,
                                                  WORK, LDWORK,
                             CUBLAS_SADDR(zone),  A1,   LDA1 );

                /* A2 = A2 - Wc * W */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(MorseNoTrans), morse_cublas_const(MorseNoTrans),
                             M2, N2, K,
                             CUBLAS_SADDR(mzone), WORKC, LDWORKC,
                                                  WORK,  LDWORK,
                             CUBLAS_SADDR(zone),  A2,    LDA2 );
            }
        }
        else {
            /*
             * Column or Rowwise / Forward / Right
             * -----------------------------------
             *
             * Form  H * A  or  H' * A  where A  = ( A1 A2 )
             *
             */

            /*
             * W = A1 + A2 * V':
             *      W = A1
             *      W = W + A2 * V'
             *
             */
            cudaMemcpy2DAsync( WORK, LDWORK * sizeof(cuDoubleComplex),
                               A1,   LDA1   * sizeof(cuDoubleComplex),
                               M1 * sizeof(cuDoubleComplex), K,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == MorseColumnwise ? MorseNoTrans : MorseConjTrans;
            transA2 = storev == MorseColumnwise ? MorseConjTrans : MorseNoTrans;

            cublasZgemm(CUBLAS_HANDLE
                        morse_cublas_const(MorseNoTrans), morse_cublas_const(transW),
                        M1, K, N2,
                        CUBLAS_SADDR(zone), A2   /* M1*N2 */, LDA2,
                                            V    /* N2*K  */, LDV,
                        CUBLAS_SADDR(zone), WORK /* M1*K  */, LDWORK);

            if (WORKC == NULL) {
                /* W = W * op(T) */
                CUDA_ztrmm( MorseRight, MorseUpper, trans, MorseNonUnit,
                            M2, K,
                            CUBLAS_SADDR(zone), T,    LDT,
                                                WORK, LDWORK,
                            CUBLAS_STREAM_VALUE );

                /* A1 = A1 - W = A1 - W * op(T) */
                for(j = 0; j < K; j++) {
                    cublasZaxpy(CUBLAS_HANDLE
                                M1, CUBLAS_SADDR(mzone),
                                (WORK + LDWORK*j), 1,
                                (A1 + LDA1*j), 1);
                }

                /* A2 = A2 - W * op(V)  */
                cublasZgemm(CUBLAS_HANDLE
                            morse_cublas_const(MorseNoTrans), morse_cublas_const(transA2),
                            M2, N2, K,
                            CUBLAS_SADDR(mzone), WORK /* M2*K  */, LDWORK,
                                                 V    /* K*N2  */, LDV,
                            CUBLAS_SADDR(zone),  A2   /* M2*N2 */, LDA2);

            } else {
                /* A1 = A1 - W * opt(T) */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(MorseNoTrans), morse_cublas_const(trans),
                             M1, K, K,
                             CUBLAS_SADDR(mzone), WORK, LDWORK,
                                                  T,    LDT,
                             CUBLAS_SADDR(zone),  A1,   LDA1 );

                /* Wc = op(T) * V */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(trans), morse_cublas_const(transA2),
                             K, N2, K,
                             CUBLAS_SADDR(zone),  T,     LDT,
                                                  V,     LDV,
                             CUBLAS_SADDR(zzero), WORKC, LDWORKC );

                /* A2 = A2 - W * Wc */
                cublasZgemm( CUBLAS_HANDLE
                             morse_cublas_const(MorseNoTrans), morse_cublas_const(MorseNoTrans),
                             M2, N2, K,
                             CUBLAS_SADDR(mzone), WORK,  LDWORK,
                                                  WORKC, LDWORKC,
                             CUBLAS_SADDR(zone),  A2,    LDA2 );
            }
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return MORSE_ERR_NOT_SUPPORTED;
    }

    (void)L;
    return MORSE_SUCCESS;
}
