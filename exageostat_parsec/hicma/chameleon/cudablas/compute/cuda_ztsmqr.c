/**
 *
 * @file cuda_ztsmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztsmqr GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_ztsmqr(
        MORSE_enum side, MORSE_enum trans,
        int M1, int N1,
        int M2, int N2,
        int K, int IB,
              cuDoubleComplex *A1,    int LDA1,
              cuDoubleComplex *A2,    int LDA2,
        const cuDoubleComplex *V,     int LDV,
        const cuDoubleComplex *T,     int LDT,
              cuDoubleComplex *WORK,  int LDWORK,
              cuDoubleComplex *WORKC, int LDWORKC,
        CUBLAS_STREAM_PARAM)
{
    int i, i1, i3;
    int NQ, NW;
    int kb;
    int ic = 0;
    int jc = 0;
    int mi = M1;
    int ni = N1;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        return -1;
    }

    /* NQ is the order of Q */
    if (side == MorseLeft) {
        NQ = M2;
        NW = IB;
    }
    else {
        NQ = N2;
        NW = M1;
    }

    if ((trans != MorseNoTrans) && (trans != MorseConjTrans)) {
        return -2;
    }
    if (M1 < 0) {
        return -3;
    }
    if (N1 < 0) {
        return -4;
    }
    if ( (M2 < 0) ||
         ( (M2 != M1) && (side == MorseRight) ) ){
        return -5;
    }
    if ( (N2 < 0) ||
         ( (N2 != N1) && (side == MorseLeft) ) ){
        return -6;
    }
    if ((K < 0) ||
        ( (side == MorseLeft)  && (K > M1) ) ||
        ( (side == MorseRight) && (K > N1) ) ) {
        return -7;
    }
    if (IB < 0) {
        return -8;
    }
    if (LDA1 < chameleon_max(1,M1)){
        return -10;
    }
    if (LDA2 < chameleon_max(1,M2)){
        return -12;
    }
    if (LDV < chameleon_max(1,NQ)){
        return -14;
    }
    if (LDT < chameleon_max(1,IB)){
        return -16;
    }
    if (LDWORK < chameleon_max(1,NW)){
        return -18;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
        return MORSE_SUCCESS;

    if (((side == MorseLeft)  && (trans != MorseNoTrans))
        || ((side == MorseRight) && (trans == MorseNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ((K-1) / IB)*IB;
        i3 = -IB;
    }

    for(i = i1; (i > -1) && (i < K); i += i3) {
        kb = chameleon_min(IB, K-i);

        if (side == MorseLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M1 - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N1 - i;
            jc = i;
        }
        /*
         * Apply H or H' (NOTE: CORE_zparfb used to be CORE_ztsrfb)
         */
        CUDA_zparfb(
                side, trans, MorseForward, MorseColumnwise,
                mi, ni, M2, N2, kb, 0,
                A1 + LDA1*jc+ic, LDA1,
                A2, LDA2,
                V + LDV*i, LDV,
                T + LDT*i, LDT,
                WORK, LDWORK, WORKC, LDWORKC, CUBLAS_STREAM_VALUE );
    }
    return MORSE_SUCCESS;
}
