/**
 *
 * @file cuda_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherfb GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int
CUDA_zherfb( MORSE_enum uplo, int n,
             int k, int ib, int nb,
             const cuDoubleComplex *A, int lda,
             const cuDoubleComplex *T, int ldt,
             cuDoubleComplex *C, int ldc,
             cuDoubleComplex *WORK, int ldwork,
             CUBLAS_STREAM_PARAM )
{
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        cudablas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        cudablas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        cudablas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        cudablas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        cudablas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < chameleon_max(1,n)) && (n > 0) ) {
        cudablas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < chameleon_max(1,ib)) && (ib > 0) ) {
        cudablas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < chameleon_max(1,n)) && (n > 0) ) {
        cudablas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == MorseLower) {
        /* Left */
        CUDA_zunmqrt( MorseLeft, MorseConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
        /* Right */
        CUDA_zunmqrt( MorseRight, MorseNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
    }
    else {
        /* Right */
        CUDA_zunmlqt( MorseRight, MorseConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
        /* Left */
        CUDA_zunmlqt( MorseLeft, MorseNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
    }
    return 0;
}
