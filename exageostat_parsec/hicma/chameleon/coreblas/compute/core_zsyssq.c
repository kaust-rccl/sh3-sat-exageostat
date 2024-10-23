/**
 *
 * @file core_zsyssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zsyssq CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "coreblas/lapacke.h"
#include "coreblas.h"

#define UPDATE( __nb, __value )                                         \
    if (__value != 0. ){                                                \
        if ( *scale < __value ) {                                       \
            *sumsq = __nb + (*sumsq) * ( *scale / __value ) * ( *scale / __value ); \
            *scale = __value;                                           \
        } else {                                                        \
            *sumsq = *sumsq + __nb * ( __value / *scale ) *  ( __value / *scale ); \
        }                                                               \
    }

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zsyssq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to N-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*n*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( A( i, j ) ) ), abs( aimag( A( i, j ) ) ) ),
 *          i,j
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile triangular part of the
 * symmetric tile A defined by uplo.
 * See also LAPACK zlassq.f
 *
 *******************************************************************************
 *
 *  @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = MorseLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = MorseUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 *  @param[in] N
 *          The number of columns and rows in the tile A.
 *
 *  @param[in] A
 *          The N-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,N).
 *
 *  @param[in,out] scale
 *          On entry, the value  scale  in the equation above.
 *          On exit, scale is overwritten with the value scl.
 *
 *  @param[in,out] sumsq
 *          On entry, the value  sumsq  in the equation above.
 *          On exit, SUMSQ is overwritten with the value ssq.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *
 */

int CORE_zsyssq(MORSE_enum uplo, int N,
                const MORSE_Complex64_t *A, int LDA,
                double *scale, double *sumsq)
{
    int i, j;
    double tmp;
    double *ptr;

    if ( uplo == MorseUpper ) {
        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA );

            for(i=0; i<j; i++, ptr++) {

                tmp = fabs(*ptr);
                UPDATE( 2., tmp );

#if defined(PRECISION_z) || defined(PRECISION_c)
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 2., tmp );
#endif
            }

            /* Diagonal */
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );

#if defined(PRECISION_z) || defined(PRECISION_c)
            ptr++;
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
#endif
        }
    } else {

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA + j);

            /* Diagonal */
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
            ptr++;

#if defined(PRECISION_z) || defined(PRECISION_c)
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
            ptr++;
#endif

            for(i=j+1; i<N; i++, ptr++) {

                tmp = fabs(*ptr);
                UPDATE( 2., tmp );

#if defined(PRECISION_z) || defined(PRECISION_c)
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 2., tmp );
#endif
            }
        }
    }
    return MORSE_SUCCESS;
}
