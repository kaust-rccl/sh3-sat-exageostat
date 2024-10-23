/**
 *
 * @file core_zlantr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlantr CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include <math.h>
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_PLASMA_zlantr returns the value
 *
 *     zlantr = ( max(abs(A(i,j))), NORM = MorseMaxNorm
 *              (
 *              ( norm1(A),         NORM = MorseOneNorm
 *              (
 *              ( normI(A),         NORM = MorseInfNorm
 *              (
 *              ( normF(A),         NORM = MorseFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = MorseMaxNorm: Max norm
 *          = MorseOneNorm: One norm
 *          = MorseInfNorm: Infinity norm
 *          = MorseFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = MorseNonUnit: A is non unit;
 *          = MorseUnit:    A us unit.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *          If uplo == MorseUpper, M <= N. When M = 0, CORE_zlantr returns 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *          If uplo == MorseLower, N <= M. When N = 0, CORE_zlantr returns 0.
 *
 * @param[in] A
 *          The LDA-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension (MAX(1,LWORK)), where LWORK >= M when norm =
 *          MorseInfNorm, or LWORK >= N when norm = MorseOneNorm; otherwise,
 *          work is not referenced.
 *
 * @param[out] normA
 *          On exit, normA is the norm of matrix A.
 *
 */

void CORE_zlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag,
                 int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
#if defined(LAPACKE_CORRECT_DLANTR)
    *normA = LAPACKE_zlantr_work(
        LAPACK_COL_MAJOR,
        morse_lapack_const(norm),
        morse_lapack_const(uplo),
        morse_lapack_const(diag),
        M, N, A, LDA, work);
#else
    const MORSE_Complex64_t *tmpA;
    double value;
    int i, j, imax;
    int idiag = (diag == MorseUnit) ? 1 : 0;

    if ( chameleon_min(M, N) == 0 ) {
        *normA = 0;
        return;
    }

    switch ( norm ) {
    case MorseMaxNorm:
        if ( diag == MorseUnit ) {
            *normA = 1.;
        } else {
            *normA = 0.;
        }

        if ( uplo == MorseUpper ) {
            M = chameleon_min(M, N);
            for (j = 0; j < N; j++) {
                tmpA = A+(j*LDA);
                imax = chameleon_min(j+1-idiag, M);

                for (i = 0; i < imax; i++) {
                    value = cabs( *tmpA );
                    *normA = ( value > *normA ) ? value : *normA;
                    tmpA++;
                }
            }
        } else {
            N = chameleon_min(M, N);
            for (j = 0; j < N; j++) {
                tmpA = A + j * (LDA+1) + idiag;

                for (i = j+idiag; i < M; i++) {
                    value = cabs( *tmpA );
                    *normA = ( value > *normA ) ? value : *normA;
                    tmpA++;
                }
            }
        }
        break;

    case MorseOneNorm:
        CORE_ztrasm( MorseColumnwise, uplo, diag, M, N,
                     A, LDA, work );
        if ( uplo == MorseLower )
            N = chameleon_min(M,N);

        *normA = 0;
        for (i = 0; i < N; i++) {
            *normA = ( work[i] > *normA ) ? work[i] : *normA;
        }
        break;

    case MorseInfNorm:
        CORE_ztrasm( MorseRowwise, uplo, diag, M, N,
                     A, LDA, work );
        if ( uplo == MorseUpper )
            M = chameleon_min(M,N);

        *normA = 0;
        for (i = 0; i < M; i++) {
            *normA = ( work[i] > *normA ) ? work[i] : *normA;
        }
        break;

    case MorseFrobeniusNorm:
    {
        double scale = 0.;
        double sumsq = 1.;
        CORE_ztrssq( uplo, diag, M, N,
                     A, LDA, &scale, &sumsq );

        *normA = scale * sqrt( sumsq );
    }
    break;
    default:
        coreblas_error(1, "Illegal value of norm");
        return;
    }
#endif
}
