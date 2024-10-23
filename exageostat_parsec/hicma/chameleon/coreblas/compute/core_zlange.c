/**
 *
 * @file core_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlange CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_PLASMA_zlange returns the value
 *
 *     zlange = ( max(abs(A(i,j))), NORM = MorseMaxNorm
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
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension (MAX(1,LWORK)), where LWORK >= M when NORM =
 *          MorseInfNorm; otherwise, WORK is not referenced.
 *
 * @param[out] normA
 *          On exit, normA is the norm of matrix A.
 *
 */

void CORE_zlange(int norm, int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
    *normA = LAPACKE_zlange_work(
        LAPACK_COL_MAJOR,
        morse_lapack_const(norm),
        M, N, A, LDA, work);
}
