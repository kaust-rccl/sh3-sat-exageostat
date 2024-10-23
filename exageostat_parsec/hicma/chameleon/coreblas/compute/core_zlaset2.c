/**
 *
 * @file core_zlaset2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlaset2 CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
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
 *  CORE_zlaset2 - Sets the elements of the matrix A to alpha.
 *  Not LAPACK compliant! Read below.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = MorseUpper: STRICT Upper part of A is set to alpha;
 *          = MorseLower: STRICT Lower part of A is set to alpha;
 *          = MorseUpperLower: ALL elements of A are set to alpha.
 *          Not LAPACK Compliant.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set to alpha accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
void CORE_zlaset2(MORSE_enum uplo, int M, int N,
                  MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA)
{
    if (uplo == MorseUpper) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            morse_lapack_const(uplo),
            M, N-1, alpha, alpha, A+LDA, LDA);
    }
    else if (uplo == MorseLower) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            morse_lapack_const(uplo),
            M-1, N, alpha, alpha, A+1, LDA);
    }
    else {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            morse_lapack_const(uplo),
            M, N, alpha, alpha, A, LDA);
    }
}


