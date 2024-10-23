/**
 *
 * @file codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset Quark codelet
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
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlaset_quark(Quark *quark)
{
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *A;
    int LDA;

    quark_unpack_args_7(quark, uplo, M, N, alpha, beta, A, LDA);
    CORE_zlaset(uplo, M, N, alpha, beta, A, LDA);
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zlaset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = MorseUpper: Upper part of A is set;
 *          = MorseLower: Lower part of A is set;
 *          = MorseUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
void MORSE_TASK_zlaset(const MORSE_option_t *options,
                       MORSE_enum uplo, int M, int N,
                       MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                       const MORSE_desc_t *A, int Am, int An, int LDA)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASET;
    QUARK_Insert_Task(opt->quark, CORE_zlaset_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(MORSE_Complex64_t),         &alpha, VALUE,
        sizeof(MORSE_Complex64_t),         &beta,  VALUE,
        sizeof(MORSE_Complex64_t)*LDA*N,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}
