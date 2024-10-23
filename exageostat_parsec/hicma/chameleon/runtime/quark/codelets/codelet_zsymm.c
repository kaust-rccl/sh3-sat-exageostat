/**
 *
 * @file codelet_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
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

void CORE_zsymm_quark(Quark *quark)
{
    MORSE_enum side;
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    CORE_zsymm(side, uplo,
        M, N,
        alpha, A, LDA,
        B, LDB,
        beta, C, LDC);
}

void MORSE_TASK_zsymm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_SYMM;
    QUARK_Insert_Task(opt->quark, CORE_zsymm_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &side,    VALUE,
        sizeof(MORSE_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(MORSE_Complex64_t),         &alpha,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(MORSE_Complex64_t),         &beta,    VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}
