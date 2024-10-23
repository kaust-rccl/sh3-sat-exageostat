/**
 *
 * @file codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro Quark codelet
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlatro_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int M;
    int N;
    const MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    quark_unpack_args_8(quark, uplo, trans, M, N, A, LDA, B, LDB);
    CORE_zlatro(uplo, trans, M, N, A, LDA, B, LDB);
}

void MORSE_TASK_zlatro(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zlatro_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &uplo,  VALUE,
        sizeof(MORSE_enum),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,   VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), OUTPUT,
        sizeof(int),                     &ldb,   VALUE,
        0);
}
