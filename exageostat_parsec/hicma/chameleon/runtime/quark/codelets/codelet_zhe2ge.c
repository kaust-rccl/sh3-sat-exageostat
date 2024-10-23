/**
 *
 * @file codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge Quark codelet
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
static inline void CORE_zhe2ge_quark(Quark *quark)
{
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    CORE_zhe2ge(uplo, M, N, A, LDA, B, LDB);
}


void MORSE_TASK_zhe2ge(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zhe2ge_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), OUTPUT,
        sizeof(int),                     &ldb,    VALUE,
        0);
}
