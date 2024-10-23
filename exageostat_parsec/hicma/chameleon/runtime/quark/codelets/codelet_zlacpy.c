/**
 *
 * @file codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
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

static inline void CORE_zlacpy_quark(Quark *quark)
{
    MORSE_enum uplo;
    int M;
    int N;
    int displA;
    MORSE_Complex64_t *A;
    int LDA;
    int displB;
    MORSE_Complex64_t *B;
    int LDB;

    quark_unpack_args_9(quark, uplo, M, N, displA, A, LDA, displB, B, LDB);
    CORE_zlacpy(uplo, M, N, A + displA, LDA, B + displB, LDB);
}

void MORSE_TASK_zlacpyx(const MORSE_option_t *options,
                        MORSE_enum uplo, int m, int n, int nb,
                        int displA, const MORSE_desc_t *A, int Am, int An, int lda,
                        int displB, const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zlacpy_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(int),                     &displA, VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(A, MORSE_Complex64_t, Am, An),             INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(int),                     &displB, VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),             OUTPUT,
        sizeof(int),                     &ldb,    VALUE,
        0);
}

void MORSE_TASK_zlacpy(const MORSE_option_t *options,
                       MORSE_enum uplo, int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    MORSE_TASK_zlacpyx( options, uplo, m, n, nb,
                        0, A, Am, An, lda,
                        0, B, Bm, Bn, ldb );
}
