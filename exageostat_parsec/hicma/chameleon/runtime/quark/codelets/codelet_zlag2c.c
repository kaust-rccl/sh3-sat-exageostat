/**
 *
 * @file codelet_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlag2c_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex32_t *B;
    int ldb;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    int info;

    quark_unpack_args_8(quark, m, n, A, lda, B, ldb, sequence, request);
    CORE_zlag2c( m, n, A, lda, B, ldb);
    if (sequence->status == MORSE_SUCCESS && info != 0)
        RUNTIME_sequence_flush(quark, sequence, request, info);
}

void MORSE_TASK_zlag2c(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAG2C;
    QUARK_Insert_Task(opt->quark, CORE_zlag2c_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INPUT,
                      sizeof(int),                        &lda,       VALUE,
                      sizeof(MORSE_Complex32_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex32_t, Bm, Bn),                 OUTPUT,
                      sizeof(int),                        &ldb,       VALUE,
                      sizeof(MORSE_sequence_t*),           &(options->sequence),  VALUE,
                      sizeof(MORSE_request_t*),            &(options->request),   VALUE,
                      0);
}

void CORE_clag2z_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex32_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;

    quark_unpack_args_6(quark, m, n, A, lda, B, ldb);
    CORE_clag2z( m, n, A, lda, B, ldb);
}

void MORSE_TASK_clag2z(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    QUARK_Insert_Task(opt->quark, CORE_clag2z_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                        &m,     VALUE,
                      sizeof(int),                        &n,     VALUE,
                      sizeof(MORSE_Complex32_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex32_t, Am, An),             INPUT,
                      sizeof(int),                        &lda,   VALUE,
                      sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),             INOUT,
                      sizeof(int),                        &ldb,   VALUE,
                      0);
}
