/**
 *
 * @file codelet_ztrssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrssq Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_ztrssq_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum diag;
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    quark_unpack_args_7( quark, uplo, diag, m, n, A, lda, SCALESUMSQ );
    CORE_ztrssq( uplo, diag, m, n, A, lda, &SCALESUMSQ[0], &SCALESUMSQ[1]);
}

void MORSE_TASK_ztrssq( const MORSE_option_t *options,
                        MORSE_enum uplo, MORSE_enum diag,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(opt->quark, CORE_ztrssq_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &uplo, VALUE,
        sizeof(MORSE_enum),              &diag, VALUE,
        sizeof(int),                     &m,    VALUE,
        sizeof(int),                     &n,    VALUE,
        sizeof(MORSE_Complex64_t)*lda*n, RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,  VALUE,
        sizeof(double)*2,                RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INOUT,
        0);
}
