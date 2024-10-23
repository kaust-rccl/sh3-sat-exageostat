/**
 *
 * @file codelet_ztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri Quark codelet
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

void CORE_ztrtri_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum diag;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    CORE_ztrtri(uplo, diag, N, A, LDA, &info);
    if ( (sequence->status == MORSE_SUCCESS) && (info > 0) ) {
        RUNTIME_sequence_flush( (MORSE_context_t*)quark, sequence, request, iinfo+info );
    }
}

void MORSE_TASK_ztrtri(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag,
                       int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(
        opt->quark, CORE_ztrtri_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo,      VALUE,
        sizeof(MORSE_enum),                &diag,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(MORSE_sequence_t*),           &(options->sequence),  VALUE,
        sizeof(MORSE_request_t*),            &(options->request),   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}
