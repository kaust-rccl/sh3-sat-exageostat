/**
 *
 * @file codelet_zlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal PaRSEC codelet
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
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
static inline int
CORE_zlascal_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &alpha, &A, &LDA);

    CORE_zlascal( uplo, M, N, alpha, A, LDA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zlascal(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int m, int n, int nb,
                        MORSE_Complex64_t alpha,
                        const MORSE_desc_t *A, int Am, int An, int lda)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlascal_parsec, options->priority, "lascal",
        sizeof(MORSE_enum),        &uplo,  VALUE,
        sizeof(int),               &m,     VALUE,
        sizeof(int),               &n,     VALUE,
        sizeof(MORSE_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,              RTBLKADDR(A, MORSE_Complex64_t, Am, An), INOUT | AFFINITY,
        sizeof(int),               &lda,   VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}


