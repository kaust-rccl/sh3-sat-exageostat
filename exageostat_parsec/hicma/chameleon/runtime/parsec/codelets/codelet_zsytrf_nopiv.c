/**
 *
 * @file codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zsytrf_nopiv_parsec( parsec_execution_stream_t *context,
                          parsec_task_t             *this_task )
{
    MORSE_enum uplo;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int iinfo;

    parsec_dtd_unpack_args(
        this_task, &uplo, &n, &A, &lda, &iinfo );

    CORE_zsytf2_nopiv( uplo, n, A, lda );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zsytrf_nopiv(const MORSE_option_t *options,
                             MORSE_enum uplo, int n, int nb,
                             const MORSE_desc_t *A, int Am, int An, int lda,
                             int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zsytrf_nopiv_parsec, options->priority, "sytrf_nopiv",
        sizeof(MORSE_enum),              &uplo,                VALUE,
        sizeof(int),                     &n,                   VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int),                     &lda,                 VALUE,
        sizeof(int),                     &iinfo,               VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
