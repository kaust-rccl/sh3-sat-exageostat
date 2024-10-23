/**
 *
 * @file codelet_zgetrf.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgetrf_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int *IPIV;
    MORSE_bool *check_info;
    int iinfo;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &IPIV, &check_info, &iinfo );

    CORE_zgetrf( m, n, A, lda, IPIV, &info );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zgetrf(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgetrf_parsec, options->priority, "getrf",
        sizeof(int),        &m,                          VALUE,
        sizeof(int),        &n,                          VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int),        &lda,                        VALUE,
        sizeof(int)*nb,      IPIV,                        SCRATCH,
        sizeof(MORSE_bool), &check_info,                 VALUE,
        sizeof(int),        &iinfo,                      VALUE,
        PARSEC_DTD_ARG_END );
}
