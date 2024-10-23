/**
 *
 * @file codelet_zasum.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zasum PaRSEC codelet
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
CORE_dzasum_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum storev;
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int lda;
    double *work;

    parsec_dtd_unpack_args(
        this_task, &storev, &uplo, &M, &N, &A, &lda, &work );

    CORE_dzasum( storev, uplo, M, N, A, lda, work );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_dzasum(const MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, int M, int N,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_dzasum_parsec, options->priority, "zasum",
        sizeof(MORSE_enum),    &storev,                           VALUE,
        sizeof(MORSE_enum),    &uplo,                             VALUE,
        sizeof(int),           &M,                                VALUE,
        sizeof(int),           &N,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),           &lda,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),     INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
