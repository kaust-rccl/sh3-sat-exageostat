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
 * @brief Chameleon zhe2ge PaRSEC codelet
 *
 * @version 1.0.0
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
CORE_zhe2ge_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum uplo;
    int M;
    int N;
    const MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &A, &LDA, &B, &LDB);

    CORE_zhe2ge( uplo, M, N, A, LDA, B, LDB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}


void MORSE_TASK_zhe2ge(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zhe2ge_parsec, options->priority, "he2ge",
        sizeof(MORSE_enum), &uplo,   VALUE,
        sizeof(int),        &m,      VALUE,
        sizeof(int),        &n,      VALUE,
        PASSED_BY_REF,       RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT ,
        sizeof(int),        &lda,    VALUE,
        PASSED_BY_REF,       RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), OUTPUT | AFFINITY,
        sizeof(int),        &ldb,    VALUE,
        PARSEC_DTD_ARG_END );

    (void)mb;
}
