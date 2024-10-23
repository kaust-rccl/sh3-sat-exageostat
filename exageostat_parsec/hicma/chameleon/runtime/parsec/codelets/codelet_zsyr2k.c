/**
 *
 * @file codelet_zsyr2k.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k PaRSEC codelet
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
CORE_zsyr2k_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int n;
    int k;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int ldc;

    parsec_dtd_unpack_args(
        this_task, &uplo, &trans, &n, &k, &alpha, &A, &lda, &B, &ldb, &beta, &C, &ldc );

    CORE_zsyr2k( uplo, trans, n, k,
                alpha, A, lda,
                        B, ldb,
                beta,  C, ldc);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zsyr2k(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zsyr2k_parsec, options->priority, "syr2k",
        sizeof(MORSE_enum),     &uplo,                  VALUE,
        sizeof(MORSE_enum),     &trans,                 VALUE,
        sizeof(int),            &n,                     VALUE,
        sizeof(int),            &k,                     VALUE,
        sizeof(MORSE_Complex64_t), &alpha,              VALUE,
        PASSED_BY_REF,          RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),            &lda,                   VALUE,
        PASSED_BY_REF,          RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | INPUT,
        sizeof(int),            &ldb,                   VALUE,
        sizeof(MORSE_Complex64_t), &beta,               VALUE,
        PASSED_BY_REF,          RTBLKADDR( C, MORSE_Complex64_t, Cm, Cn ), morse_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int),            &ldc,                   VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
