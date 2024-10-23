/**
 *
 * @file codelet_zhemm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c
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
CORE_zhemm_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    MORSE_enum side;
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int LDC;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &M, &N, &alpha, &A, &LDA, &B, &LDB, &beta, &C, &LDC );

    CORE_zhemm( side, uplo, M, N,
               alpha, A, LDA,
                      B, LDB,
               beta,  C, LDC);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zhemm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zhemm_parsec, options->priority, "hemm",
        sizeof(MORSE_enum),                &side,    VALUE,
        sizeof(MORSE_enum),                &uplo,    VALUE,
        sizeof(int),                       &m,       VALUE,
        sizeof(int),                       &n,       VALUE,
        sizeof(MORSE_Complex64_t),         &alpha,   VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),           &lda,                 VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | INPUT,
        sizeof(int),           &ldb,                 VALUE,
        sizeof(MORSE_Complex64_t),         &beta,    VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, MORSE_Complex64_t, Cm, Cn ), morse_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int),           &ldc,                 VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
