/**
 *
 * @file codelet_zlag2c.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
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
CORE_zlag2c_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex32_t *B;
    int ldb;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &B, &ldb );

    CORE_zlag2c( m, n, A, lda, B, ldb );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zlag2c(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(PARSEC_dtd_taskpool, CORE_zlag2c_parsec, "lag2c",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),                        &lda,       VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex32_t, Bm, Bn ),     OUTPUT | AFFINITY,
        sizeof(int),                        &ldb,       VALUE,
        PARSEC_DTD_ARG_END );
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
static inline int
CORE_clag2z_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
{
    int m;
    int n;
    MORSE_Complex32_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &B, &ldb );

    CORE_clag2z( m, n, A, lda, B, ldb );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_clag2z(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_clag2z_parsec, options->priority, "lag2z",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex32_t, Am, An ),     INPUT,
        sizeof(int),                        &lda,       VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | OUTPUT | AFFINITY,
        sizeof(int),                        &ldb,       VALUE,
        PARSEC_DTD_ARG_END );
}
