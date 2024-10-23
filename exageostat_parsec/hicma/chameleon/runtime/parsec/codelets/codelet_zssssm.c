/**
 *
 * @file codelet_zssssm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm PaRSEC codelet
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
CORE_zssssm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *L1;
    int ldl1;
    MORSE_Complex64_t *L2;
    int ldl2;
    int *IPIV;

    parsec_dtd_unpack_args(
        this_task, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &L1, &ldl1, &L2, &ldl2, &IPIV );

    CORE_zssssm( m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zssssm(const MORSE_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                       const MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zssssm_parsec, options->priority, "ssssm",
        sizeof(int),           &m1,                                VALUE,
        sizeof(int),           &n1,                                VALUE,
        sizeof(int),           &m2,                                VALUE,
        sizeof(int),           &n2,                                VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A1, MORSE_Complex64_t, A1m, A1n ), morse_parsec_get_arena_index( A1 ) | INOUT,
        sizeof(int),           &lda1,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A2, MORSE_Complex64_t, A2m, A2n ), morse_parsec_get_arena_index( A2 ) | INOUT | AFFINITY,
        sizeof(int),           &lda2,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( L1, MORSE_Complex64_t, L1m, L1n ), morse_parsec_get_arena_index( L1 ) | INPUT,
        sizeof(int),           &ldl1,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( L2, MORSE_Complex64_t, L2m, L2n ), morse_parsec_get_arena_index( L2 ) | INPUT,
        sizeof(int),           &ldl2,                              VALUE,
        sizeof(int*),          &IPIV,                              VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
