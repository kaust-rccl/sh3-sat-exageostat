/**
 *
 * @file codelet_zttqrt.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zttqrt PaRSEC codelet
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
CORE_zttqrt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *TAU;
    MORSE_Complex64_t *WORK;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &ib, &A1, &lda1, &A2, &lda2, &T, &ldt, &TAU, &WORK );

    CORE_zttqrt( m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zttqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zttqrt_parsec, options->priority,  "ttqrt",
        sizeof(MORSE_enum),    &m,                                 VALUE,
        sizeof(int),           &n,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A1, MORSE_Complex64_t, A1m, A1n ), morse_parsec_get_arena_index( A1 ) | INOUT,
        sizeof(int),           &lda1,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A2, MORSE_Complex64_t, A2m, A2n ), morse_parsec_get_arena_index( A2 ) | INOUT | AFFINITY,
        sizeof(int),           &lda2,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ), morse_parsec_get_arena_index( T ) | OUTPUT,
        sizeof(int),           &ldt,                               VALUE,
        sizeof(MORSE_Complex64_t)*nb,       NULL,                  SCRATCH,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,                  SCRATCH,
        PARSEC_DTD_ARG_END );
}
