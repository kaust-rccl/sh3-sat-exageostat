/**
 *
 * @file codelet_zaxpy.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy PaRSEC codelet
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
CORE_zaxpy_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    int M;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int incA;
    MORSE_Complex64_t *B;
    int incB;

    parsec_dtd_unpack_args(
        this_task, &M, &alpha, &A, &incA, &B, &incB );

    CORE_zaxpy( M, alpha, A, incA, B, incB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zaxpy(const MORSE_option_t *options,
                      int M, MORSE_Complex64_t alpha,
                      const MORSE_desc_t *A, int Am, int An, int incA,
                      const MORSE_desc_t *B, int Bm, int Bn, int incB)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zaxpy_parsec, options->priority,  "axpy",
        sizeof(int),               &M,     VALUE,
        sizeof(MORSE_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),               &incA, VALUE,
        PASSED_BY_REF,  RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int),               &incB, VALUE,
        PARSEC_DTD_ARG_END );
}
