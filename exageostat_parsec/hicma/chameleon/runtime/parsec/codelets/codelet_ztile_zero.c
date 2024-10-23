/**
 *
 * @file codelet_ztile_zero.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztile_zero PaRSEC codelet
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
CORE_ztile_zero_parsec( parsec_execution_stream_t *context,
                        parsec_task_t             *this_task )
{
    int X1;
    int X2;
    int Y1;
    int Y2;
    MORSE_Complex64_t *A;
    int lda;
    int x, y;

    parsec_dtd_unpack_args(
        this_task, &X1, &X2, &Y1, &Y2, &A, &lda );

    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            A[lda * x + y] = 0.0;

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_ztile_zero( const MORSE_option_t *options,
                            int X1, int X2, int Y1, int Y2,
                            const MORSE_desc_t *A, int Am, int An, int lda )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztile_zero_parsec, options->priority, "tile zero",
        sizeof(int),       &X1,                       VALUE,
        sizeof(int),       &X2,                       VALUE,
        sizeof(int),       &Y1,                       VALUE,
        sizeof(int),       &Y2,                       VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | OUTPUT | AFFINITY,
        sizeof(int),       &lda,                      VALUE,
        PARSEC_DTD_ARG_END );
}
