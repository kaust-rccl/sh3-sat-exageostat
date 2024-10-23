/**
 *
 * @file codelet_ztile_zero.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztile_zero StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 */
void MORSE_TASK_ztile_zero( const MORSE_option_t *options,
                            int X1, int X2, int Y1, int Y2,
                            const MORSE_desc_t *A, int Am, int An, int lda )
{
    struct starpu_codelet *codelet;
    codelet = &cl_ztile_zero;
    void (*callback)(void*) = options->profiling ? cl_zlacpy_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE, &X1,  sizeof(int),
        STARPU_VALUE, &X2,  sizeof(int),
        STARPU_VALUE, &Y1,  sizeof(int),
        STARPU_VALUE, &Y2,  sizeof(int),
        STARPU_W,     RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE, &lda, sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback, NULL,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztile_zero",
#endif
        0);
}

/**
 *
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_ztile_zero_cpu_func(void *descr[], void *cl_arg)
{
    int X1;
    int X2;
    int Y1;
    int Y2;
    MORSE_Complex64_t *A;
    int lda;

    int x, y;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &X1, &X2, &Y1, &Y2, &lda);

    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            A[lda*x+y] = 0.0;

}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztile_zero, 1, cl_ztile_zero_cpu_func)
