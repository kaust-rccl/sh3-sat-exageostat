/**
 *
 * @file codelet_zaxpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy StarPU codelet
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-07-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_zaxpy(const MORSE_option_t *options,
                      int M, MORSE_Complex64_t alpha,
                      const MORSE_desc_t *A, int Am, int An, int incA,
                      const MORSE_desc_t *B, int Bm, int Bn, int incB)
{
    struct starpu_codelet *codelet = &cl_zaxpy;
    void (*callback)(void*) = options->profiling ? cl_zaxpy_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_RW(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &M,                           sizeof(int),
            STARPU_VALUE,    alpha,                       sizeof(MORSE_Complex64_t),
            STARPU_R,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,    &incA,                        sizeof(int),
            STARPU_RW,       RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
            STARPU_VALUE,    &incB,                        sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zaxpy",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zaxpy_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int incA;
    MORSE_Complex64_t *B;
    int incB;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &M, &alpha, &incA, &incB);
    CORE_zaxpy(M, alpha, A, incA, B, incB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zaxpy, 2, cl_zaxpy_cpu_func)
