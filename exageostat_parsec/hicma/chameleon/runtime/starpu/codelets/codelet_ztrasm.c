/**
 *
 * @file codelet_ztrasm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrasm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_ztrasm(const MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn)
{
    struct starpu_codelet *codelet = &cl_ztrasm;
    void (*callback)(void*) = options->profiling ? cl_ztrasm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_RW(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &storev,                    sizeof(MORSE_enum),
        STARPU_VALUE,    &uplo,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &diag,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &M,                         sizeof(int),
        STARPU_VALUE,    &N,                         sizeof(int),
        STARPU_R,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_RW,       RTBLKADDR(B, double, Bm, Bn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "ztrasm",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrasm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum storev;
    MORSE_enum uplo;
    MORSE_enum diag;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int lda;
    double *work;

    A     = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    work  = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &storev, &uplo, &diag, &M, &N, &lda);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztrasm, 2, cl_ztrasm_cpu_func)
