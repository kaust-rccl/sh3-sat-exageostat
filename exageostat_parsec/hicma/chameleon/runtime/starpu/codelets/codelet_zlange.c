/**
 *
 * @file codelet_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_zlange(const MORSE_option_t *options,
                       MORSE_enum norm, int M, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn)
{
    (void)NB;
    struct starpu_codelet *codelet = &cl_zlange;
    void (*callback)(void*) = options->profiling ? cl_zlange_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_W(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &norm,              sizeof(MORSE_enum),
        STARPU_VALUE,    &M,                 sizeof(int),
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_R,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &LDA,               sizeof(int),
        STARPU_SCRATCH,  options->ws_worker,
        STARPU_W,        RTBLKADDR(B, double, Bm, Bn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlange",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlange_cpu_func(void *descr[], void *cl_arg)
{
    double *normA;
    MORSE_enum norm;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    double *work;

    A     = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    work  = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    normA = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &norm, &M, &N, &LDA);
    CORE_zlange( norm, M, N, A, LDA, work, normA);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlange, 3, cl_zlange_cpu_func)

void MORSE_TASK_zlange_max(const MORSE_option_t *options,
                           const MORSE_desc_t *A, int Am, int An,
                           const MORSE_desc_t *B, int Bm, int Bn)
{
    struct starpu_codelet *codelet = &cl_zlange_max;
    void (*callback)(void*) = options->profiling ? cl_zlange_callback : NULL;

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( B, Bm, Bn ) )
    {
        starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_R,        RTBLKADDR(A, double, Am, An),
            STARPU_RW,       RTBLKADDR(B, double, Bm, Bn),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zlange_max",
#endif
            0);
    }
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlange_max_cpu_func(void *descr[], void *cl_arg)
{
    double *A;
    double *normA;

    A     = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    normA = (double *)STARPU_MATRIX_GET_PTR(descr[1]);

    if ( *A > *normA )
        *normA = *A;

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlange_max, 2, cl_zlange_max_cpu_func)
