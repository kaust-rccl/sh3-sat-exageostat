/**
 *
 * @file codelet_ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt StarPU codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztpqrt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int L;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *WORK;

    A    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    T    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &M, &N, &L, &ib,
                                &lda, &ldb, &ldt );

    CORE_ztpqrt( M, N, L, ib,
                 A, lda, B, ldb, T, ldt, WORK );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztpqrt, 4, cl_ztpqrt_cpu_func)

void
MORSE_TASK_ztpqrt( const MORSE_option_t *options,
                   int M, int N, int L, int ib, int nb,
                   const MORSE_desc_t *A, int Am, int An, int lda,
                   const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                   const MORSE_desc_t *T, int Tm, int Tn, int ldt )
{
    struct starpu_codelet *codelet = &cl_ztpqrt;
    void (*callback)(void*) = options->profiling ? cl_ztpqrt_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_ACCESS_RW(B, Bm, Bn);
    MORSE_ACCESS_W(T, Tm, Tn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,    sizeof(int),
        STARPU_RW,     RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE, &lda,   sizeof(int),
        STARPU_RW,     RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE, &ldb,   sizeof(int),
        STARPU_W,      RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),
        STARPU_VALUE, &ldt,   sizeof(int),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztpqrt",
#endif
        0);

    (void)ib; (void)nb;
}
