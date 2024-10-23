/**
 *
 * @file codelet_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmlqt StarPU codelet
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
static void cl_ztpmlqt_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    const MORSE_Complex64_t *V;
    int ldv;
    const MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t *WORK;

    V    = (const MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T    = (const MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    A    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    B    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[4]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib,
                                &ldv, &ldt, &lda, &ldb );

    CORE_ztpmlqt( side, trans, M, N, K, L, ib,
                  V, ldv, T, ldt, A, lda, B, ldb, WORK );
}

#if defined(CHAMELEON_USE_CUDA) && 0
static void cl_ztpmlqt_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    const cuDoubleComplex *V;
    int ldv;
    const cuDoubleComplex *T;
    int ldt;
    cuDoubleComplex *A;
    int lda;
    cuDoubleComplex *B;
    int ldb;
    cuDoubleComplex *W;

    V = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    T = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]);
    W = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[4]); /* 2*ib*nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib,
                                &ldv, &ldt, &lda, &ldb );

    RUNTIME_getStream(stream);

    CUDA_ztpmlqt(
            side, trans, M, N, K, L, ib,
            V, ldv, T, ldt, A, lda, B, ldb,
            W, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztpmlqt, 5, cl_ztpmlqt_cpu_func)
//CODELETS(ztpmlqt, 5, cl_ztpmlqt_cpu_func, cl_ztpmlqt_cuda_func, STARPU_CUDA_ASYNC)

void
MORSE_TASK_ztpmlqt( const MORSE_option_t *options,
                    MORSE_enum side, MORSE_enum trans,
                    int M, int N, int K, int L, int ib, int nb,
                    const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                    const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                    const MORSE_desc_t *A, int Am, int An, int lda,
                    const MORSE_desc_t *B, int Bm, int Bn, int ldb )
{
    struct starpu_codelet *codelet = &cl_ztpmlqt;
    void (*callback)(void*) = options->profiling ? cl_ztpmlqt_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(V, Vm, Vn);
    MORSE_ACCESS_R(T, Tm, Tn);
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_ACCESS_RW(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE, &side,  sizeof(MORSE_enum),
        STARPU_VALUE, &trans, sizeof(MORSE_enum),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &K,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,     sizeof(int),
        STARPU_R,      RTBLKADDR(V, MORSE_Complex64_t, Vm, Vn),
        STARPU_VALUE, &ldv,   sizeof(int),
        STARPU_R,      RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),
        STARPU_VALUE, &ldt,   sizeof(int),
        STARPU_RW,     RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE, &lda,   sizeof(int),
        STARPU_RW,     RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE, &ldb,   sizeof(int),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, (( L == 0 ) ? "ztsmlq" : "ztpmlqt"),
#endif
        0);

    (void)ib; (void)nb;
}
