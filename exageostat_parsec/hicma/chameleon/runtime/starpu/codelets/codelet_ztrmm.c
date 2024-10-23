/**
 *
 * @file codelet_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void MORSE_TASK_ztrmm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztrmm;
    void (*callback)(void*) = options->profiling ? cl_ztrmm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_RW(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,      &side,                sizeof(MORSE_enum),
        STARPU_VALUE,      &uplo,                sizeof(MORSE_enum),
        STARPU_VALUE,    &transA,                sizeof(MORSE_enum),
        STARPU_VALUE,      &diag,                sizeof(MORSE_enum),
        STARPU_VALUE,         &m,                        sizeof(int),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(MORSE_Complex64_t),
        STARPU_R,                 RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,       &lda,                        sizeof(int),
        STARPU_RW,                 RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE,       &ldb,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztrmm",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrmm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum uplo;
    MORSE_enum transA;
    MORSE_enum diag;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha, &LDA, &LDB);
    CORE_ztrmm(side, uplo,
        transA, diag,
        M, N,
        alpha, A, LDA,
        B, LDB);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_ztrmm_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum uplo;
    MORSE_enum transA;
    MORSE_enum diag;
    int M;
    int N;
    cuDoubleComplex alpha;
    const cuDoubleComplex *A;
    int LDA;
    cuDoubleComplex *B;
    int LDB;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha, &LDA, &LDB);

    RUNTIME_getStream(stream);

    CUDA_ztrmm(
        side, uplo,
        transA, diag,
        M, N,
        &alpha, A, LDA,
        B, LDB,
        stream);

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */


/*
 * Codelet definition
 */
CODELETS(ztrmm, 2, cl_ztrmm_cpu_func, cl_ztrmm_cuda_func, STARPU_CUDA_ASYNC)
