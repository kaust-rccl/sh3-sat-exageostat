/**
 *
 * @file codelet_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
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
void MORSE_TASK_zsymm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zsymm;
    void (*callback)(void*) = options->profiling ? cl_zsymm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_R(B, Bm, Bn);
    MORSE_ACCESS_RW(C, Cm, Cn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &side,                sizeof(MORSE_enum),
        STARPU_VALUE,    &uplo,                sizeof(MORSE_enum),
        STARPU_VALUE,       &m,                        sizeof(int),
        STARPU_VALUE,       &n,                        sizeof(int),
        STARPU_VALUE,   &alpha,         sizeof(MORSE_Complex64_t),
        STARPU_R,               RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,     &lda,                        sizeof(int),
        STARPU_R,               RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE,     &ldb,                        sizeof(int),
        STARPU_VALUE,    &beta,         sizeof(MORSE_Complex64_t),
        STARPU_RW,               RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),
        STARPU_VALUE,     &ldc,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zsymm",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zsymm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int LDC;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &LDA, &LDB, &beta, &LDC);
    CORE_zsymm(side, uplo,
        M, N,
        alpha, A, LDA,
        B, LDB,
        beta, C, LDC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zsymm_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum uplo;
    int M;
    int N;
    cuDoubleComplex alpha;
    const cuDoubleComplex *A;
    int LDA;
    const cuDoubleComplex *B;
    int LDB;
    cuDoubleComplex beta;
    cuDoubleComplex *C;
    int LDC;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &LDA, &LDB, &beta, &LDC);

    RUNTIME_getStream(stream);

    CUDA_zsymm(
        side, uplo,
        M, N,
        &alpha, A, LDA,
        B, LDB,
        &beta, C, LDC,
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
CODELETS(zsymm, 3, cl_zsymm_cpu_func, cl_zsymm_cuda_func, STARPU_CUDA_ASYNC)
