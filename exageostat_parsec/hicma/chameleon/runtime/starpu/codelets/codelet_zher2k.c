/**
 *
 * @file codelet_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k StarPU codelet
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
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void MORSE_TASK_zher2k(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zher2k;
    void (*callback)(void*) = options->profiling ? cl_zher2k_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_R(B, Bm, Bn);
    MORSE_ACCESS_RW(C, Cm, Cn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,      &uplo,                sizeof(MORSE_enum),
        STARPU_VALUE,     &trans,                sizeof(MORSE_enum),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,         &k,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(MORSE_Complex64_t),
        STARPU_R,                 RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,       &lda,                        sizeof(int),
        STARPU_R,                 RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE,       &ldb,                        sizeof(int),
        STARPU_VALUE,      &beta,                     sizeof(double),
        STARPU_RW,                 RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),
        STARPU_VALUE,       &ldc,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zher2k",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zher2k_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int n;
    int k;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    double beta;
    MORSE_Complex64_t *C;
    int ldc;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);
    CORE_zher2k(uplo, trans,
                n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zher2k_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int n;
    int k;
    cuDoubleComplex alpha;
    const cuDoubleComplex *A;
    int lda;
    const cuDoubleComplex *B;
    int ldb;
    double beta;
    cuDoubleComplex *C;
    int ldc;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);

    RUNTIME_getStream(stream);

    CUDA_zher2k( uplo, trans,
                 n, k, &alpha, A, lda, B, ldb, &beta, C, ldc,
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
CODELETS(zher2k, 3, cl_zher2k_cpu_func, cl_zher2k_cuda_func, STARPU_CUDA_ASYNC)
