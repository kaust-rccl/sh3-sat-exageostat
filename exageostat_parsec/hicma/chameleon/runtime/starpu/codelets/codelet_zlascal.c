/**
 *
 * @file codelet_zlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Dalal Sukkari
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
 *  CORE_zlascal adds to matrices together.
 *
 *       A <- alpha * A
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

void MORSE_TASK_zlascal(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int m, int n, int nb,
                        MORSE_Complex64_t alpha,
                        const MORSE_desc_t *A, int Am, int An, int lda)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zlascal;
    void (*callback)(void*) = options->profiling ? cl_zlascal_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,              sizeof(MORSE_enum),
        STARPU_VALUE,    &m,                  sizeof(int),
        STARPU_VALUE,    &n,                  sizeof(int),
        STARPU_VALUE,    &alpha,              sizeof(MORSE_Complex64_t),
        STARPU_RW,         RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlascal",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zlascal_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &alpha, &LDA);
    CORE_zlascal(uplo, M, N, alpha, A, LDA);
    return;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlascal, 1, cl_zlascal_cpu_func)
