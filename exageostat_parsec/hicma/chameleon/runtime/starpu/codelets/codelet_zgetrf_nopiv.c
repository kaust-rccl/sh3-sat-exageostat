/**
 *
 * @file codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv StarPU codelet
 *
 * @version 1.0.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2013-02-01
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgetrf_nopiv computes an LU factorization of a general diagonal
 *  dominant M-by-N matrix A witout pivoting.
 *
 *  The factorization has the form
 *     A = L * U
 *  where L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 3 BLAS version of the algorithm.
 *  WARNING: Your matrix need to be diagonal dominant if you want to call this
 *  routine safely.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  @param[in] IB
 *          The block size to switch between blocked and unblocked code.
 *
 *  @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  @param[in] LDA
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 */

void MORSE_TASK_zgetrf_nopiv(const MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_nopiv;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_nopiv_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_nopiv",
#endif
        0);
}

/*
 * Codelet CPU
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_nopiv_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    int iinfo;
    int info = 0;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &iinfo);
    CORE_zgetrf_nopiv(m, n, ib, A, lda, &info);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_nopiv, 1, cl_zgetrf_nopiv_cpu_func)
