/**
 *
 * @file codelet_zgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv StarPU codelet
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
 *  CORE_zgetrf_incpiv computes an LU factorization of a general M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factors L and U from the factorization
 *         A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile was interchanged with row IPIV(i).
 *
 * @param[out] INFO
 *         See returned value.
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

void MORSE_TASK_zgetrf_incpiv(const MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_incpiv;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_incpiv_callback : NULL;

    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_ACCESS_W(L, Lm, Ln);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_W,         RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),
        STARPU_VALUE,    &ldl,               sizeof(int),
        STARPU_VALUE,    &IPIV,              sizeof(int*),
        STARPU_VALUE,    &check_info,        sizeof(MORSE_bool),
        STARPU_VALUE,    &iinfo,             sizeof(int),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_incpiv",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_incpiv_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A;
    int lda, ldl;
    int *IPIV;
    MORSE_bool check_info;
    int iinfo;

    int info = 0;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldl, &IPIV, &check_info, &iinfo, &h_work);
    CORE_zgetrf_incpiv(m, n, ib, A, lda, IPIV, &info);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_incpiv, 3, cl_zgetrf_incpiv_cpu_func)
