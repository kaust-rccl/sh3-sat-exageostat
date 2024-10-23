/**
 *
 * @file codelet_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf StarPU codelet
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
 *  CORE_ztstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular NB-by-N tile U on top of a M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] NB
 *
 * @param[in,out] U
 *         On entry, the NB-by-N upper triangular tile.
 *         On exit, the new factor U from the factorization
 *
 * @param[in] LDU
 *         The leading dimension of the array U.  LDU >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factor L from the factorization
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[in,out] L
 *         On entry, the IB-by-N lower triangular tile.
 *         On exit, the interchanged rows form the tile A in case of pivoting.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,IB).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile U was interchanged with row IPIV(i) of the tile A.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *
 * @param[out] INFO
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

void MORSE_TASK_ztstrf(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *U, int Um, int Un, int ldu,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztstrf;
    void (*callback)(void*) = options->profiling ? cl_ztstrf_callback : NULL;
    MORSE_starpu_ws_t *d_work = (MORSE_starpu_ws_t*)(options->ws_host);

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(U, Um, Un);
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_ACCESS_W(L, Lm, Ln);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_VALUE,    &nb,                        sizeof(int),
        STARPU_RW,        RTBLKADDR(U, MORSE_Complex64_t, Um, Un),
        STARPU_VALUE,    &ldu,                       sizeof(int),
        STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_W,         RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),
        STARPU_VALUE,    &ldl,                       sizeof(int),
        STARPU_VALUE,    &IPIV,                      sizeof(int*),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &d_work,                    sizeof(MORSE_starpu_ws_t *),
        STARPU_VALUE,    &nb,                        sizeof(int),
        STARPU_VALUE,    &check_info,                sizeof(MORSE_bool),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztstrf",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztstrf_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *d_work;
    int m;
    int n;
    int ib;
    int nb;
    MORSE_Complex64_t *U;
    int ldu;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *L;
    int ldl;
    int *IPIV;
    MORSE_Complex64_t *WORK;
    int ldwork;
    MORSE_bool check_info;
    int iinfo;

    int info = 0;

    U = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    L = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &nb, &ldu, &lda, &ldl, &IPIV, &d_work, &ldwork, &check_info, &iinfo);

    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztstrf, 4, cl_ztstrf_cpu_func)

