/**
 *
 * @file codelet_zgessm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm StarPU codelet
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
 *  CORE_zgessm applies the factors L computed by CORE_zgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *         The number of columns of the tile L. K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_zgetrf_incpiv.
 *
 * @param[in] L
 *         The M-by-K lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,M).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */

void MORSE_TASK_zgessm(const MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       const MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       const MORSE_desc_t *A, int Am, int An, int lda)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgessm;
    void (*callback)(void*) = options->profiling ? cl_zgessm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(L, Lm, Ln);
    MORSE_ACCESS_R(D, Dm, Dn);
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,     &m,                        sizeof(int),
        STARPU_VALUE,     &n,                        sizeof(int),
        STARPU_VALUE,     &k,                        sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_VALUE,          &IPIV,                      sizeof(int*),
        STARPU_R,             RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),
        STARPU_VALUE,   &ldl,                        sizeof(int),
        STARPU_R,             RTBLKADDR(D, MORSE_Complex64_t, Dm, Dn),
        STARPU_VALUE,   &ldd,                        sizeof(int),
        STARPU_RW,            RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,   &lda,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgessm",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zgessm_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    int ldl;
    MORSE_Complex64_t *D;
    int ldd;
    MORSE_Complex64_t *A;
    int lda;

    D = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &k, &ib, &IPIV, &ldl, &ldd, &lda);
    CORE_zgessm(m, n, k, ib, IPIV, D, ldd, A, lda);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgessm, 3, cl_zgessm_cpu_func)
