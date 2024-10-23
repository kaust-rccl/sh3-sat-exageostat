/**
 *
 * @file codelet_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm StarPU codelet
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
 *  CORE_zssssm applies the LU factorization update from a complex
 *  matrix formed by a lower triangular IB-by-K tile L1 on top of a
 *  M2-by-K tile L2 to a second complex matrix formed by a M1-by-N1
 *  tile A1 on top of a M2-by-N2 tile A2 (N1 == N2).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M1
 *         The number of rows of the tile A1.  M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1.  N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2 and of the tile L2.
 *         M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2.  N2 >= 0.
 *
 * @param[in] K
 *         The number of columns of the tiles L1 and L2.  K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M2).
 *
 * @param[in] L1
 *         The IB-by-K lower triangular tile as returned by
 *         CORE_ztstrf.
 *
 * @param[in] LDL1
 *         The leading dimension of the array L1.  LDL1 >= max(1,IB).
 *
 * @param[in] L2
 *         The M2-by-K tile as returned by CORE_ztstrf.
 *
 * @param[in] LDL2
 *         The leading dimension of the array L2.  LDL2 >= max(1,M2).
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_ztstrf.
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */

void MORSE_TASK_zssssm(const MORSE_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                       const MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zssssm;
    void (*callback)(void*) = options->profiling ? cl_zssssm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A1, A1m, A1n);
    MORSE_ACCESS_RW(A2, A2m, A2n);
    MORSE_ACCESS_R(L1, L1m, L1n);
    MORSE_ACCESS_R(L2, L2m, L2n);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m1,                        sizeof(int),
        STARPU_VALUE,    &n1,                        sizeof(int),
        STARPU_VALUE,    &m2,                        sizeof(int),
        STARPU_VALUE,    &n2,                        sizeof(int),
        STARPU_VALUE,     &k,                        sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_RW,            RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n),
        STARPU_VALUE,  &lda1,                        sizeof(int),
        STARPU_RW,            RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n),
        STARPU_VALUE,  &lda2,                        sizeof(int),
        STARPU_R,             RTBLKADDR(L1, MORSE_Complex64_t, L1m, L1n),
        STARPU_VALUE,  &ldl1,                        sizeof(int),
        STARPU_R,             RTBLKADDR(L2, MORSE_Complex64_t, L2m, L2n),
        STARPU_VALUE,  &ldl2,                        sizeof(int),
        STARPU_VALUE,          &IPIV,                      sizeof(int*),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zssssm",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zssssm_cpu_func(void *descr[], void *cl_arg)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *L1;
    int ldl1;
    MORSE_Complex64_t *L2;
    int ldl2;
    int *IPIV;

    A1 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    L1 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    L2 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    starpu_codelet_unpack_args(cl_arg, &m1, &n1, &m2, &n2, &k, &ib, &lda1, &lda2, &ldl1, &ldl2, &IPIV);
    CORE_zssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zssssm, 4, cl_zssssm_cpu_func)

