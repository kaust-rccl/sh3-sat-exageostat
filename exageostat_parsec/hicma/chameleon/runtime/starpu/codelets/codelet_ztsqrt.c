/**
 *
 * @file codelet_ztsqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsqrt StarPU codelet
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
 * CORE_ztsqrt computes a QR factorization of a rectangular matrix
 * formed by coupling a complex N-by-N upper triangular tile A1
 * on top of a complex M-by-N tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of columns of the tile A2. M >= 0.
 *
 * @param[in] N
 *         The number of rows of the tile A1.
 *         The number of columns of the tiles A1 and A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, all the elements with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

void MORSE_TASK_ztsqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztsqrt;
    void (*callback)(void*) = options->profiling ? cl_ztsqrt_callback : NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A1, A1m, A1n);
    MORSE_ACCESS_RW(A2, A2m, A2n);
    MORSE_ACCESS_W(T, Tm, Tn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n),
        STARPU_VALUE,    &lda1,              sizeof(int),
        STARPU_RW,        RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n),
        STARPU_VALUE,    &lda2,              sizeof(int),
        STARPU_W,         RTBLKADDR(T,  MORSE_Complex64_t, Tm,  Tn ),
        STARPU_VALUE,    &ldt,               sizeof(int),
        /* max( nb * (ib+1), ib * (ib+nb) ) */
        STARPU_SCRATCH,   options->ws_worker,
        /* 2 * ib * (nb+ib) + nb */
        STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, A2->get_rankof(A2, A2m, A2n),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztsqrt",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztsqrt_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *TAU, *WORK;

    A1 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    T  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    TAU= (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* nb + ib*nb */

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda1, &lda2, &ldt, &h_work);

    WORK = TAU + chameleon_max( m, n );
    CORE_ztsqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztsqrt, 4, cl_ztsqrt_cpu_func)
