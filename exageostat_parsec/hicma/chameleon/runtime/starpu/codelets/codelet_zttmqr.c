/**
 *
 * @file codelet_zttmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zttmqr StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Dulceneia Becker
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
 *  CORE_zttmqr overwrites the general complex M1-by-N1 tile A1 and
 *  M2-by-N2 tile A2 (N1 == N2) with
 *
 *                        SIDE = 'L'        SIDE = 'R'
 *    TRANS = 'N':         Q * | A1 |       | A1 | * Q
 *                             | A2 |       | A2 |
 *
 *    TRANS = 'C':      Q**H * | A1 |       | A1 | * Q**H
 *                             | A2 |       | A2 |
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zttqrt.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg MorseLeft  : apply Q or Q**H from the Left;
 *         @arg MorseRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg MorseNoTrans   :  No transpose, apply Q;
 *         @arg MorseConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] M1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2. M2 >= 0.
 *         M2 = M1 if side == MorseRight.
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
 *         N2 = N1 if side == MorseLeft.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTTQRT in the first k columns of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-N1 if side == MorseLeft
 *             LDWORK-by-IB if side == MorseRight
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == MorseLeft
 *             LDWORK >= max(1,M1) if side == MorseRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

void MORSE_TASK_zttmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zttmqr;
    void (*callback)(void*) = options->profiling ? cl_zttmqr_callback : NULL;
    int ldwork = side == MorseLeft ? ib : nb;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A1, A1m, A1n);
    MORSE_ACCESS_RW(A2, A2m, A2n);
    MORSE_ACCESS_R(V, Vm, Vn);
    MORSE_ACCESS_R(T, Tm, Tn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &side,              sizeof(MORSE_enum),
        STARPU_VALUE,    &trans,             sizeof(MORSE_enum),
        STARPU_VALUE,    &m1,                sizeof(int),
        STARPU_VALUE,    &n1,                sizeof(int),
        STARPU_VALUE,    &m2,                sizeof(int),
        STARPU_VALUE,    &n2,                sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n),
        STARPU_VALUE,    &lda1,              sizeof(int),
        STARPU_RW,        RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n),
        STARPU_VALUE,    &lda2,              sizeof(int),
        STARPU_R,         RTBLKADDR(V, MORSE_Complex64_t, Vm, Vn),
        STARPU_VALUE,    &ldv,               sizeof(int),
        STARPU_R,         RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),
        STARPU_VALUE,    &ldt,               sizeof(int),
        /* max( ib*nb, 2*ib*nb ) */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &ldwork,            sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, A2->get_rankof(A2, A2m, A2n),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zttmqr",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zttmqr_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum trans;
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
    MORSE_Complex64_t *V;
    int ldv;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *WORK;
    int ldwork;

    A1   = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2   = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    V    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    T    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[4]); /* ib * nb */

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib,
                               &lda1, &lda2, &ldv, &ldt, &ldwork);

    CORE_zttmqr(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zttmqr_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum side;
    MORSE_enum trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    cuDoubleComplex *A1;
    int lda1;
    cuDoubleComplex *A2;
    int lda2;
    cuDoubleComplex *V;
    int ldv;
    cuDoubleComplex *T;
    int ldt;
    cuDoubleComplex *W, *WC;
    int ldwork;
    int ldworkc;

    A1 = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    V  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    T  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]);
    W  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[4]); /* 2*ib*nb */

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib,
                               &lda1, &lda2, &ldv, &ldt, &ldwork);

    WC = W + ib * (side == MorseLeft ? m1 : n1);
    ldworkc = (side == MorseLeft) ? m2 : ib;

    RUNTIME_getStream(stream);

    CUDA_zttmqr(
            side, trans, m1, n1, m2, n2, k, ib,
            A1, lda1, A2, lda2, V, ldv, T, ldt,
            W, ldwork, WC, ldworkc, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zttmqr, 5, cl_zttmqr_cpu_func, cl_zttmqr_cuda_func, STARPU_CUDA_ASYNC)
