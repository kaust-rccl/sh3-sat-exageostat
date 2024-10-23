/**
 *
 * @file codelet_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmqr_hetra1 StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
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
void MORSE_TASK_ztsmqr_hetra1(const MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              const MORSE_desc_t *V,  int Vm,  int Vn,  int ldv,
                              const MORSE_desc_t *T,  int Tm,  int Tn,  int ldt)
{
    struct starpu_codelet *codelet = &cl_ztsmqr_hetra1;
    void (*callback)(void*) = options->profiling ? cl_ztsmqr_hetra1_callback : NULL;

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
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &ldwork,            sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztsmqr_hetra1",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztsmqr_hetra1_cpu_func(void *descr[], void *cl_arg)
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

    /* TODO: manage workspace */
    MORSE_Complex64_t *WORK;
    int ldwork;

    A1    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2    = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    V     = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    T     = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[4]);

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k,
                               &ib, &lda1, &lda2, &ldv, &ldt, &ldwork);
    CORE_ztsmqr_hetra1(side, trans, m1, n1, m2, n2, k,
                       ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztsmqr_hetra1, 5, cl_ztsmqr_hetra1_cpu_func)
