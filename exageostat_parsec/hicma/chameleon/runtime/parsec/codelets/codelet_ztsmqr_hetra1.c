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
 * @brief Chameleon ztsmqr_hetra1 PaRSEC codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztsmqr_hetra1_parsec( parsec_execution_stream_t *context,
                           parsec_task_t             *this_task )
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

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &V, &ldv, &T, &ldt, &WORK, &ldwork);

    CORE_ztsmqr_hetra1( side, trans, m1, n1, m2, n2, k, ib,
                        A1, lda1, A2, lda2,
                        V, ldv, T, ldt,
                        WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_ztsmqr_hetra1(const MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    int ldwork = side == MorseLeft ? ib : nb;

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztsmqr_hetra1_parsec, options->priority, "tsmqr_hetra1",
        sizeof(MORSE_enum), &side,   VALUE,
        sizeof(MORSE_enum), &trans,  VALUE,
        sizeof(int),        &m1,     VALUE,
        sizeof(int),        &n1,     VALUE,
        sizeof(int),        &m2,     VALUE,
        sizeof(int),        &n2,     VALUE,
        sizeof(int),        &k,      VALUE,
        sizeof(int),        &ib,     VALUE,
        PASSED_BY_REF,       RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n), INOUT,
        sizeof(int),        &lda1,   VALUE,
        PASSED_BY_REF,       RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n), INOUT | AFFINITY,
        sizeof(int),        &lda2,   VALUE,
        PASSED_BY_REF,       RTBLKADDR(V,  MORSE_Complex64_t, Vm,  Vn),  INPUT,
        sizeof(int),        &ldv,    VALUE,
        PASSED_BY_REF,       RTBLKADDR(T,  MORSE_Complex64_t, Tm,  Tn),  INPUT,
        sizeof(int),        &ldt,    VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb, NULL, SCRATCH,
        sizeof(int),        &ldwork, VALUE,
        PARSEC_DTD_ARG_END );
}
