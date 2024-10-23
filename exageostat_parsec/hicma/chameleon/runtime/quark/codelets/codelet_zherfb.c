/**
 *
 * @file codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb Quark codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zherfb_quark(Quark *quark)
{
    MORSE_enum uplo;
    int n;
    int k;
    int ib;
    int nb;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *C;
    int ldc;
    MORSE_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_13(quark, uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zherfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
}

void MORSE_TASK_zherfb(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zherfb_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo, VALUE,
        sizeof(int),                       &n,    VALUE,
        sizeof(int),                       &k,    VALUE,
        sizeof(int),                       &ib,   VALUE,
        sizeof(int),                       &nb,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An), (uplo == MorseUpper) ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(int),                       &lda,  VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,    RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn), INPUT,
        sizeof(int),                       &ldt,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn), (uplo == MorseUpper) ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                       &ldc,  VALUE,
        sizeof(MORSE_Complex64_t)*2*nb*nb,  NULL, SCRATCH,
        sizeof(int),                       &nb,   VALUE,
        0);
}
