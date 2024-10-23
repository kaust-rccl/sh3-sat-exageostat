/**
 *
 * @file codelet_zunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr Quark codelet
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
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zunmqr_quark(Quark *quark)
{
    MORSE_enum side;
    MORSE_enum trans;
    int m;
    int n;
    int k;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *C;
    int ldc;
    MORSE_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                         A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zunmqr(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zunmqr overwrites the general complex M-by-N tile C with
 *
 *                    SIDE = 'L'     SIDE = 'R'
 *    TRANS = 'N':      Q * C          C * Q
 *    TRANS = 'C':      Q**H * C       C * Q**H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zgeqrt. Q is of order M if SIDE = 'L' and of order N
 *  if SIDE = 'R'.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg MorseLeft  : apply Q or Q**H from the Left;
 *         @arg MorseRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg MorseNoTrans   :  No transpose, apply Q;
 *         @arg MorseConjTrans :  Transpose, apply Q**H.
 *
 * @param[in] M
 *         The number of rows of the tile C.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile C.  N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *         If SIDE = MorseLeft,  M >= K >= 0;
 *         if SIDE = MorseRight, N >= K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] A
 *         Dimension:  (LDA,K)
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgeqrt in the first k columns of its array argument A.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         If SIDE = MorseLeft,  LDA >= max(1,M);
 *         if SIDE = MorseRight, LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the M-by-N tile C.
 *         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
 *
 * @param[in] LDC
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *         If SIDE = MorseLeft,  LDWORK >= max(1,N);
 *         if SIDE = MorseRight, LDWORK >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
void MORSE_TASK_zunmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_UNMQR;
    QUARK_Insert_Task(opt->quark, CORE_zunmqr_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &side,  VALUE,
        sizeof(MORSE_enum),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(int),                     &k,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT | QUARK_REGION_L,
        sizeof(int),                     &lda,   VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,  RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn), INPUT,
        sizeof(int),                     &ldt,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn), INOUT,
        sizeof(int),                     &ldc,   VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,  NULL,      SCRATCH,
        sizeof(int),                     &nb,    VALUE,
        0);
}
