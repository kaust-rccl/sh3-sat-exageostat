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
 * @brief Chameleon zttmqr Quark codelet
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
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static void
CORE_zttmqr_quark( Quark *quark )
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

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib,
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_zttmqr(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

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
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
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
 *         CORE_ZTTQRT in the first k rows of its array argument V.
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
 *         Workspace array of size LDWORK-by-N1.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK. LDWORK >= max(1,IB).
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
    int ldwork = side == MorseLeft ? ib : nb;

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TTMQR;
    QUARK_Insert_Task(opt->quark, CORE_zttmqr_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &side,  VALUE,
        sizeof(MORSE_enum),              &trans, VALUE,
        sizeof(int),                     &m1,    VALUE,
        sizeof(int),                     &n1,    VALUE,
        sizeof(int),                     &m2,    VALUE,
        sizeof(int),                     &n2,    VALUE,
        sizeof(int),                     &k,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb, RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n), INOUT,
        sizeof(int),                     &lda1,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb, RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n), INOUT | LOCALITY,
        sizeof(int),                     &lda2,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb, RTBLKADDR(V, MORSE_Complex64_t, Vm, Vn),    INPUT | QUARK_REGION_U | QUARK_REGION_D,
        sizeof(int),                     &ldv,   VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb, RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),    INPUT,
        sizeof(int),                     &ldt,   VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb, NULL,          SCRATCH,
        sizeof(int),                     &ldwork, VALUE,
        0);
}
