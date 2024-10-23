/**
 *
 * @file codelet_ztslqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztslqt Quark codelet
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

void CORE_ztslqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *TAU;
    MORSE_Complex64_t *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_ztslqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_ztslqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M
 *  lower triangular tile A1 and a complex M-by-N tile A2:
 *
 *    | A1 A2 | = L * Q
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
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

void MORSE_TASK_ztslqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TSLQT;
    QUARK_Insert_Task(opt->quark, CORE_ztslqt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n), INOUT | QUARK_REGION_L | QUARK_REGION_D,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n), INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,    RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),    OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(MORSE_Complex64_t)*nb,       NULL,          SCRATCH,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        0);
}
