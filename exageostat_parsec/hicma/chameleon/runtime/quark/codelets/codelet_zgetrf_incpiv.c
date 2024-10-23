/**
 *
 * @file codelet_zgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv Quark codelet
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

void CORE_zgetrf_incpiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    int *IPIV;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    MORSE_bool check_info;
    int iinfo;

    int info;

    quark_unpack_args_10(quark, m, n, ib, A, lda, IPIV, sequence, request, check_info, iinfo);
    CORE_zgetrf_incpiv(m, n, ib, A, lda, IPIV, &info);
    if ( (info != MORSE_SUCCESS) && check_info ) {
        RUNTIME_sequence_flush( (MORSE_context_t*)quark, sequence, request, iinfo+info );
    }
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgetrf_incpiv computes an LU factorization of a general M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factors L and U from the factorization
 *         A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile was interchanged with row IPIV(i).
 *
 * @param[out] INFO
 *         See returned value.
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 */
void MORSE_TASK_zgetrf_incpiv(const MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(opt->quark, CORE_zgetrf_incpiv_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                     INOUT,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(MORSE_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(MORSE_request_t*),            &(options->request),       VALUE,
        sizeof(MORSE_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);

    (void)L; (void)Lm; (void)Ln; (void)ldl;
}
