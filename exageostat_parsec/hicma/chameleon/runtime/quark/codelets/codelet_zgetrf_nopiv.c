/**
 *
 * @file codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv Quark codelet
 *
 * @version 1.0.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2013-02-01
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zgetrf_nopiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    int iinfo;
    int info;

    quark_unpack_args_8(quark, m, n, ib, A, lda, sequence, request, iinfo);
    CORE_zgetrf_nopiv(m, n, ib, A, lda, &info);
    if ( info != MORSE_SUCCESS ) {
        RUNTIME_sequence_flush( (MORSE_context_t*)quark, sequence, request, iinfo+info );
    }
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgetrf_nopiv computes an LU factorization of a general diagonal
 *  dominant M-by-N matrix A witout pivoting.
 *
 *  The factorization has the form
 *     A = L * U
 *  where L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 3 BLAS version of the algorithm.
 *  WARNING: Your matrix need to be diagonal dominant if you want to call this
 *  routine safely.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  @param[in] IB
 *          The block size to switch between blocked and unblocked code.
 *
 *  @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  @param[in] LDA
 *          The leading dimension of the array A.  LDA >= max(1,M).
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
void MORSE_TASK_zgetrf_nopiv(const MORSE_option_t *options,
                             int m, int n, int ib, int nb,
                             const MORSE_desc_t *A, int Am, int An, int lda,
                             int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(
        opt->quark, CORE_zgetrf_nopiv_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                     INOUT,
        sizeof(int),                        &lda,           VALUE,
        sizeof(MORSE_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(MORSE_request_t*),            &(options->request),       VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
