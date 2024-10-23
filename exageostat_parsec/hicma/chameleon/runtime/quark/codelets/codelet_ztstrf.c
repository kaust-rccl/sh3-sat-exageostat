/**
 *
 * @file codelet_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf Quark codelet
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
#include "coreblas/cblas.h"
#include <math.h>

void CORE_ztstrf_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    int nb;
    MORSE_Complex64_t *U;
    int ldu;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *L;
    int ldl;
    int *IPIV;
    MORSE_Complex64_t *WORK;
    int ldwork;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    MORSE_bool check_info;
    int iinfo;

    int info;

    quark_unpack_args_17(quark, m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, sequence, request, check_info, iinfo);
    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);
    if ( (info != MORSE_SUCCESS) && check_info ) {
        RUNTIME_sequence_flush( (MORSE_context_t*)quark, sequence, request, iinfo+info );
    }
}

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_ztstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular NB-by-N tile U on top of a M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] NB
 *
 * @param[in,out] U
 *         On entry, the NB-by-N upper triangular tile.
 *         On exit, the new factor U from the factorization
 *
 * @param[in] LDU
 *         The leading dimension of the array U.  LDU >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factor L from the factorization
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[in,out] L
 *         On entry, the IB-by-N lower triangular tile.
 *         On exit, the interchanged rows form the tile A in case of pivoting.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,IB).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile U was interchanged with row IPIV(i) of the tile A.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *
 * @param[out] INFO
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
void MORSE_TASK_ztstrf(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *U, int Um, int Un, int ldu,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TSTRF;
    QUARK_Insert_Task(opt->quark, CORE_ztstrf_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(int),                        &nb,            VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(U, MORSE_Complex64_t, Um, Un),                     INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(int),                        &ldu,           VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,    RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),                     OUTPUT,
        sizeof(int),                        &ldl,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,                  SCRATCH,
        sizeof(int),                        &nb,            VALUE,
        sizeof(MORSE_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(MORSE_request_t*),            &(options->request),       VALUE,
        sizeof(MORSE_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
