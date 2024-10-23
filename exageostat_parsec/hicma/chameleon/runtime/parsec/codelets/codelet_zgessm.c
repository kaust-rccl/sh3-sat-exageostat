/**
 *
 * @file codelet_zgessm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgessm applies the factors L computed by CORE_zgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *         The number of columns of the tile L. K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_zgetrf_incpiv.
 *
 * @param[in] L
 *         The M-by-K lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,M).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */
static inline int
CORE_zgessm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    MORSE_Complex64_t *L;
    int ldl;
    MORSE_Complex64_t *D;
    int ldd;
    MORSE_Complex64_t *A;
    int lda;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &k, &ib, &IPIV, &L, &ldl, &D, &ldd, &A, &lda );

    CORE_zgessm( m, n, k, ib, IPIV, D, ldd, A, lda );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zgessm(const MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb, int *IPIV,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       const MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       const MORSE_desc_t *A, int Am, int An, int lda)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgessm_parsec, options->priority, "gessm",
        sizeof(int),           &m,                                VALUE,
        sizeof(int),           &n,                                VALUE,
        sizeof(int),           &k,                                VALUE,
        sizeof(int),           &ib,                               VALUE,
        sizeof(int*),          &IPIV,                             VALUE,
        PASSED_BY_REF,         RTBLKADDR( L, MORSE_Complex64_t, Lm, Ln ), morse_parsec_get_arena_index( L ) | INPUT,
        sizeof(int),           &ldl,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( D, MORSE_Complex64_t, Dm, Dn ), morse_parsec_get_arena_index( D ) | INPUT,
        sizeof(int),           &ldd,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int),           &lda,                              VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
