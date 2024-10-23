/**
 *
 * @file codelet_ztplqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztplqt PaRSEC codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztplqt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int M;
    int N;
    int L;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *WORK;

    parsec_dtd_unpack_args(
        this_task, &M, &N, &L, &ib, &A, &lda, &B, &ldb, &T, &ldt, &WORK );

    CORE_ztplqt( M, N, L, ib,
                 A, lda, B, ldb, T, ldt, WORK );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_ztplqt( const MORSE_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztplqt_parsec, options->priority, "tplqt",
        sizeof(int),   &M,   VALUE,
        sizeof(int),   &N,   VALUE,
        sizeof(int),   &L,   VALUE,
        sizeof(int),   &ib,  VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INOUT,
        sizeof(int),   &lda, VALUE,
        PASSED_BY_REF,  RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int),   &ldb, VALUE,
        PASSED_BY_REF,  RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ), morse_parsec_get_arena_index( T ) | OUTPUT,
        sizeof(int),   &ldt, VALUE,
        sizeof(MORSE_Complex64_t)*(ib+1)*nb, NULL, SCRATCH,
        PARSEC_DTD_ARG_END );
}
