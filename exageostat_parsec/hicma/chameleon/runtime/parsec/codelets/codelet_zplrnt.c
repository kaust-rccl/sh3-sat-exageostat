/**
 *
 * @file codelet_zplrnt.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zplrnt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &bigM, &m0, &n0, &seed );

    CORE_zplrnt( m, n, A, lda, bigM, m0, n0, seed );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zplrnt( const MORSE_option_t *options,
                        int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplrnt_parsec, options->priority, "zplrnt",
        sizeof(int),       &m,                          VALUE,
        sizeof(int),       &n,                          VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | OUTPUT | AFFINITY,
        sizeof(int),       &lda,                        VALUE,
        sizeof(int),       &bigM,                       VALUE,
        sizeof(int),       &m0,                         VALUE,
        sizeof(int),       &n0,                         VALUE,
        sizeof(unsigned long long int),       &seed,    VALUE,
        PARSEC_DTD_ARG_END );
}
