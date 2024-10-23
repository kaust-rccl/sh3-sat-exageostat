/**
 *
 * @file pztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define T(m,n) T,  m,  n

/**
 *  Parallel tile QR factorization - dynamic scheduling
 */
void morse_pztpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T,
                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldbm;
    int tempkm, tempkn, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = chameleon_max( B->m - L, 1 );
    int maxmt = (maxm % B->mb == 0) ? (maxm / B->mb) : (maxm / B->mb + 1);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * ztsqrt  = A->nb * (ib+1)
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * ztpmqrt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < A->nt; k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);

        for (m = 0; m < maxmt; m++) {
            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
            templm = m == maxmt-1 ? tempmm       : 0;
            ldbm = BLKLDD(B, m);
            /* TT kernel */
            MORSE_TASK_ztpqrt(
                &options,
                tempmm, tempkn, templm, ib, T->nb,
                A(k, k), ldak,
                B(m, k), ldbm,
                T(m, k), T->mb );

            for (n = k+1; n < B->nt; n++) {
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkm, templm, ib, T->nb,
                    B(m, k), ldbm,
                    T(m, k), T->mb,
                    A(k, n), ldak,
                    B(m, n), ldbm );
            }
        }

        maxmt = chameleon_min( B->mt, maxmt+1 );

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
