/**
 *
 * @file pztpgqrt.c
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

#define V1(m,n) V1,  m,  n
#define T1(m,n) T1,  m,  n
#define V2(m,n) V2,  m,  n
#define T2(m,n) T2,  m,  n
#define Q1(m,n) Q1,  m,  n
#define Q2(m,n) Q2,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)    D, k, 0
#else
#define D(k)    V1, k, k
#endif

/**
 *  Parallel tile QR factorization - dynamic scheduling
 */
void morse_pztpgqrt( int L,
                     MORSE_desc_t *V1, MORSE_desc_t *T1,
                     MORSE_desc_t *V2, MORSE_desc_t *T2,
                     MORSE_desc_t *Q1, MORSE_desc_t *Q2,
                     MORSE_desc_t *D,
                     MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldvk, ldvm;
    int ldqk, ldqm;
    int tempkm, tempkn, tempkk, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = chameleon_max( Q2->m - L, 1 );
    int maxmt = (maxm % Q2->mb == 0) ? (maxm / Q2->mb) : (maxm / Q2->mb + 1);
    int maxmtk;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;
    /*
     * ztpmqrt = Q1->nb * ib
     */
    ws_worker = Q1->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * ztpmqrt = 2 * Q1->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * Q1->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = V1->nt-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == V1->mt-1 ? V1->m-k*V1->mb : V1->mb;
        tempkk = k == V1->nt-1 ? V1->n-k*V1->nb : V1->nb;
        tempkn = k == Q1->nt-1 ? Q1->n-k*Q1->nb : Q1->nb;
        ldvk = BLKLDD(V1, k);
        ldqk = BLKLDD(Q1, k);

        /* Equivalent to the tsmqr step on Q1,Q2 */
        maxmtk = chameleon_min( Q2->mt, maxmt+k ) - 1;
        for (m = maxmtk; m > -1; m--) {
            tempmm = m == Q2->mt-1 ? Q2->m-m*Q2->mb : Q2->mb;
            templm = m == maxmtk ? tempmm : 0;
            ldvm = BLKLDD(V2, m);
            ldqm = BLKLDD(Q2, m);

            for (n = k; n < Q2->nt; n++) {
                tempnn = n == Q2->nt-1 ? Q2->n-n*Q2->nb : Q2->nb;
                /* TT kernel */
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseNoTrans,
                    tempmm, tempnn, tempkn, templm, ib, T2->nb,
                    V2(m, k), ldvm,
                    T2(m, k), T2->mb,
                    Q1(k, n), ldqk,
                    Q2(m, n), ldqm );
            }
        }

        for (m = Q1->mt - 1; m > k; m--) {
            tempmm = m == Q1->mt-1 ? Q1->m-m*Q1->mb : Q1->mb;
            ldvm = BLKLDD(V1, m);
            ldqm = BLKLDD(Q1, m);
            for (n = k; n < Q1->nt; n++) {
                tempnn = n == Q1->nt-1 ? Q1->n-n*Q1->nb : Q1->nb;
                /* TS kernel */
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseNoTrans,
                    tempmm, tempnn, tempkn, 0, ib, T1->nb,
                    V1(m, k), ldvm,
                    T1(m, k), T1->mb,
                    Q1(k, n), ldqk,
                    Q1(m, n), ldqm );
            }
        }

#if defined(CHAMELEON_COPY_DIAG)
        MORSE_TASK_zlacpy(
            &options,
            MorseLower, tempkm, tempkk, V1->nb,
            V1(k, k), ldvk,
            D(k), ldvk );
#if defined(CHAMELEON_USE_CUDA)
        MORSE_TASK_zlaset(
            &options,
            MorseUpper, tempkm, tempkk,
            0., 1.,
            D(k), ldvk );
#endif
#endif
        for (n = k; n < Q1->nt; n++) {
            tempnn = n == Q1->nt-1 ? Q1->n-n*Q1->nb : Q1->nb;
            MORSE_TASK_zunmqr(
                &options,
                MorseLeft, MorseNoTrans,
                tempkm, tempnn, tempkk, ib, T1->nb,
                D(k), ldvk,
                T1(k, k), T1->mb,
                Q1(k, n), ldqk);
        }

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
