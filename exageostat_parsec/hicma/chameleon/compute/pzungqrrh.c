/**
 *
 * @file pzungqrrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqrrh parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  (m),  (n)
#define Q(m,n) Q,  (m),  (n)
#define T(m,n) T,  (m),  (n)
#define T2(m,n) T,  (m),  (n)+(A->nt)
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n) D, ((m)/BS), 0
#else
#define D(m,n) A, (m), (n)
#endif

/**
 *  Parallel construction of Q using tile V (application to identity;
 *  reduction Householder) - dynamic scheduling
 */
void morse_pzungqrrh(MORSE_desc_t *A, MORSE_desc_t *Q,
                     MORSE_desc_t *T,  MORSE_desc_t *D, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, M, RD, lastRD;
    int ldaM, ldam, ldaMRD;
    int ldqM, ldqm, ldqMRD;
    int tempkn, tempMm, tempnn, tempmm, tempMRDm, tempkmin;
    int ib;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);
    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        lastRD = 0;
        for (RD = BS; RD < A->mt-k; RD *= 2)
            lastRD = RD;
        for (RD = lastRD; RD >= BS; RD /= 2) {
            for (M = k; M+RD < A->mt; M += 2*RD) {
                tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                ldqM   = BLKLDD(Q, M   );
                ldqMRD = BLKLDD(Q, M+RD);
                ldaMRD = BLKLDD(A, M+RD);
                for (n = k; n < Q->nt; n++) {
                    tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                    RUNTIME_data_migrate( sequence, Q(M, n),
                                          Q->get_rankof( Q, M+RD, n ) );
                    RUNTIME_data_migrate( sequence, Q(M+RD, n),
                                          Q->get_rankof( Q, M+RD, n ) );

                    /* TT kernel */
                    MORSE_TASK_ztpmqrt(
                        &options,
                        MorseLeft, MorseNoTrans,
                        tempMRDm, tempnn, tempkn, tempMRDm, ib, T->nb,
                        A (M+RD, k), ldaMRD,
                        T2(M+RD, k), T->mb,
                        Q (M,    n), ldqM,
                        Q (M+RD, n), ldqMRD);
                }

                RUNTIME_data_flush( sequence, A (M+RD, k) );
                RUNTIME_data_flush( sequence, T2(M+RD, k) );
            }
        }
        for (M = k; M < A->mt; M += BS) {
            tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
            tempkmin = chameleon_min(tempMm, tempkn);
            ldaM = BLKLDD(A, M);
            ldqM = BLKLDD(Q, M);
            for (m = chameleon_min(M+BS, A->mt)-1; m > M; m--) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldqm = BLKLDD(Q, m);
                ldam = BLKLDD(A, m);

                for (n = k; n < Q->nt; n++) {
                    tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                    RUNTIME_data_migrate( sequence, Q(M, n),
                                          Q->get_rankof( Q, m, n ) );
                    RUNTIME_data_migrate( sequence, Q(m, n),
                                          Q->get_rankof( Q, m, n ) );

                    /* TS kernel */
                    MORSE_TASK_ztpmqrt(
                        &options,
                        MorseLeft, MorseNoTrans,
                        tempmm, tempnn, tempkn, 0, ib, T->nb,
                        A(m, k), ldam,
                        T(m, k), T->mb,
                        Q(M, n), ldqM,
                        Q(m, n), ldqm);
                }
                RUNTIME_data_flush( sequence, A(m, k) );
                RUNTIME_data_flush( sequence, T(m, k) );
            }

#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseLower, tempMm, tempkmin, A->nb,
                A(M, k), ldaM,
                D(M, k), ldaM );
#if defined(CHAMELEON_USE_CUDA)
            MORSE_TASK_zlaset(
                &options,
                MorseUpper, tempMm, tempkmin,
                0., 1.,
                D(M, k), ldaM );
#endif
#endif
            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                /* Restore the original location of the tiles */
                RUNTIME_data_migrate( sequence, Q(M, n),
                                      Q->get_rankof( Q, M, n ) );

                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseNoTrans,
                    tempMm, tempnn,
                    tempkmin, ib, T->nb,
                    D(M, k), ldaM,
                    T(M, k), T->mb,
                    Q(M, n), ldqM);
            }
            RUNTIME_data_flush( sequence, D(M, k) );
            RUNTIME_data_flush( sequence, T(M, k) );
        }
        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
