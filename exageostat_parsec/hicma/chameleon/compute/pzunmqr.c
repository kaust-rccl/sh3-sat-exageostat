/**
 *
 * @file pzunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)   D,  k,  0
#else
#define D(k)   D,  k,  k
#endif

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void morse_pzunmqr(MORSE_enum side, MORSE_enum trans,
                   MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldbk, ldam, ldan, ldbm;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib, minMT, minM;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    if (A->m > A->n) {
        minM  = A->n;
        minMT = A->nt;
    } else {
        minM  = A->m;
        minMT = A->mt;
    }

    if (D == NULL) {
        D = A;
    }

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
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

    if (side == MorseLeft ) {
        if (trans == MorseConjTrans) {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = 0; k < minMT; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm   = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempkm, tempkmin, A->nb,
                    A(k, k), ldak,
                    D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempkm, tempkmin,
                    0., 1.,
                    D(k), ldak );
#endif
#endif
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    MORSE_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        D(k), ldak,
                        T(k, k), T->mb,
                        B(k, n), ldbk);
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                for (m = k+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(k, n),
                                              B->get_rankof( B, m, n ) );

                        /* TS kernel */
                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(k, n), ldbk,
                            B(m, n), ldbm);
                    }

                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < B->nt; n++) {
                    RUNTIME_data_migrate( sequence, B(k, n),
                                          B->get_rankof( B, k, n ) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
        /*
         *  MorseLeft / MorseNoTrans
         */
        else {
            for (k = minMT-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkm   = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (m = B->mt-1; m > k; m--) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(k, n),
                                              B->get_rankof( B, m, n ) );

                        /* TS kernel */
                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(k, n), ldbk,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempkm, tempkmin, A->nb,
                    A(k, k), ldak,
                    D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempkm, tempkmin,
                    0., 1.,
                    D(k), ldak );
#endif
#endif
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    RUNTIME_data_migrate( sequence, B(k, n),
                                          B->get_rankof( B, k, n ) );

                    MORSE_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        D(k), ldak,
                        T(k, k), T->mb,
                        B(k, n), ldbk);
                }
                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );
                RUNTIME_iteration_pop(morse);
            }
        }
    }
    /*
     *  MorseRight / MorseConjTrans
     */
    else {
        if (trans == MorseConjTrans) {
            for (k = minMT-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkn   = k == B->nt - 1 ? B->n - k * B->nb : B->nb;
                tempkmin = k == minMT - 1 ? minM - k * A->nb : A->nb;
                ldak = BLKLDD(A, k);
                for (n = B->nt-1; n > k; n--) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, k),
                                              B->get_rankof( B, m, n ) );

                        /* TS kernel */
                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, k), ldbm,
                            B(m, n), ldbm);
                    }

                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempkn, tempkmin, A->nb,
                    A(k, k), ldak,
                    D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempkn, tempkmin,
                    0., 1.,
                    D(k), ldak );
#endif
#endif
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);

                    RUNTIME_data_migrate( sequence, B(m, k),
                                          B->get_rankof( B, m, k ) );

                    MORSE_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        D(k), ldak,
                        T(k, k), T->mb,
                        B(m, k), ldbm);
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                RUNTIME_iteration_pop(morse);
            }
        }
        /*
         *  MorseRight / MorseNoTrans
         */
        else {
            for (k = 0; k < minMT; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkn   = k == B->nt-1 ? B->n-k*B->nb : B->nb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempkn, tempkmin, A->nb,
                    A(k, k), ldak,
                    D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempkn, tempkmin,
                    0., 1.,
                    D(k), ldak );
#endif
#endif
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    MORSE_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        D(k), ldak,
                        T(k, k), T->mb,
                        B(m, k), ldbm);
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                for (n = k+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, k),
                                              B->get_rankof( B, m, n ) );

                        /* TS kernel */
                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, k), ldbm,
                            B(m, n), ldbm);
                    }

                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < B->mt; m++) {
                    RUNTIME_data_migrate( sequence, B(m, k),
                                          B->get_rankof( B, m, k ) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
