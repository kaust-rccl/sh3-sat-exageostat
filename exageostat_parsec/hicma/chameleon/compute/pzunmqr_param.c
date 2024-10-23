/**
 *
 * @file pzunmqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr_param parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-17
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define T(m,n) T,  m,  n
#define D(m,n) D,  m,  n

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void morse_pzunmqr_param(const libhqr_tree_t *qrtree,
                         MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int ldam, ldan, ldbm, ldbp;
    int tempnn, tempkmin, tempmm, tempkn;
    int ib, K, L;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    K = chameleon_min(A->mt, A->nt);

    if (D == NULL) {
        D = A;
    }

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

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    if (side == MorseLeft ) {
        if (trans == MorseConjTrans) {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempmm, tempkmin, A->nb,
                        A(m, k), ldam,
                        D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempmm, tempkmin,
                        0., 1.,
                        D(m, k), ldam );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(m, k), ldam,
                            T(m, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = k+1; i < A->mt; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    if(qrtree->gettype(qrtree, k, m) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(p, n),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(p, n), ldbp,
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
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = A->mt-1; i > k; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    if(qrtree->gettype(qrtree, k, m) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = k; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(p, n),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(p, n), ldbp,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempmm, tempkmin, A->nb,
                        A(m, k), ldam,
                        D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempmm, tempkmin,
                        0., 1.,
                        D(m, k), ldam );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(m, k), ldam,
                            T(m, k), T->mb,
                            B(m, n), ldbm);
                    }

                    RUNTIME_data_flush( sequence, D(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }
    /*
     *  MorseRight / MorseConjTrans
     */
    else {
        if (trans == MorseConjTrans) {
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == A->nt-1 ? A->n - k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = A->nt-1; i > k; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    if( qrtree->gettype(qrtree, k, n) == 0 ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, p),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempnn, tempkmin, A->nb,
                        A(n, k), ldan,
                        D(n, k), ldan );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempnn, tempkmin,
                        0., 1.,
                        D(n, k), ldan );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }
                RUNTIME_iteration_pop(morse);
            }
        }
        /*
         *  MorseRight / MorseNoTrans
         */
        else {
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempnn, tempkmin, A->nb,
                        A(n, k), ldan,
                        D(n, k), ldan );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempnn, tempkmin,
                        0., 1.,
                        D(n, k), ldan );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                /* Setting the order of tiles */
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = k+1; i < A->nt; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    if( qrtree->gettype(qrtree, k, n) == 0 ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->mb;
                        T = TT;
                    }

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, p),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
