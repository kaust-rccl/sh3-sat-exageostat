/**
 *
 * @file pzunmlq_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq_param parallel algorithm
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
 *  Parallel application of Q using tile V - LQ factorization - dynamic scheduling
 */
void morse_pzunmlq_param(const libhqr_tree_t *qrtree,
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
    int ldbm, ldak, ldbp;
    int tempnn, temppn, tempkmin, tempmm, tempkm;
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
     * zunmlq = A->nb * ib
     * ztsmlq = A->nb * ib
     * zttmlq = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    if (side == MorseLeft ) {
        if (trans == MorseNoTrans) {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            temppn, tempnn, tempkmin, ib, T->nb,
                            D(k, p), ldak,
                            T(k, p), T->mb,
                            B(p, n), ldbp);
                    }

                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = k+1; i < A->nt; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbp = BLKLDD(B, p);
                    ldbm = BLKLDD(B, m);

                    if(qrtree->gettype(qrtree, k, m) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->nb;
                        T = TT;
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(p, n),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmlqt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkm, chameleon_min( L, tempnn ), ib, T->nb,
                            A(k, m), ldak,
                            T(k, m), T->mb,
                            B(p, n), ldbp,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(k, m) );
                    RUNTIME_data_flush( sequence, T(k, m) );
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
         *  MorseLeft / MorseConjTrans
         */
        else {
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = A->nt-1; i > k; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbp = BLKLDD(B, p);
                    ldbm = BLKLDD(B, m);

                    if(qrtree->gettype(qrtree, k, m) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->nb;
                        T = TT;
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(p, n),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmlqt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkm, chameleon_min(L, tempnn), ib, T->nb,
                            A(k, m), ldak,
                            T(k, m), T->mb,
                            B(p, n), ldbp,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(k, m) );
                    RUNTIME_data_flush( sequence, T(k, m) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(p, n),
                                              B->get_rankof( B, p, n ) );

                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            temppn, tempnn, tempkmin, ib, T->nb,
                            D(k, p), ldak,
                            T(k, p), T->mb,
                            B(p, n), ldbp);
                    }

                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }
    /*
     *  MorseRight / MorseNoTrans
     */
    else {
        if (trans == MorseNoTrans) {
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                /* Setting the order of the tiles*/
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = A->nt-1; i > k; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    if(qrtree->gettype(qrtree, k, n) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempnn;
                        T = TT;
                    }
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, p),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmlqt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkm, L, ib, T->nb,
                            A(k, n), ldak,
                            T(k, n), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(k, n) );
                    RUNTIME_data_flush( sequence, T(k, n) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                        RUNTIME_data_migrate( sequence, B(m, p),
                                              B->get_rankof( B, m, p ) );

                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p), ldak,
                            T(k, p), T->mb,
                            B(m, p), ldbm);
                    }

                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
        /*
         *  MorseRight / MorseConjTrans
         */
        else {
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p), ldak,
                            T(k, p), TS->mb,
                            B(m, p), ldbm);
                    }

                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                /* Setting the order of tiles */
                libhqr_walk_stepk(qrtree, k, tiles + (k+1));

                for (i = k+1; i < A->nt; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    if(qrtree->gettype(qrtree, k, n) == 0){
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempnn;
                        T = TT;
                    }

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        RUNTIME_data_migrate( sequence, B(m, p),
                                              B->get_rankof( B, m, n ) );
                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        MORSE_TASK_ztpmlqt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkm, L, ib, T->nb,
                            A(k, n), ldak,
                            T(k, n), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(k, n) );
                    RUNTIME_data_flush( sequence, T(k, n) );
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
