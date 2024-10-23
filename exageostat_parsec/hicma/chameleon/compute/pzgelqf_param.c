/**
 *
 * @file pzgelqf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf_param parallel algorithm
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
#include "libhqr.h"

#define A(m,n) A, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)


/*
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 */
void morse_pzgelqf_param( const libhqr_tree_t *qrtree, MORSE_desc_t *A,
                          MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int K, L;
    int ldak, ldam;
    int tempkmin, tempkm, tempnn, tempmm, temppn;
    int ib;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    if ( D == NULL ) {
        D = A;
    }

    /*
     * zgelqt = A->nb * (ib+1)
     * zunmlq = A->nb * ib
     * ztslqt = A->nb * (ib+1)
     * zttlqt = A->nb * (ib+1)
     * ztsmlq = A->nb * ib
     * zttmlq = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */

    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        ldak = BLKLDD(A, k);

        T = TS;
        /* The number of geqrt to apply */
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            p = qrtree->getm(qrtree, k, i);
            temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
            tempkmin = chameleon_min(tempkm, temppn);

            MORSE_TASK_zgelqt(
                &options,
                tempkm, temppn, ib, T->nb,
                A( k, p), ldak,
                T(k, p), T->mb);
            if ( k < (A->mt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkm, temppn, A->nb,
                    A(k, p), ldak,
                    D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkm, temppn,
                    0., 1.,
                    D(k, p), ldak );
#endif
#endif
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseConjTrans,
                   tempmm, temppn, tempkmin, ib, T->nb,
                    D(k, p), ldak,
                    T(k, p), T->mb,
                    A(m, p), ldam);
            }
            RUNTIME_data_flush( sequence, D(k, p) );
            RUNTIME_data_flush( sequence, T(k, p) );
        }

        /* Setting the order of the tiles */
        libhqr_walk_stepk( qrtree, k, tiles + (k+1) );

        for (i = k+1; i < A->nt; i++) {
            n = tiles[i];
            p = qrtree->currpiv(qrtree, k, n);

            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

            if (qrtree->gettype(qrtree, k, n) == 0) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempnn;
            }

            RUNTIME_data_migrate( sequence, A(k, p),
                                  A->get_rankof( A, k, n ) );
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );

            MORSE_TASK_ztplqt(
                &options,
                tempkm, tempnn, chameleon_min(L, tempkm), ib, T->nb,
                A(k, p), ldak,
                A(k, n), ldak,
                T(k, n), T->mb);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                RUNTIME_data_migrate( sequence, A(m, p),
                                      A->get_rankof( A, m, n ) );
                RUNTIME_data_migrate( sequence, A(m, n),
                                      A->get_rankof( A, m, n ) );

                MORSE_TASK_ztpmlqt(
                    &options,
                    MorseRight, MorseConjTrans,
                    tempmm, tempnn, tempkm, L, ib, T->nb,
                    A(k, n), ldak,
                    T(k, n), T->mb,
                    A(m, p), ldam,
                    A(m, n), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, n) );
            RUNTIME_data_flush( sequence, T(k, n) );
        }

        /* Restore the original location of the tiles */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );
        }

        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
