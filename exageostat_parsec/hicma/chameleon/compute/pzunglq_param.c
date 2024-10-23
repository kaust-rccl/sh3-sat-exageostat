/**
 *
 * @file pzunglq_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq_param parallel algorithm
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

#define A(m,n) A, (m), (n)
#define Q(m,n) Q, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)

/**
 *  Parallel construction of Q using tile V - dynamic scheduling
 */
void morse_pzunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
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
    int ldak, ldqm;
    int tempkm, tempkmin, temppn, tempnn, tempmm;
    int ib;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

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

    tiles = (int*)calloc( qrtree->mt, sizeof(int));

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        ldak = BLKLDD(A, k);

        /* Setting the order of the tiles*/
        libhqr_walk_stepk(qrtree, k, tiles + (k+1));

        for (i = A->nt-1; i > k; i--) {
            n = tiles[i];
            p = qrtree->currpiv(qrtree, k, n);

            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

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
            for (m = k; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                ldqm = BLKLDD(Q, m);

                RUNTIME_data_migrate( sequence, Q(m, p),
                                      Q->get_rankof( Q, m, n ) );
                RUNTIME_data_migrate( sequence, Q(m, n),
                                      Q->get_rankof( Q, m, n ) );

                MORSE_TASK_ztpmlqt(
                    &options,
                    MorseRight, MorseNoTrans,
                    tempmm, tempnn, tempkm, L, ib, T->nb,
                    A(k, n), ldak,
                    T(k, n), T->mb,
                    Q(m, p), ldqm,
                    Q(m, n), ldqm);
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
            for (m = k; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                ldqm = BLKLDD(Q, m);

                RUNTIME_data_migrate( sequence, Q(m, p),
                                      Q->get_rankof( Q, m, p ) );

                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseNoTrans,
                    tempmm, temppn, tempkmin, ib, T->nb,
                    D(k, p), ldak,
                    T(k, p), T->mb,
                    Q(m, p), ldqm);
            }
            RUNTIME_data_flush( sequence, D(k, p) );
            RUNTIME_data_flush( sequence, T(k, p) );
        }

        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
