/**
 *
 * @file pzungqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr_param parallel algorithm
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

#define A(m,n) A, m, n
#define Q(m,n) Q, m, n
#define T(m,n) T, m, n
#define D(m,n) D, m, n

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void morse_pzungqr_param(const libhqr_tree_t *qrtree,
                         MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p, L;
    int ldam, ldqm, ldqp;
    int tempmm, tempnn, tempkmin, tempkn;
    int ib, minMT;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    if (A->m > A->n) {
        minMT = A->nt;
    } else {
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

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */

    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = minMT-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* Setting the order of tiles */
        libhqr_walk_stepk(qrtree, k, tiles + (k+1));

        for (i = Q->mt-1; i > k; i--) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
            ldam = BLKLDD(A, m);
            ldqm = BLKLDD(Q, m);
            ldqp = BLKLDD(Q, p);

            if(qrtree->gettype(qrtree, k , m) == 0) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempmm;
            }

            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                RUNTIME_data_migrate( sequence, Q(p, n),
                                      Q->get_rankof( Q, m, n ) );
                RUNTIME_data_migrate( sequence, Q(m, n),
                                      Q->get_rankof( Q, m, n ) );

                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseNoTrans,
                    tempmm, tempnn, tempkn, L, ib, T->nb,
                    A(m, k), ldam,
                    T(m, k), T->mb,
                    Q(p, n), ldqp,
                    Q(m, n), ldqm);
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
            ldqm = BLKLDD(Q, m);

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

            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                /* Restore the original location of the tiles */
                RUNTIME_data_migrate( sequence, Q(m, n),
                                      Q->get_rankof( Q, m, n ) );

                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseNoTrans,
                    tempmm, tempnn, tempkmin, ib, T->nb,
                    D(m, k), ldam,
                    T(m, k), T->mb,
                    Q(m, n), ldqm);
            }
            RUNTIME_data_flush( sequence, D(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
