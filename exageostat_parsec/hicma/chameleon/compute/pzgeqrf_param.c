/**
 *
 * @file pzgeqrf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf_param parallel algorithm
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

#define A(m,n)  A, (m), (n)
#define T(m,n)  T, (m), (n)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 */
void morse_pzgeqrf_param( const libhqr_tree_t *qrtree, MORSE_desc_t *A,
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
    int ldap, ldam;
    int tempkmin, tempkn, tempnn, tempmm;
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
     * zgeqrt = A->nb * (ib+1)
     * zunmqr = A->nb * ib
     * ztsqrt = A->nb * (ib+1)
     * zttqrt = A->nb * (ib+1)
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(morse, k);
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* The number of geqrt to apply */
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            m = qrtree->getm(qrtree, k, i);
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            tempkmin = chameleon_min(tempmm, tempkn);
            ldam = BLKLDD(A, m);

            T = TS;

            MORSE_TASK_zgeqrt(
                &options,
                tempmm, tempkn, ib, T->nb,
                A(m, k), ldam,
                T(m, k), T->mb);
            if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempmm, A->nb, A->nb,
                    A(m, k), ldam,
                    D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempmm, A->nb,
                    0., 1.,
                    D(m, k), ldam );
#endif
#endif
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkmin, ib, T->nb,
                    D(m, k), ldam,
                    T(m, k), T->mb,
                    A(m, n), ldam);
            }
            RUNTIME_data_flush( sequence, D(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Setting the order of the tiles */
        libhqr_walk_stepk( qrtree, k, tiles + (k+1) );

        for (i = k+1; i < A->mt; i++) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldap = BLKLDD(A, p);
            ldam = BLKLDD(A, m);

            if (qrtree->gettype(qrtree, k, m) == 0) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempmm;
            }

            RUNTIME_data_migrate( sequence, A(p, k),
                                  A->get_rankof( A, m, k ) );
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );

            MORSE_TASK_ztpqrt(
                &options,
                tempmm, tempkn, chameleon_min(L, tempkn), ib, T->nb,
                A(p, k), ldap,
                A(m, k), ldam,
                T(m, k), T->mb);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                RUNTIME_data_migrate( sequence, A(p, n),
                                      A->get_rankof( A, m, n ) );
                RUNTIME_data_migrate( sequence, A(m, n),
                                      A->get_rankof( A, m, n ) );

                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, A->nb, L, ib, T->nb,
                    A(m, k), ldam,
                    T(m, k), T->mb,
                    A(p, n), ldap,
                    A(m, n), ldam);
            }
            RUNTIME_data_flush( sequence, A(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Restore the original location of the tiles */
        for (n = k; n < A->nt; n++) {
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );
        }

        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
