/**
 *
 * @file pzgeqrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)   D,  k,  0
#else
#define D(k)   D,  k,  k
#endif

/**
 *  Parallel tile QR factorization - dynamic scheduling
 */
void morse_pzgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempnn, tempmm;
    int ib;
    int minMNT = chameleon_min(A->mt, A->nt);

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
     * ztsmqr = A->nb * ib
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

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < minMNT; k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);
        MORSE_TASK_zgeqrt(
            &options,
            tempkm, tempkn, ib, T->nb,
            A(k, k), ldak,
            T(k, k), T->mb);
        if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseLower, A->mb, A->nb, A->nb,
                A(k, k), ldak,
                D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
            MORSE_TASK_zlaset(
                &options,
                MorseUpper, A->mb, A->nb,
                0., 1.,
                D(k), ldak );
#endif
#endif
        }
        for (n = k+1; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zunmqr(
                &options,
                MorseLeft, MorseConjTrans,
                tempkm, tempnn, tempkm, ib, T->nb,
                D(k), ldak,
                T(k, k), T->mb,
                A(k, n), ldak);
        }
        RUNTIME_data_flush( sequence, D(k)    );
        RUNTIME_data_flush( sequence, T(k, k) );

        for (m = k+1; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            RUNTIME_data_migrate( sequence, A(k, k),
                                  A->get_rankof( A, m, k ) );

            /* TS kernel */
            MORSE_TASK_ztpqrt(
                &options,
                tempmm, tempkn, 0, ib, T->nb,
                A(k, k), ldak,
                A(m, k), ldam,
                T(m, k), T->mb);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                RUNTIME_data_migrate( sequence, A(k, n),
                                      A->get_rankof( A, m, n ) );

                /* TS kernel */
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, A->nb, 0, ib, T->nb,
                    A(m, k), ldam,
                    T(m, k), T->mb,
                    A(k, n), ldak,
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

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
