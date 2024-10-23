/**
 *
 * @file pzgelqfrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqfrh parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A,  (m),  (n)
#define T(m,n)  T,  (m),  (n)
#define T2(m,n) T,  (m),  (n)+A->nt
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n) D, ((n)/BS), 0
#else
#define D(m,n) A,  (m),  (n)
#endif

/*
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 */
void morse_pzgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, N, RD;
    int ldak, ldam;
    int tempkmin, tempkm, tempNn, tempnn, tempmm, tempNRDn;
    int ib;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

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
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        ldak = BLKLDD(A, k);
        for (N = k; N < A->nt; N += BS) {
            tempNn = N == A->nt-1 ? A->n-N*A->nb : A->nb;
            tempkmin = chameleon_min(tempkm, tempNn);
            MORSE_TASK_zgelqt(
                &options,
                tempkm, tempNn, ib, T->nb,
                A(k, N), ldak,
                T(k, N), T->mb);
            if ( k < (A->mt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkm, tempNn, A->nb,
                    A(k, N), ldak,
                    D(k, N), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkm, tempNn,
                    0., 1.,
                    D(k, N), ldak );
#endif
#endif
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseConjTrans,
                    tempmm, tempNn, tempkmin, ib, T->nb,
                    D(k, N), ldak,
                    T(k, N), T->mb,
                    A(m, N), ldam);
            }
            RUNTIME_data_flush( sequence, D(k, N) );
            RUNTIME_data_flush( sequence, T(k, N) );

            for (n = N+1; n < chameleon_min(N+BS, A->nt); n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                RUNTIME_data_migrate( sequence, A(k, N),
                                      A->get_rankof( A, k, n ) );

                /* TS kernel */
                MORSE_TASK_ztplqt(
                    &options,
                    tempkm, tempnn, 0, ib, T->nb,
                    A(k, N), ldak,
                    A(k, n), ldak,
                    T(k, n), T->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);

                    RUNTIME_data_migrate( sequence, A(m, N),
                                          A->get_rankof( A, m, n ) );

                    MORSE_TASK_ztpmlqt(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, tempnn, tempkm, 0, ib, T->nb,
                        A(k, n), ldak,
                        T(k, n), T->mb,
                        A(m, N), ldam,
                        A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, n) );
                RUNTIME_data_flush( sequence, T(k, n) );
            }
        }
        for (RD = BS; RD < A->nt-k; RD *= 2) {
            for (N = k; N+RD < A->nt; N += 2*RD) {
                tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;

                RUNTIME_data_migrate( sequence, A(k, N),
                                      A->get_rankof( A, k, N+RD ) );
                RUNTIME_data_migrate( sequence, A(k, N+RD),
                                      A->get_rankof( A, k, N+RD ) );

                /* TT kernel */
                MORSE_TASK_ztplqt(
                    &options,
                    tempkm, tempNRDn, chameleon_min(tempNRDn, tempkm), ib, T->nb,
                    A (k, N   ), ldak,
                    A (k, N+RD), ldak,
                    T2(k, N+RD), T->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam   = BLKLDD(A, m );

                    RUNTIME_data_migrate( sequence, A(m, N),
                                          A->get_rankof( A, m, N+RD ) );
                    RUNTIME_data_migrate( sequence, A(m, N+RD),
                                          A->get_rankof( A, m, N+RD ) );

                    MORSE_TASK_ztpmlqt(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, tempNRDn, tempkm, tempNRDn, ib, T->nb,
                        A (k, N+RD), ldak,
                        T2(k, N+RD), T->mb,
                        A (m, N   ), ldam,
                        A (m, N+RD), ldam);
                }
                RUNTIME_data_flush( sequence, A (k, N+RD) );
                RUNTIME_data_flush( sequence, T2(k, N+RD) );
            }
        }

        /* Restore the original location of the tiles */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );
        }

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
