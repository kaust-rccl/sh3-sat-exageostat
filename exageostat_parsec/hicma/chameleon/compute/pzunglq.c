/**
 *
 * @file pzunglq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define Q(m,n) Q,  m,  n
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)   D,  k,  0
#else
#define D(k)   D,  k,  k
#endif

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void morse_pzunglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldqm;
    int tempnn, tempmm, tempkmin, tempkn;
    int tempAkm, tempAkn;
    int ib, minMT;

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
     * zunmlq = A->nb * ib
     * ztsmlq = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = minMT-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempAkm  = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempAkn  = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempkmin = chameleon_min( tempAkn, tempAkm );
        tempkn   = k == Q->nt-1 ? Q->n-k*Q->nb : Q->nb;
        ldak = BLKLDD(A, k);
        for (n = Q->nt-1; n > k; n--) {
            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;
            for (m = 0; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                ldqm = BLKLDD(Q, m);

                RUNTIME_data_migrate( sequence, Q(m, k),
                                      Q->get_rankof( Q, m, n ) );

                /* TS kernel */
                MORSE_TASK_ztpmlqt(
                    &options,
                    MorseRight, MorseNoTrans,
                    tempmm, tempnn, tempAkm, 0, ib, T->nb,
                    A(k, n), ldak,
                    T(k, n), T->mb,
                    Q(m, k), ldqm,
                    Q(m, n), ldqm);
            }
            RUNTIME_data_flush( sequence, A(k, n) );
            RUNTIME_data_flush( sequence, T(k, n) );
        }
#if defined(CHAMELEON_COPY_DIAG)
        MORSE_TASK_zlacpy(
            &options,
            MorseUpper, tempkmin, tempkn, A->nb,
            A(k, k), ldak,
            D(k), ldak );
#if defined(CHAMELEON_USE_CUDA)
        MORSE_TASK_zlaset(
            &options,
            MorseLower, tempkmin, tempkn,
            0., 1.,
            D(k), ldak );
#endif
#endif
        for (m = k; m < Q->mt; m++) {
            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
            ldqm = BLKLDD(Q, m);

            RUNTIME_data_migrate( sequence, Q(m, k),
                                  Q->get_rankof( Q, m, k ) );

            MORSE_TASK_zunmlq(
                &options,
                MorseRight, MorseNoTrans,
                tempmm, tempkn, tempkmin, ib, T->nb,
                D(k), ldak,
                T(k, k), T->mb,
                Q(m, k), ldqm);
        }
        RUNTIME_data_flush( sequence, D(k)    );
        RUNTIME_data_flush( sequence, T(k, k) );

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
