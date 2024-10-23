/**
 *
 * @file pzgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv parallel algorithm
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
//ALLOC_WS :  ib*L->nb
//WS_ADD :  ib*L->nb
#include "control/common.h"
#include <stdlib.h>

#define A(_m_,_n_) A, _m_, _n_
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)   D, k, 0
#else
#define D(k)   A, k, k
#endif
#define L(_m_,_n_) L,  _m_,  _n_
#define IPIV(_m_,_n_) &(IPIV[(int64_t)A->mb*((int64_t)(_m_)+(int64_t)A->mt*(int64_t)(_n_))])

/**
 *  Parallel tile LU factorization - dynamic scheduling
 */
void morse_pzgetrf_incpiv( MORSE_desc_t *A, MORSE_desc_t *L, MORSE_desc_t *D, int *IPIV,
                           MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int minMNT = chameleon_min(A->mt, A->nt);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zgetrf_incpiv = 0
     * zgessm        = 0
     * ztstrf        = A->mb * ib
     * zssssm        = 0
     */
    ws_worker = A->mb * ib;

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < minMNT; k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);
        MORSE_TASK_zgetrf_incpiv(
            &options,
            tempkm, tempkn, ib, L->nb,
            A(k, k), ldak,
            L(k, k), L->mb,
            IPIV(k, k),
            k == A->mt-1, A->nb*k);

        if ( k < (minMNT-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower, tempkm, tempkn, A->nb,
                A(k, k), ldak,
                D(k), ldak);
#endif
        }

        for (n = k+1; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zgessm(
                &options,
                tempkm, tempnn, tempkm, ib, L->nb,
                IPIV(k, k),
                L(k, k), L->mb,
                D(k), ldak,
                A(k, n), ldak);
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            MORSE_TASK_ztstrf(
                &options,
                tempmm, tempkn, ib, L->nb,
                A(k, k), ldak,
                A(m, k), ldam,
                L(m, k), L->mb,
                IPIV(m, k),
                m == A->mt-1, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zssssm(
                    &options,
                    A->nb, tempnn, tempmm, tempnn, A->nb, ib, L->nb,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    L(m, k), L->mb,
                    A(m, k), ldam,
                    IPIV(m, k));
            }
        }

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
