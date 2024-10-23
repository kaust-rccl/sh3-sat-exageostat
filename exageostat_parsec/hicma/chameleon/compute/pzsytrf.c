/**
 *
 * @file pzsytrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf parallel algorithm
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @date 2014-10-09
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void morse_pzsytrf(MORSE_enum uplo, MORSE_desc_t *A,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn;
    size_t ws_host   = 0;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    RUNTIME_options_ws_alloc( &options, 0, ws_host );

    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);

            MORSE_TASK_zsytrf_nopiv(
                &options,
                MorseLower, tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, MorseLower, MorseTrans, MorseNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                MORSE_TASK_zsyrk(
                    &options,
                    MorseLower, MorseNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k), ldan,
                     1.0, A(n, n), ldan);

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k), ldam,
                               A(n, k), ldan,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(n, k) );
            }

            RUNTIME_iteration_pop(morse);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            MORSE_TASK_zsytrf_nopiv(
                &options,
                MorseUpper,
                tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, MorseUpper, MorseTrans, MorseNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zsyrk(
                    &options,
                    MorseUpper, MorseTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    MORSE_TASK_zgemm(
                        &options,
                        MorseTrans, MorseNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m), ldak,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, m) );
            }

            RUNTIME_iteration_pop(morse);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
