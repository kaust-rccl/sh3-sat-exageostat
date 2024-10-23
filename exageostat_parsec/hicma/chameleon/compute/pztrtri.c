/**
 *
 * @file pztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 */
void morse_pztrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldam, ldak;
    int tempkn, tempkm, tempmm, tempnn;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, diag,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k), ldak,
                           A(m, k), ldam);
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                for (n = 0; n < k; n++) {
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(m, k) );
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, diag,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, diag,
                tempkn, A->mb,
                A(k, k), ldak, A->nb*k);

            RUNTIME_iteration_pop(morse);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, diag,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k), ldak,
                           A(k, n), ldak);
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                for (m = 0; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, n) );
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, diag,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, diag,
                tempkm, A->mb,
                A(k, k), ldak, A->mb*k);

            RUNTIME_iteration_pop(morse);
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
