/**
 *
 * @file pzlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum parallel algorithm
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
 *  Parallel UU' or L'L operation - dynamic scheduling
 */
void morse_pzlauum(MORSE_enum uplo, MORSE_desc_t *A,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldak, ldam, ldan;
    int tempkm, tempkn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for(n = 0; n < k; n++) {
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n), ldak,
                    1.0, A(n, n), ldan);

                for(m = n+1; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m), ldak,
                             A(k, n), ldak,
                        1.0, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseLeft, uplo, MorseConjTrans, MorseNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k), ldak,
                         A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k), ldak);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (k = 0; k < A->mt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);

            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k), ldam,
                    1.0, A(m, m), ldam);

                for (n = m+1; n < k; n++){
                    ldan = BLKLDD(A, n);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k), ldam,
                             A(n, k), ldan,
                        1.0, A(m, n), ldam);
                }
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseRight, uplo, MorseConjTrans, MorseNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k), ldak);
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
