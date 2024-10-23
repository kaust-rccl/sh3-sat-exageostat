/**
 *
 * @file pzlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal parallel algorithm
 *
 * @version 1.0.0
 * @author Dalal Sukkari
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  Parallel scale of a matrix A
 */
void morse_pzlascal(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A,
                    MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int ldam, ldan;
    int minmnt = chameleon_min(A->mt, A->nt);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;

    RUNTIME_options_init(&options, morse, sequence, request);

    switch(uplo) {
    case MorseLower:
        for (n = 0; n < minmnt; n++) {
            tempnm = n == A->mt-1 ? A->m-n*A->mb : A->mb;
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            ldan = BLKLDD(A, n);

            MORSE_TASK_zlascal(
                &options,
                MorseLower, tempnm, tempnn, A->mb,
                alpha, A(n, n), ldan);

            for (m = n+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zlascal(
                    &options,
                    MorseUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case MorseUpper:
        for (m = 0; m < minmnt; m++) {
            tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
            tempmn = m == A->nt-1 ? A->n-m*A->nb : A->nb;
            ldam = BLKLDD(A, m);

            MORSE_TASK_zlascal(
                &options,
                MorseUpper, tempmm, tempmn, A->mb,
                alpha, A(m, m), ldam);

            for (n = m+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                MORSE_TASK_zlascal(
                    &options,
                    MorseUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case MorseUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
            ldam = BLKLDD(A, m);

            for (n = 0; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                MORSE_TASK_zlascal(
                    &options,
                    MorseUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
