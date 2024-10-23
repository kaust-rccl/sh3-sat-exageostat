/**
 *
 * @file pzlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
/**
 *
 */
/**
 *
 */
void morse_pzlacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    switch (uplo) {
    /*
     *  MorseUpper
     */
    case MorseUpper:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            if (m < A->nt) {
                Y = m == A->nt-1 ? A->n-m*A->nb : A->nb;
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper,
                    X, Y, A->mb,
                    A(m, m), ldam,
                    B(m, m), ldbm);
            }
            for (n = m+1; n < A->nt; n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpperLower,
                    X, Y, A->mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
        break;
    /*
     *  MorseLower
     */
    case MorseLower:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            if (m < A->nt) {
                Y = m == A->nt-1 ? A->n-m*A->nb : A->nb;
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower,
                    X, Y, A->mb,
                    A(m, m), ldam,
                    B(m, m), ldbm);
            }
            for (n = 0; n < chameleon_min(m, A->nt); n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpperLower,
                    X, Y, A->mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
        break;
    /*
     *  MorseUpperLower
     */
    case MorseUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            for (n = 0; n < A->nt; n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpperLower,
                    X, Y, A->mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
