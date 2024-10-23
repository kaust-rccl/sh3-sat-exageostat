/**
 *
 * @file pztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2011-11-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n

/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 */
void morse_pztradd(MORSE_enum uplo, MORSE_enum trans,
                   MORSE_Complex64_t alpha, MORSE_desc_t *A,
                   MORSE_Complex64_t beta,  MORSE_desc_t *B,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int ldam, ldan, ldbm, ldbn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    switch(uplo){
    case MorseLower:
        if (trans == MorseNoTrans) {
            for (n = 0; n < chameleon_min(B->mt,B->nt); n++) {
                tempnm = n == B->mt-1 ? B->m-n*B->mb : B->mb;
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                ldan = BLKLDD(A, n);
                ldbn = BLKLDD(B, n);

                MORSE_TASK_ztradd(
                    &options,
                    uplo, trans, tempnm, tempnn, B->mb,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn);

                for (m = n+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (n = 0; n < chameleon_min(B->mt,B->nt); n++) {
                tempnm = n == B->mt-1 ? B->m-n*B->mb : B->mb;
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                ldan = BLKLDD(A, n);
                ldbn = BLKLDD(B, n);

                MORSE_TASK_ztradd(
                    &options,
                    uplo, trans, tempnm, tempnn, B->mb,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn);

                for (m = n+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                    ldbm = BLKLDD(B, m);

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        break;
    case MorseUpper:
        if (trans == MorseNoTrans) {
            for (m = 0; m < chameleon_min(B->mt,B->nt); m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                tempmn = m == B->nt-1 ? B->n-m*B->nb : B->nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                MORSE_TASK_ztradd(
                    &options,
                    uplo, trans, tempmm, tempmn, B->mb,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm);

                for (n = m+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (m = 0; m < chameleon_min(B->mt,B->nt); m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                tempmn = m == B->nt-1 ? B->n-m*B->nb : B->nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                MORSE_TASK_ztradd(
                    &options,
                    uplo, trans, tempmm, tempmn, B->mb,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm);

                for (n = m+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        break;
    case MorseUpperLower:
    default:
        if (trans == MorseNoTrans) {
            for (m = 0; m < B->mt; m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (m = 0; m < B->mt; m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                ldbm = BLKLDD(B, m);

                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    MORSE_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
    }

    RUNTIME_options_finalize(&options, morse);
}
