/**
 *
 * @file pzlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
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
/**
 *  Parallel initialization a 2-D array A to BETA on the diagonal and
 *  ALPHA on the offdiagonals.
 */
void morse_pzlaset(MORSE_enum uplo,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                          MORSE_desc_t *A,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int i, j;
    int ldai, ldaj;
    int tempim;
    int tempjm, tempjn;
    int minmn = chameleon_min(A->mt, A->nt);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;

    RUNTIME_options_init(&options, morse, sequence, request);

    if (uplo == MorseLower) {
       for (j = 0; j < minmn; j++){
           tempjm = j == A->mt-1 ? A->m-j*A->mb : A->mb;
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           ldaj = BLKLDD(A, j);
           MORSE_TASK_zlaset(
               &options,
               MorseLower, tempjm, tempjn, alpha, beta,
               A(j, j), ldaj);

           for (i = j+1; i < A->mt; i++){
               tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
               ldai = BLKLDD(A, i);
               MORSE_TASK_zlaset(
                   &options,
                   MorseUpperLower, tempim, tempjn, alpha, alpha,
                   A(i, j), ldai);
           }
       }
    }
    else if (uplo == MorseUpper) {
       for (j = 1; j < A->nt; j++){
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           for (i = 0; i < chameleon_min(j, A->mt); i++){
               tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
               ldai = BLKLDD(A, i);
               MORSE_TASK_zlaset(
                   &options,
                   MorseUpperLower, tempim, tempjn, alpha, alpha,
                   A(i, j), ldai);
           }
       }
       for (j = 0; j < minmn; j++){
           tempjm = j == A->mt-1 ? A->m-j*A->mb : A->mb;
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           ldaj = BLKLDD(A, j);
           MORSE_TASK_zlaset(
               &options,
               MorseUpper, tempjm, tempjn, alpha, beta,
               A(j, j), ldaj);
       }
    }
    else {
       for (i = 0; i < A->mt; i++){
           tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
           ldai = BLKLDD(A, i);
           for (j = 0; j < A->nt; j++){
               tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
               MORSE_TASK_zlaset(
                   &options,
                   MorseUpperLower, tempim, tempjn, alpha, alpha,
                   A(i, j), ldai);
           }
       }
       for (j = 0; j < minmn; j++){
           tempjm = j == A->mt-1 ? A->m-j*A->mb : A->mb;
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           ldaj = BLKLDD(A, j);
           MORSE_TASK_zlaset(
               &options,
               MorseUpperLower, tempjm, tempjn, alpha, beta,
               A(j, j), ldaj);
       }
    }
    RUNTIME_options_finalize(&options, morse);
}
