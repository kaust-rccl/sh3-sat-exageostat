/**
 *
 * @file pzplgsy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgsy parallel algorithm
 *
 * @version 1.0.0
 * @comment This file is a copy of pzplgsy.c,
            wich has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Rade Mathis
 * @author Florent Pruvost
 * @date 2016-08-01
 * @precisions normal z -> c d s
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  morse_pzplgsy - Generate a random symmetric (positive definite if 'bump' is large enough) half-matrix by tiles.
 */
void morse_pzplgsy( MORSE_Complex64_t bump, MORSE_enum uplo, MORSE_desc_t *A,
                    unsigned long long int seed,
                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);

        /*
         * MorseLower
         */
        if (uplo == MorseLower) {
            for (n = 0; n <= m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                MORSE_TASK_zplgsy(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        /*
         * MorseUpper
         */
        else if (uplo == MorseUpper) {
            for (n = m; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                MORSE_TASK_zplgsy(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        /*
         * MorseUpperLower
         */
        else {
            for (n = 0; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                MORSE_TASK_zplgsy(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
