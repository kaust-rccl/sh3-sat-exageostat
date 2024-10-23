/**
 *
 * @file pzlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define SA(m,n) SA,  m,  n
#define SB(m,n) SB,  m,  n
/**
 *
 */
/**
 *
 */
void morse_pclag2z(MORSE_desc_t *SA, MORSE_desc_t *B,
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

    for(m = 0; m < SA->mt; m++) {
        X = m == SA->mt-1 ? SA->m-m*SA->mb : SA->mb;
        ldam = BLKLDD(SA, m);
        ldbm = BLKLDD(B, m);
        for(n = 0; n < SA->nt; n++) {
            Y = n == SA->nt-1 ? SA->n-n*SA->nb : SA->nb;
            MORSE_TASK_clag2z(
                &options,
                X, Y, SA->mb,
                SA(m, n), ldam,
                B(m, n), ldbm);
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
