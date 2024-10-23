/**
 *
 * @file pzgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm parallel algorithm
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

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n
#define C(m, n) C,  m,  n
/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 */
void morse_pzgemm(MORSE_enum transA, MORSE_enum transB,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_Complex64_t beta,  MORSE_desc_t *C,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    MORSE_Complex64_t zbeta;
    MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < C->mt; m++) {
        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
        ldcm = BLKLDD(C, m);
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            /*
             *  A: MorseNoTrans / B: MorseNoTrans
             */
            if (transA == MorseNoTrans) {
                ldam = BLKLDD(A, m);
                if (transB == MorseNoTrans) {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: MorseNoTrans / B: Morse[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Morse[Conj]Trans / B: MorseNoTrans
             */
            else {
                if (transB == MorseNoTrans) {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: Morse[Conj]Trans / B: Morse[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
            RUNTIME_data_flush( sequence, C(m, n) );
        }
        if (transA == MorseNoTrans) {
            for (k = 0; k < A->nt; k++) {
                RUNTIME_data_flush( sequence, A(m, k) );
            }
        } else {
            for (k = 0; k < A->mt; k++) {
                RUNTIME_data_flush( sequence, A(k, m) );
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
