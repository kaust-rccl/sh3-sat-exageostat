/**
 *
 * @file pzhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile Hermitian matrix-matrix multiplication - dynamic scheduling
 */
void morse_pzhemm(MORSE_enum side, MORSE_enum uplo,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_Complex64_t beta, MORSE_desc_t *C,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldam, ldan, ldak, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    MORSE_Complex64_t zbeta;
    MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for(m = 0; m < C->mt; m++) {
        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
        ldcm = BLKLDD(C, m);
        for(n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            /*
             *  MorseLeft / MorseLower
             */
            if (side == MorseLeft) {
                ldam = BLKLDD(A, m);
                if (uplo == MorseLower) {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = k == C->mt-1 ? C->m-k*C->mb : C->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k), ldam,  /* lda * K */
                                       B(k, n), ldbk,  /* ldb * Y */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == m) {
                                MORSE_TASK_zhemm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                MORSE_TASK_zgemm(
                                    &options,
                                    MorseConjTrans, MorseNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(k, m), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  MorseLeft / MorseUpper
                 */
                else {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = k == C->mt-1 ? C->m-k*C->mb : C->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseConjTrans, MorseNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m), ldak,  /* ldak * X */
                                       B(k, n), ldbk,  /* ldb  * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == m) {
                                MORSE_TASK_zhemm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * K */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                MORSE_TASK_zgemm(
                                    &options,
                                    MorseNoTrans, MorseNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(m, k), ldam,  /* lda * K */
                                           B(k, n), ldbk,  /* ldb * Y */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
            /*
             *  MorseRight / MorseLower
             */
            else {
                ldan = BLKLDD(A, n);
                ldbm = BLKLDD(B, m);
                if (uplo == MorseLower) {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseConjTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldb * K */
                                       A(n, k), ldan,  /* lda * K */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == n) {
                                MORSE_TASK_zhemm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                MORSE_TASK_zgemm(
                                    &options,
                                    MorseNoTrans, MorseNoTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,  /* ldb  * K */
                                           A(k, n), ldak,  /* ldak * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  MorseRight / MorseUpper
                 */
                else {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldb  * K */
                                       A(k, n), ldak,  /* ldak * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == n) {
                                MORSE_TASK_zhemm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                MORSE_TASK_zgemm(
                                    &options,
                                    MorseNoTrans, MorseConjTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,  /* ldb * K */
                                           A(n, k), ldan,  /* lda * K */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
