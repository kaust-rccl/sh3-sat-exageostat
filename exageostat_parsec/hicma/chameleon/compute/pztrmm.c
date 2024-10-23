/**
 *
 * @file pztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm parallel algorithm
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
 *  Parallel tile triangular matrix-matrix multiplication - dynamic scheduling
 */
void morse_pztrmm(MORSE_enum side, MORSE_enum uplo,
                         MORSE_enum trans, MORSE_enum diag,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLeft / MorseUpper / MorseNoTrans
     */
    if (side == MorseLeft) {
        if (uplo == MorseUpper) {
            if (trans == MorseNoTrans) {
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = m+1; k < A->mt; k++) {
                            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                            ldbk = BLKLDD(B, k);
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  MorseLeft / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (m = B->mt-1; m > -1; m--) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            ldbk = BLKLDD(B, k);
                            ldak = BLKLDD(A, k);
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
        /*
         *  MorseLeft / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (m = B->mt-1; m > -1; m--) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            ldbk = BLKLDD(B, k);
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  MorseLeft / MorseLower / Morse[Conj]Trans
             */
            else {
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldam = BLKLDD(A, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m), ldam,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = m+1; k < A->mt; k++) {
                            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                            ldak = BLKLDD(A, k);
                            ldbk = BLKLDD(B, k);
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
    }
    /*
     *  MorseRight / MorseUpper / MorseNoTrans
     */
    else {
        if (uplo == MorseUpper) {
            if (trans == MorseNoTrans) {
                for (n = B->nt-1; n > -1; n--) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            ldak = BLKLDD(A, k);
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  MorseRight / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = n+1; k < A->mt; k++) {
                            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
        /*
         *  MorseRight / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = n+1; k < A->mt; k++) {
                            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                            ldak = BLKLDD(A, k);
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  MorseRight / MorseLower / Morse[Conj]Trans
             */
            else {
                for (n = B->nt-1; n > -1; n--) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n), ldan,  /* lda * tempkm */
                                   B(m, n), ldbm); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zone,  B(m, n), ldbm);
                        }
                    }
                }
            }
        }
    }

    RUNTIME_options_finalize(&options, morse);
}
