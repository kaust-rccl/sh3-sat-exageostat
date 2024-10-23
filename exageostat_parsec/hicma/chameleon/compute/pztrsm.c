/**
 *
 * @file pztrsm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
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
#define B(m,n) B,  m,  n
/**
 *  Parallel tile triangular solve - dynamic scheduling
 */
void morse_pztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum trans, MORSE_enum diag,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    MORSE_Complex64_t zone       = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone      = (MORSE_Complex64_t)-1.0;
    MORSE_Complex64_t minvalpha  = (MORSE_Complex64_t)-1.0 / alpha;
    MORSE_Complex64_t lalpha;

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
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m-(B->mt-1)*B->mb : B->mb;
                    ldak = BLKLDD(A, B->mt-1-k);
                    ldbk = BLKLDD(B, B->mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(B->mt-1-k, B->mt-1-k), ldak,  /* lda * tempkm */
                                    B(B->mt-1-k,        n), ldbk); /* ldb * tempnn */
                    }
                    RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-k) );
                    for (m = k+1; m < B->mt; m++) {
                        ldam = BLKLDD(A, B->mt-1-m);
                        ldbm = BLKLDD(B, B->mt-1-m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                B->mb, tempnn, tempkm, A->mb,
                                mzone,  A(B->mt-1-m, B->mt-1-k), ldam,
                                        B(B->mt-1-k, n       ), ldbk,
                                lalpha, B(B->mt-1-m, n       ), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(B->mt-1-m, B->mt-1-k) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, B(B->mt-1-k, n) );
                    }
                }
            }
            /*
             *  MorseLeft / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(k, k), ldak,
                                    B(k, n), ldbk);
                    }
                    RUNTIME_data_flush( sequence, A(k, k) );
                    for (m = k+1; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  A(k, m), ldak,
                                        B(k, n), ldbk,
                                lalpha, B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(k, m) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, B(k, n) );
                    }

                }
            }
        }
        /*
         *  MorseLeft / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(k, k), ldak,
                                    B(k, n), ldbk);
                    }
                    RUNTIME_data_flush( sequence, A(k, k) );
                    for (m = k+1; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldam = BLKLDD(A, m);
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  A(m, k), ldam,
                                        B(k, n), ldbk,
                                lalpha, B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(m, k) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, B(k, n) );
                    }
                }
            }
            /*
             *  MorseLeft / MorseLower / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m-(B->mt-1)*B->mb : B->mb;
                    ldak = BLKLDD(A, B->mt-1-k);
                    ldbk = BLKLDD(B, B->mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(B->mt-1-k, B->mt-1-k), ldak,
                                    B(B->mt-1-k,        n), ldbk);
                    }
                    RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-k) );
                    for (m = k+1; m < B->mt; m++) {
                        ldbm = BLKLDD(B, B->mt-1-m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                B->mb, tempnn, tempkm, A->mb,
                                mzone,  A(B->mt-1-k, B->mt-1-m), ldak,
                                        B(B->mt-1-k, n       ), ldbk,
                                lalpha, B(B->mt-1-m, n       ), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-m) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, B(B->mt-1-k, n) );
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
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            lalpha, A(k, k), ldak,  /* lda * tempkn */
                                    B(m, k), ldbm); /* ldb * tempkn */
                    }
                    RUNTIME_data_flush( sequence, A(k, k) );
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = k+1; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  B(m, k), ldbm,  /* ldb * B->mb   */
                                        A(k, n), ldak,  /* lda * tempnn */
                                lalpha, B(m, n), ldbm); /* ldb * tempnn */
                        }
                        RUNTIME_data_flush( sequence, B(m, k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, A(k, n) );
                    }
                }
            }
            /*
             *  MorseRight / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n-(B->nt-1)*B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt-1-k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            alpha, A(B->nt-1-k, B->nt-1-k), ldak,  /* lda * tempkn */
                                   B(       m, B->nt-1-k), ldbm); /* ldb * tempkn */
                        RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-k) );

                        for (n = k+1; n < B->nt; n++) {
                            ldan = BLKLDD(A, B->nt-1-n);
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, B->nb, tempkn, A->mb,
                                minvalpha, B(m,        B->nt-1-k), ldbm,  /* ldb  * tempkn */
                                           A(B->nt-1-n, B->nt-1-k), ldan, /* A->mb * tempkn (Never last row) */
                                zone,      B(m,        B->nt-1-n), ldbm); /* ldb  * B->nb   */
                        }
                        RUNTIME_data_flush( sequence, B(m,        B->nt-1-k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, A(B->nt-1-n, B->nt-1-k) );
                    }
                }
            }
        }
        /*
         *  MorseRight / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n-(B->nt-1)*B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            lalpha, A(B->nt-1-k, B->nt-1-k), ldak,  /* lda * tempkn */
                                    B(       m, B->nt-1-k), ldbm); /* ldb * tempkn */
                        RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-k) );

                        for (n = k+1; n < B->nt; n++) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, B->nb, tempkn, A->mb,
                                mzone,  B(m,        B->nt-1-k), ldbm,  /* ldb * tempkn */
                                        A(B->nt-1-k, B->nt-1-n), ldak,  /* lda * B->nb   */
                                lalpha, B(m,        B->nt-1-n), ldbm); /* ldb * B->nb   */
                        }
                        RUNTIME_data_flush( sequence, B(m,        B->nt-1-k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-n) );
                    }
                }
            }
            /*
             *  MorseRight / MorseLower / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            alpha, A(k, k), ldak,  /* lda * tempkn */
                                   B(m, k), ldbm); /* ldb * tempkn */
                        RUNTIME_data_flush( sequence, A(k, k) );

                        for (n = k+1; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            ldan = BLKLDD(A, n);
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                minvalpha, B(m, k), ldbm,  /* ldb  * tempkn */
                                           A(n, k), ldan, /* ldan * tempkn */
                                zone,      B(m, n), ldbm); /* ldb  * tempnn */
                        }
                        RUNTIME_data_flush( sequence, B(m, k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        RUNTIME_data_flush( sequence, A(n, k) );
                    }

                }
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
