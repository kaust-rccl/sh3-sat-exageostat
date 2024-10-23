/**
 *
 * @file pzpotrimm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrimm parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Ali M Charara
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void morse_pzpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldbm, ldcm;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn, tempkn;

    MORSE_Complex64_t alpha = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t beta  = (MORSE_Complex64_t) 0.0;
    MORSE_Complex64_t zbeta;
    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;


    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);

            MORSE_TASK_zpotrf(
                &options,
                MorseLower, tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, MorseLower, MorseConjTrans, MorseNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    MorseLower, MorseNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k), ldan,
                     1.0, A(n, n), ldan);

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k), ldam,
                               A(n, k), ldan,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(n, k) );
            }

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(morse, A->nt + k);

            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k), ldak,
                           A(m, k), ldam);
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                for (n = 0; n < k; n++) {
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(m, k) );
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, MorseNonUnit,
                tempkn, A->mb,
                A(k, k), ldak, A->nb*k);

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, 2*A->nt + k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for(n = 0; n < k; n++) {
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n), ldak,
                    1.0, A(n, n), ldan);

                for(m = n+1; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m), ldak,
                             A(k, n), ldak,
                        1.0, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseLeft, uplo, MorseConjTrans, MorseNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k), ldak,
                         A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k), ldak);

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZSYMM Right / Lower
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(morse, 3*A->nt + k);

            tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
            ldak = BLKLDD(A, k);
            zbeta = k == 0 ? beta : zone;

            for (m = 0; m < C->mt; m++) {
                tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                ldbm = BLKLDD(B, m);
                ldcm = BLKLDD(C, m);

                for (n = 0; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    ldan = BLKLDD(A, n);

                    if (k < n) {
                       MORSE_TASK_zgemm(
                           &options,
                           MorseNoTrans, MorseTrans,
                           tempmm, tempnn, tempkn, A->mb,
                           alpha, B(m, k), ldbm,  /* ldbm * K */
                                  A(n, k), ldan,  /* ldan * K */
                           zbeta, C(m, n), ldcm); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                           MORSE_TASK_zsymm(
                               &options,
                               MorseRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k), ldak,  /* ldak * Y */
                                      B(m, k), ldbm,  /* ldbm * Y */
                               zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                        else {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldbm * K */
                                       A(k, n), ldak,  /* ldak * Y */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(m, k) );
            }
            for (n = 0; n <= k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
            }

            RUNTIME_iteration_pop(morse);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(morse, k);

            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            MORSE_TASK_zpotrf(
                &options,
                MorseUpper,
                tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, MorseUpper, MorseConjTrans, MorseNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zherk(
                    &options,
                    MorseUpper, MorseConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m), ldak,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, m) );
            }

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, A->nt + k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k), ldak,
                           A(k, n), ldak);
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                for (m = 0; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, n) );
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, MorseNonUnit,
                tempkm, A->mb,
                A(k, k), ldak, A->mb*k);

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(morse, 2*A->nt + k);

            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);

            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k), ldam,
                    1.0, A(m, m), ldam);

                for (n = m+1; n < k; n++){
                    ldan = BLKLDD(A, n);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k), ldam,
                             A(n, k), ldan,
                        1.0, A(m, n), ldam);
                }
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseRight, uplo, MorseConjTrans, MorseNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k), ldak);

            RUNTIME_iteration_pop(morse);
        }
        /*
         *  ZSYMM Right / Upper
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(morse, 3*A->nt + k);

            tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
            ldak = BLKLDD(A, k);
            zbeta = k == 0 ? beta : zone;

            for (m = 0; m < C->mt; m++) {
                tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                ldbm = BLKLDD(B, m);
                ldcm = BLKLDD(C, m);

                for (n = 0; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    ldan = BLKLDD(A, n);

                    if (k < n) {
                        MORSE_TASK_zgemm(
                            &options,
                            MorseNoTrans, MorseNoTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, B(m, k), ldbm,  /* ldbm * K */
                                   A(k, n), ldak,  /* ldak * Y */
                            zbeta, C(m, n), ldcm); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                            MORSE_TASK_zsymm(
                                &options,
                                MorseRight, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k), ldak,  /* ldak * Y */
                                       B(m, k), ldbm,  /* ldbm * Y */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                        else {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldbm * K */
                                       A(n, k), ldan,  /* ldan * K */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(m, k) );
            }
            for (m = 0; m <= k; m++) {
                RUNTIME_data_flush( sequence, A(m, k) );
            }

            RUNTIME_iteration_pop(morse);
        }
    }

    RUNTIME_options_finalize(&options, morse);
}
