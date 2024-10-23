/**
 *
 * @file pzlanhe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
//ALLOC_WS :  A->mb
//#include <stdlib.h>
//#include <math.h>
//WS_ADD :  A->mb
#include "control/common.h"

#define A(m, n) A, m, n
#define VECNORMS_STEP1(m, n) VECNORMS_STEP1, m, n
#define VECNORMS_STEP2(m, n) VECNORMS_STEP2, m, n
#define RESULT(m, n) RESULT, m, n
/**
 *
 */
/**
 *
 */
void morse_pzlanhe(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_desc_t *VECNORMS_STEP1 = NULL;
    MORSE_desc_t *VECNORMS_STEP2 = NULL;
    MORSE_desc_t *RESULT         = NULL;
    MORSE_context_t *morse;
    MORSE_option_t options;

    int workm, workn;
    int tempkm, tempkn;
    int ldam;
    int m, n;
    /* int part_p, part_q; */

    /* part_p = A->myrank / A->q; */
    /* part_q = A->myrank % A->q; */

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    *result = 0.0;
    switch ( norm ) {
    /*
     *  MorseOneNorm / MorseInfNorm
     */
    case MorseOneNorm:
    case MorseInfNorm:
        /* Init workspace handle for the call to zlanhe */
        RUNTIME_options_ws_alloc( &options, A->mb, 0 );

        workm = A->m;
        workn = chameleon_max( A->nt, A->q );
        MORSE_Desc_Create(&(VECNORMS_STEP1), NULL, MorseRealDouble, A->mb, 1, A->mb,
                          workm, workn, 0, 0, workm, workn, A->p, A->q);

        MORSE_Desc_Create(&(VECNORMS_STEP2), NULL, MorseRealDouble, A->mb, 1, A->mb,
                          workm, 1, 0, 0, workm, 1, A->p, A->q);

        MORSE_Desc_Create(&(RESULT), NULL, MorseRealDouble, 1, 1, 1,
                          1, 1, 0, 0, 1, 1, 1, 1);

        /* Zeroes my intermediate vectors */
        for(m = 0; m < A->mt; m++) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            for(n = 0; n < workn; n++) {
                MORSE_TASK_dlaset(
                    &options,
                    MorseUpperLower, tempkm, 1,
                    0., 0.,
                    VECNORMS_STEP1(m, n), 1);
            }
        }

        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            /* compute sums of absolute values on diagonal tile m */
            MORSE_TASK_dzasum(
                &options,
                MorseRowwise, uplo, tempkm, tempkm,
                A(m, m), ldam, VECNORMS_STEP1(m, m));

            /*
             *  MorseLower
             */
            if (uplo == MorseLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n < m; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute sums of absolute values on rows of tile m */
                    MORSE_TASK_dzasum(
                        &options,
                        MorseRowwise, MorseUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    MORSE_TASK_dzasum(
                        &options,
                        MorseColumnwise, MorseUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(n, m));
                }
            }
            /*
             *  MorseUpper
             */
            else {
//                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
//                    n < A->mt; n+=A->q) {
                for(n = m+1; n < A->mt; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute sums of absolute values on rows of tile m */
                    MORSE_TASK_dzasum(
                        &options,
                        MorseRowwise, MorseUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    MORSE_TASK_dzasum(
                        &options,
                        MorseColumnwise, MorseUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(n, m));
                }
            }
        }

        /* compute vector sum between tiles in rows */
        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            MORSE_TASK_dlaset(
                &options,
                MorseUpperLower, tempkm, 1,
                0., 0.,
                VECNORMS_STEP2(m, 0), 1);
            for(n = 0; n < A->nt; n++) {
                MORSE_TASK_dgeadd(
                    &options,
                    MorseNoTrans, tempkm, 1, A->mb,
                    1.0, VECNORMS_STEP1(m, n), tempkm,
                    1.0, VECNORMS_STEP2(m, 0), tempkm);
            }
            /*
             * Compute max norm of each segment of the final vector in the
             * previous workspace
             */
            MORSE_TASK_dlange(
                &options,
                MorseMaxNorm, tempkm, 1, A->nb,
                VECNORMS_STEP2(m, 0), tempkm,
                VECNORMS_STEP1(m, 0));
        }

        /* Initialize RESULT array */
        MORSE_TASK_dlaset(
            &options,
            MorseUpperLower, 1, 1,
            0., 0.,
            RESULT(0,0), 1);

        /* compute max norm between tiles in the column */
        if (A->myrank % A->q == 0) {
            for(m = 0; m < A->mt; m++) {
                MORSE_TASK_dlange_max(
                    &options,
                    VECNORMS_STEP1(m, 0),
                    RESULT(0,0));
            }
        }

        /* Scatter norm over processus */
        for(m = 0; m < A->p; m++) {
            for(n = 0; n < A->q; n++) {
                MORSE_TASK_dlacpy(
                    &options,
                    MorseUpperLower, 1, 1, 1,
                    RESULT(0,0), 1,
                    VECNORMS_STEP1(m, n), 1 );
            }
        }
        MORSE_Desc_Flush( VECNORMS_STEP2, sequence );
        MORSE_Desc_Flush( VECNORMS_STEP1, sequence );
        MORSE_Desc_Flush( RESULT, sequence );
        RUNTIME_sequence_wait(morse, sequence);
        *result = *(double *)VECNORMS_STEP1->get_blkaddr(VECNORMS_STEP1, A->myrank / A->q, A->myrank % A->q );
        MORSE_Desc_Destroy( &(VECNORMS_STEP1) );
        MORSE_Desc_Destroy( &(VECNORMS_STEP2) );
        MORSE_Desc_Destroy( &(RESULT) );
        break;
    /*
     *  MorseFrobeniusNorm
     */
    case MorseFrobeniusNorm:
        workm = chameleon_max( A->mt, A->p );
        workn = chameleon_max( A->nt, A->q );

        MORSE_Desc_Create(&(VECNORMS_STEP1), NULL, MorseRealDouble, 1, 2, 2,
                          workm, 2*workn, 0, 0, workm, 2*workn, A->p, A->q);
        MORSE_Desc_Create(&(RESULT), NULL, MorseRealDouble, 1, 2, 2,
                          1, 2, 0, 0, 1, 2, 1, 1);

        /* Compute local norm to each tile */
        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            /* Zeroes my intermediate vectors */
            for(n = A->myrank % A->q; n < workn; n+=A->q) {
                MORSE_TASK_dlaset(
                    &options,
                    MorseUpperLower, 1, 2,
                    1., 0.,
                    VECNORMS_STEP1(m,n), 1);
            }

            /* compute norm on diagonal tile m */
            MORSE_TASK_zhessq(
                &options,
                uplo, tempkm,
                A(m, m), ldam,
                VECNORMS_STEP1(m, m));

            /*
             *  MorseLower
             */
            if (uplo == MorseLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n < m; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute norm on the lower part */
                    MORSE_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    MORSE_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
            /*
             *  MorseUpper
             */
            else {
//                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
//                    n < A->mt; n+=A->q) {
                for(n = m+1; n < A->mt; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute norm on the lower part */
                    MORSE_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    MORSE_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
        }

        /* Initialize arrays */
        MORSE_TASK_dlaset(
            &options,
            MorseUpperLower, 1, 2,
            1., 0.,
            RESULT(0,0), 1);

        /* Compute accumulation of scl and ssq */
        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            /*
             *  MorseLower
             */
            if (uplo == MorseLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n <= m; n++) {
                    MORSE_TASK_dplssq(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
            /*
             *  MorseUpper
             */
            else {
//                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
//                    n < A->mt; n+=A->q) {
                for(n = m; n < A->mt; n++) {
                    MORSE_TASK_dplssq(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
        }

        /* Compute scl * sqrt(ssq) */
        MORSE_TASK_dplssq2(
            &options,
            RESULT(0,0));

        /* Copy max norm in tiles to dispatch on every nodes */
        for(m = 0; m < A->p; m++) {
            for(n = 0; n < A->q; n++) {
                MORSE_TASK_dlacpy(
                    &options,
                    MorseUpperLower, 1, 1, 1,
                    RESULT(0,0), 1,
                    VECNORMS_STEP1(m, n), 1 );
            }
        }

        MORSE_Desc_Flush( VECNORMS_STEP1, sequence );
        MORSE_Desc_Flush( RESULT, sequence );
        RUNTIME_sequence_wait(morse, sequence);
        *result = *(double *)VECNORMS_STEP1->get_blkaddr(VECNORMS_STEP1, A->myrank / A->q, A->myrank % A->q );
        MORSE_Desc_Destroy( &(VECNORMS_STEP1) );
        MORSE_Desc_Destroy( &(RESULT) );
        break;

        /*
         *  MorseMaxNorm
         */
        case MorseMaxNorm:
        default:
            /* Init workspace handle for the call to zlange but unused */
            RUNTIME_options_ws_alloc( &options, 1, 0 );

            workm = chameleon_max( A->mt, A->p );
            workn = chameleon_max( A->nt, A->q );

            MORSE_Desc_Create(&(VECNORMS_STEP1), NULL, MorseRealDouble, 1, 1, 1,
                              workm, workn, 0, 0, workm, workn, A->p, A->q);
            MORSE_Desc_Create(&(RESULT), NULL, MorseRealDouble, 1, 1, 1,
                              1, 1, 0, 0, 1, 1, 1, 1);

            /* Compute local maximum to each tile */
            for(m = 0; m < A->mt; m++) {
                tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zlanhe(
                    &options,
                    MorseMaxNorm, uplo, tempkm, A->nb,
                    A(m, m), ldam,
                    VECNORMS_STEP1(m, m));

                /*
                 *  MorseLower
                 */
                if (uplo == MorseLower) {
                    //for(n = A->myrank % A->q; n < m; n+=A->q) {
                    for(n = 0; n < m; n++) {
                        tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                        MORSE_TASK_zlange(
                            &options,
                            MorseMaxNorm, tempkm, tempkn, A->nb,
                            A(m, n), ldam,
                            VECNORMS_STEP1(m, n));
                    }
                }
                /*
                 *  MorseUpper
                 */
                else {
                    //for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                    //  n < A->mt; n+=A->q) {
                    for(n = m+1; n < A->mt; n++) {
                        tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                        MORSE_TASK_zlange(
                            &options,
                            MorseMaxNorm, tempkm, tempkn, A->nb,
                            A(m, n), ldam,
                            VECNORMS_STEP1(m, n));
                    }
                }
            }

            /* Initialize RESULT array */
            MORSE_TASK_dlaset(
                &options,
                MorseUpperLower, 1, 1,
                0., 0.,
                RESULT(0,0), 1);

            /* Compute max norm between tiles */
            for(m = 0; m < A->mt; m++) {
                /*
                 *  MorseLower
                 */
                if (uplo == MorseLower) {
                    //for(n = A->myrank % A->q; n < m; n+=A->q) {
                    for(n = 0; n <= m; n++) {
                        MORSE_TASK_dlange_max(
                            &options,
                            VECNORMS_STEP1(m, n),
                            RESULT(0,0));
                    }
                }
                /*
                 *  MorseUpper
                 */
                else {
                    //for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                    //  n < A->mt; n+=A->q) {
                    for(n = m; n < A->mt; n++) {
                        MORSE_TASK_dlange_max(
                            &options,
                            VECNORMS_STEP1(m, n),
                            RESULT(0,0));
                    }
                }
            }

            /* Copy max norm in tiles to dispatch on every nodes */
            for(m = 0; m < A->p; m++) {
                for(n = 0; n < A->q; n++) {
                    MORSE_TASK_dlacpy(
                        &options,
                        MorseUpperLower, 1, 1, 1,
                        RESULT(0,0), 1,
                        VECNORMS_STEP1(m, n), 1 );
                }
            }

            MORSE_Desc_Flush( VECNORMS_STEP1, sequence );
            MORSE_Desc_Flush( RESULT, sequence );
            RUNTIME_sequence_wait(morse, sequence);
            *result = *(double *)VECNORMS_STEP1->get_blkaddr(VECNORMS_STEP1, A->myrank / A->q, A->myrank % A->q );
            MORSE_Desc_Destroy( &(VECNORMS_STEP1) );
            MORSE_Desc_Destroy( &(RESULT) );
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
