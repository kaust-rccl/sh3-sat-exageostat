/**
 *
 * @file time_zgetri_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  MORSE_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MORSE_zgetri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF(M, N) + FMULS_GETRI( N ))
#define _FADDS (FADDS_GETRF(M, N) + FADDS_GETRI( N ))

//#define GETRI_SYNC

#include "./timing.c"

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
#if 0
static int check_getri_factorization(MORSE_desc_t *descA1, MORSE_desc_t *descA2, int *IPIV)
{
    int info_factorization;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double *work = (double *)malloc((descA1->m)*sizeof(double));
    double eps = LAPACKE_dlamch_work('e');
    MORSE_desc_t        *descB, *descX;
    MORSE_Complex64_t *b = (MORSE_Complex64_t *)malloc((descA1->m)*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *x = (MORSE_Complex64_t *)malloc((descA1->m)*sizeof(MORSE_Complex64_t));

    MORSE_Desc_Create(&descB, b, MorseComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                      descA1->m, 1, 0, 0, descA1->m, 1, 1, 1);
    MORSE_Desc_Create(&descX, x, MorseComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                      descA1->m, 1, 0, 0, descA1->m, 1, 1, 1);

    MORSE_zplrnt_Tile( descX, 537 );
    MORSE_zlacpy_Tile( MorseUpperLower, descX, descB);

    MORSE_zgetrs_Tile( MorseNoTrans, descA2, IPIV, descX );

    Xnorm = MORSE_zlange_Tile(MorseInfNorm, descX,  work);
    Anorm = MORSE_zlange_Tile(MorseInfNorm, descA1, work);
    Bnorm = MORSE_zlange_Tile(MorseInfNorm, descB,  work);

    MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans,
                       (MORSE_Complex64_t)1.,  descA1, descX,
                       (MORSE_Complex64_t)-1., descB);

    Rnorm = MORSE_zlange_Tile(MorseInfNorm, descB, work);

    if (getenv("MORSE_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*(descA1->m)*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The factorization is suspicious ! \n");
        info_factorization = 1;
     }
    else{
        printf("-- The factorization is CORRECT ! \n");
        info_factorization = 0;
    }
    free(x); free(b); free(work);
    MORSE_Desc_Destroy(&descB);
    MORSE_Desc_Destroy(&descX);

    return info_factorization;
}
#endif

/*------------------------------------------------------------------------
 *  Check the accuracy of the computed inverse
 */

static int check_getri_inverse(MORSE_desc_t *descA1, MORSE_desc_t *descA2, int *IPIV, double *dparam )
{
    double Rnorm, Anorm, Ainvnorm, result;
    double *W = (double *)malloc(descA1->n*sizeof(double));
    MORSE_Complex64_t *work = (MORSE_Complex64_t *)malloc(descA1->n*descA1->n*sizeof(MORSE_Complex64_t));
    double eps = LAPACKE_dlamch_work('e');
    MORSE_desc_t        *descW;

    MORSE_Desc_Create(&descW, work, MorseComplexDouble,  descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, descA1->n, 0, 0, descA1->m, descA1->n);

    MORSE_zlaset_Tile( MorseUpperLower, (MORSE_Complex64_t)0., (MORSE_Complex64_t)1., descW);
    MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans,
                       (MORSE_Complex64_t)-1., descA2, descA1,
                       (MORSE_Complex64_t)1.,  descW);

    Anorm    = MORSE_zlange_Tile(MorseInfNorm, descA1, W);
    Ainvnorm = MORSE_zlange_Tile(MorseInfNorm, descA2, W);
    Rnorm    = MORSE_zlange_Tile(MorseInfNorm, descW,  W);

    dparam[IPARAM_ANORM] = Anorm;
    dparam[IPARAM_BNORM] = Ainvnorm;

    result = Rnorm / ( (Anorm*Ainvnorm)*descA1->m*eps ) ;
    dparam[IPARAM_RES] = Rnorm;

    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        dparam[IPARAM_XNORM] = -1.;
    }
    else{
        dparam[IPARAM_XNORM] = 0.;
    }

    MORSE_Desc_Destroy(&descW);
    free(W);
    free(work);

    return MORSE_SUCCESS;
}

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    MORSE_desc_t descW;
    int ret = 0;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,      1, MORSE_Complex64_t, MorseComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA2, check, MORSE_Complex64_t, MorseComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );

    MORSE_Alloc_Workspace_zgetri_Tile_Async(descA, &descW);
    MORSE_zplrnt_Tile( descA, 3453 );

    if ( check ) {
        MORSE_zlacpy_Tile( MorseUpperLower, descA, descA2 );
    }

    /* MORSE ZGETRF / ZTRTRI / ZTRSMRV  */
    {
#if defined(TRACE_BY_SEQUENCE)
        MORSE_sequence_t *sequence;
        MORSE_request_t request[4] = { MORSE_REQUEST_INITIALIZER,
                                       MORSE_REQUEST_INITIALIZER,
                                       MORSE_REQUEST_INITIALIZER,
                                       MORSE_REQUEST_INITIALIZER };

        MORSE_Sequence_Create(&sequence);

        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            MORSE_zgetrf_Tile_Async( descA, piv, sequence, &request[0] );
            MORSE_Sequence_Wait(sequence);

            MORSE_ztrtri_Tile_Async( MorseUpper, MorseNonUnit, descA, sequence, &request[1] );
            MORSE_Sequence_Wait(sequence);

            MORSE_ztrsmrv_Tile_Async( MorseRight, MorseLower, MorseNoTrans, MorseUnit,
                                      (MORSE_Complex64_t) 1.0, descA, &descW,
                                      sequence, &request[2] );
            MORSE_Sequence_Wait(sequence);

            MORSE_zlaswpc_Tile_Async( descA, 1, descA->m, piv, -1,
                                      sequence, &request[3] );
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);
            STOP_TIMING();

        } else {

            START_TIMING();
            MORSE_zgetrf_Tile_Async( descA, piv, sequence, &request[0]);
            MORSE_ztrtri_Tile_Async( MorseUpper, MorseNonUnit,
                                     descA, sequence, &request[1] );
            MORSE_ztrsmrv_Tile_Async( MorseRight, MorseLower, MorseNoTrans, MorseUnit,
                                      (MORSE_Complex64_t) 1.0,
                                      descA, &descW, sequence, &request[2] );
            MORSE_zlaswpc_Tile_Async( descA, 1, descA->m, piv, -1,
                                      sequence, &request[3] );

            /* Wait for everything */
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait( sequence );
            STOP_TIMING();

        }

        MORSE_Sequence_Destroy(sequence[0]);
        MORSE_Sequence_Destroy(sequence[1]);
        MORSE_Sequence_Destroy(sequence[2]);
        MORSE_Sequence_Destroy(sequence[3]);

#else
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            MORSE_zgetrf_Tile(descA, piv);
            MORSE_ztrtri_Tile(MorseUpper, MorseNonUnit, descA);
            MORSE_ztrsmrv_Tile(MorseRight, MorseLower, MorseNoTrans, MorseUnit,
                                (MORSE_Complex64_t) 1.0, descA, &descW);
            MORSE_zlaswpc_Tile(descA, 1, descA->m, piv, -1);
            STOP_TIMING();

        } else {

            MORSE_sequence_t *sequence;
            MORSE_request_t request[2] = { MORSE_REQUEST_INITIALIZER,
                                          MORSE_REQUEST_INITIALIZER };

            MORSE_Sequence_Create(&sequence);

            START_TIMING();
            MORSE_zgetrf_Tile_Async(descA, piv, sequence, &request[0]);
            MORSE_zgetri_Tile_Async(descA, piv, &descW, sequence, &request[1]);
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);
            STOP_TIMING();

            MORSE_Sequence_Destroy(sequence);
        }
#endif
    }

    /* Check the solution */
    if ( check )
    {
        ret = check_getri_inverse(descA2, descA, piv, dparam);

        PASTE_CODE_FREE_MATRIX( descA2 );
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free(descW.mat);
    free( piv );

    return ret;
}
