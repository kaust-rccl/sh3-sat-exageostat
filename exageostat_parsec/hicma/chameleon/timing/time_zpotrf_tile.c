/**
 *
 * @file time_zpotrf_tile.c
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

#define _NAME  "MORSE_zpotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = MorseUpper;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,  1,     MORSE_Complex64_t, MorseComplexDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  check, MORSE_Complex64_t, MorseComplexDouble, LDB, N, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, check, MORSE_Complex64_t, MorseComplexDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  check, MORSE_Complex64_t, MorseComplexDouble, LDB, N, NRHS );
    MORSE_zplghe_Tile( (double)N, MorseUpperLower, descA, 51 );

    /* Save A for check */
    if (check == 1){
        MORSE_zlacpy_Tile(MorseUpperLower, descA, descAC);
    }

    //RUNTIME_zlocality_allrestrict( STARPU_CUDA );

    /* MORSE ZPOTRF */
    START_TIMING();
    MORSE_zpotrf_Tile(uplo, descA);
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        /* Initialize and save B */
        MORSE_zplrnt_Tile( descB, 7672 );
        MORSE_zlacpy_Tile(MorseUpperLower, descB, descX);

        /* Compute the solution */
        MORSE_zpotrs_Tile( uplo, descA, descX );

        /* Check solution */
        dparam[IPARAM_ANORM] = MORSE_zlange_Tile(MorseInfNorm, descAC);
        dparam[IPARAM_BNORM] = MORSE_zlange_Tile(MorseInfNorm, descB);
        dparam[IPARAM_XNORM] = MORSE_zlange_Tile(MorseInfNorm, descX);
        MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = MORSE_zlange_Tile(MorseInfNorm, descB);

        PASTE_CODE_FREE_MATRIX( descB  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descX  );

    }
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
