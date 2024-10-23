/**
 *
 * @file time_zposv_tile.c
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

#define _NAME  "MORSE_zposv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    MORSE_enum uplo = MorseUpper;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1,      MORSE_Complex64_t, MorseComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX, 1,      MORSE_Complex64_t, MorseComplexDouble, LDB, N, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, check, MORSE_Complex64_t, MorseComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  check, MORSE_Complex64_t, MorseComplexDouble, LDB, N, NRHS );

    /* Initialize AT and bT for Symmetric Positif Matrix */
    MORSE_zplghe_Tile((double)N, MorseUpperLower, descA, 51 );
    MORSE_zplrnt_Tile( descX, 7732 );

    /* Save AT and bT for check */
    if (check == 1){
        MORSE_zlacpy_Tile(MorseUpperLower, descA, descAC);
        MORSE_zlacpy_Tile(MorseUpperLower, descX, descB);
    }
    /* MORSE ZPOSV */
    START_TIMING();
    MORSE_zposv_Tile(uplo, descA, descX);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        dparam[IPARAM_ANORM] = MORSE_zlange_Tile(MorseInfNorm, descAC);
        dparam[IPARAM_BNORM] = MORSE_zlange_Tile(MorseInfNorm, descB);
        dparam[IPARAM_XNORM] = MORSE_zlange_Tile(MorseInfNorm, descX);
        MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = MORSE_zlange_Tile(MorseInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  );
      }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descX );

    return 0;
}
