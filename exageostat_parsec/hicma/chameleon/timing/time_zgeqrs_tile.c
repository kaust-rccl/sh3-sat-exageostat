/**
 *
 * @file time_zgeqrs_tile.c
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

#define _NAME  "MORSE_zgeqrs_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRS( M, N, NRHS )
#define _FADDS FADDS_GEQRS( M, N, NRHS )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    MORSE_desc_t *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    check = 1;
    M = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,  1, MORSE_Complex64_t, MorseComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );

    MORSE_zplrnt_Tile( descA, 5373 );

    /* Save A for check */
    if (check == 1 && M == N){
        MORSE_zlacpy_Tile(MorseUpperLower, descA, descAC);
    }

    /* Allocate Workspace */
    MORSE_Alloc_Workspace_zgels_Tile(M, N, &descT, P, Q);
    memset(descT->mat, 0, (descT->llm*descT->lln)*sizeof(MorseComplexDouble));

    /* MORSE ZGEQRF */
    MORSE_zgeqrf_Tile( descA, descT );

    /* Check the solution */
    if ( check && M == N )
    {
        /* Initialize and save B */
        MORSE_zplrnt_Tile( descX, 2264 );
        MORSE_zlacpy_Tile(MorseUpperLower, descX, descB);

        /* Compute the solution */
        START_TIMING();
        MORSE_zgeqrs_Tile( descA, descT, descX );
        STOP_TIMING();

        /* Check solution */
        dparam[IPARAM_ANORM] = MORSE_zlange_Tile(MorseInfNorm, descAC);
        dparam[IPARAM_BNORM] = MORSE_zlange_Tile(MorseInfNorm, descB);
        dparam[IPARAM_XNORM] = MORSE_zlange_Tile(MorseInfNorm, descX);
        MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = MORSE_zlange_Tile(MorseInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  )
    }

    /* Free data */
    MORSE_Dealloc_Workspace(&descT);
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
