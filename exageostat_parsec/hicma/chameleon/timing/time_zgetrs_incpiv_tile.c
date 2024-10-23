/**
 *
 * @file time_zgetrs_incpiv_tile.c
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

#define _NAME  "MORSE_zgetrs_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    MORSE_desc_t *descL;
    int *piv;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    check = 1;

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, MORSE_Complex64_t, MorseComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  check, MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, check, MORSE_Complex64_t, MorseComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  check, MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );

    MORSE_zplrnt_Tile(descA, 3456);

    /* Allocate Workspace */
    MORSE_Alloc_Workspace_zgesv_incpiv_Tile(chameleon_min(M,N), &descL, &piv, P, Q);

    /* Save A for check */
    if (check == 1){
        MORSE_zlacpy_Tile(MorseUpperLower, descA, descAC);
    }

    /* MORSE ZGETRF_NOPIV */
    MORSE_zgetrf_incpiv_Tile( descA, descL, piv );

    /* Check the solution */
    if ( check )
    {
        /* Initialize and save B */
        MORSE_zplrnt_Tile( descX, 7732 );
        MORSE_zlacpy_Tile(MorseUpperLower, descX, descB);

        /* Compute the solution */
        START_TIMING();
        MORSE_zgetrs_incpiv_Tile( descA, descL, piv, descX );
        STOP_TIMING();

        /* Check solution */
        dparam[IPARAM_ANORM] = MORSE_zlange_Tile(MorseInfNorm, descAC);
        dparam[IPARAM_BNORM] = MORSE_zlange_Tile(MorseInfNorm, descB);
        dparam[IPARAM_XNORM] = MORSE_zlange_Tile(MorseInfNorm, descX);
        MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = MORSE_zlange_Tile(MorseInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  );
    }

    /* Deallocate Workspace */
    MORSE_Dealloc_Workspace(&descL);

    PASTE_CODE_FREE_MATRIX( descA );
    free( piv );

    return 0;
}
