/**
 *
 * @file time_zgemm_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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

#define _NAME  "MORSE_zgemm_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(M, N, K)
#define _FADDS FADDS_GEMM(M, N, K)

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    MORSE_Complex64_t alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );


    LDB = chameleon_max(K, iparam[IPARAM_LDB]);
    LDC = chameleon_max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, MORSE_Complex64_t, MorseComplexDouble, LDA, M, K );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, MORSE_Complex64_t, MorseComplexDouble, LDB, K, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, MORSE_Complex64_t, MorseComplexDouble, LDC, M, N );

    /* Initialiaze Data */
    MORSE_zplrnt_Tile( descA, 5373 );
    MORSE_zplrnt_Tile( descB, 7672 );
    MORSE_zplrnt_Tile( descC, 6387 );

#if !defined(CHAMELEON_SIMULATION)
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta);
#else
    alpha = 1.5;
    beta = -2.3;
#endif

    /* Save C for check */
    PASTE_TILE_TO_LAPACK( descC, C2, check, MORSE_Complex64_t, LDC, N );

    START_TIMING();
    MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, alpha, descA, descB, beta, descC );
    STOP_TIMING();

#if !defined(CHAMELEON_SIMULATION)
    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descA, A, check, MORSE_Complex64_t, LDA, K );
        PASTE_TILE_TO_LAPACK( descB, B, check, MORSE_Complex64_t, LDB, N );
        PASTE_TILE_TO_LAPACK( descC, C, check, MORSE_Complex64_t, LDC, N );

        dparam[IPARAM_RES] = z_check_gemm( MorseNoTrans, MorseNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]),
                                           &(dparam[IPARAM_BNORM]),
                                           &(dparam[IPARAM_XNORM]));

        free(A); free(B); free(C); free(C2);
    }
#endif

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descC );
    return 0;
}
