/**
 *
 * @file time_zgeqrf_hqr_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-06-08
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  MORSE_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MORSE_zgeqrf_param"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF(M, N)
#define _FADDS FADDS_GEQRF(M, N)

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    MORSE_desc_t *TS;
    MORSE_desc_t *TT;
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;
    int hlvl, llvl, qr_a, domino;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, MORSE_Complex64_t,  MorseComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA0, ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  ( check && M == N ), MORSE_Complex64_t, MorseComplexDouble, LDB, M, NRHS );

    MORSE_zplrnt_Tile( descA, 5373 );

    /* Save A for check */
    if (check == 1 && M == N){
        MORSE_zlacpy_Tile(MorseUpperLower, descA, descA0);
    }

    /* Allocate Workspace */
    MORSE_Alloc_Workspace_zgels(M, N, &TS, P, Q);
    memset(TS->mat, 0, (TS->llm*TS->lln)*sizeof(MorseComplexDouble));
    MORSE_Alloc_Workspace_zgels(M, N, &TT, P, Q);
    memset(TT->mat, 0, (TT->llm*TT->lln)*sizeof(MorseComplexDouble));

    /* Initialize matrix */
    matrix.mt = TS->mt;
    matrix.nt = TS->nt;
    matrix.nodes = 1;
    matrix.p = 1;

    /* Initialize qrtree  */
    hlvl = iparam[IPARAM_HIGHLVL_TREE];
    llvl = iparam[IPARAM_LOWLVL_TREE];
    qr_a = iparam[IPARAM_RHBLK];
    domino =  iparam[IPARAM_QR_DOMINO];

    libhqr_init_hqr( &qrtree,
                     ( M >= N ) ? LIBHQR_QR : LIBHQR_LQ,
                     &matrix, llvl, hlvl, qr_a, P, domino, 0);

    START_TIMING();
    MORSE_zgeqrf_param_Tile(&qrtree, descA, TS, TT );
    STOP_TIMING();

    /* Check the solution */
    if ( check && M == N)
    {
         /* Initialize and save B */
        MORSE_zplrnt_Tile( descX, 2264 );
        MORSE_zlacpy_Tile(MorseUpperLower, descX, descB);

        /* Compute the solution */
        MORSE_zgeqrs_param_Tile(&qrtree, descA, TS, TT, descX );

        /* Check solution */
        dparam[IPARAM_ANORM] = MORSE_zlange_Tile(MorseInfNorm, descA0);
        dparam[IPARAM_BNORM] = MORSE_zlange_Tile(MorseInfNorm, descB);
        dparam[IPARAM_XNORM] = MORSE_zlange_Tile(MorseInfNorm, descX);
        MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descA0, descX, -1.0, descB );
        dparam[IPARAM_RES] = MORSE_zlange_Tile(MorseInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  )
        PASTE_CODE_FREE_MATRIX( descA0 )
        PASTE_CODE_FREE_MATRIX( descB  )
      }

    /* Free Workspace */
    libhqr_finalize( &qrtree );
    MORSE_Dealloc_Workspace( &TS );
    MORSE_Dealloc_Workspace( &TT );
    free( descA );

    return 0;
}
