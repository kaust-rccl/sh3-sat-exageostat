/**
 *
 * @file time_zgesvd_tile.c
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

#define _NAME  "MORSE_zheev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEBRD( M, N )
#define _FADDS FADDS_GEBRD( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    MORSE_desc_t *descT;
    int jobu  = MorseVec;
    int jobvt = MorseVec;
    int INFO;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, MORSE_Complex64_t, MorseComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( VT, (jobvt == MorseVec), MORSE_Complex64_t, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( U,  (jobu  == MorseVec), MORSE_Complex64_t, M, M );
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, double, N, 1 );

    /* Initialiaze Data */
    MORSE_zplrnt_Tile(descA, 51 );

    /* Allocate Workspace */
    MORSE_Alloc_Workspace_zgesvd(N, N, &descT, 1, 1);

    if ( jobu == MorseVec ) {
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U,  M);
    }
    if ( jobvt == MorseVec ) {
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }

    START_TIMING();
    INFO = MORSE_zgesvd_Tile(jobu, jobvt, descA, S, descT, U, M, VT, N);
    STOP_TIMING();

    if( INFO != 0 ) {
        printf(" ERROR OCCURED INFO %d\n",INFO);
    }

    /* DeAllocate Workspace */
    MORSE_Dealloc_Workspace(&descT);

    if ( U != NULL ) {
        free( U );
    }
    if ( VT != NULL) {
        free( VT );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( S );

    (void)dparam;
    return 0;
}
