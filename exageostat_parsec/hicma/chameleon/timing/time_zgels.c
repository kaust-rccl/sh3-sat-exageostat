/**
 *
 * @file time_zgels.c
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

#define _NAME  "MORSE_zgels"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    MORSE_desc_t *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,    1,     MORSE_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( x,    1,     MORSE_Complex64_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( Acpy, check, MORSE_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( b,    check, MORSE_Complex64_t, LDB, NRHS);

     /* Initialiaze Data */
    MORSE_zplrnt( M, N,    A, LDA,  453 );
    MORSE_zplrnt( M, NRHS, x, LDB, 5673 );

    MORSE_Alloc_Workspace_zgels(M, N, &T, P, Q);
    memset(T->mat, 0, (T->llm*T->lln)*sizeof(MorseComplexDouble));

    /* Save A and b  */
    if (check) {
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N,    A, LDA, Acpy, LDA);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, x, LDB, b,    LDB);
    }

    START_TIMING();
    MORSE_zgels( MorseNoTrans, M, N, NRHS, A, LDA, T, x, LDB );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, Acpy, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(b);
    }

    MORSE_Dealloc_Workspace( &T );
    free( A );
    free( x );

    return 0;
}
