/**
 *
 * @file time_zpotrf.c
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

#define _NAME  "MORSE_zpotrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = MorseLower;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, MORSE_Complex64_t, LDA, N );

    /* Initialiaze Data */
    MORSE_zplghe( (double)N, MorseUpperLower, N, A, LDA, 51 );

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( A2, check, MORSE_Complex64_t, A, LDA, N    );

    /* MORSE ZPOSV */
    START_TIMING();
    MORSE_zpotrf(uplo, N, A, LDA);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        PASTE_CODE_ALLOCATE_MATRIX( B, check, MORSE_Complex64_t, LDB, NRHS );
        MORSE_zplrnt( N, NRHS, B, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( X,  check, MORSE_Complex64_t, B, LDB, NRHS );

        MORSE_zpotrs(uplo, N, NRHS, A, LDA, X, LDB);

        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, A2, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

        free(A2); free(B); free(X);
      }

    free(A);

    return 0;
}
