/**
 *
 * @file testing_zposv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zposv testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Bilel Hadri, Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <morse.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
#include "testing_zauxiliary.h"

static int check_factorization(int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, int , double);
static int check_solution(int, int, MORSE_Complex64_t*, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, double);

int testing_zposv(int argc, char **argv)
{
    int hres = 0;

    /* Check for number of arguments*/
    if (argc != 4){
        USAGE("POSV", "N LDA NRHS LDB",
              "   - N    : the size of the matrix\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the RHS B\n");
        return -1;
    }

    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);
    int NRHS  = atoi(argv[2]);
    int LDB   = atoi(argv[3]);
    double eps;
    int uplo;
    int info_solution, info_factorization;
    int trans1, trans2;

    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B1 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B2 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A1) || (!A2)|| (!B1) || (!B2) )
    {
        free(A1); free(A2);
        free(B1); free(B2);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work( 'e' );

    uplo = MorseUpper;
    trans1 = uplo == MorseUpper ? MorseConjTrans : MorseNoTrans;
    trans2 = uplo == MorseUpper ? MorseNoTrans : MorseConjTrans;

    /*-------------------------------------------------------------
    *  TESTING ZPOSV
    */

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    MORSE_zplghe( (double)N, MorseUpperLower, N, A1, LDA, 51 );
    MORSE_zlacpy( MorseUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    MORSE_zplrnt( N, NRHS, B1, LDB, 371 );
    MORSE_zlacpy( MorseUpperLower, N, NRHS, B1, LDB, B2, LDB );

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZPOSV ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* MORSE ZPOSV */
    MORSE_zposv(uplo, N, NRHS, A2, LDA, B2, LDB);

    /* Check the factorization and the solution */
    info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ( (info_solution == 0) && (info_factorization == 0) ) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZPOSV ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("***************************************************\n");
        printf(" - TESTING ZPOSV ... FAILED !\n");    hres++;
        printf("***************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZPOTRF + ZPOTRS
    */

    /* Initialize A1 and A2 for Symmetric Positif Matrix */
    MORSE_zplghe( (double)N, MorseUpperLower, N, A1, LDA, 51 );
    MORSE_zlacpy( MorseUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    MORSE_zplrnt( N, NRHS, B1, LDB, 371 );
    MORSE_zlacpy( MorseUpperLower, N, NRHS, B1, LDB, B2, LDB );

    /* Morse routines */
    MORSE_zpotrf(uplo, N, A2, LDA);
    MORSE_zpotrs(uplo, N, NRHS, A2, LDA, B2, LDB);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZPOTRF + ZPOTRS ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the factorization and the solution */
    info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)&(info_factorization == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING ZPOTRF + ZPOTRS ............ PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("****************************************************\n");
        printf(" - TESTING ZPOTRF + ZPOTRS ... FAILED !\n");
        printf("****************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZPOTRF + ZPTRSM + ZTRSM
    */

    /* Initialize A1 and A2 for Symmetric Positif Matrix */
    MORSE_zplghe( (double)N, MorseUpperLower, N, A1, LDA, 51 );
    MORSE_zlacpy( MorseUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    MORSE_zplrnt( N, NRHS, B1, LDB, 371 );
    MORSE_zlacpy( MorseUpperLower, N, NRHS, B1, LDB, B2, LDB );

    /* MORSE routines */
    MORSE_zpotrf(uplo, N, A2, LDA);
    MORSE_ztrsm(MorseLeft, uplo, trans1, MorseNonUnit,
                 N, NRHS, 1.0, A2, LDA, B2, LDB);
    MORSE_ztrsm(MorseLeft, uplo, trans2, MorseNonUnit,
                 N, NRHS, 1.0, A2, LDA, B2, LDB);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZPOTRF + ZTRSM + ZTRSM  ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the factorization and the solution */
    info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)&(info_factorization == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING ZPOTRF + ZTRSM + ZTRSM ..... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" - TESTING ZPOTRF + ZTRSM + ZTRSM ... FAILED !\n");
        printf("***************************************************\n");
    }

    free(A1); free(A2); free(B1); free(B2);

    return hres;
}


/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
static int check_factorization(int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, int uplo, double eps)
{
    double Anorm, Rnorm;
    MORSE_Complex64_t alpha;
    int info_factorization;
    int i,j;

    MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *L1       = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *L2       = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
    double *work              = (double *)malloc(N*sizeof(double));

    memset((void*)L1, 0, N*N*sizeof(MORSE_Complex64_t));
    memset((void*)L2, 0, N*N*sizeof(MORSE_Complex64_t));

    alpha= 1.0;

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == MorseUpper){
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_ztrmm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, N, N, CBLAS_SADDR(alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, N, Residual, N );
    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, N, A1, LDA );

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

static int check_solution(int N, int NRHS, MORSE_Complex64_t *A1, int LDA, MORSE_Complex64_t *B1, MORSE_Complex64_t *B2, int LDB, double eps )
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    MORSE_Complex64_t alpha, beta;
    double *work = (double *)malloc(N*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    Xnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, NRHS, B2, LDB );
    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, N,    A1, LDA );
    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, NRHS, B1, LDB );

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B2, LDB, CBLAS_SADDR(beta), B1, LDB);
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'i', N, NRHS, B1, LDB );

    if (getenv("MORSE_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
