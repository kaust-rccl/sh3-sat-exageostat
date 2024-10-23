/**
 *
 * @file testing_zgesv_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesv_incpiv testing
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

static int check_solution(int, int , MORSE_Complex64_t *, int, MORSE_Complex64_t *, MORSE_Complex64_t *, int, double);

int testing_zgesv_incpiv(int argc, char **argv)
{
    int hres = 0;
    /* Check for valid arguments*/
    if (argc != 4){
        USAGE("GESV_INCPIV", "N LDA NRHS LDB",
              "   - N    : the size of the matrix\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the matrix B\n");
        return -1;
    }

    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);
    int NRHS  = atoi(argv[2]);
    int LDB   = atoi(argv[3]);
    double eps;
    int info_solution;
    int i,j;
    int LDAxN = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(LDA*N*(sizeof*A1));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(LDA*N*(sizeof*A2));
    MORSE_Complex64_t *B1 = (MORSE_Complex64_t *)malloc(LDB*NRHS*(sizeof*B1));
    MORSE_Complex64_t *B2 = (MORSE_Complex64_t *)malloc(LDB*NRHS*(sizeof*B2));
    MORSE_desc_t *L;
    int *IPIV;

    /* Check if unable to allocate memory */
    if ( (!A1) || (!A2)|| (!B1) || (!B2) )
    {
        free(A1); free(A2);
        free(B1); free(B2);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work( 'e' );

    /*----------------------------------------------------------
    *  TESTING ZGESV
    */

    /* Initialize A1 and A2 Matrix */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for ( i = 0; i < N; i++)
        for (  j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for ( i = 0; i < N; i++)
        for ( j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* MORSE ZGESV */
    MORSE_Alloc_Workspace_zgesv_incpiv(N, &L, &IPIV, 1, 1);
    MORSE_zgesv_incpiv(N, NRHS, A2, LDA, L, IPIV, B2, LDB);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON INCPIV ZGESV ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the factorization and the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING INCPIV ZGESV ............... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("************************************************\n");
        printf(" - TESTING INCPIV ZGESV ... FAILED !\n");    hres++;
        printf("************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZGETRF + ZGETRS
    */

    /* Initialize A1 and A2  */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for ( i = 0; i < N; i++)
        for (  j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for ( i = 0; i < N; i++)
        for ( j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* Morse routines */
    MORSE_zgetrf_incpiv(N, N, A2, LDA, L, IPIV);
    MORSE_zgetrs_incpiv(MorseNoTrans, N, NRHS, A2, LDA, L, IPIV, B2, LDB);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGETRF + ZGETRS ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING INCPIV ZGETRF + ZGETRS ..... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" - TESTING INCPIV ZGETRF + ZGETRS ... FAILED !\n");
        printf("***************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZGETRF + ZTRSMPL + ZTRSM
    */

    /* Initialize A1 and A2  */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for ( i = 0; i < N; i++)
        for (  j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for ( i = 0; i < N; i++)
        for ( j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* MORSE routines */
    MORSE_zgetrf_incpiv(N, N, A2, LDA, L, IPIV);
    MORSE_ztrsmpl(N, NRHS, A2, LDA, L, IPIV, B2, LDB);
    MORSE_ztrsm(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit,
                 N, NRHS, 1.0, A2, LDA, B2, LDB);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON INCPIV ZGETRF + ZTRSMPL + ZTRSM  ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING INCPIV ZGETRF + ZTRSMPL + ZTRSM ... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("**************************************************\n");
        printf(" - TESTING INCPIV ZGETRF + ZTRSMPL + ZTRSM ... FAILED !\n");
        printf("**************************************************\n");
    }

    free(A1); free(A2); free(B1); free(B2); free(IPIV);
    MORSE_Dealloc_Workspace( &L );

    return hres;
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
