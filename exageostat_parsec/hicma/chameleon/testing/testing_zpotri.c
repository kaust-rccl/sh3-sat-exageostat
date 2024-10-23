/**
 *
 * @file testing_zpotri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotri testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
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

static int check_factorization(int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, int, double);
static int check_inverse(int, MORSE_Complex64_t *, MORSE_Complex64_t *, int, int, double);

int testing_zpotri(int argc, char **argv)
{
    int hres = 0;

    /* Check for number of arguments*/
    if (argc != 2){
        USAGE("POTRI", "N LDA",
              "   - N    : the size of the matrix\n"
              "   - LDA  : leading dimension of the matrix A\n");
        return -1;
    }

    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);
    double eps;
    int uplo;
    int info_inverse, info_factorization;

    MORSE_Complex64_t *A1   = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2   = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *WORK = (MORSE_Complex64_t *)malloc(2*LDA*sizeof(MORSE_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A1) || (!A2) || (!WORK) )
    {
        free(A1); free(A2);
        free(WORK);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work( 'e' );

    uplo = MorseUpper;

    /*-------------------------------------------------------------
     *  TESTING ZPOTRI
     */

    /* Initialize A1 and A2 for Symmetric Positif Matrix */
    MORSE_zplghe( (double)N, MorseUpperLower, N, A1, LDA, 51 );
    MORSE_zlacpy( MorseUpperLower, N, N, A1, LDA, A2, LDA );

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZPOTRI ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* MORSE ZPOTRF */
    MORSE_zpotrf(uplo, N, A2, LDA);

    /* Check the factorization */
    info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);

    /* MORSE ZPOTRI */
    MORSE_zpotri(uplo, N, A2, LDA);

    /* Check the inverse */
    info_inverse = check_inverse(N, A1, A2, LDA, uplo, eps);

    if ( (info_inverse == 0) && (info_factorization == 0) ) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZPOTRI ..................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("***************************************************\n");
        printf(" - TESTING ZPOTRI ... FAILED !\n");    hres++;
        printf("***************************************************\n");
    }

    free(A1); free(A2); free(WORK);

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

    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', N, N, Residual, N, work );
    Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', N, N, A1, LDA, work );

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
 *  Check the accuracy of the computed inverse
 */

static int check_inverse(int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, int uplo, double eps )
{
    int info_inverse;
    int i, j;
    double Rnorm, Anorm, Ainvnorm, result;
    MORSE_Complex64_t alpha, beta, zone;
    MORSE_Complex64_t *work = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));

    alpha = -1.0;
    beta  = 0.0;
    zone = 1.0;

    /* Rebuild the other part of the inverse matrix */
    if(uplo == MorseUpper){
        for(j=0; j<N; j++) {
            for(i=0; i<j; i++) {
                *(A2+j+i*LDA) = *(A2+i+j*LDA);
            }
        }
        cblas_zhemm(CblasColMajor, CblasLeft, CblasUpper, N, N, CBLAS_SADDR(alpha), A2, LDA, A1, LDA, CBLAS_SADDR(beta), work, N);

    }
    else {
        for(j=0; j<N; j++) {
            for(i=j; i<N; i++) {
                *(A2+j+i*LDA) = *(A2+i+j*LDA);
            }
        }
        cblas_zhemm(CblasColMajor, CblasLeft, CblasLower, N, N, CBLAS_SADDR(alpha), A2, LDA, A1, LDA, CBLAS_SADDR(beta), work, N);
    }

    /* Add the identity matrix to work */
    for(i=0; i<N; i++)
        *(work+i+i*N) = *(work+i+i*N) + zone;

    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'o', N, N, work, N );
    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'o', N, N, A1, LDA );
    Ainvnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'o', N, N, A2, LDA );

    if (getenv("MORSE_TESTING_VERBOSE")) {
        printf( "||A||_1=%f\n||Ainv||_1=%f\n||Id - A*Ainv||_1=%e\n", Anorm, Ainvnorm, Rnorm );
    }

    result = Rnorm / ( (Anorm*Ainvnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the inverse \n");
    printf("-- ||Id - A*Ainv||_1/((||A||_1||Ainv||_1).N.eps) = %e \n", result);

    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The inverse is suspicious ! \n");
        info_inverse = 1;
    }
    else{
        printf("-- The inverse is CORRECT ! \n");
        info_inverse = 0;
    }

    free(work);

    return info_inverse;
}
