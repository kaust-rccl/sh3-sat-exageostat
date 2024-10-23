/**
 *
 * @file testing_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
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

static int check_solution(MORSE_enum uplo, MORSE_enum trans, int N, int K,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                          MORSE_Complex64_t *B, int LDB,
                          double beta,  MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC);


int testing_zher2k(int argc, char **argv)
{
    int hres = 0;
    /* Check for number of arguments*/
    if ( argc != 7 ){
        USAGE("HER2K", "alpha beta M N LDA LDB LDC",
              "   - alpha : alpha coefficient\n"
              "   - beta : beta coefficient\n"
              "   - N : number of columns and rows of matrix C and number of row of matrix A and B\n"
              "   - K : number of columns of matrix A and B\n"
              "   - LDA : leading dimension of matrix A\n"
              "   - LDB : leading dimension of matrix B\n"
              "   - LDC : leading dimension of matrix C\n");
        return -1;
    }

    MORSE_Complex64_t alpha = (MORSE_Complex64_t) atol(argv[0]);
    double beta  = (double) atol(argv[1]);
    int N     = atoi(argv[2]);
    int K     = atoi(argv[3]);
    int LDA   = atoi(argv[4]);
    int LDB   = atoi(argv[5]);
    int LDC   = atoi(argv[6]);
    int NKmax = max(N, K);

    double eps;
    int info_solution;
    int u, t;
    size_t LDAxK = LDA*NKmax;
    size_t LDBxK = LDB*NKmax;
    size_t LDCxN = LDC*N;

    MORSE_Complex64_t *A      = (MORSE_Complex64_t *)malloc(LDAxK*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B      = (MORSE_Complex64_t *)malloc(LDBxK*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *C      = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cinit  = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cfinal = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!B) || (!C) || (!Cinit) || (!Cfinal) ){
        free(A); free(B); free(C);
        free(Cinit); free(Cfinal);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZHER2K ROUTINE -------  \n");
    printf("            Size of the Matrix C %d by %d\n", N, K);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
    *  TESTING ZHER2K
    */

    /* Initialize A,B */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxK, A);
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxK, B);

    /* Initialize C */
    MORSE_zplghe( (double)0., MorseUpperLower, N, C, LDC, 51 );

    for (u=0; u<2; u++) {
        for (t=0; t<3; t++) {
            if (trans[t] == MorseTrans) continue;

            memcpy(Cinit,  C, LDCxN*sizeof(MORSE_Complex64_t));
            memcpy(Cfinal, C, LDCxN*sizeof(MORSE_Complex64_t));

            /* MORSE ZHER2K */
            MORSE_zher2k(uplo[u], trans[t], N, K, alpha, A, LDA, B, LDB, beta, Cfinal, LDC);

            /* Check the solution */
            info_solution = check_solution(uplo[u], trans[t], N, K,
                                           alpha, A, LDA, B, LDB, beta, Cinit, Cfinal, LDC);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING ZHER2K (%5s, %s) ........... PASSED !\n", uplostr[u], transstr[t]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING ZHER2K (%5s, %s) ... FAILED !\n", uplostr[u], transstr[t]);    hres++;
                printf("************************************************\n");
            }
        }
    }

    free(A); free(B); free(C);
    free(Cinit); free(Cfinal);

    return hres;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(MORSE_enum uplo, MORSE_enum trans, int N, int K,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                          MORSE_Complex64_t *B, int LDB,
                          double beta, MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC)
{
    int info_solution;
    double Anorm, Bnorm, Cinitnorm, Cmorsenorm, Clapacknorm, Rnorm, result;
    double eps;
    MORSE_Complex64_t beta_const;

    double *work = (double *)malloc(max(N, K)* sizeof(double));

    beta_const  = -1.0;
    Anorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                      (trans == MorseNoTrans) ? N : K,
                                      (trans == MorseNoTrans) ? K : N, A, LDA, work);
    Bnorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                      (trans == MorseNoTrans) ? N : K,
                                      (trans == MorseNoTrans) ? K : N, B, LDB, work);
    Cinitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref,    LDC, work);
    Cmorsenorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cmorse, LDC, work);

    cblas_zher2k(CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 N, K, CBLAS_SADDR(alpha), A, LDA, B, LDB, (beta), Cref, LDC);

    Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref, LDC, work);

    cblas_zaxpy(LDC*N, CBLAS_SADDR(beta_const), Cmorse, 1, Cref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref, LDC, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Cinitnorm %e, Cmorsenorm %e, Clapacknorm %e\n",
           Rnorm, Anorm, Cinitnorm, Cmorsenorm, Clapacknorm);

    result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * N * eps);
    printf("============\n");
    printf("Checking the norm of the difference against reference ZHER2K \n");
    printf("-- ||Cmorse - Clapack||_oo/((||A||_oo+||C||_oo).N.eps) = %e \n", result);

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
         printf("-- The solution is suspicious ! \n");
         info_solution = 1;
    }
    else {
         printf("-- The solution is CORRECT ! \n");
         info_solution= 0 ;
    }

    free(work);

    return info_solution;
}
