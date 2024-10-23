/**
 *
 * @file testing_zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
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

static int check_solution(MORSE_enum uplo, MORSE_enum trans, int N, int K,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                          MORSE_Complex64_t beta,  MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC);


int testing_zsyrk(int argc, char **argv)
{
    int hres = 0;
    /* Check for number of arguments*/
    if ( argc != 6){
        USAGE("SYRK", "alpha beta M N LDA LDC",
              "   - alpha : alpha coefficient\n"
              "   - beta : beta coefficient\n"
              "   - N : number of columns and rows of matrix C and number of row of matrix A\n"
              "   - K : number of columns of matrix A\n"
              "   - LDA : leading dimension of matrix A\n"
              "   - LDC : leading dimension of matrix C\n");
        return -1;
    }

    MORSE_Complex64_t alpha = (MORSE_Complex64_t) atol(argv[0]);
    MORSE_Complex64_t beta  = (MORSE_Complex64_t) atol(argv[1]);
    int N     = atoi(argv[2]);
    int K     = atoi(argv[3]);
    int LDA   = atoi(argv[4]);
    int LDC   = atoi(argv[5]);
    int NKmax = max(N, K);

    double eps;
    int info_solution;
    int u, t;
    size_t LDAxK = LDA*NKmax;
    size_t LDCxN = LDC*N;

    MORSE_Complex64_t *A      = (MORSE_Complex64_t *)malloc(LDAxK*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *C      = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cinit  = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cfinal = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!C) || (!Cinit) || (!Cfinal) ){
        free(A); free(C);
        free(Cinit); free(Cfinal);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZSYRK ROUTINE -------  \n");
    printf("            Size of the Matrix A %d by %d\n", N, K);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
    *  TESTING ZSYRK
    */

    /* Initialize A */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxK, A);

    /* Initialize C */
    MORSE_zplgsy( (double)0., MorseUpperLower, N, C, LDC, 51 );

    for (u=0; u<2; u++) {
        for (t=0; t<2; t++) {
            memcpy(Cinit,  C, LDCxN*sizeof(MORSE_Complex64_t));
            memcpy(Cfinal, C, LDCxN*sizeof(MORSE_Complex64_t));

            /* MORSE ZSYRK */
            MORSE_zsyrk(uplo[u], trans[t], N, K, alpha, A, LDA, beta, Cfinal, LDC);

            /* Check the solution */
            info_solution = check_solution(uplo[u], trans[t], N, K,
                                           alpha, A, LDA, beta, Cinit, Cfinal, LDC);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING ZSYRK (%5s, %s) ........... PASSED !\n", uplostr[u], transstr[t]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING ZSYRK (%5s, %s) ... FAILED !\n", uplostr[u], transstr[t]);    hres++;
                printf("************************************************\n");
            }
        }
    }

    free(A); free(C);
    free(Cinit); free(Cfinal);

    return hres;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(MORSE_enum uplo, MORSE_enum trans, int N, int K,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                          MORSE_Complex64_t beta,  MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC)
{
    int info_solution;
    double Anorm, Cinitnorm, Cmorsenorm, Clapacknorm, Rnorm;
    double eps;
    MORSE_Complex64_t beta_const;
    double result;
    double *work = (double *)malloc(max(N, K)* sizeof(double));

    beta_const  = -1.0;
    Anorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                (trans == MorseNoTrans) ? N : K,
                                (trans == MorseNoTrans) ? K : N, A, LDA, work);
    Cinitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref,    LDC, work);
    Cmorsenorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cmorse, LDC, work);

    cblas_zsyrk(CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                N, K, CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), Cref, LDC);

    Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref, LDC, work);

    cblas_zaxpy(LDC*N, CBLAS_SADDR(beta_const), Cmorse, 1, Cref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', N, N, Cref, LDC, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Cinitnorm %e, Cmorsenorm %e, Clapacknorm %e\n",
           Rnorm, Anorm, Cinitnorm, Cmorsenorm, Clapacknorm);

    result = Rnorm / ((Anorm + Cinitnorm) * N * eps);

    printf("============\n");
    printf("Checking the norm of the difference against reference ZSYRK \n");
    printf("-- ||Cmorse - Clapack||_oo/((||A||_oo+||C||_oo).N.eps) = %e \n", result);

    if ( isinf(Clapacknorm) || isinf(Cmorsenorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
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
