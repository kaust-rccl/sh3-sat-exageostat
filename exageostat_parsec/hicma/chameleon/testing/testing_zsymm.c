/**
 *
 * @file testing_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm testing
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

static int check_solution(MORSE_enum transA, MORSE_enum transB, int M, int N,
                          MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                          MORSE_Complex64_t *B, int LDB,
                          MORSE_Complex64_t beta, MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC);

int testing_zsymm(int argc, char **argv)
{
    int hres = 0;
    /* Check for number of arguments*/
    if ( argc != 7 ){
        USAGE("SYMM", "alpha beta M N K LDA LDB LDC",
              "   - alpha : alpha coefficient \n"
              "   - beta : beta coefficient \n"
              "   - M : number of rows of matrices A and C \n"
              "   - N : number of columns of matrices B and C \n"
              "   - LDA : leading dimension of matrix A \n"
              "   - LDB : leading dimension of matrix B \n"
              "   - LDC : leading dimension of matrix C\n");
        return -1;
    }

    MORSE_Complex64_t alpha = (MORSE_Complex64_t) atol(argv[0]);
    MORSE_Complex64_t beta  = (MORSE_Complex64_t) atol(argv[1]);
    int M     = atoi(argv[2]);
    int N     = atoi(argv[3]);
    int LDA   = atoi(argv[4]);
    int LDB   = atoi(argv[5]);
    int LDC   = atoi(argv[6]);
    int MNmax = max(M, N);

    double eps;
    int info_solution;
    int i, j, s, u;
    int LDAxM = LDA*MNmax;
    int LDBxN = LDB*N;
    int LDCxN = LDC*N;

    MORSE_Complex64_t *A      = (MORSE_Complex64_t *)malloc(LDAxM*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B      = (MORSE_Complex64_t *)malloc(LDBxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *C      = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cinit  = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Cfinal = (MORSE_Complex64_t *)malloc(LDCxN*sizeof(MORSE_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!B) || (!C) || (!Cinit) || (!Cfinal) )
    {
        free(A); free(B); free(C);
        free(Cinit); free(Cfinal);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZSYMM ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
    *  TESTING ZSYMM
    */

    /* Initialize A */
    MORSE_zplgsy( (double)0., MorseUpperLower, MNmax, A, LDA, 51 );

    /* Initialize B */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxN, B);

    /* Initialize C */
    LAPACKE_zlarnv_work(IONE, ISEED, LDCxN, C);

    for (s=0; s<2; s++) {
        for (u=0; u<2; u++) {

            /* Initialize  Cinit / Cfinal */
            for ( i = 0; i < M; i++)
                for (  j = 0; j < N; j++)
                    Cinit[LDC*j+i] = C[LDC*j+i];
            for ( i = 0; i < M; i++)
                for (  j = 0; j < N; j++)
                    Cfinal[LDC*j+i] = C[LDC*j+i];

            /* MORSE ZSYMM */
            MORSE_zsymm(side[s], uplo[u], M, N, alpha, A, LDA, B, LDB, beta, Cfinal, LDC);

            /* Check the solution */
            info_solution = check_solution(side[s], uplo[u], M, N, alpha, A, LDA, B, LDB, beta, Cinit, Cfinal, LDC);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING ZSYMM (%5s, %5s) ....... PASSED !\n", sidestr[s], uplostr[u]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING ZSYMM (%s, %s) ... FAILED !\n", sidestr[s], uplostr[u]);    hres++;
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

static int check_solution(MORSE_enum side, MORSE_enum uplo, int M, int N,
                   MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                   MORSE_Complex64_t *B, int LDB,
                   MORSE_Complex64_t beta, MORSE_Complex64_t *Cref, MORSE_Complex64_t *Cmorse, int LDC)
{
    int info_solution, NrowA;
    double Anorm, Bnorm, Cinitnorm, Cmorsenorm, Clapacknorm, Rnorm;
    double eps;
    MORSE_Complex64_t beta_const;
    double result;
    double *work = (double *)malloc(max(M, N)* sizeof(double));

    beta_const  = (MORSE_Complex64_t)-1.0;

    NrowA = (side == MorseLeft) ? M : N;
    Anorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', NrowA, NrowA, A,       LDA, work);
    Bnorm       = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M,     N,     B,       LDB, work);
    Cinitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M,     N,     Cref,    LDC, work);
    Cmorsenorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M,     N,     Cmorse, LDC, work);

    cblas_zsymm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC);

    Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    cblas_zaxpy(LDC * N, CBLAS_SADDR(beta_const), Cmorse, 1, Cref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Bnorm %e, Cinitnorm %e, Cmorsenorm %e, Clapacknorm %e\n",Rnorm,Anorm,Bnorm,Cinitnorm,Cmorsenorm,Clapacknorm);

    result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * N * eps);

    printf("============\n");
    printf("Checking the norm of the difference against reference ZSYMM \n");
    printf("-- ||Cmorse - Clapack||_oo/((||A||_oo+||B||_oo+||C||_oo).N.eps) = %e \n", result );

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
