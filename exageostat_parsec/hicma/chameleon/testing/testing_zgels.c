/**
 *
 * @file testing_zgels.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Bilel Hadri
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

static int check_orthogonality(int, int, int, MORSE_Complex64_t*, double);
static int check_factorization(int, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, MORSE_Complex64_t*, double);
static int check_solution(int, int, int, MORSE_Complex64_t*, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, double);

int testing_zgels(int argc, char **argv)
{
    int hres = 0;
    int mode = 0;

    if ( argc < 1 ){
        goto usage;
    } else {
        mode = atoi(argv[0]);
    }

    /* Check for number of arguments*/
    if ( ((mode == 0) && (argc != 6)) ||
         ((mode != 0) && (argc != 7)) ){
      usage:
        USAGE("GELS", "MODE M N LDA NRHS LDB [RH]",
              "   - MODE : 0: flat, 1: tree (RH needed)\n"
              "   - M    : number of rows of the matrix A\n"
              "   - N    : number of columns of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the matrix B\n"
              "   - RH   : Size of each subdomains\n");
        return -1;
    }

    int M     = atoi(argv[1]);
    int N     = atoi(argv[2]);
    int LDA   = max( atoi(argv[3]), M );
    int NRHS  = atoi(argv[4]);
    int LDB   = max( max( atoi(argv[5]), M ), N );
    int rh;

    int K = min(M, N);
    double eps;
    int info_ortho, info_solution, info_factorization;
    int i,j;
    int LDAxN = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B1 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B2 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q  = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_desc_t *T;

    /* Check if unable to allocate memory */
    if ( (!A1) || (!A2) || (!B1) || (!B2) || (!Q) )
    {
        free(A1); free(A2);
        free(B1); free(B2);
        free(Q);
        printf("Out of Memory \n ");
        return -2;
    }

    if ( mode ) {
        rh = atoi(argv[6]);

        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_TREE_HOUSEHOLDER);
        MORSE_Set(MORSE_HOUSEHOLDER_SIZE, rh);
    }

    MORSE_Alloc_Workspace_zgels(M, N, &T, 1, 1);
    memset(T->mat, 0, (T->llm*T->lln)*sizeof(MORSE_Complex64_t));
    eps = LAPACKE_dlamch_work('e');

    /*----------------------------------------------------------
    *  TESTING ZGELS
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, B1, LDB, B2, LDB );

    /* MORSE ZGELS */
    MORSE_zgels(MorseNoTrans, M, N, NRHS, A2, LDA, T, B2, LDB);

    /* MORSE ZGELS */
    if (M >= N)
       /* Building the economy-size Q */
       MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
    else
       /* Building the economy-size Q */
       MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGELS ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality(M, N, LDA, Q, eps);
    info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)&(info_factorization == 0)&(info_ortho == 0)) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZGELS ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZGELS ... FAILED !\n");    hres++;
        printf("************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZGEQRF + ZGEQRS or ZGELQF + ZGELQS
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
             B2[LDB*j+i] = B1[LDB*j+i];

    if (M >= N) {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGEQRF + ZGEQRS ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n", eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        /* Morse routines */
        MORSE_zgeqrf(M, N, A2, LDA, T);
        MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zgeqrs(M, N, NRHS, A2, LDA, T, B2, LDB);

        /* Check the orthogonality, factorization and the solution */
        info_ortho = check_orthogonality(M, N, LDA, Q, eps);
        info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
        info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

        if ((info_solution == 0)&(info_factorization == 0)&(info_ortho == 0)) {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGEQRF + ZGEQRS ............ PASSED !\n");
            printf("***************************************************\n");
        }
        else{
            printf("***************************************************\n");
            printf(" - TESTING ZGEQRF + ZGEQRS ... FAILED !\n");
            printf("***************************************************\n");
        }
    }
    else  {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGELQF + ZGELQS ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n", eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        /* Morse routines */
        MORSE_zgelqf(M, N, A2, LDA, T);
        MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zgelqs(M, N, NRHS, A2, LDA, T, B2, LDB);

       /* Check the orthogonality, factorization and the solution */
       info_ortho = check_orthogonality(M, N, LDA, Q, eps);
       info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
       info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

       if ( (info_solution == 0) & (info_factorization == 0) & (info_ortho == 0) ) {
          printf("***************************************************\n");
          printf(" ---- TESTING ZGELQF + ZGELQS ............ PASSED !\n");
          printf("***************************************************\n");
       }
       else {
          printf("***************************************************\n");
          printf(" - TESTING ZGELQF + ZGELQS ... FAILED !\n");
          printf("***************************************************\n");
        }
    }

    /*----------------------------------------------------------
    *  TESTING ZGEQRF + ZORMQR + ZTRSM
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* MORSE ZGEQRF+ ZUNMQR + ZTRSM */
    memset((void*)Q, 0, LDA*N*sizeof(MORSE_Complex64_t));
    for (i = 0; i < K; i++)
        Q[LDA*i+i] = 1.0;

    if (M >= N) {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGEQRF + ZUNMQR + ZTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        MORSE_zgeqrf(M, N, A2, LDA, T);
        MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zunmqr(MorseLeft, MorseConjTrans, M, NRHS, N, A2, LDA, T, B2, LDB);
        MORSE_ztrsm(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, N, NRHS, 1.0, A2, LDA, B2, LDB);
    }
    else {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGELQF + ZUNMLQ + ZTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        MORSE_zgelqf(M, N, A2, LDA, T);
        MORSE_ztrsm(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, M, NRHS, 1.0, A2, LDA, B2, LDB);
        MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zunmlq(MorseLeft, MorseConjTrans, N, NRHS, M, A2, LDA, T, B2, LDB);
    }

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality(M, N, LDA, Q, eps);
    info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ( (info_solution == 0) & (info_factorization == 0) & (info_ortho == 0) ) {
        if (M >= N) {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGEQRF + ZUNMQR + ZTRSM .... PASSED !\n");
            printf("***************************************************\n");
        }
        else {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGELQF + ZTRSM + ZUNMLQ .... PASSED !\n");
            printf("***************************************************\n");
        }
    }
    else {
        if (M >= N) {
            printf("***************************************************\n");
            printf(" - TESTING ZGEQRF + ZUNMQR + ZTRSM ... FAILED !\n");
            printf("***************************************************\n");
        }
        else {
            printf("***************************************************\n");
            printf(" - TESTING ZGELQF + ZTRSM + ZUNMLQ ... FAILED !\n");
            printf("***************************************************\n");
        }
    }

    free(A1); free(A2); free(B1); free(B2); free(Q);
    MORSE_Dealloc_Workspace( &T );

    return hres;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

static int check_orthogonality(int M, int N, int LDQ, MORSE_Complex64_t *Q, double eps)
{
    double alpha, beta;
    double normQ;
    int info_ortho;
    int i;
    int minMN = min(M, N);

    double *work = (double *)malloc(minMN*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix USE DLASET?*/
    MORSE_Complex64_t *Id = (MORSE_Complex64_t *) malloc(minMN*minMN*sizeof(MORSE_Complex64_t));
    memset((void*)Id, 0, minMN*minMN*sizeof(MORSE_Complex64_t));
    for (i = 0; i < minMN; i++) {
        Id[i*minMN+i] = (MORSE_Complex64_t)1.0;
    }

    /* Perform Id - Q'Q */
    if (M >= N) {
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    }
    else {
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);
    }
    normQ = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'i', 'u', minMN, Id, minMN, work );

    printf("============\n");
    printf("Checking the orthogonality of Q \n");
    printf("||Id-Q'*Q||_oo / (N*eps) = %e \n", normQ/(minMN*eps));

    if ( isnan(normQ / (minMN * eps)) || isinf(normQ / (minMN * eps)) || (normQ / (minMN * eps) > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    free(work); free(Id);

    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the factorization QR
 */

static int check_factorization(int M, int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, MORSE_Complex64_t *Q, double eps )
{
    double Anorm, Rnorm;
    MORSE_Complex64_t alpha, beta;
    int info_factorization;
    int i,j;

    MORSE_Complex64_t *Ql       = (MORSE_Complex64_t *)malloc(M*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(M*N*sizeof(MORSE_Complex64_t));
    double *work              = (double *)malloc(max(M,N)*sizeof(double));

    alpha=1.0;
    beta=0.0;

    if (M >= N) {
        /* Extract the R */
        MORSE_Complex64_t *R = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
        memset((void*)R, 0, N*N*sizeof(MORSE_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', M, N, A2, LDA, R, N);

        /* Perform Ql=Q*R */
        memset((void*)Ql, 0, M*N*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, CBLAS_SADDR(alpha), Q, LDA, R, N, CBLAS_SADDR(beta), Ql, M);
        free(R);
    }
    else {
        /* Extract the L */
        MORSE_Complex64_t *L = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
        memset((void*)L, 0, M*M*sizeof(MORSE_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', M, N, A2, LDA, L, M);

    /* Perform Ql=LQ */
        memset((void*)Ql, 0, M*N*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, CBLAS_SADDR(alpha), L, M, Q, LDA, CBLAS_SADDR(beta), Ql, M);
        free(L);
    }

    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i]-Ql[j*M+i];

    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', M, N, Residual, M, work );
    Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', M, N, A2, LDA, work );

    if (M >= N) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }
    else {
        printf("============\n");
        printf("Checking the LQ Factorization \n");
        printf("-- ||A-LQ||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }

    if (isnan(Rnorm / (Anorm * N *eps)) || isinf(Rnorm / (Anorm * N *eps)) || (Rnorm / (Anorm * N * eps) > 60.0) ) {
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(work); free(Ql); free(Residual);

    return info_factorization;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, MORSE_Complex64_t *X, int LDB, double eps)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    MORSE_Complex64_t alpha, beta;
    double result;
    double *work = (double *)malloc(max(M, N)* sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', M, N,    A, LDA, work );
    Bnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', N, NRHS, B, LDB, work );
    Xnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', M, NRHS, X, LDB, work );

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, CBLAS_SADDR(alpha), A, LDA, X, LDB, CBLAS_SADDR(beta), B, LDB);

    if (M >= N) {
       MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(M*NRHS*sizeof(MORSE_Complex64_t));
       memset((void*)Residual, 0, M*NRHS*sizeof(MORSE_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Residual, M);
       Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', M, NRHS, Residual, M, work );
       free(Residual);
    }
    else {
       MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(N*NRHS*sizeof(MORSE_Complex64_t));
       memset((void*)Residual, 0, N*NRHS*sizeof(MORSE_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Residual, N);
       Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'i', N, NRHS, Residual, N, work );
       free(Residual);
    }

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
