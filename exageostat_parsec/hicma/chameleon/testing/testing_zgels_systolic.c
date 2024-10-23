/**
 *
 * @file testing_zgels_systolic.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels_systolic testing
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Boucherie Raphael
 * @date 2017-05-17
 * @precisions normal z -> c d s
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <morse.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
#include "testing_zauxiliary.h"

static int check_orthogonality(int, int, int, MORSE_Complex64_t*, double);
static int check_factorization(int, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, MORSE_Complex64_t*, double);
static int check_solution(int, int, int, MORSE_Complex64_t*, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, double);

int testing_zgels_systolic(int argc, char **argv)
{
    int hres = 0;

    if ( argc < 1 ){
        goto usage;
    }

    /* Check for number of arguments*/
    if ( argc != 7 ) {
      usage:
        USAGE("GELS_SYSTOLIC", "M N LDA NRHS LDB",
              "   - M      : number of rows of the matrix A\n"
              "   - N      : number of columns of the matrix A\n"
              "   - LDA    : leading dimension of the matrix A\n"
              "   - NRHS   : number of RHS\n"
              "   - LDB    : leading dimension of the matrix B\n"
              "   - P      : size of the highest level reduction tree\n"
              "   - Q      : size of the middle reduction trees\n"
              );
        return -1;
    }

    int M      = atoi(argv[0]);
    int N      = atoi(argv[1]);
    int LDA    = max( atoi(argv[2]), M );
    int NRHS   = atoi(argv[3]);
    int LDB    = max( max( atoi(argv[4]), M ), N );
    int p      = atoi(argv[5]);
    int q      = atoi(argv[6]);
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    int K = min(M, N);
    double eps;
    int info_ortho, info_solution, info_factorization;
    int LDAxN    = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B1 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B2 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q  = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_desc_t *TS;
    MORSE_desc_t *TT = NULL;

    /* Check if unable to allocate memory */
    if ( (!A1) || (!A2) || (!B1) || (!B2) || (!Q) )
    {
        free(A1); free(A2);
        free(B1); free(B2);
        free(Q);
        printf("Out of Memory \n ");
        return -2;
    }

    MORSE_Alloc_Workspace_zgels(M, N, &TS, 1, 1);
    MORSE_Alloc_Workspace_zgels(M, N, &TT, 1, 1);
    memset(TS->mat, 0, (TS->llm*TS->lln)*sizeof(MORSE_Complex64_t));
    memset(TT->mat, 0, (TT->llm*TT->lln)*sizeof(MORSE_Complex64_t));

    eps = LAPACKE_dlamch_work( 'e' );

    /*----------------------------------------------------------
     *  TESTING ZGEQRF_PARAM
     */

    /* Initialize matrix */
    matrix.mt = TS->mt;
    matrix.nt = TS->nt;
    matrix.nodes = 1;
    matrix.p = 1;

    libhqr_init_sys( &qrtree,
                     ( M >= N ) ? LIBHQR_QR : LIBHQR_LQ,
                     &matrix, p, q );

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, B1, LDB, B2, LDB );

    /* MORSE ZGELS */
    MORSE_zgels_param(&qrtree, MorseNoTrans, M, N, NRHS, A2, LDA, TS, TT, B2, LDB);
    //MORSE_zgels(MorseNoTrans, M, N, NRHS, A2, LDA, TS, B2, LDB);

    /* MORSE ZGELS */
    if (M >= N)
        /* Building the economy-size Q */
        MORSE_zungqr_param(&qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
    else
        /* Building the economy-size Q */
        MORSE_zunglq_param(&qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
        //MORSE_zunglq(M, N, K, A2, LDA, TS, Q, LDA);


    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGELS_SYSTOLIC ROUTINE -------  \n");
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
        printf(" ---- TESTING ZGELS_SYSTOLIC ............... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZGELS_SYSTOLIC ... FAILED !\n");    hres++;
        printf("************************************************\n");
    }

    /*-------------------------------------------------------------
     *  TESTING ZGEQRF + ZGEQRS or ZGELQF + ZGELQS
     */
    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, B1, LDB, B2, LDB );

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
        MORSE_zgeqrf_param( &qrtree, M, N, A2, LDA, TS, TT );
        MORSE_zungqr_param( &qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
        MORSE_zgeqrs_param( &qrtree, M, N, NRHS, A2, LDA, TS,TT, B2, LDB);

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
        MORSE_zgelqf_param(&qrtree, M, N, A2, LDA, TS, TT);
        //MORSE_zunglq(M, N, K, A2, LDA, TS, Q, LDA);
        MORSE_zunglq_param(&qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
        MORSE_zgelqs_param(&qrtree, M, N, NRHS, A2, LDA, TS, TT, B2, LDB);
        //MORSE_zgelqs(M, N, NRHS, A2, LDA, TS, B2, LDB);

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
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, B1, LDB, B2, LDB );

    /* Initialize Q */
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', LDA, N, 0., 1., Q, LDA );

    /* MORSE ZGEQRF+ ZUNMQR + ZTRSM */
    if (M >= N) {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGEQRF + ZUNMQR + ZTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        MORSE_zgeqrf_param( &qrtree, M, N, A2, LDA, TS, TT );
        MORSE_zungqr_param( &qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
        MORSE_zunmqr_param( &qrtree, MorseLeft, MorseConjTrans, M, NRHS, N, A2, LDA, TS, TT, B2, LDB);
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

        MORSE_zgelqf_param(&qrtree, M, N, A2, LDA, TS, TT);
        MORSE_ztrsm(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, M, NRHS, 1.0, A2, LDA, B2, LDB);
        MORSE_zunglq_param(&qrtree, M, N, K, A2, LDA, TS, TT, Q, LDA);
        //MORSE_zunglq(M, N, K, A2, LDA, TS, Q, LDA);
        MORSE_zunmlq_param(&qrtree, MorseLeft, MorseConjTrans, N, NRHS, M, A2, LDA, TS, TT, B2, LDB);
        //MORSE_zunmlq(MorseLeft, MorseConjTrans, N, NRHS, M, A2, LDA, TS, B2, LDB);
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

    libhqr_finalize( &qrtree );

    free(A1); free(A2); free(B1); free(B2); free(Q);
    MORSE_Dealloc_Workspace( &TS );
    MORSE_Dealloc_Workspace( &TT );

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
    int minMN = min(M, N);

    double *work = (double *)malloc(minMN*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix USE DLASET?*/
    MORSE_Complex64_t *Id = (MORSE_Complex64_t *) malloc(minMN*minMN*sizeof(MORSE_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN );

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_zlansy( LAPACK_COL_MAJOR, 'I', 'U', minMN, Id, minMN );

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

    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, N, Residual, M );
    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, N, A2, LDA );

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

    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, N,    A, LDA );
    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, NRHS, B, LDB );
    Xnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, NRHS, X, LDB );

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, CBLAS_SADDR(alpha), A, LDA, X, LDB, CBLAS_SADDR(beta), B, LDB);

    if (M >= N) {
        MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(M*NRHS*sizeof(MORSE_Complex64_t));
        memset((void*)Residual, 0, M*NRHS*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Residual, M);
        Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, NRHS, Residual, M );
        free(Residual);
    }
    else {
        MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(N*NRHS*sizeof(MORSE_Complex64_t));
        memset((void*)Residual, 0, N*NRHS*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Residual, N);
        Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, NRHS, Residual, N );
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
