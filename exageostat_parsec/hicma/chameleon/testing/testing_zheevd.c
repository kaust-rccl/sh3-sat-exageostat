/**
 *
 * @file testing_zheevd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zheevd testing
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
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

static int check_orthogonality(int, int, MORSE_Complex64_t*, int, double);
static int check_reduction(int, int, int, MORSE_Complex64_t*, double*, int, MORSE_Complex64_t*, double);
static int check_solution(int, double*, double*, double);

int testing_zheevd(int argc, char **argv)
{
    /* Check for number of arguments*/
    if (argc != 3) {
        USAGE("HEEVD", "MODE N LDA",
              "   - MODE : mode used to generate the matrix\n"
              "   - N    : size of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n");
        return -1;
    }

    int    mode  = atoi(argv[0]);
    int    N     = atoi(argv[1]);
    int    LDA   = atoi(argv[2]);
    double eps   = LAPACKE_dlamch_work('e');
    double dmax  = 1.0;
    double rcond = 1.0e6;
    int INFO     = -1;

    MORSE_enum uplo    = MorseUpper;
    MORSE_enum vec     = MorseVec;
    int info_ortho     = 0;
    int info_solution  = 0;
    int info_reduction = 0;
    int LDAxN          = LDA * N;

    MORSE_Complex64_t *A1   = NULL;
    MORSE_Complex64_t *A2   = (MORSE_Complex64_t *)malloc(LDAxN * sizeof(MORSE_Complex64_t));
    double            *W1   = (double *)malloc(N * sizeof(double));
    double            *W2   = (double *)malloc(N * sizeof(double));
    MORSE_Complex64_t *work = (MORSE_Complex64_t *)malloc(3* N * sizeof(MORSE_Complex64_t));
    MORSE_desc_t      *T;

    /* Check if unable to allocate memory */
    if ( (!A2) || (!W1) || (!W2) || !(work) )
    {
        free(A2);
        free(W1); free(W2);
        free(work);
        printf("Out of Memory \n ");
        return -2;
    }

    MORSE_Alloc_Workspace_zheevd(N, N, &T, 1, 1);

    /*----------------------------------------------------------
    *  TESTING ZHEEVD
    */
    /* Initialize A1 */
    if (mode == 0){
        int i;
        for (i=0; i<N; i++){
            W1[i] = (double )i+1;
        }
    }
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, N, N,
                         morse_lapack_const(MorseDistSymmetric), ISEED,
                         morse_lapack_const(MorseHermGeev), W1, mode, rcond,
                         dmax, N, N,
                         morse_lapack_const(MorseNoPacking), A2, LDA, work );

    /*
     * Sort the eigenvalue because when computing the tridiag
     * and then the eigenvalue of the DSTQR are sorted.
     * So to avoid testing fail when having good results W1 should be sorted
     */
    LAPACKE_dlasrt_work( 'I', N, W1 );

    if ( vec == MorseVec ) {
        A1 = (MORSE_Complex64_t *)malloc(LDAxN*sizeof(MORSE_Complex64_t));

        /* Copy A2 into A1 */
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, A2, LDA, A1, LDA);
    }

    /*
     * MORSE ZHEEVD
     */
    INFO = MORSE_zheevd(vec, uplo, N, A2, LDA, W2, T);

    if (INFO != 0) {
        printf(" ERROR OCCURED INFO %d\n", INFO);
        goto fin;
    }

    printf("\n");
    printf("------ TESTS FOR MORSE ZHEEVD ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the orthogonality, reduction and the eigen solutions */
    if (vec == MorseVec) {
        info_ortho = check_orthogonality(N, N, A2, LDA, eps);
        info_reduction = check_reduction(uplo, N, 1, A1, W2, LDA, A2, eps);
    }
    info_solution = check_solution(N, W1, W2, eps);

    if ( (info_solution == 0) & (info_ortho == 0) & (info_reduction == 0) ) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZHEEVD ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZHEEVD ... FAILED !\n");
        printf("************************************************\n");
    }

 fin:
    MORSE_Dealloc_Workspace(&T);
    free(A2);
    free(W1);
    free(W2);
    free(work);
    if (A1 != NULL) free(A1);

    return 0;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(int M, int N, MORSE_Complex64_t *Q, int LDQ, double eps)
{
    double  done  =  1.0;
    double  mdone = -1.0;
    double  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    double *work = (double *)malloc(minMN*sizeof(double));

    /* Build the idendity matrix */
    MORSE_Complex64_t *Id = (MORSE_Complex64_t *) malloc(minMN*minMN*sizeof(MORSE_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, done, Q, LDQ, mdone, Id, minMN);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, done, Q, LDQ, mdone, Id, minMN);

    normQ = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), 'U', minMN, Id, minMN, work);

    result = normQ / (minMN * eps);
    printf(" ======================================================\n");
    printf(" ||Id-Q'*Q||_oo / (minMN*eps)          : %15.3E \n",  result );
    printf(" ======================================================\n");

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
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
 *  Check the reduction
 */
static int check_reduction(int uplo, int N, int bw, MORSE_Complex64_t *A, double *D, int LDA, MORSE_Complex64_t *Q, double eps )
{
    (void) bw;
    MORSE_Complex64_t zone  =  1.0;
    MORSE_Complex64_t mzone = -1.0;
    MORSE_Complex64_t *TEMP     = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
    double *work = (double *)malloc(N*sizeof(double));
    double Anorm, Rnorm, result;
    int info_reduction;
    int i;

    /* Compute TEMP =  Q * LAMBDA */
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseUpperLower), N, N, Q, LDA, TEMP, N);

    for (i = 0; i < N; i++){
        cblas_zdscal(N, D[i], &(TEMP[i*N]), 1);
    }
    /* Compute Residual = A - Q * LAMBDA * Q^H */
    /* A is Hermetian but both upper and lower
     * are assumed valable here for checking
     * otherwise it need to be symetrized before
     * checking.
     */
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseUpperLower), N, N, A, LDA, Residual, N);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N, CBLAS_SADDR(mzone), TEMP, N,  Q, LDA, CBLAS_SADDR(zone), Residual,     N);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseOneNorm), N, N, Residual, N,   work);
    Anorm = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseOneNorm), morse_lapack_const(uplo), N, A, LDA, work);

    result = Rnorm / ( Anorm * N * eps);
    if ( uplo == MorseLower ){
        printf(" ======================================================\n");
        printf(" ||A-Q*LAMBDA*Q'||_oo/(||A||_oo.N.eps) : %15.3E \n",  result );
        printf(" ======================================================\n");
    }else{
        printf(" ======================================================\n");
        printf(" ||A-Q'*LAMBDA*Q||_oo/(||A||_oo.N.eps) : %15.3E \n",  result );
        printf(" ======================================================\n");
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }

    free(TEMP); free(Residual);
    free(work);

    return info_reduction;
}
/*------------------------------------------------------------
 *  Check the eigenvalues
 */
static int check_solution(int N, double *E1, double *E2, double eps)
{
    int info_solution, i;
    double resid;
    double maxtmp;
    double maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    double maxeig = max( fabs(E1[0]), fabs(E2[0]) );

    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = max(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = max(maxtmp, maxeig);
        maxel  = max(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" | D - eigcomputed | / (|D| * N * eps) : %15.3E \n",  maxel );
    printf(" ======================================================\n");

    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The eigenvalues are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The eigenvalues are CORRECT ! \n");
        info_solution = 0;
    }
    return info_solution;
}
