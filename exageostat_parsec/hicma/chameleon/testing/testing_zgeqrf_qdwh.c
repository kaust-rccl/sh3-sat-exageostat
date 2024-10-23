/**
 *
 * @file testing_zgeqrf_qdwh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf_qdwh testing
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
#include <coreblas/coreblas_z.h>
#include "testing_zauxiliary.h"

static int check_orthogonality(int, int, const MORSE_Complex64_t*, int, double);
static int check_factorization(int, int, const MORSE_Complex64_t*, int, const MORSE_Complex64_t*, int, MORSE_Complex64_t*, int, double);

int testing_zgeqrf_qdwh(int argc, char **argv)
{
    int hres = 0;

    if ( argc != 4 ) {
        USAGE("GEQRF_QDWH", "optid M NB LDA",
              "   - optid: Take into account the fact that A2 is Id or not\n"
              "   - M    : number of rows of the matrix A1 and A2\n"
              "   - NB   : tile size\n"
              "   - IB   : inner tile size\n");
        return -1;
    }

    int optid = atoi(argv[0]) ? 1: 0;
    int M  = atoi(argv[1]);
    int NB = atoi(argv[2]);
    int IB = atoi(argv[3]);
    int MxM = M * M;
    int LDA = 2*M;
    double eps;
    int info_ortho, info_factorization;

    /**
     * Compute A = QR with
     *
     * A = [ A1 ]  and Q = [ Q1 ]
     *     [ A2 ]        = [ Q2 ]
     *
     * and where A1 is the same size as A2
     *
     */
    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q1 = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q2 = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A  = (MORSE_Complex64_t *)malloc(2*M*M*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q;
    MORSE_desc_t *T1, *T2;

    /* Check if unable to allocate memory */
    if ( (!A) || (!A1) || (!A2) || (!Q1) || (!Q2) ){
        free(A); free(A1); free(A2);
        free(Q1); free(Q2);
        printf("Out of Memory \n ");
        return -2;
    }

    MORSE_Disable(MORSE_AUTOTUNING);
    MORSE_Set(MORSE_TILE_SIZE, NB);
    MORSE_Set(MORSE_INNER_BLOCK_SIZE, IB);

    MORSE_Alloc_Workspace_zgels(M, M, &T1, 1, 1);
    MORSE_Alloc_Workspace_zgels(M, M, &T2, 1, 1);

    eps = LAPACKE_dlamch('e');

    /* Initialize A1, A2, and A */
    LAPACKE_zlarnv_work(IONE, ISEED, MxM, A1);
    LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, M, 0., 1., A2, M );

    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, M, A1, M, A,     LDA );
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, M, A2, M, A + M, LDA );

    /* Factorize A */
    MORSE_zgeqrf( M, M, A1, M, T1 );
    MORSE_ztpqrt( M, M, optid ? M : 0,
                  A1, M,
                  A2, M, T2 );

    /* Generate the Q */
    MORSE_ztpgqrt( M, M, M, (optid) ? M : 0,
                   A1, M, T1, A2, M, T2, Q1, M, Q2, M );

    /* Copy Q in a single matrix */
    Q = (MORSE_Complex64_t *)malloc(2*M*M*sizeof(MORSE_Complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, M, Q1, M, Q,     LDA );
    free(Q1);
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, M, Q2, M, Q + M, LDA );
    free(Q2);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGELS ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, M);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality( 2*M, M, Q, LDA, eps );
    info_factorization = check_factorization( 2*M, M, A, LDA, A1, M, Q, LDA, eps );

    if ((info_factorization == 0) & (info_ortho == 0)) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZGELS ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZGELS ... FAILED !\n");    hres++;
        printf("************************************************\n");
    }

    free(A); free(A1); free(A2); free(Q);
    MORSE_Dealloc_Workspace( &T1 );
    MORSE_Dealloc_Workspace( &T2 );

    return hres;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

static int
check_orthogonality( int M, int N,
                     const MORSE_Complex64_t *Q, int LDQ,
                     double eps )
{
    MORSE_Complex64_t *Id;
    double alpha, beta;
    double normQ;
    int info_ortho;
    int minMN = min(M, N);

    double *work = (double *)malloc(minMN*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix */
    Id = (MORSE_Complex64_t *) malloc(minMN*minMN*sizeof(MORSE_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN );

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, beta, Q, LDQ, alpha, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, beta, Q, LDQ, alpha, Id, M);

    normQ = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', 'U', minMN, Id, minMN, work );

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

static int
check_factorization(int M, int N,
                    const MORSE_Complex64_t *A, int LDA,
                    const MORSE_Complex64_t *R, int LDR,
                          MORSE_Complex64_t *Q, int LDQ,
                    double eps )
{
    double Anorm, Rnorm;
    MORSE_Complex64_t alpha;
    int info_factorization;
    double *work = (double *)malloc(max(M,N)*sizeof(double));

    alpha = 1.0;

    if (M >= N) {
        /* Perform Q = Q * R */
        cblas_ztrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, CBLAS_SADDR(alpha), R, LDR, Q, LDQ);
    }
    else {
        /* Perform Q = L * Q */
        cblas_ztrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, N, CBLAS_SADDR(alpha), R, LDR, Q, LDQ);
    }

    /* Compute the Residual */
    CORE_zgeadd( MorseNoTrans, M, N, -1., A, LDA, 1., Q, LDQ );

    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'I', M, N, Q, LDQ, work );
    Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'I', M, N, A, LDA, work );

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

    free(work);

    return info_factorization;
}
