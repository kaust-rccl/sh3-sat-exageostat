/**
 *
 * @file step1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step1 example
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-29
 *
 */
#include "step1.h"

/*
 * @brief step1 introduces the simple MORSE interface
 * @details This program solves a linear system AX=B with matrix A symmetric
 * positive definite.
 * The matrix A is first factorized using the Cholesky factorization, A = LL^T.
 * Then the solution X is calculated thanks to forward and back substitutions.
 * We use MORSE routines with the LAPACK interface which means the routines
 * accepts the same matrix format as LAPACK.
 * Note that we copy the matrix to get it in our own tile structures.
 * This means you can get an overhead coming from copies.
 * The resulting code is very close to step0 so that users of CBLAS/LAPACKE
 * should not be lost. The only things new are the MORSE_Init/Morse_Finalize
 * calls, necessary to initialize/finalize the runtime system and mpi, and
 * MORSE_Set to give some specific parameters.
 * This code allows you to expoit parallelism coming from all the cores of your
 * computer and from gpus if you have properly linked with pthread and CUDA
 * ( + CUBLAS optionnaly ).
 * The precision is: double
 */
int main(int argc, char *argv[]) {
    size_t N;    // matrix order
    size_t NRHS; // number of RHS vectors
    int NCPU; // number of cores to use
    int NGPU; // number of gpus (cuda devices) to use
    int UPLO = MorseUpper; // where is stored L

    /* declarations to time the program and evaluate performances */
    double fmuls, fadds, flops, gflops, cpu_time;

    /* variable to check the numerical results */
    double anorm, bnorm, xnorm, eps, res;
    int hres;

    /* initialize some parameters with default values */
    int iparam[IPARAM_SIZEOF];
    memset(iparam, 0, IPARAM_SIZEOF*sizeof(int));
    init_iparam(iparam);

    /* read arguments */
    read_args(argc, argv, iparam);
    N    = iparam[IPARAM_N];
    NRHS = iparam[IPARAM_NRHS];

    /* compute the algorithm complexity to evaluate performances */
    fadds = (double)( FADDS_POTRF(N) + 2 * FADDS_TRSM(N,NRHS) );
    fmuls = (double)( FMULS_POTRF(N) + 2 * FMULS_TRSM(N,NRHS) );
    flops = 1e-9 * (fmuls + fadds);

    /* initialize the number of thread if not given by the user in argv
     * It makes sense only if this program is linked with pthread and
     * multithreaded BLAS and LAPACK */
    if ( iparam[IPARAM_THRDNBR] == -1 ) {
        get_thread_count( &(iparam[IPARAM_THRDNBR]) );
    }
    NCPU = iparam[IPARAM_THRDNBR];
    NGPU = 0;

    /* print informations to user */
    print_header( argv[0], iparam);

    /* Initialize MORSE with main parameters */
    if ( MORSE_Init( NCPU, NGPU ) != MORSE_SUCCESS ) {
        fprintf(stderr, "Error initializing MORSE library\n");
        return EXIT_FAILURE;
    }

    /*
     * Allocate memory for our data using a C macro (see step1.h)
     *     - matrix A                   : size N x N
     *     - set of RHS vectors B       : size N x NRHS
     *     - set of solutions vectors X : size N x NRHS
     */
    double *A    = malloc( N * N    * sizeof(double) );
    double *Acpy = malloc( N * N    * sizeof(double) );
    double *B    = malloc( N * NRHS * sizeof(double) );
    double *X    = malloc( N * NRHS * sizeof(double) );

    /* generate A matrix with random values such that it is spd */
    MORSE_dplgsy( (double)N, MorseUpperLower, N, A, N, 51 );

    /* generate RHS */
    MORSE_dplrnt( N, NRHS, B, N, 5673 );

    /* copy A before facto. in order to check the result */
    memcpy(Acpy, A, N*N*sizeof(double));

    /* copy B in X before solving */
    memcpy(X, B, N*NRHS*sizeof(double));

    /************************************************************/
    /* solve the system AX = B using the Cholesky factorization */
    /************************************************************/

    cpu_time = -CHAMELEON_timer();

    /* Cholesky factorization:
     * A is replaced by its factorization L or L^T depending on uplo */
    MORSE_dpotrf( UPLO, N, A, N );

    /* Solve:
     * B is stored in X on entry, X contains the result on exit.
     * Forward and back substitutions
     */
    MORSE_dpotrs(UPLO, N, NRHS, A, N, X, N);

    cpu_time += CHAMELEON_timer();

    /* print informations to user */
    gflops = flops / cpu_time;
    printf( "%9.3f %9.2f\n", cpu_time, gflops);
    fflush( stdout );

    /************************************************************/
    /* check if solve is correct i.e. AX-B = 0                  */
    /************************************************************/

    /* compute norms to check the result */
    anorm = MORSE_dlange( MorseInfNorm, N, N, Acpy, N);
    bnorm = MORSE_dlange( MorseInfNorm, N, NRHS, B, N);
    xnorm = MORSE_dlange( MorseInfNorm, N, NRHS, X, N);

    /* compute A*X-B, store the result in B */
    MORSE_dgemm(MorseNoTrans, MorseNoTrans,
                N, NRHS, N, 1.0, Acpy, N, X, N, -1.0, B, N);
    res = MORSE_dlange( MorseInfNorm, N, NRHS, B, N);

    /* check residual and print a message */
    eps = LAPACKE_dlamch_work( 'e' );

    /*
     * if hres = 0 then the test succeed
     * else the test failed
     */
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 );
    printf( "   ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN\n");
    if (hres) {
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e FAILURE \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));
    }
    else {
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e SUCCESS \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));
    }

    /* deallocate data */
    free(A);
    free(Acpy);
    free(B);
    free(X);

    /* Finalize MORSE */
    MORSE_Finalize();

    return EXIT_SUCCESS;
}
