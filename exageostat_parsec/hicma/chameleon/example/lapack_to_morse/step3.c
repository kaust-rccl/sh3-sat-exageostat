/**
 *
 * @file step3.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step3 example
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-29
 *
 */
#include "step3.h"

/*
 * @brief step3 indicates how to give your own tile matrix to MORSE.
 * @details This program is a copy of step2 but instead of using a predefined
 * way for accessing tile data (i.e with MORSE_Desc_Create), we will indicate
 * how to create a MORSE descriptor with an arbitrary tile matrix structure
 * by calling MORSE_Desc_Create_User function.
 * During this step we do not use classical LAPACK matrices (1D array) anymore.
 */
int main(int argc, char *argv[]) {
    size_t N; // matrix order
    int NB;   // number of rows and columns in tiles
    int NRHS; // number of RHS vectors
    int NCPU; // number of cores to use
    int NGPU; // number of gpus (cuda devices) to use
    int UPLO = MorseUpper; // where is stored L

    /* descriptors necessary for calling MORSE tile interface  */
    MORSE_desc_t *descA = NULL, *descAC = NULL, *descB = NULL, *descX = NULL;

    /* Array of pointers to double arrays (representing tiles) */
    double **matA = NULL;

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

    /* Question morse to get the block (tile) size (number of columns) */
    MORSE_Get( MORSE_TILE_SIZE, &NB );

    /* allocate tile data */
    matA = allocate_tile_matrix(N, N, NB);

    /*
     * This function is very similar to MORSE_Desc_Create but the way to
     * access matrix tiles can be controled by the user.
     * To do so, three functions with a precise prototype must be given:
     *     - void* get_blkaddr(const MORSE_desc_t *, int, int)
     *     returns a pointer to the tile m, n
     *     - int   get_blkldd (const MORSE_desc_t *, int)
     *     returns the leading dimension of the tile m, n
     *     - int   get_rankof (const MORSE_desc_t *, int, int)
     *     returns the MPI rank of the tile m, n (0 here because we do not
     *     intend to use this program with MPI)
     */
    MORSE_Desc_Create_User(&descA, matA, MorseRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1,
                           user_getaddr_arrayofpointers,
                           user_getblkldd_arrayofpointers,
                           user_getrankof_zero);

    /*
     * We use the classical MORSE way for accessing tiles for descripotrs
     * B, X and AC. This to show you can define different way to consider tiles
     * in your matrix. The only thing important is to have well defined
     * functions get_blkaddr, get_blkldd, get_rankof corresponding to your data.
     * Note that this call of MORSE_Desc_Create_User routine with
     * morse_getaddr_ccrb, morse_getblkldd_ccrband morse_getrankof_2d functions
     * is equivalent to a call to MORSE_Desc_Create (morse_get... are the
     * functions used inside MORSE_Desc_Create).
     */
    MORSE_Desc_Create(&descB, NULL, MorseRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    MORSE_Desc_Create(&descX, NULL, MorseRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    MORSE_Desc_Create(&descAC, NULL, MorseRealDouble,
                      NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1);

    /* generate A matrix with random values such that it is spd */
    MORSE_dplgsy_Tile( (double)N, MorseUpperLower, descA, 51 );

    /* generate RHS */
    MORSE_dplrnt_Tile( descB, 5673 );

    /* copy A before facto. in order to check the result */
    MORSE_dlacpy_Tile(MorseUpperLower, descA, descAC);

    /* copy B in X before solving
     * same sense as memcpy(X, B, N*NRHS*sizeof(double)) but for descriptors */
    MORSE_dlacpy_Tile(MorseUpperLower, descB, descX);

    /************************************************************/
    /* solve the system AX = B using the Cholesky factorization */
    /************************************************************/

    cpu_time = -CHAMELEON_timer();

    /* Cholesky factorization:
     * A is replaced by its factorization L or L^T depending on uplo */
    MORSE_dpotrf_Tile( UPLO, descA );

    /* Solve:
     * B is stored in X on entry, X contains the result on exit.
     * Forward and back substitutions
     */
    MORSE_dpotrs_Tile( UPLO, descA, descX );

    cpu_time += CHAMELEON_timer();

    /* print informations to user */
    gflops = flops / cpu_time;
    printf( "%9.3f %9.2f\n", cpu_time, gflops);
    fflush( stdout );

    /************************************************************/
    /* check if solve is correct i.e. AX-B = 0                  */
    /************************************************************/

    /* compute norms to check the result */
    anorm = MORSE_dlange_Tile( MorseInfNorm, descAC);
    bnorm = MORSE_dlange_Tile( MorseInfNorm, descB);
    xnorm = MORSE_dlange_Tile( MorseInfNorm, descX);

    /* compute A*X-B, store the result in B */
    MORSE_dgemm_Tile( MorseNoTrans, MorseNoTrans,
                      1.0, descAC, descX, -1.0, descB );
    res = MORSE_dlange_Tile( MorseInfNorm, descB );

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

    /* free the matrix of tiles */
    deallocate_tile_matrix(matA, N, N, NB);
    descA->mat = NULL;

    /* free descriptors descA, descB, descX, descAC */
    MORSE_Desc_Destroy( &descA );
    MORSE_Desc_Destroy( &descB );
    MORSE_Desc_Destroy( &descX );
    MORSE_Desc_Destroy( &descAC );

    /* Finalize MORSE */
    MORSE_Finalize();

    return EXIT_SUCCESS;
}
