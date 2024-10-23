/**
 *
 * @file out_of_core.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon out_of_core example
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-29
 *
 */
#include "out_of_core.h"

/*
 * @brief ooc is driver example routine to test the out-of-core feature with StarPU
 * @details TODO: write some details
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
    NB   = iparam[IPARAM_NB];
    NRHS = iparam[IPARAM_NRHS];

    /* compute the algorithm complexity to evaluate performances */
    fadds = (double)( FADDS_POTRF(N) + 2 * FADDS_TRSM(N,NRHS) );
    fmuls = (double)( FMULS_POTRF(N) + 2 * FMULS_TRSM(N,NRHS) );
    flops = 1e-9 * (fmuls + fadds);

    /* initialize the number of thread if not given by the user in argv */
    if ( iparam[IPARAM_THRDNBR] == -1 ) {
        get_thread_count( &(iparam[IPARAM_THRDNBR]) );
    }
    NCPU = iparam[IPARAM_THRDNBR];
    NGPU = 0;

    /* print informations to user */
    print_header( argv[0], iparam);

    /* check that o direct will work */
    if (iparam[IPARAM_OUTOFCORE] > 0) {
        if (! will_o_direct_work(NB)) {
            print_o_direct_wont_work();
            return EXIT_FAILURE;
        }
        char maxMemoryAllowed[32];
        sprintf (maxMemoryAllowed, "%d", iparam[IPARAM_OUTOFCORE]);
        setenv ("STARPU_LIMIT_CPU_MEM", maxMemoryAllowed, 1);
    }

     /* Initialize MORSE with main parameters */
    if ( MORSE_Init( NCPU, NGPU ) != MORSE_SUCCESS ) {
        fprintf(stderr, "Error initializing MORSE library\n");
        return EXIT_FAILURE;
    }
    MORSE_Set(MORSE_TILE_SIZE, NB);

    /* limit ram memory */
    if (iparam[IPARAM_OUTOFCORE] > 0) {
        int new_dd = starpu_disk_register (&starpu_disk_unistd_o_direct_ops,
                                           (void*) "./ooc/", 1024*1024*10);
        if (new_dd == -ENOENT){
            fprintf(stderr, "Can't write on ./ooc/\n");
        	return EXIT_FAILURE;
        }
    }

    MORSE_Desc_Create_User(&descA, NULL, MorseRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1,
                           morse_getaddr_null, // specific function
                           morse_getblkldd_ccrb,
                           morse_getrankof_2d);
    MORSE_Desc_Create(&descB,  NULL, MorseRealDouble,
                      NB, NB,  NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    MORSE_Desc_Create(&descX,  NULL, MorseRealDouble,
                      NB, NB,  NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    MORSE_Desc_Create(&descAC, NULL, MorseRealDouble,
                      NB, NB,  NB*NB, N, N, 0, 0, N, N, 1, 1);

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
    hres = ( (res / N / eps / (anorm * xnorm + bnorm )) > 100.0 );
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

    /* free descriptors descA, descB, descX, descAC */
    MORSE_Desc_Destroy( &descA );
    MORSE_Desc_Destroy( &descB );
    MORSE_Desc_Destroy( &descX );
    MORSE_Desc_Destroy( &descAC );

    /* Finalize MORSE */
    MORSE_Finalize();

    return EXIT_SUCCESS;
}
