/**
 *
 * @file zposv_users_functions.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zposv_users_functions example
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-13
 * @precisions normal z -> s d c
 *
 */
#include "posv_users_functions.h"

/**
 *  Function that allocate an array of pointers to square tiles (allocated to 0)
 */
MORSE_Complex64_t **allocate_tile_matrix(int m, int n, int nb){
    int i;
    int mt, nt;
    MORSE_Complex64_t **mat;

    /* compute number of tiles in rows and columns */
    mt = (m%nb==0) ? (m/nb) : (m/nb+1);
    nt = (n%nb==0) ? (n/nb) : (n/nb+1);
    mat = malloc( mt*nt*sizeof(MORSE_Complex64_t*) );
    if (!mat){
        printf ("\nIn allocate_tile_matrix, memory Allocation Failure of mat !\n\n");
        exit (EXIT_FAILURE);
    }
    for (i=0; i<mt*nt; i++){
        *(mat+i) = calloc( nb*nb, sizeof(MORSE_Complex64_t) );
        if (!*(mat+i)){
            printf ("\nIn allocate_tile_matrix, memory Allocation Failure of *(mat+i) !\n\n");
            exit (EXIT_FAILURE);
        }
    }
    return mat;
}

int main(int argc, char *argv[]) {

    int N, NB, NCPU, NCUDA;
    MORSE_desc_t *descA = NULL, *descAC = NULL, *descB = NULL, *descX = NULL;
    MORSE_Complex64_t **matA = NULL;

    /* initialize some parameters with default values */
    int iparam[IPARAM_SIZEOF];
    memset(iparam, 0, IPARAM_SIZEOF*sizeof(int));
    init_iparam(iparam);
    iparam[IPARAM_UPLO] = MorseUpper;

    /* read arguments */
    read_args(argc, argv, iparam);
    N  = iparam[IPARAM_N];
    NB = iparam[IPARAM_NB];

    /* allocate my tile data */
    matA = allocate_tile_matrix(N, N, NB);

    /* initialize the number of thread if not given by the user in argv*/
    if ( iparam[IPARAM_THRDNBR] == -1 ) {
      get_thread_count( &(iparam[IPARAM_THRDNBR]) );
      iparam[IPARAM_THRDNBR] -= iparam[IPARAM_NCUDAS];
    }
    NCPU  = iparam[IPARAM_THRDNBR];
    NCUDA = iparam[IPARAM_NCUDAS];

    /* initialize MORSE */
    MORSE_Init( NCPU, NCUDA);
    MORSE_Set(MORSE_TILE_SIZE,        NB );
    MORSE_Set(MORSE_INNER_BLOCK_SIZE, iparam[IPARAM_IB] );

    print_header( argv[0], iparam);

    /*
     * initialize the structure required for MORSE sequential task-based algorithms
     * MORSE_desc_t is a structure wrapping your data allowing MORSE to get
     * pointers to tiles
     */
    MORSE_Desc_Create_User(&descA, matA, MorseComplexDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1,
                           user_getaddr_arrayofpointers,
                           user_getblkldd_arrayofpointers,
                           user_getrankof_2d);

    /* generate A matrix with random values such that it is hermitian */
    MORSE_zplghe_Tile( (double)N, MorseUpperLower, descA, 51 );

    /* generate RHS */
    MORSE_Desc_Create(&descB, NULL, MorseComplexDouble,
                      NB, NB, NB*NB, N, 1, 0, 0, N, 1, 1, 1);
    MORSE_zplrnt_Tile( descB, 7672 );

    /* copy A before facto. in order to check the result */
    MORSE_Desc_Create(&descAC, NULL, MorseComplexDouble,
                      NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1);
    MORSE_zlacpy_Tile(MorseUpperLower, descA, descAC);

    /* copy B before solving in order to check the result */
    MORSE_Desc_Create(&descX, NULL, MorseComplexDouble,
                      NB, NB, NB*NB, N, 1, 0, 0, N, 1, 1, 1);
    MORSE_zlacpy_Tile(MorseUpperLower, descB, descX);

    /* solve the system AX = B using the Cholesky factorization,
     * A is replaced by its factorization L or L^T depending on uplo
     * B is stored in X on entry, X contains the result on exit */
    MORSE_zposv_Tile(iparam[IPARAM_UPLO], descA, descX);
    /* note that this call is equivalent to
     * MORSE_zpotrf_Tile(uplo, descA) followed by
     * MORSE_zpotrs_Tile(uplo, descA, descX )
     * */

    /* check if solve is correct i.e. AX-B = 0 */
    /* compute norms to check the result */
    double anorm = MORSE_zlange_Tile(MorseInfNorm, descAC);
    double bnorm = MORSE_zlange_Tile(MorseInfNorm, descB);
    double xnorm = MORSE_zlange_Tile(MorseInfNorm, descX);

    /* compute A*X-B, store the result in B */
    MORSE_zgemm_Tile( MorseNoTrans, MorseNoTrans, 1.0, descAC, descX, -1.0, descB );
    double res = MORSE_zlange_Tile(MorseInfNorm, descB);

    /* check residual and print a message */
    #if defined(CHAMELEON_SIMULATION)
    double    eps = 0.;
    #else
    double    eps = LAPACKE_dlamch_work( 'e' );
    #endif
    /*
     * if hres = 0 then the test succeed
     * else the test failed
     */
    int hres = 0;
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 );
    printf( "   ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN\n");
    if (hres)
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e FAILURE \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));
    else
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e SUCCESS \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));

    /* destroy MORSE specific structures */
    MORSE_Desc_Destroy( &descX  );
    MORSE_Desc_Destroy( &descB  );
    MORSE_Desc_Destroy( &descAC );
    MORSE_Desc_Destroy( &descA  );

    /* Finalize MORSE */
    MORSE_Finalize();

    return EXIT_SUCCESS;
}

