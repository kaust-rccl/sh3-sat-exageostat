#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif


/* Common functions for all steps of the tutorial */

static void get_thread_count(int *thrdnbr) {
#if defined WIN32 || defined WIN64
    sscanf( getenv( "NUMBER_OF_PROCESSORS" ), "%d", thrdnbr );
#else
    *thrdnbr = sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

static int startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}


/* define complexity of algorithms - see Lawn 41 page 120 */
#define FMULS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n) + 0.5) * (double)(__n) + (1. / 3.)))
#define FADDS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n)      ) * (double)(__n) - (1. / 6.)))
#define FMULS_TRSM(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)+1.))
#define FADDS_TRSM(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)-1.))

#include <coreblas/lapacke.h>
#include <morse.h>

/* Integer parameters for step1 */
enum iparam_step1 {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_NRHS,           /* Number of RHS                              */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in step1.c main program */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_NRHS          ] = 1;
 }

/**
 * Print how to use the program
 */
static void show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --n=X            dimension (N). (default: 500)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "\n");
}

/**
 * Read arguments following step1 program call
 */
static void read_args(int argc, char *argv[], int *iparam){
    int i;
    for (i = 1; i < argc && argv[i]; ++i) {
        if ( startswith( argv[i], "--help") || startswith( argv[i], "-help") ||
             startswith( argv[i], "--h") || startswith( argv[i], "-h") ) {
            show_help( argv[0] );
            exit(0);
        } else if (startswith( argv[i], "--n=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_N]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }
}

/**
 * Print a header message to summarize main parameters
 */
static void print_header(char *prog_name, int * iparam) {
#if defined(CHAMELEON_SIMULATION)
    double    eps = 0.;
#else
    double    eps = LAPACKE_dlamch_work( 'e' );
#endif

    printf( "#\n"
            "# MORSE %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb gpus:    %d\n"
            "# N:          %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            0,
            iparam[IPARAM_N],
            128,
            32,
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

/*
 * test external application link with chameleon
 */
int main(int argc, char *argv[]) {

    size_t i, j;
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

    int major, minor, patch;
    MORSE_Version(&major, &minor, &patch);

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
    gflops = 0.0;
    cpu_time = 0.0;

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
    MORSE_Init( NCPU, NGPU );

    /*
     * allocate memory for our data using a C macro (see step1.h)
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

    /* Cholesky facorization:
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
    hres = 0;
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

    /* deallocate data */
    free(A);
    free(Acpy);
    free(B);
    free(X);

    /* Finalize MORSE */
    MORSE_Finalize();

    return EXIT_SUCCESS;
}
