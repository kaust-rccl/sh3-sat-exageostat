/**
 *
 * @file basic_posv.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon basic_posv example header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-13
 *
 */
#ifndef BASIC_POSV_H
#define BASIC_POSV_H

#if defined( _WIN32 ) || defined( _WIN64 )
#define int64_t __int64
#endif

/* Define these so that the Microsoft VC compiler stops complaining
   about scanf and friends */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

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

#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <morse.h>
#include <coreblas.h>

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

static void
show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "  --gpus=X         Number of GPU workers (default: 0)\n"
            "\n"
            "  --n=X            dimension (N) (default: 500)\n"
            "  --nb=X           Nb size. (default: 128)\n"
            "  --ib=X           IB size. (default: 32)\n"
            "\n");
}

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

enum iparam_examples {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_THRDNBR_SUBGRP, /* Number of cores in a subgroup (NUMA node)  */
    IPARAM_SCHEDULER,      /* What scheduler do we choose (dyn, stat)    */
    IPARAM_M,              /* Number of rows of the matrix               */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_K,              /* RHS or K                                   */
    IPARAM_LDA,            /* Leading dimension of A                     */
    IPARAM_LDB,            /* Leading dimension of B                     */
    IPARAM_LDC,            /* Leading dimension of C                     */
    IPARAM_IB,             /* Inner-blocking size                        */
    IPARAM_NB,             /* Number of columns in a tile                */
    IPARAM_MB,             /* Number of rows in a tile                   */
    IPARAM_UPLO,           /* Consider upper or lower part of sym. mat.  */
    IPARAM_NITER,          /* Number of iteration of each test           */
    IPARAM_WARMUP,         /* Run one test to load dynamic libraries     */
    IPARAM_CHECK,          /* Checking activated or not                  */
    IPARAM_VERBOSE,        /* How much noise do we want?                 */
    IPARAM_AUTOTUNING,     /* Disable/enable autotuning                  */
    IPARAM_INPUTFMT,       /* Input format (Use only for getmi/gecfi)    */
    IPARAM_OUTPUTFMT,      /* Output format (Use only for getmi/gecfi)   */
    IPARAM_TRACE,          /* Generate trace on the first non warmup run */
    IPARAM_DAG,            /* Do we require to output the DOT file?      */
    IPARAM_ASYNC,          /* Asynchronous calls                         */
    IPARAM_MX,             /* */
    IPARAM_NX,             /* */
    IPARAM_RHBLK,          /* Householder reduction parameter for QR/LQ  */
    IPARAM_INPLACE,        /* InPlace/OutOfPlace translation mode        */

    IPARAM_INVERSE,
    IPARAM_NCUDAS,
    IPARAM_NMPI,
    IPARAM_P,              /* Parameter for 2D cyclic distribution       */
    IPARAM_Q,              /* Parameter for 2D cyclic distribution       */
    /* Added for StarPU version */
    IPARAM_PROFILE,
    IPARAM_PRINT_ERRORS,
    IPARAM_PARALLEL_TASKS,
    IPARAM_NO_CPU,
    IPARAM_BOUND,
    /* End */
    IPARAM_SIZEOF
};

enum dparam_examples {
  IPARAM_TIME,
  IPARAM_ANORM,
  IPARAM_BNORM,
  IPARAM_XNORM,
  IPARAM_RNORM,
  IPARAM_AinvNORM,
  IPARAM_RES,
  /* Begin section for hydra integration tool */
  IPARAM_THRESHOLD_CHECK, /* Maximum value accepted for: |Ax-b||/N/eps/(||A||||x||+||b||) */
  /* End section for hydra integration tool  */
  IPARAM_DNBPARAM
};

static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_THRDNBR_SUBGRP] = 1;
    iparam[IPARAM_SCHEDULER     ] = 0;
    iparam[IPARAM_M             ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_K             ] = 1;
    iparam[IPARAM_LDA           ] = -1;
    iparam[IPARAM_LDB           ] = -1;
    iparam[IPARAM_LDC           ] = -1;
    iparam[IPARAM_MB            ] = 128;
    iparam[IPARAM_UPLO          ] = MorseUpper;
    iparam[IPARAM_NB            ] = 128;
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NITER         ] = 1;
    iparam[IPARAM_WARMUP        ] = 1;
    iparam[IPARAM_CHECK         ] = 0;
    iparam[IPARAM_VERBOSE       ] = 0;
    iparam[IPARAM_AUTOTUNING    ] = 0;
    iparam[IPARAM_INPUTFMT      ] = 0;
    iparam[IPARAM_OUTPUTFMT     ] = 0;
    iparam[IPARAM_TRACE         ] = 0;
    iparam[IPARAM_DAG           ] = 0;
    iparam[IPARAM_ASYNC         ] = 1;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_RHBLK         ] = 0;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_RHBLK         ] = 0;
    iparam[IPARAM_INPLACE       ] = MORSE_OUTOFPLACE;

    iparam[IPARAM_INVERSE       ] = 0;
    iparam[IPARAM_NCUDAS        ] = 0;
    iparam[IPARAM_NMPI          ] = 1;
    iparam[IPARAM_P             ] = 1;
    iparam[IPARAM_Q             ] = 1;
    iparam[IPARAM_PROFILE       ] = 0;
    iparam[IPARAM_PRINT_ERRORS  ] = 0;
    iparam[IPARAM_PARALLEL_TASKS] = 0;
    iparam[IPARAM_NO_CPU        ] = 0;
    iparam[IPARAM_BOUND         ] = 0;
}

static void read_args(int argc, char *argv[], int *iparam){
    int i;
    for (i = 1; i < argc && argv[i]; ++i) {
        if ( startswith( argv[i], "--help") || startswith( argv[i], "-help") ||
             startswith( argv[i], "--h") || startswith( argv[i], "-h") ) {
            show_help( argv[0] );
            exit(0);
        } else if (startswith( argv[i], "--n=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_N]) );
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NB]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else if (startswith( argv[i], "--gpus=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NCUDAS]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }
}

static void print_header(char *prog_name, int * iparam) {
    const char *bound_header   = iparam[IPARAM_BOUND]   ? "   thGflop/s" : "";
    const char *check_header   = iparam[IPARAM_CHECK]   ? "     ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN" : "";
    const char *inverse_header = iparam[IPARAM_INVERSE] ? " ||I-A*Ainv||       ||A||    ||Ainv||       ||Id - A*Ainv||/((||A|| ||Ainv||).N.eps)" : "";
#if defined(CHAMELEON_SIMULATION)
    double    eps = 0.;
#else
    double    eps = LAPACKE_dlamch_work( 'e' );
#endif

    printf( "#\n"
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb GPUs:    %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            iparam[IPARAM_NCUDAS],
            iparam[IPARAM_NB],
            iparam[IPARAM_IB],
            eps );

    printf( "#     M       N  K/NRHS   seconds   Gflop/s Deviation%s%s\n",
            bound_header, iparam[IPARAM_INVERSE] ? inverse_header : check_header);
    printf( "# %5.0d   %5.0d   %5.0d\n", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_K]);
    return;
}

#endif /* BASIC_POSV_H */
