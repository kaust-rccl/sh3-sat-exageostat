/**
 *
 * @file out_of_core.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon out_of_core example header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2016-08-23
 *
 */
#ifndef OOC_H
#define OOC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#define int64_t __int64
#endif

/* Define these so that the Microsoft VC compiler stops complaining
   about scanf and friends */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif

#include <starpu.h>
#include "coreblas/lapacke.h"
#include "morse.h"
#include "control/common.h"

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

/* define some tools to time the program */
#include <chameleon/chameleon_timer.h>

/* Integer parameters */
enum iparam_ooc {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_NB,             /* Number of columns in a tile                */
    IPARAM_NRHS,           /* Number of RHS                              */
    IPARAM_OUTOFCORE,      /* if > 0 --> how many memory accepted incore */
                           /* else --> do not use ooc.                   */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_NB            ] = 128;
    iparam[IPARAM_NRHS          ] = 1;
    iparam[IPARAM_OUTOFCORE     ] = 2000;
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
            "  --nb=X           NB size. (default: 128)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "  --ooc=N          Allow to store N MiB in main memory. (default: )\n"
            "\n");
}

/**
 * Read arguments following ooc program call
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
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NB]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else if (startswith( argv[i], "--ooc=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_OUTOFCORE]) );
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
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb gpus:    %d\n"
            "# N:          %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "# ooc:        %d\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            0,
            iparam[IPARAM_N],
            iparam[IPARAM_NB],
            32,
            eps,
            iparam[IPARAM_OUTOFCORE]);

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

// Checking if all block size is a multiple of 4096 Bytes
static int
will_o_direct_work(int nb) {
    if ((nb * nb * sizeof(float)) % 4096 != 0)
        return 0;
    return 1;
}

static void
print_o_direct_wont_work(void) {
    fprintf(stderr, "\n[chameleon] Using out-of-core in o_direct force your blocks' size to be\n"
                    "multiples of 4096. Tip : chose 'n' and 'nb' as both multiples of 32.\n");
}

#endif /* OOC_H */
