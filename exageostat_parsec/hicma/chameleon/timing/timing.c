/**
 *
 * @file timing.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Dulceneia Becker
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
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
#include <limits.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif

#include <morse.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#include <coreblas.h>
#endif
#include "flops.h"
#include "timing.h"
#include "control/auxiliary.h"

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif /* defined(CHAMELEON_USE_MPI */

#if defined(CHAMELEON_SCHED_STARPU)
#include <starpu.h>
#endif /* defined(CHAMELEON_SCHED_STARPU) */


#if defined(CHAMELEON_HAVE_GETOPT_H)
#include <getopt.h>
#endif /* defined(CHAMELEON_HAVE_GETOPT_H) */

static int RunTest(int *iparam, _PREC *dparam, double *t_);
static inline void* morse_getaddr_null(const MORSE_desc_t *A, int m, int n)
{
    (void)A;(void)m;(void)n;
    return (void*)( NULL );
}

int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

static int
Test(int64_t n, int *iparam) {
    int      i, j, iter;
    int      thrdnbr, niter;
    int64_t  M, N, K, NRHS;
    double  *t;
#if defined(CHAMELEON_SIMULATION)
    _PREC    eps = 0.;
#else
    _PREC    eps = _LAMCH( 'e' );
#endif
    _PREC    dparam[IPARAM_DNBPARAM];
    double   fmuls, fadds, fp_per_mul, fp_per_add;
    double   sumgf, sumgf2, sumt, sd, flops, gflops;
    char    *s;
    char    *env[] = {
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "ACML_NUM_THREADS",
        "ATLAS_NUM_THREADS",
        "BLAS_NUM_THREADS", ""
    };
    int gnuplot = 0;

/*
 * if hres = 0 then the test succeed
 * if hres = n then the test failed n times
 */
    int hres = 0;

    memset( &dparam, 0, IPARAM_DNBPARAM * sizeof(_PREC) );
    dparam[IPARAM_THRESHOLD_CHECK] = 100.0;

    thrdnbr = iparam[IPARAM_THRDNBR];
    niter   = iparam[IPARAM_NITER];

    M    = iparam[IPARAM_M];
    N    = iparam[IPARAM_N];
    K    = iparam[IPARAM_K];
    NRHS = K;
    (void)M;(void)N;(void)K;(void)NRHS;

    if ( (n < 0) || (thrdnbr < 0 ) ) {
        if (gnuplot && (MORSE_My_Mpi_Rank() == 0) ) {
            printf( "set title '%d_NUM_THREADS: ", thrdnbr );
            for (i = 0; env[i][0]; ++i) {
                s = getenv( env[i] );

                if (i) printf( " " ); /* separating space */

                for (j = 0; j < 5 && env[i][j] && env[i][j] != '_'; ++j)
                    printf( "%c", env[i][j] );

                if (s)
                    printf( "=%s", s );
                else
                    printf( "->%s", "?" );
            }
            printf( "'\n" );
            printf( "%s\n%s\n%s\n%s\n%s%s%s\n",
                    "set xlabel 'Matrix size'",
                    "set ylabel 'Gflop/s'",
                    "set key bottom",
                    gnuplot > 1 ? "set terminal png giant\nset output 'timeplot.png'" : "",
                    "plot '-' using 1:5 title '", _NAME, "' with linespoints" );
        }
        return 0;
    }

    if ( MORSE_My_Mpi_Rank() == 0)
        printf( "%7d %7d %7d ", iparam[IPARAM_M], iparam[IPARAM_N], iparam[IPARAM_K] );
    fflush( stdout );

    t = (double*)malloc(niter*sizeof(double));
    memset(t, 0, niter*sizeof(double));

    if (sizeof(_TYPE) == sizeof(_PREC)) {
        fp_per_mul = 1;
        fp_per_add = 1;
    } else {
        fp_per_mul = 6;
        fp_per_add = 2;
    }

    fadds = (double)(_FADDS);
    fmuls = (double)(_FMULS);
    flops = 1e-9 * (fmuls * fp_per_mul + fadds * fp_per_add);

    if ( iparam[IPARAM_WARMUP] ) {
      int status = RunTest( iparam, dparam, &(t[0]));
      if (status != MORSE_SUCCESS) {
          free(t);
          return status;
      }
    }

    sumgf  = 0.0;
    double sumgf_upper  = 0.0;
    sumgf2 = 0.0;
    sumt   = 0.0;

    for (iter = 0; iter < niter; iter++)
    {
        if( iter == 0 ) {
          if ( iparam[IPARAM_TRACE] )
            iparam[IPARAM_TRACE] = 2;
          if ( iparam[IPARAM_DAG] )
            iparam[IPARAM_DAG] = 2;
          if ( iparam[IPARAM_PROFILE] )
            iparam[IPARAM_PROFILE] = 2;

          int status = RunTest( iparam, dparam, &(t[iter]));
          if (status != MORSE_SUCCESS) return status;

          iparam[IPARAM_TRACE] = 0;
          iparam[IPARAM_DAG] = 0;
          iparam[IPARAM_PROFILE] = 0;
        }
        else {
          int status = RunTest( iparam, dparam, &(t[iter]));
          if (status != MORSE_SUCCESS) return status;
        }
        gflops = flops / t[iter];

#if defined (CHAMELEON_SCHED_STARPU)
        /* TODO: create chameleon interface encapsulating this instead */
        if (iparam[IPARAM_BOUND])
        {
            double upper_gflops = 0.0;
            double tmin = 0.0;
            double integer_tmin = 0.0;
            starpu_bound_compute(&tmin, &integer_tmin, 0);
            upper_gflops  = (flops / (tmin / 1000.0));
            sumgf_upper += upper_gflops;
        }
#endif
        sumt   += t[iter];
        sumgf  += gflops;
        sumgf2 += gflops*gflops;
    }

    gflops = sumgf / niter;
    sd = sqrt((sumgf2 - (sumgf*sumgf)/niter)/niter);

    if ( MORSE_My_Mpi_Rank() == 0) {
        printf( "%9.3f %9.2f +-%7.2f  ", sumt/niter, gflops, sd);

        if (iparam[IPARAM_BOUND]) {
            printf(" %9.2f",  sumgf_upper/niter);
        }

        if ( iparam[IPARAM_CHECK] ){
            hres = ( dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ) > dparam[IPARAM_THRESHOLD_CHECK] );

            if (hres) {
                printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e FAILURE",
                    dparam[IPARAM_RES], dparam[IPARAM_ANORM], dparam[IPARAM_XNORM], dparam[IPARAM_BNORM],
                    dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ));
            }
            else {
                printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e SUCCESS",
                    dparam[IPARAM_RES], dparam[IPARAM_ANORM], dparam[IPARAM_XNORM], dparam[IPARAM_BNORM],
                    dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ));
            }
        }

        if ( iparam[IPARAM_INVERSE] ) {
            printf( " %8.5e %8.5e %8.5e     %8.5e",
                    dparam[IPARAM_RNORM], dparam[IPARAM_ANORM], dparam[IPARAM_AinvNORM],
                    dparam[IPARAM_RNORM] /((dparam[IPARAM_ANORM] * dparam[IPARAM_AinvNORM])*n*eps));
        }

        printf("\n");

        fflush( stdout );
    }
    free(t);

    return hres;
}

static inline int
startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

static int
get_range(char *range, int *start_p, int *stop_p, int *step_p) {
    char *s, *s1, buf[21];
    int colon_count, copy_len, nbuf=20, n;
    int start=1000, stop=10000, step=1000;

    colon_count = 0;
    for (s = strchr( range, ':'); s; s = strchr( s+1, ':'))
        colon_count++;

    if (colon_count == 0) { /* No colon in range. */
        if (sscanf( range, "%d", &start ) < 1 || start < 1)
            return -1;
        step = start / 10;
        if (step < 1) step = 1;
        stop = start + 10 * step;

    } else if (colon_count == 1) { /* One colon in range.*/
        /* First, get the second number (after colon): the stop value. */
        s = strchr( range, ':' );
        if (sscanf( s+1, "%d", &stop ) < 1 || stop < 1)
            return -1;

        /* Next, get the first number (before colon): the start value. */
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start > stop || start < 1)
            return -1;

        /* Let's have 10 steps or less. */
        step = (stop - start) / 10;
        if (step < 1)
            step = 1;
    } else if (colon_count == 2) { /* Two colons in range. */
        /* First, get the first number (before the first colon): the start value. */
        s = strchr( range, ':' );
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start < 1)
            return -1;

        /* Next, get the second number (after the first colon): the stop value. */
        s1 = strchr( s+1, ':' );
        n = s1 - (s + 1);
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, s+1, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &stop ) < 1 || stop < start)
            return -1;

        /* Finally, get the third number (after the second colon): the step value. */
        if (sscanf( s1+1, "%d", &step ) < 1 || step < 1)
            return -1;
    } else

        return -1;

    *start_p = start;
    *stop_p = stop;
    *step_p = step;

    return 0;
}

static void
show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  -h  --help               Show this help\n"
            "\n"
            "  Machine parameters:\n"
            "    -t, --threads=x        Number of CPU workers (default: automatic detection through runtime)\n"
            "    -g, --gpus=x           Number of GPU workers (default: 0)\n"
            "    -P, --P=x              Rows (P) in the PxQ process grid (deafult: 1)\n"
            "        --nocpu            All GPU kernels are exclusively executed on GPUs (default: 0)\n"
            "\n"
            "  Matrix parameters:\n"
            "    -m, --m, --M=x         Dimension (M) of the matrices (default: N)\n"
            "    -n, --n, --N=x         Dimension (N) of the matrices\n"
            "    -N, --n_range=R        Range of N values\n"
            "                           with R=Start:Stop:Step (default: 500:5000:500)\n"
            "    -k, --k, --K, --nrhs=x Dimension (K) of the matrices or number of right-hand size (default: 1)\n"
            "    -b, --nb=x             NB size. (default: 320)\n"
            "    -i, --ib=x             IB size. (default: 32)\n"
            //"    -x, --mx=x             ?\n" todo
            //"    -X, --nx=x             ?\n" todo
            "\n"
            "  Check/prints:\n"
            "        --niter=x          Number of iterations performed for each test (default: 1)\n"
            "    -W, --nowarnings       Do not show warnings\n"
            "    -w, --nowarmup         Cancel the warmup run to pre-load libraries\n"
            "    -c, --check            Check result\n"
            "    -C, --inv              Check on inverse\n"
            "        --mode=x           Change xLATMS matrix mode generation for SVD/EVD (default: 4)\n"
            "                           Must be between 0 and 20 included\n"
            "\n"
            "  Profiling:\n"
            "    -T, --trace            Enable trace generation\n"
            "        --progress         Display progress indicator\n"
            "    -d, --dag              Enable DAG generation\n"
            "                           Generates a dot_dag_file.dot.\n"
            "    -p, --profile          Print profiling informations\n"
            "\n"
            "  HQR options:\n"
            "    -a, --qr_a, --rhblk=N  If N > 0, enable Householder mode for QR and LQ factorization\n"
            "                           N is the size of each subdomain (default: -1)\n"
            "    -l, --llvl=x           Tree used for low level reduction inside nodes (default: -1)\n"
            "    -L, --hlvl=x           Tree used for high level reduction between nodes, only if P > 1 (default: -1).\n"
            "                           (-1: Automatic, 0: Flat, 1: Greedy, 2: Fibonacci, 3: Binary, 4: Replicated greedy)\n"
            "    -D, --domino           Enable the domino between upper and lower trees.\n"
            "\n"
            "  Advanced options\n"
            "        --nobigmat         Disable single large matrix allocation for multiple tiled allocations\n"
            "    -s, --sync             Enable synchronous calls in wrapper function such as POTRI\n"
            "    -o, --ooc              Enable out-of-core (available only with StarPU)\n"
            "    -G, --gemm3m           Use gemm3m complex method\n"
            "        --bound            Compare result to area bound\n"
            "\n");
}


static void
print_header(char *prog_name, int * iparam) {
    const char *bound_header   = iparam[IPARAM_BOUND]   ? "   thGflop/s" : "";
    const char *check_header   = iparam[IPARAM_CHECK]   ? "     ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN" : "";
    const char *inverse_header = iparam[IPARAM_INVERSE] ? " ||I-A*Ainv||       ||A||    ||Ainv||       ||Id - A*Ainv||/((||A|| ||Ainv||).N.eps)" : "";
#if defined(CHAMELEON_SIMULATION)
    _PREC    eps = 0.;
#else
    _PREC    eps = _LAMCH( 'e' );
#endif

    printf( "#\n"
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb GPUs:    %d\n"
#if defined(CHAMELEON_USE_MPI)
            "# Nb mpi:     %d\n"
            "# PxQ:        %dx%d\n"
#endif
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
#if defined(CHAMELEON_USE_MPI)
            iparam[IPARAM_NMPI],
            iparam[IPARAM_P], iparam[IPARAM_Q],
#endif
            iparam[IPARAM_NB],
            iparam[IPARAM_IB],
            eps );

    printf( "#     M       N  K/NRHS   seconds   Gflop/s Deviation%s%s\n",
            bound_header, iparam[IPARAM_INVERSE] ? inverse_header : check_header);
    return;
}

#define GETOPT_STRING "ht:g:P:8M:m:N:n:K:k:b:i:x:X:1:WwcCT2dpa:l:L:D9:3soG45"
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
static struct option long_options[] =
{
    {"help",          no_argument,       0,      'h'},
    // Configuration
    {"threads",       required_argument, 0,      't'},
    {"gpus",          required_argument, 0,      'g'},
    {"P",             required_argument, 0,      'P'},
    {"nocpu",         no_argument,       0,      '8'},
    // Matrix parameters
    {"M",             required_argument, 0,      'm'},
    {"m",             required_argument, 0,      'm'},
    {"N",             required_argument, 0,      'n'},
    {"n",             required_argument, 0,      'n'},
    {"n_range",       required_argument, 0,      'N'},
    {"K",             required_argument, 0,      'K'},
    {"k",             required_argument, 0,      'k'},
    {"nrhs",          required_argument, 0,      'k'},
    {"nb",            required_argument, 0,      'b'},
    {"ib",            required_argument, 0,      'i'},
    {"mx",            required_argument, 0,      'x'},
    {"nx",            required_argument, 0,      'X'},
    // Check/prints
    {"niter",         required_argument, 0,      '1'},
    {"nowarnings",    no_argument,       0,      'W'},
    {"nowarmup",      no_argument,       0,      'w'},
    {"check",         no_argument,       0,      'c'},
    {"inv",           no_argument,       0,      'C'},
    // Profiling
    {"trace",         no_argument,       0,      'T'},
    {"progress",      no_argument,       0,      '2'},
    {"dag",           no_argument,       0,      'd'},
    {"profile",       no_argument,       0,      'p'},
    // HQR options
    {"rhblk",         required_argument, 0,      'a'},
    {"qr_a",          required_argument, 0,      'a'},
    {"llvl",          required_argument, 0,      'l'},
    {"hlvl",          required_argument, 0,      'L'},
    {"domino",        no_argument,       0,      'D'},
    // Other
    {"mode",          required_argument, 0,      '9'},
    {"nobigmat",      no_argument,       0,      '3'},
    {"sync",          no_argument,       0,      's'},
    {"ooc",           no_argument,       0,      'o'},
    {"gemm3m",        no_argument,       0,      'G'},
    {"bound",         no_argument,       0,      '5'},
    {0, 0, 0, 0}
};
#endif  /* defined(CHAMELEON_HAVE_GETOPT_LONG) */

static void
set_iparam_default(int *iparam){

    memset(iparam, 0, IPARAM_SIZEOF*sizeof(int));

    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_THRDNBR_SUBGRP] = 1;
    iparam[IPARAM_M             ] = -1;
    iparam[IPARAM_N             ] = -1;
    iparam[IPARAM_K             ] = 1;
    iparam[IPARAM_LDA           ] = -1;
    iparam[IPARAM_LDB           ] = -1;
    iparam[IPARAM_LDC           ] = -1;
    iparam[IPARAM_MB            ] = 320;
    iparam[IPARAM_NB            ] = 320;
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NITER         ] = 1;
    iparam[IPARAM_WARMUP        ] = 1;
    iparam[IPARAM_BIGMAT        ] = 1;
    iparam[IPARAM_ASYNC         ] = 1;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_INPLACE       ] = MORSE_OUTOFPLACE;
    iparam[IPARAM_NMPI          ] = 1;
    iparam[IPARAM_P             ] = 1;
    iparam[IPARAM_Q             ] = 1;
    iparam[IPARAM_PRINT_WARNINGS] = 1;
    iparam[IPARAM_LOWLVL_TREE   ] = -1;
    iparam[IPARAM_HIGHLVL_TREE  ] = -1;
    iparam[IPARAM_QR_TS_SZE     ] = -1;
    iparam[IPARAM_QR_HLVL_SZE   ] = -1;
    iparam[IPARAM_QR_DOMINO     ] = -1;
}

static inline int
read_integer_from_options(int long_index, int opt_char)
{
    char *endptr;
    long int value;
    (void) long_index;

    value = strtol(optarg, &endptr, 10);
    if ( *optarg == '\0' || *endptr != '\0' ) {
#ifdef CHAMELEON_HAVE_GETOPT_LONG
        if ( long_index < 0 ) {
#endif
            fprintf(stderr, "Invalid numeric value '%s' for '-%c' parameter\n", optarg, opt_char);
#ifdef CHAMELEON_HAVE_GETOPT_LONG
        } else {
            fprintf(stderr, "Invalid numeric value '%s' for '--%s' parameter\n", optarg, long_options[long_index].name);
        }
#endif
        exit(EXIT_FAILURE);
    }
    if ( value > INT_MAX || value < INT_MIN ) {
        fprintf(stderr, "Out of range integer '%ld'\n", value);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

void
parse_arguments(int *_argc, char ***_argv, int *iparam, int *start, int *stop, int*step)
{
    int opt = -1;
    int c;
    int argc = *_argc;
    char **argv = *_argv;

    optind = 1;
    do {
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
        opt = -1;
        c = getopt_long(argc, argv, GETOPT_STRING,
                             long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(CHAMELEON_HAVE_GETOPT_LONG) */

        switch(c)
        {
        // Configuration
        case 't' : iparam[IPARAM_THRDNBR       ] = read_integer_from_options(opt, c); break;
        case 'g' : iparam[IPARAM_NCUDAS        ] = read_integer_from_options(opt, c); break;
        case 'P' : iparam[IPARAM_P             ] = read_integer_from_options(opt, c); break;
        case '8' : iparam[IPARAM_NO_CPU        ] = 1; break;
        // Matrix parameters
        case 'M' :
        case 'm' : iparam[IPARAM_M             ] = read_integer_from_options(opt, c); break;
        case 'n' : iparam[IPARAM_N             ] = read_integer_from_options(opt, c); break;
        case 'N' : get_range(optarg, start, stop, step); break;
        case 'K' :
        case 'k' : iparam[IPARAM_K             ] = read_integer_from_options(opt, c); break;
        case 'b' : iparam[IPARAM_NB            ] = read_integer_from_options(opt, c);
                   iparam[IPARAM_MB            ] = read_integer_from_options(opt, c); break;
        case 'i' : iparam[IPARAM_IB            ] = read_integer_from_options(opt, c); break;
        case 'x' : iparam[IPARAM_MX            ] = read_integer_from_options(opt, c); break;
        case 'X' : iparam[IPARAM_NX            ] = read_integer_from_options(opt, c); break;
        // Check/prints
        case '1' : iparam[IPARAM_NITER         ] = read_integer_from_options(opt, c); break;
        case 'W' : iparam[IPARAM_PRINT_WARNINGS] = 0; break;
        case 'w' : iparam[IPARAM_WARMUP        ] = 0; break;
        case 'c' : iparam[IPARAM_CHECK         ] = 1; break;
        case 'C' : iparam[IPARAM_INVERSE       ] = 1; break;
        // Profiling
        case 'T' : iparam[IPARAM_TRACE         ] = 1; break;
        case '2' : iparam[IPARAM_PROGRESS      ] = 1; break;
        case 'd' : iparam[IPARAM_DAG           ] = 1; break;
        case 'p' : iparam[IPARAM_PROFILE       ] = 1; break;
        // HQR options
        case 'a' : iparam[IPARAM_RHBLK         ] = read_integer_from_options(opt, c); break;
        case 'l' : iparam[IPARAM_LOWLVL_TREE   ] = read_integer_from_options(opt, c); break;
        case 'L' : iparam[IPARAM_HIGHLVL_TREE  ] = read_integer_from_options(opt, c); break;
        case 'D' : iparam[IPARAM_QR_DOMINO     ] = 1; break;
        //Other
        case '9' : iparam[IPARAM_MODE          ] = read_integer_from_options(opt, c); break;
        case '3' : iparam[IPARAM_BIGMAT        ] = 0; break;
        case 's' : iparam[IPARAM_ASYNC         ] = 0; break;
        case 'o' : iparam[IPARAM_OOC           ] = 1; break;
        case 'G' : iparam[IPARAM_GEMM3M        ] = 1; break;
        case '5' : iparam[IPARAM_BOUND         ] = 1; break;
        case 'h' :
        case '?' :
            show_help(argv[0]); exit(EXIT_FAILURE);
        default:
            break;
        }
    } while(-1 != c);
}

int
main(int argc, char *argv[]) {
    int i, m, n, mx, nx;
    int status;
    int nbnode = 1;
    int start =  500;
    int stop  = 5000;
    int step  =  500;
    int iparam[IPARAM_SIZEOF];
    int success = 0;

    set_iparam_default(iparam);

    parse_arguments(&argc, &argv, iparam, &start, &stop, &step);

#if !defined(CHAMELEON_USE_CUDA)
    if (iparam[IPARAM_NCUDAS] != 0){
    	fprintf(stderr, "ERROR: CHAMELEON_USE_CUDA is not defined. "
    			"The number of CUDA devices must be set to 0 (--gpus=0).\n");
    	return EXIT_FAILURE;
    }
#endif

    n  = iparam[IPARAM_N];
    m  = iparam[IPARAM_M];
    mx = iparam[IPARAM_MX];
    nx = iparam[IPARAM_NX];

    /* Initialize MORSE */
    MORSE_Init( iparam[IPARAM_THRDNBR],
                iparam[IPARAM_NCUDAS] );

    /* Get the number of threads set by the runtime */
    iparam[IPARAM_THRDNBR] = MORSE_GetThreadNbr();

    /* Stops profiling here to avoid profiling uninteresting routines.
       It will be reactivated in the time_*.c routines with the macro START_TIMING() */
    RUNTIME_stop_profiling();

    MORSE_Disable(MORSE_AUTOTUNING);
    MORSE_Set(MORSE_TILE_SIZE,        iparam[IPARAM_NB] );
    MORSE_Set(MORSE_INNER_BLOCK_SIZE, iparam[IPARAM_IB] );

    /* Householder mode */
    if (iparam[IPARAM_RHBLK] < 1) {
        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_FLAT_HOUSEHOLDER);
    } else {
        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_TREE_HOUSEHOLDER);
        MORSE_Set(MORSE_HOUSEHOLDER_SIZE, iparam[IPARAM_RHBLK]);
    }

    if (iparam[IPARAM_PROFILE] == 1) {
        MORSE_Enable(MORSE_PROFILING_MODE);
    }

    if (iparam[IPARAM_PROGRESS] == 1) {
        MORSE_Enable(MORSE_PROGRESS);
    }

    if (iparam[IPARAM_PRINT_WARNINGS] == 0) {
        MORSE_Disable(MORSE_WARNINGS);
    }

    if (iparam[IPARAM_GEMM3M] == 1) {
        MORSE_Enable(MORSE_GEMM3M);
    }

#if defined(CHAMELEON_USE_MPI)
    nbnode = MORSE_Comm_size();
    iparam[IPARAM_NMPI] = nbnode;
    /* Check P */
    if ( (iparam[IPARAM_P] > 1) &&
         (nbnode % iparam[IPARAM_P] != 0) ) {
      fprintf(stderr, "ERROR: %d doesn't divide the number of node %d\n",
              iparam[IPARAM_P], nbnode );
      return EXIT_FAILURE;
    }
#endif
    iparam[IPARAM_Q] = nbnode / iparam[IPARAM_P];

    /* Layout conversion */
    MORSE_Set(MORSE_TRANSLATION_MODE, iparam[IPARAM_INPLACE]);

    if ( MORSE_My_Mpi_Rank() == 0 ) {
        print_header( argv[0], iparam);
    }

    if (step < 1) step = 1;

    status = Test( -1, iparam ); /* print header */
    if (status != MORSE_SUCCESS) return status;
    if ( n == -1 ){
        for (i = start; i <= stop; i += step)
        {
            if ( nx > 0 ) {
                iparam[IPARAM_M] = i;
                iparam[IPARAM_N] = chameleon_max(1, i/nx);
            }
            else if ( mx > 0 ) {
                iparam[IPARAM_M] = chameleon_max(1, i/mx);
                iparam[IPARAM_N] = i;
            }
            else {
                if ( m == -1 ) {
                    iparam[IPARAM_M] = i;
                }
                iparam[IPARAM_N] = i;
            }
            status = Test( iparam[IPARAM_N], iparam );
            if (status != MORSE_SUCCESS) {
                return status;
            }
            success += status;
        }
    }
    else {
        if ( m == -1 ) {
            iparam[IPARAM_M] = n;
        }
        iparam[IPARAM_N] = n;
        status = Test( iparam[IPARAM_N], iparam );
        if (status != MORSE_SUCCESS) return status;
        success += status;
    }
    MORSE_Finalize();
    return success;
}

