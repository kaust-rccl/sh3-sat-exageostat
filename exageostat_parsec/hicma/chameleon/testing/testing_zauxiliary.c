/**
 *
 * @file testing_zauxiliary.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon MORSE_Complex64_t auxiliary testings routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cédric Castagnède
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif
#include <morse.h>
#include "testing_zauxiliary.h"

int   IONE     = 1;
int   ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int   format[6]= { MorseCM, MorseRM, MorseCCRB, MorseCRRB, MorseRCRB, MorseRRRB };
int   side[2]  = { MorseLeft,    MorseRight };
int   uplo[2]  = { MorseUpper,   MorseLower };
int   diag[2]  = { MorseNonUnit, MorseUnit  };
int   trans[3] = { MorseNoTrans, MorseTrans, MorseConjTrans };
int   itype[3]  = { 1, 2, 3 };
int   storev[2] = { MorseRowwise, MorseColumnwise };
int   norm[4]   = { MorseMaxNorm, MorseOneNorm, MorseInfNorm, MorseFrobeniusNorm };

char *formatstr[6]= { "CM", "RM", "CCRB", "CRRB", "RCRB", "RRRB"};
char *sidestr[2]  = { "Left ", "Right" };
char *uplostr[2]  = { "Upper", "Lower" };
char *diagstr[2]  = { "NonUnit", "Unit   " };
char *transstr[3] = { "N", "T", "H" };
char *itypestr[3]  = { "inv(U')xAxinv(U) or inv(L)xAxinv(L')", "UxAxU' or L'xAxL", "UxAxU' or L'xAxL" };
char *storevstr[2] = { "Rowwise", "Columnwise" };
char *normstr[4]  = { "Max", "One", "Inf", "Fro" };

#define map_cm(m, n, i, j)   ((i) + (j) * (m))
#define map_rm(m, n, i, j)   ((i) * (n) + (j))

int map_CM(int m, int n, int mb, int nb, int i, int j)
{
    int hres = map_cm(m, n, i, j);
    (void)mb;
    (void)nb;
    (void)n;
    return hres;
}

int map_RM(int m, int n, int mb, int nb, int i, int j)
{
    int hres = map_rm(m, n, i, j);
    (void)mb;
    (void)nb;
    (void)m;
    return hres;
}

int map_CCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_CRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

int (*formatmap[6])(int, int, int, int, int, int) = {  map_CM, map_RM, map_CCRB, map_CRRB, map_RCRB, map_RRRB };

int main (int argc, char **argv)
{
    int ncores, ngpus;
    int info = 0;
    char func[32];

    /* Check for number of arguments*/
    if ( argc < 4) {
        printf(" Proper Usage is : ./ztesting ncores ngpus FUNC ...\n"
               "   - ncores : number of cores\n"
               "   - ngpus  : number of GPUs\n"
               "   - FUNC   : name of function to test\n"
               "   - ... plus arguments depending on the testing function \n");
        exit(1);
    }

    sscanf( argv[1], "%d",   &ncores );
    sscanf( argv[2], "%d",   &ngpus  );
    sscanf( argv[3], "%31s",  func   );

    /* Initialize MORSE */
    /*if(nthreads_per_worker)
     MORSE_InitPar(ncores/nthreads_per_worker, ncudas, nthreads_per_worker);
     else*/
    MORSE_Init( ncores, ngpus);
    MORSE_Disable(MORSE_AUTOTUNING);
    MORSE_Set(MORSE_TILE_SIZE,         32 );
    MORSE_Set(MORSE_INNER_BLOCK_SIZE,   5 );

    argc -= 4;
    argv += 4;
    info  = 0;

    /*
     * Norms
     */
    if ( strcmp(func, "LANGE") == 0 ) {
        info += testing_zlange( argc, argv );
    }
    /*
     * Blas Level 3
     */
    else if ( strcmp(func, "GEMM") == 0 ) {
        info += testing_zgemm( argc, argv );
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    else if ( strcmp(func, "HEMM") == 0 ) {
        info += testing_zhemm( argc, argv );
    }
    else if ( strcmp(func, "HERK") == 0 ) {
        info += testing_zherk( argc, argv );
    }
    else if ( strcmp(func, "HER2K") == 0 ) {
        info += testing_zher2k( argc, argv );
    }
#endif
    else if ( strcmp(func, "SYMM") == 0 ) {
        info += testing_zsymm( argc, argv );
    }
    else if ( strcmp(func, "SYRK") == 0 ) {
        info += testing_zsyrk( argc, argv );
    }
    else if ( strcmp(func, "SYR2K") == 0 ) {
        info += testing_zsyr2k( argc, argv );
    }
    else if ( strcmp(func, "TRMM") == 0 ) {
        info += testing_ztrmm( argc, argv );
    }
    else if ( strcmp(func, "TRSM") == 0 ) {
        info += testing_ztrsm( argc, argv );
    }
    else if ( strcmp(func, "PEMV") == 0 ) {
        info += testing_zpemv( argc, argv );
    }
    else if ( strcmp(func, "GEADD") == 0 ) {
        info = testing_zgeadd( argc, argv );
    }
    /*
     * Linear system
     */
    else if ( strcmp(func, "POSV") == 0 ) {
        info += testing_zposv( argc, argv );
    }
    else if ( strcmp(func, "GELS") == 0 ) {
        info += testing_zgels( argc, argv );
    }
    else if ( strcmp(func, "GESV_INCPIV") == 0 ) {
        info += testing_zgesv_incpiv( argc, argv );
    }
    else if ( strcmp(func, "GELS_HQR") == 0 ) {
        info += testing_zgels_hqr( argc, argv );
    }
    else if ( strcmp(func, "GELS_SYSTOLIC") == 0 ) {
        info += testing_zgels_systolic( argc, argv );
    }
    /* else if ( strcmp(func, "GESV") == 0 ) { */
    /*     info += testing_zgesv( argc, argv ); */
    /* } */
    /*
     * Matrix inversion
     */
    else if ( strcmp(func, "POTRI") == 0 ) {
        info += testing_zpotri( argc, argv );
    }
    /* else if ( strcmp(func, "GETRI") == 0 ) { */
    /*      info += testing_zgetri( argc, argv ); */
    /* } */
    /*
     * Eigenvalue Problems
     */
    /* else if ( strcmp(func, "HEEV") == 0 ) { */
    /*     info += testing_zheev( argc, argv ); */
    /* } */
    else if ( strcmp(func, "HEEVD") == 0 ) {
        info += testing_zheevd( argc, argv );
    }
    /* else if ( strcmp(func, "HEGV") == 0 ) { */
    /*     info += testing_zhegv( argc, argv ); */
    /* } */
    /* else if ( strcmp(func, "HEGVD") == 0 ) { */
    /*     info += testing_zhegv( argc, argv ); */
    /* } */
    /* else if ( strcmp(func, "HEGST") == 0 ) { */
    /*     info += testing_zhegst( argc, argv ); */
    /* } */
    /*
     * Singular Value Decomposition
     */
    else if ( strcmp(func, "GESVD") == 0 ) {
        info += testing_zgesvd( argc, argv );
    }
#ifdef DOUBLE
    /*
     * Mixed precision
     */
    /* else if ( strcmp(func, "CPOSV") == 0 ) { */
    /*     info += testing_zcposv( argc, argv ); */
    /* } */
    /* else if ( strcmp(func, "CGESV") == 0 ) { */
    /*     info += testing_zcgesv( argc, argv ); */
    /* } */
    /* else if ( strcmp(func, "CUNGESV") == 0 ) { */
    /*     info += testing_zcungesv( argc, argv ); */
    /* } */
#endif
    /*
     * Layout Transformation
     */
    /* else if ( strcmp(func, "GECFI") == 0 ) { */
    /*     info += testing_zgecfi( argc, argv ); */
    /* } */
    /* else if ( strcmp(func, "GETMI") == 0 ) { */
    /*     info += testing_zgetmi( argc, argv ); */
    /* } */
    else if ( strcmp(func, "GEQRF_QDWH") == 0 ) {
        info += testing_zgeqrf_qdwh( argc, argv );
    }
    else {
        fprintf(stderr, "Function unknown\n");
    }

    if ( info == -1 ) {
        printf( "TESTING %s FAILED : incorrect number of arguments\n", func);
    } else if ( info == -2 ) {
        printf( "TESTING %s FAILED : not enough memory\n", func);
    }

    MORSE_Finalize();

    return info;
}
