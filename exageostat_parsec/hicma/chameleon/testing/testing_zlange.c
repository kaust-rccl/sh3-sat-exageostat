/**
 *
 * @file testing_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
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
#include "testing_zauxiliary.h"

int testing_zlange(int argc, char **argv)
{
    int hres = 0;
    /* Check for number of arguments*/
    if ( argc < 3) {
        USAGE("LANGE", "M N LDA",
              "   - M      : number of rows of matrices A and C\n"
              "   - N      : number of columns of matrices B and C\n"
              "   - LDA    : leading dimension of matrix A\n");
        return -1;
    }
    int M = atoi(argv[0]);
    int N = atoi(argv[1]);
    int LDA = atoi(argv[2]);
    int LDAxN = LDA*N;
    int n, u;
    double eps;

    /* Allocate Data */
    MORSE_Complex64_t *A    = (MORSE_Complex64_t *)malloc(LDAxN*sizeof(MORSE_Complex64_t));
    double            *work = (double*) malloc(max(M,N)*sizeof(double));
    double normmorse, normlapack, result;

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZLANGE ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZLANGE
     */

    /* Initialize A, B, C */
    MORSE_zplrnt( M, N, A, LDA, 3436 );

    /* MORSE ZLANGE */
    for(n=0; n<4; n++) {
        normmorse  = MORSE_zlange(norm[n], M, N, A, LDA);
        normlapack = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(norm[n]), M, N, A, LDA, work);
        printf("LAPACK %e, CHAMELEON %e\n", normlapack, normmorse);
        result = fabs(normmorse - normlapack) / (normlapack * eps);

        switch(norm[n]) {
        case MorseMaxNorm:
            /* result should be perfectly equal */
            break;
        case MorseInfNorm:
            /* Sum order on the line can differ */
            result = result / (double)N;
            break;
        case MorseOneNorm:
            /* Sum order on the column can differ */
            result = result / (double)M;
            break;
        case MorseFrobeniusNorm:
            /* Sum order on every element can differ */
            result = result / ((double)M * (double)N);
            break;
        }

        printf("***************************************************\n");
        if ( result < 1. ) {
            printf(" ---- TESTING ZLANGE (%s)............... PASSED !\n", normstr[n]);
        }
        else {
            printf(" - TESTING ZLANGE (%s)... FAILED !\n", normstr[n]);    hres++;
        }
        printf("***************************************************\n");

    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* MORSE ZLANTR */
    for(n=1; n<3; n++) {
        for(u=0; u<2; u++) {
            int d;
            for(d=0; d<2; d++) {
                normmorse = MORSE_zlantr(norm[n], uplo[u], diag[d], M, N, A, LDA);
                normlapack = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, morse_lapack_const(norm[n]), morse_lapack_const(uplo[u]),
                                                 morse_lapack_const(diag[d]), M, N, A, LDA, work);
                printf("LAPACK %e, CHAMELEON %e\n", normlapack, normmorse);

                result = fabs(normmorse - normlapack) / (normlapack * eps);
                switch(norm[n]) {
                case MorseMaxNorm:
                    /* result should be perfectly equal */
                    break;
                case MorseInfNorm:
                    /* Sum order on the line can differ */
                    result = result / (double)N;
                    break;
                case MorseOneNorm:
                    /* Sum order on the column can differ */
                    result = result / (double)M;
                    break;
                case MorseFrobeniusNorm:
                    /* Sum oreder on every element can differ */
                    result = result / ((double)M * (double)N);
                    break;
                }

                printf("***************************************************\n");
                if ( result < 1. ) {
                    printf(" ---- TESTING ZLANTR (%s, %s, %s)......... PASSED !\n",
                           normstr[n], uplostr[u], diagstr[d]);
                }
                else {
                    printf(" - TESTING ZLANTR (%s, %s, %s)... FAILED !\n",
                           normstr[n], uplostr[u], diagstr[d]);
                }
                printf("***************************************************\n");
            }
        }
    }
#endif

    /* MORSE ZLANSY */
    for(n=0; n<4; n++) {
        for(u=0; u<2; u++) {
            normmorse = MORSE_zlansy(norm[n], uplo[u], min(M,N), A, LDA);
            normlapack = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, morse_lapack_const(norm[n]), morse_lapack_const(uplo[u]), min(M,N), A, LDA, work);
            printf("LAPACK %e, CHAMELEON %e\n", normlapack, normmorse);

            result = fabs(normmorse - normlapack) / (normlapack * eps);
            switch(norm[n]) {
            case MorseMaxNorm:
                /* result should be perfectly equal */
                break;
            case MorseInfNorm:
                /* Sum order on the line can differ */
                result = result / (double)N;
                break;
            case MorseOneNorm:
                /* Sum order on the column can differ */
                result = result / (double)M;
                break;
            case MorseFrobeniusNorm:
                /* Sum oreder on every element can differ */
                result = result / ((double)M * (double)N);
                break;
            }

            printf("***************************************************\n");
            if ( result < 1. ) {
                printf(" ---- TESTING ZLANSY (%s, %s)......... PASSED !\n", normstr[n], uplostr[u]);
            }
            else {
                printf(" - TESTING ZLANSY (%s, %s)... FAILED !\n", normstr[n], uplostr[u]);
            }
            printf("***************************************************\n");
        }
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* MORSE ZLANHE */
    {
        int j;
        for (j=0; j<min(M,N); j++) {
            A[j*LDA+j] -= I*cimag(A[j*LDA+j]);
        }
    }

    for(n=0; n<4; n++) {
        for(u=0; u<2; u++) {
            normmorse = MORSE_zlanhe(norm[n], uplo[u], min(M,N), A, LDA);
            normlapack = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, morse_lapack_const(norm[n]), morse_lapack_const(uplo[u]), min(M,N), A, LDA, work);
            printf("LAPACK %e, CHAMELEON %e\n", normlapack, normmorse);

            result = fabs(normmorse - normlapack) / (normlapack * eps);
            switch(norm[n]) {
            case MorseMaxNorm:
                /* result should be perfectly equal */
                break;
            case MorseInfNorm:
                /* Sum order on the line can differ */
                result = result / (double)N;
                break;
            case MorseOneNorm:
                /* Sum order on the column can differ */
                result = result / (double)M;
                break;
            case MorseFrobeniusNorm:
                /* Sum oreder on every element can differ */
                result = result / ((double)M * (double)N);
                break;
            }

            printf("***************************************************\n");
            if ( result < 1. ) {
                printf(" ---- TESTING ZLANHE (%s, %s)......... PASSED !\n", normstr[n], uplostr[u]);
            }
            else {
                printf(" - TESTING ZLANHE (%s, %s)... FAILED !\n", normstr[n], uplostr[u]);
            }
            printf("***************************************************\n");
        }
    }
#endif

    free(A);
    free(work);
    return 0 /*hres*/;
}
