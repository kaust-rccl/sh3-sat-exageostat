/**
 *
 * @file timing_zauxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#ifndef TIMING_ZAUXILIARY_H
#define TIMING_ZAUXILIARY_H

int    z_check_orthogonality   (int M, int N, int LDQ, MORSE_Complex64_t *Q);
int    z_check_QRfactorization (int M, int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, MORSE_Complex64_t *Q);
int    z_check_LLTfactorization(int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, int uplo);
double z_check_gemm(MORSE_enum transA, MORSE_enum transB, int M, int N, int K,
                   MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                   MORSE_Complex64_t *B, int LDB,
                   MORSE_Complex64_t beta, MORSE_Complex64_t *Cmorse,
                   MORSE_Complex64_t *Cref, int LDC,
                   double *Cinitnorm, double *Cmorsenorm, double *Clapacknorm );

double z_check_trsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum trans, MORSE_enum diag,
           int M, int NRHS, MORSE_Complex64_t alpha,
           MORSE_Complex64_t *A, int LDA,
           MORSE_Complex64_t *Bmorse, MORSE_Complex64_t *Bref, int LDB,
           double *Binitnorm, double *Bmorsenorm, double *Blapacknorm );

double z_check_solution(int M, int N, int NRHS,
                      MORSE_Complex64_t *A1, int LDA,
                      MORSE_Complex64_t *B1, MORSE_Complex64_t *B2, int LDB,
                      double *anorm, double *bnorm, double *xnorm);

int zcheck_inverse(int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2,
                         int LDA, MORSE_enum uplo, double *rnorm, double *anorm, double *ainvnorm);


#endif /* TIMING_ZAUXILIARY_H */
