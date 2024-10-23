/**
 *
 * @file cudablas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU MORSE_Complex64_t kernels header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 */
#ifndef _MORSE_CUDA_ZBLAS_H_
#define _MORSE_CUDA_ZBLAS_H_

/**
 *  Declarations of cuda kernels - alphabetical order
 */
int CUDA_zgeadd( MORSE_enum trans, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *beta, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_zgemerge( MORSE_enum side, MORSE_enum diag, int M, int N, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, CUBLAS_STREAM_PARAM );
int CUDA_zgemm(  MORSE_enum transa, MORSE_enum transb, int m, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zhemm(  MORSE_enum side, MORSE_enum uplo, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zher2k( MORSE_enum uplo, MORSE_enum trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, double *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zherfb( MORSE_enum uplo, int n, int k, int ib, int nb, const cuDoubleComplex *A, int lda, const cuDoubleComplex *T, int ldt, cuDoubleComplex *C, int ldc, cuDoubleComplex *WORK, int ldwork, CUBLAS_STREAM_PARAM );
int CUDA_zherk(  MORSE_enum uplo, MORSE_enum trans, int n, int k, double *alpha, const cuDoubleComplex *A, int lda, double *beta, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_zlarfb(MORSE_enum side, MORSE_enum trans, MORSE_enum direct, MORSE_enum storev, int M, int N, int K, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *C, int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );
int CUDA_zparfb(MORSE_enum side, MORSE_enum trans, MORSE_enum direct, MORSE_enum storev, int M1, int N1, int M2, int N2, int K, int L, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LDWORK, cuDoubleComplex *WORKC, int LDWORKC, CUBLAS_STREAM_PARAM );
int CUDA_zsymm(  MORSE_enum side, MORSE_enum uplo, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zsyr2k( MORSE_enum uplo, MORSE_enum trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zsyrk(  MORSE_enum uplo, MORSE_enum trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_ztpmqrt( MORSE_enum side, MORSE_enum trans, int M, int N, int K, int L, int IB, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cuDoubleComplex *WORK, CUBLAS_STREAM_PARAM );
int CUDA_ztrmm(  MORSE_enum side, MORSE_enum uplo, MORSE_enum transa, MORSE_enum diag, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_ztrsm(  MORSE_enum side, MORSE_enum uplo, MORSE_enum transa, MORSE_enum diag, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_ztsmlq( MORSE_enum side, MORSE_enum trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LDWORK, cuDoubleComplex *WORKC, int LDWORKC, CUBLAS_STREAM_PARAM );
int CUDA_ztsmqr( MORSE_enum side, MORSE_enum trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LDWORK, cuDoubleComplex *WORKC, int LDWORKC, CUBLAS_STREAM_PARAM );
int CUDA_zttmqr( MORSE_enum side, MORSE_enum trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LDWORK, cuDoubleComplex *WORKC, int LDWORKC, CUBLAS_STREAM_PARAM );
int CUDA_zunmlqt(MORSE_enum side, MORSE_enum trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );
int CUDA_zunmqrt(MORSE_enum side, MORSE_enum trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );

#endif
