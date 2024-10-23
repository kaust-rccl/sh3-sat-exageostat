/**
 *
 * @file morse_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon MORSE_complex64_t wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#ifndef _MORSE_Z_H_
#define _MORSE_Z_H_

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
int MORSE_zgeadd(MORSE_enum trans, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t beta, MORSE_Complex64_t *B, int LDB);
//int MORSE_zgebrd(int M, int N, MORSE_Complex64_t *A, int LDA, double *D, double *E, MORSE_desc_t *descT);
//int MORSE_zgecon(MORSE_enum norm, int N, MORSE_Complex64_t *A, int LDA, double anorm, double *rcond);
//int MORSE_zpocon(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, double anorm, double *rcond);
int MORSE_zgelqf(int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT);
int MORSE_zgelqs(int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
int MORSE_zgels(MORSE_enum trans, int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
int MORSE_zgemm(MORSE_enum transA, MORSE_enum transB, int M, int N, int K, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zgeqrf(int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT);
int MORSE_zgeqrs(int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
//int MORSE_zgesv(int N, int NRHS, MORSE_Complex64_t *A, int LDA, int *IPIV, MORSE_Complex64_t *B, int LDB);
int MORSE_zgesv_incpiv(int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descL, int *IPIV, MORSE_Complex64_t *B, int LDB);
int MORSE_zgesv_nopiv(int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zgesvd(MORSE_enum jobu, MORSE_enum jobvt, int M, int N, MORSE_Complex64_t *A, int LDA, double *S, MORSE_desc_t *descT, MORSE_Complex64_t *U, int LDU, MORSE_Complex64_t *VT, int LDVT);
//int MORSE_zgetrf(int M, int N, MORSE_Complex64_t *A, int LDA, int *IPIV);
int MORSE_zgetrf_incpiv(int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descL, int *IPIV);
int MORSE_zgetrf_nopiv(int M, int N, MORSE_Complex64_t *A, int LDA);
//int MORSE_zgetri(int N, MORSE_Complex64_t *A, int LDA, int *IPIV);
//int MORSE_zgetrs(MORSE_enum trans, int N, int NRHS, MORSE_Complex64_t *A, int LDA, int *IPIV, MORSE_Complex64_t *B, int LDB);
int MORSE_zgetrs_incpiv(MORSE_enum trans, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descL, int *IPIV, MORSE_Complex64_t *B, int LDB);
int MORSE_zgetrs_nopiv(MORSE_enum trans, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zhemm(MORSE_enum side, MORSE_enum uplo, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zherk(MORSE_enum uplo, MORSE_enum trans, int N, int K, double alpha, MORSE_Complex64_t *A, int LDA, double beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zher2k(MORSE_enum uplo, MORSE_enum trans, int N, int K, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, double beta, MORSE_Complex64_t *C, int LDC);
//int MORSE_zheev(MORSE_enum jobz, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, double *W, MORSE_desc_t *descT, MORSE_Complex64_t *Q, int LDQ);
int MORSE_zheevd(MORSE_enum jobz, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, double *W, MORSE_desc_t *descT);
//int MORSE_zhegv(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, double *W, MORSE_desc_t *descT, MORSE_Complex64_t *Q, int LDQ);
//int MORSE_zhegvd(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, double *W, MORSE_desc_t *descT, MORSE_Complex64_t *Q, int LDQ);
//int MORSE_zhegst(MORSE_enum itype, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zhetrd(MORSE_enum jobz, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, double *D, double *E, MORSE_desc_t *descT, MORSE_Complex64_t *Q, int LDQ);
int MORSE_zlacpy(MORSE_enum uplo, int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
double MORSE_zlange(MORSE_enum norm, int M, int N, MORSE_Complex64_t *A, int LDA);
double MORSE_zlanhe(MORSE_enum norm, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
double MORSE_zlansy(MORSE_enum norm, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
double MORSE_zlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, int M, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zlascal(MORSE_enum uplo, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA);
int MORSE_zlaset(MORSE_enum uplo, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t beta, MORSE_Complex64_t *A, int LDA);
//int MORSE_zlaswp(int N, MORSE_Complex64_t *A, int LDA, int K1, int K2, int *IPIV, int INCX);
//int MORSE_zlaswpc(int N, MORSE_Complex64_t *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int MORSE_zlauum(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zplghe( double bump, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, unsigned long long int seed );
int MORSE_zplgsy( MORSE_Complex64_t bump, MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, unsigned long long int seed );
int MORSE_zplrnt( int M, int N, MORSE_Complex64_t *A, int LDA, unsigned long long int seed );
int MORSE_zposv(MORSE_enum uplo, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zpotrf(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zsytrf(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zpotri(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zpotrimm(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t *C, int LDC);
int MORSE_zpotrs(MORSE_enum uplo, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zsymm(MORSE_enum side, MORSE_enum uplo, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zsyrk(MORSE_enum uplo, MORSE_enum trans, int N, int K, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zsyr2k(MORSE_enum uplo, MORSE_enum trans, int N, int K, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC);
int MORSE_zsysv(MORSE_enum uplo, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_zsytrs(MORSE_enum uplo, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_ztpgqrt( int M, int N, int K, int L, MORSE_Complex64_t *V1, int LDV1, MORSE_desc_t *descT1, MORSE_Complex64_t *V2, int LDV2, MORSE_desc_t *descT2, MORSE_Complex64_t *Q1, int LDQ1, MORSE_Complex64_t *Q2, int LDQ2 );
int MORSE_ztpqrt( int M, int N, int L, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_desc_t *descT );
int MORSE_ztradd(MORSE_enum uplo, MORSE_enum trans, int M, int N, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t beta, MORSE_Complex64_t *B, int LDB);
int MORSE_ztrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_ztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_ztrsmpl(int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descL, int *IPIV, MORSE_Complex64_t *B, int LDB);
int MORSE_ztrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB);
int MORSE_ztrtri(MORSE_enum uplo, MORSE_enum diag, int N, MORSE_Complex64_t *A, int LDA);
int MORSE_zunglq(int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
int MORSE_zungqr(int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunmlq(MORSE_enum side, MORSE_enum trans, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunmqr(MORSE_enum side, MORSE_enum trans, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_Complex64_t *B, int LDB);

//int MORSE_zgecfi(int m, int n, MORSE_Complex64_t *A, MORSE_enum fin, int imb, int inb, MORSE_enum fout, int omb, int onb);
//int MORSE_zgetmi(int m, int n, MORSE_Complex64_t *A, MORSE_enum fin, int mb,  int nb);

/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
int MORSE_zgeadd_Tile(MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B);
//int MORSE_zgebrd_Tile(MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T);
//int MORSE_zgecon_Tile(MORSE_enum norm, MORSE_desc_t *A, double anorm, double *rcond);
//int MORSE_zpocon_Tile(MORSE_enum uplo, MORSE_desc_t *A, double anorm, double *rcond);
int MORSE_zgelqf_Tile(MORSE_desc_t *A, MORSE_desc_t *T);
int MORSE_zgelqs_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_zgels_Tile(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_zgemm_Tile(MORSE_enum transA, MORSE_enum transB, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C);
int MORSE_zgeqrf_Tile(MORSE_desc_t *A, MORSE_desc_t *T);
int MORSE_zgeqrs_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
//int MORSE_zgesv_Tile(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B);
int MORSE_zgesv_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_zgesv_nopiv_Tile(MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zgesvd_Tile(MORSE_enum jobu, MORSE_enum jobvt, MORSE_desc_t *A, double *S, MORSE_desc_t *T, MORSE_Complex64_t *U, int LDU, MORSE_Complex64_t *VT, int LDVT);
//int MORSE_zgetrf_Tile(MORSE_desc_t *A, int *IPIV);
int MORSE_zgetrf_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV);
int MORSE_zgetrf_nopiv_Tile(MORSE_desc_t *A);
//int MORSE_zgetri_Tile(MORSE_desc_t *A, int *IPIV);
//int MORSE_zgetrs_Tile(MORSE_enum trans, MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B);
int MORSE_zgetrs_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_zgetrs_nopiv_Tile(MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zhemm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C);
int MORSE_zherk_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C);
int MORSE_zher2k_Tile(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
//int MORSE_zheev_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, MORSE_Complex64_t *Q, int LDQ);
int MORSE_zheevd_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T);
//int MORSE_zhegv_Tile( MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q);
//int MORSE_zhegvd_Tile(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q);
//int MORSE_zhegst_Tile(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zhetrd_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_Complex64_t *Q, int LDQ);
int MORSE_zlacpy_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
double MORSE_zlange_Tile(MORSE_enum norm, MORSE_desc_t *A);
double MORSE_zlanhe_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A);
double MORSE_zlansy_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A);
double MORSE_zlantr_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A);
int MORSE_zlascal_Tile(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A);
int MORSE_zlaset_Tile(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_Complex64_t beta, MORSE_desc_t *A);
//int MORSE_zlaswp_Tile(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX);
//int MORSE_zlaswpc_Tile(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX);
int MORSE_zlauum_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_zplghe_Tile(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed );
int MORSE_zplgsy_Tile(MORSE_Complex64_t bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed );
int MORSE_zplrnt_Tile(MORSE_desc_t *A, unsigned long long int seed );
int MORSE_zposv_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zpotrf_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_zsytrf_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_zpotri_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_zpotrimm_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C);
int MORSE_zpotrs_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zsymm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C);
int MORSE_zsyrk_Tile(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *C);
int MORSE_zsyr2k_Tile(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C);
int MORSE_zsysv_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_zsytrs_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_ztpgqrt_Tile( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2 );
int MORSE_ztpqrt_Tile( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T );
int MORSE_ztradd_Tile(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B);
int MORSE_ztrmm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_ztrsm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_ztrsmpl_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_ztrsmrv_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_ztrtri_Tile(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A);
int MORSE_zunglq_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_zungqr_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_zunmlq_Tile(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_zunmqr_Tile(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);

/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
int MORSE_zgeadd_Tile_Async(MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgebrd_Tile_Async(MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgecon_Tile_Async(MORSE_enum norm, MORSE_desc_t *A, double anorm, double *rcond, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zpocon_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, double anorm, double *rcond, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgelqf_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgelqs_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgels_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgemm_Tile_Async(MORSE_enum transA, MORSE_enum transB, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgeqrf_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgeqrs_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgesv_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgesv_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgesv_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgesvd_Tile_Async(MORSE_enum jobu, MORSE_enum jobvt, MORSE_desc_t *A, double *S, MORSE_desc_t *T, MORSE_Complex64_t *U, int LDU, MORSE_Complex64_t *VT, int LDVT, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgetrf_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgetrf_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgetrf_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgetri_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgetrs_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgetrs_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgetrs_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zhemm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zherk_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zher2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zheev_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, MORSE_Complex64_t *Q, int LDQ, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zheevd_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zhegv_Tile_Async( MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zhegvd_Tile_Async(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zhegst_Tile_Async(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zhetrd_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_Complex64_t *Q, int LDQ, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlacpy_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlange_Tile_Async(MORSE_enum norm, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlanhe_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlansy_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlantr_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlascal_Tile_Async(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlaset_Tile_Async(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_Complex64_t beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zlaswp_Tile_Async(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zlaswpc_Tile_Async(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zlauum_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zplghe_Tile_Async(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_zplgsy_Tile_Async(MORSE_Complex64_t bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_zplrnt_Tile_Async(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_zposv_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zpotrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsytrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zpotri_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zpotrimm_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zpotrs_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsysv_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsytrs_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsymm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsyrk_Tile_Async(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zsyr2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztpgqrt_Tile_Async( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_ztpqrt_Tile_Async( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_ztradd_Tile_Async(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztrmm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztrsm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztrsmpl_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztrsmrv_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_ztrtri_Tile_Async(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunglq_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zungqr_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunmlq_Tile_Async(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunmqr_Tile_Async(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);

//int MORSE_zgecfi_Async(int m, int n, MORSE_Complex64_t *A, MORSE_enum f_in, int imb, int inb, MORSE_enum f_out, int omb, int onb, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zgetmi_Async(int m, int n, MORSE_Complex64_t *A, MORSE_enum f_in, int mb,  int inb, MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of libhqr dependent functions.
 */
/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
int MORSE_zgels_param(const libhqr_tree_t *qrtree, MORSE_enum trans, int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, MORSE_Complex64_t *B, int LDB);
int MORSE_zgelqf_param(const libhqr_tree_t *qrtree, int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT);
int MORSE_zgelqs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, MORSE_Complex64_t *B, int LDB);
int MORSE_zgeqrf_param(const libhqr_tree_t *qrtree, int M, int N, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT);
int MORSE_zgeqrs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunglq_param(const libhqr_tree_t *qrtree, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, MORSE_Complex64_t *B, int LDB);
int MORSE_zungqr_param(const libhqr_tree_t *qrtree, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunmlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_desc_t *TT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunmlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descT, MORSE_desc_t *TT, MORSE_Complex64_t *B, int LDB);
int MORSE_zunmqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, MORSE_Complex64_t *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, MORSE_Complex64_t *B, int LDB);
/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
int MORSE_zgels_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zgelqf_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT);
int MORSE_zgelqs_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zgeqrf_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT);
int MORSE_zgeqrs_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zunglq_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zungqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zungqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zunmlq_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_zunmqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
int MORSE_zgels_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgelqf_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgelqs_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgeqrf_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zgeqrs_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunglq_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zungqr_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunmlq_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_zunmqr_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 */
int MORSE_Alloc_Workspace_zgesv_incpiv(        int N, MORSE_desc_t **descL, int **IPIV, int p, int q);
int MORSE_Alloc_Workspace_zgetrf_incpiv(int M, int N, MORSE_desc_t **descL, int **IPIV, int p, int q);

int MORSE_Alloc_Workspace_zgebrd(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgeev( int N,        MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgehrd(int N,        MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgelqf(int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_zgels( int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_zgeqrf(int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_zgesvd(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zheev( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zheevd( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zhegv( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zhegvd( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zhetrd(int M, int N, MORSE_desc_t **descT, int p, int q);

/**
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 */
int MORSE_Alloc_Workspace_zgetri_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *W);

/* Warning: Those functions are deprecated */
int MORSE_Alloc_Workspace_zgelqf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgels_Tile( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgeqrf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_zgesv_incpiv_Tile (int N, MORSE_desc_t **descL, int **IPIV, int p, int q);
int MORSE_Alloc_Workspace_zgetrf_incpiv_Tile(int N, MORSE_desc_t **descL, int **IPIV, int p, int q);

/**
 *  Auxiliary function prototypes
 */
int MORSE_zLapack_to_Tile(MORSE_Complex64_t *Af77, int LDA, MORSE_desc_t *A);
int MORSE_zTile_to_Lapack(MORSE_desc_t *A, MORSE_Complex64_t *Af77, int LDA);

/**
 *  User Builder function prototypes
 */
int MORSE_zbuild(MORSE_enum uplo, int M, int N, MORSE_Complex64_t *A, int LDA, void *user_data, void* user_build_callback);
int MORSE_zbuild_Tile(MORSE_enum uplo,  MORSE_desc_t *A, void *user_data, void* user_build_callback );
int MORSE_zbuild_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t  *request);

END_C_DECLS

#endif
