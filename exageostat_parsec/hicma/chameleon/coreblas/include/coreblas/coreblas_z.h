/**
 *
 * @file coreblas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU MORSE_Complex64_t kernels header
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
#ifndef _MORSE_CORE_ZBLAS_H_
#define _MORSE_CORE_ZBLAS_H_

/**
 *  Declarations of serial kernels - alphabetical order
 */
void CORE_dzasum(int storev, MORSE_enum uplo, int M, int N,
                 const MORSE_Complex64_t *A, int lda, double *work);
void CORE_zbrdalg(MORSE_enum uplo, int N, int NB,
                  const MORSE_desc_t *pA, MORSE_Complex64_t *C, MORSE_Complex64_t *S,
                  int i, int j, int m, int grsiz);
int CORE_zgbelr(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgbrce(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgblrx(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zaxpy(int M, MORSE_Complex64_t alpha,
               const MORSE_Complex64_t *A, int incA,
                     MORSE_Complex64_t *B, int incB);
int CORE_zgeadd(MORSE_enum trans, int M, int N,
                      MORSE_Complex64_t alpha,
                const MORSE_Complex64_t *A, int LDA,
                      MORSE_Complex64_t beta,
                      MORSE_Complex64_t *B, int LDB);
int CORE_zlascal( MORSE_enum uplo, int m, int n,
              MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int lda );
int  CORE_zgelqt(int M, int N, int IB,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU,
                 MORSE_Complex64_t *WORK);
void CORE_zgemm(MORSE_enum transA, MORSE_enum transB,
                int M, int N, int K,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                          const MORSE_Complex64_t *B, int LDB,
                MORSE_Complex64_t beta,        MORSE_Complex64_t *C, int LDC);
int  CORE_zgeqrt(int M, int N, int IB,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU, MORSE_Complex64_t *WORK);
int CORE_zgesplit(MORSE_enum side, MORSE_enum diag,
                  int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB);
int  CORE_zgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const MORSE_Complex64_t *L, int LDL,
                 MORSE_Complex64_t *A, int LDA);
int  CORE_zgessq(int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
int CORE_zgetf2_nopiv(int M, int N,
                  MORSE_Complex64_t *A, int LDA);
int  CORE_zgetrf(int M, int N,
                 MORSE_Complex64_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_zgetrf_incpiv(int M, int N, int IB,
                        MORSE_Complex64_t *A, int LDA,
                        int *IPIV, int *INFO);
int CORE_zgetrf_nopiv(int M, int N, int IB,
                      MORSE_Complex64_t *A, int LDA,
                      int *INFO);
int  CORE_zgetrf_reclap(int M, int N,
                        MORSE_Complex64_t *A, int LDA,
                        int *IPIV, int *info);
int  CORE_zgetrf_rectil(const MORSE_desc_t A, int *IPIV, int *info);
void CORE_zgetrip(int m, int n, MORSE_Complex64_t *A,
                  MORSE_Complex64_t *work);
void CORE_zhe2ge(MORSE_enum uplo, int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *B, int LDB);
int CORE_zhbelr(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhblrx(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhbrce(MORSE_enum uplo, int N,
                MORSE_desc_t *A, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                int st, int ed, int eltsize);
void CORE_zhbtype1cb(int N, int NB,
                     MORSE_Complex64_t *A, int LDA,
                     MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     MORSE_Complex64_t *WORK);
void CORE_zhbtype2cb(int N, int NB,
                     MORSE_Complex64_t *A, int LDA,
                     MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     MORSE_Complex64_t *WORK);
void CORE_zhbtype3cb(int N, int NB,
                     MORSE_Complex64_t *A, int LDA,
                     const MORSE_Complex64_t *V, const MORSE_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     MORSE_Complex64_t *WORK);
void CORE_zhegst(int itype, MORSE_enum uplo, int N,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *B, int LDB, int *INFO);
void CORE_zhemm(MORSE_enum side, MORSE_enum uplo,
                int M, int N,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                          const MORSE_Complex64_t *B, int LDB,
                MORSE_Complex64_t beta,        MORSE_Complex64_t *C, int LDC);
void CORE_zherk(MORSE_enum uplo, MORSE_enum trans,
                int N, int K,
                double alpha, const MORSE_Complex64_t *A, int LDA,
                double beta,        MORSE_Complex64_t *C, int LDC);
void CORE_zher2k(MORSE_enum uplo, MORSE_enum trans,
                 int N, int K,
                 MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                           const MORSE_Complex64_t *B, int LDB,
                 double beta,                    MORSE_Complex64_t *C, int LDC);
int  CORE_zhessq(MORSE_enum uplo, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
int  CORE_zherfb(MORSE_enum uplo, int N, int K, int IB, int NB,
                 const MORSE_Complex64_t *A,    int LDA,
                 const MORSE_Complex64_t *T,    int LDT,
                       MORSE_Complex64_t *C,    int LDC,
                       MORSE_Complex64_t *WORK, int LDWORK);
void CORE_zlacpy(MORSE_enum uplo, int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                       MORSE_Complex64_t *B, int LDB);
void CORE_zlange(int norm, int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlanhe(int norm, MORSE_enum uplo, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlansy(int norm, MORSE_enum uplo, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlantr(int norm, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *work, double *normA);
int CORE_zlarfb_gemm(MORSE_enum side, MORSE_enum trans, MORSE_enum direct, MORSE_enum storev,
                     int M, int N, int K,
                     const MORSE_Complex64_t *V, int LDV,
                     const MORSE_Complex64_t *T, int LDT,
                           MORSE_Complex64_t *C, int LDC,
                           MORSE_Complex64_t *WORK, int LDWORK);
int CORE_zlarfx2(MORSE_enum side, int N,
                 MORSE_Complex64_t V,
                 MORSE_Complex64_t TAU,
                 MORSE_Complex64_t *C1, int LDC1,
                 MORSE_Complex64_t *C2, int LDC2);
int CORE_zlarfx2c(MORSE_enum uplo,
                  MORSE_Complex64_t V,
                  MORSE_Complex64_t TAU,
                  MORSE_Complex64_t *C1,
                  MORSE_Complex64_t *C2,
                  MORSE_Complex64_t *C3);
int CORE_zlarfx2ce(MORSE_enum uplo,
                   MORSE_Complex64_t *V,
                   MORSE_Complex64_t *TAU,
                   MORSE_Complex64_t *C1,
                   MORSE_Complex64_t *C2,
                   MORSE_Complex64_t *C3);
void CORE_zlarfy(int N,
                 MORSE_Complex64_t *A, int LDA,
                 const MORSE_Complex64_t *V,
                 const MORSE_Complex64_t *TAU,
                 MORSE_Complex64_t *WORK);
void CORE_zlaset(MORSE_enum uplo, int n1, int n2,
                 MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                 MORSE_Complex64_t *tileA, int ldtilea);
void CORE_zlaset2(MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                  MORSE_Complex64_t *tileA, int ldtilea);
void CORE_zlaswp(int N, MORSE_Complex64_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_zlaswp_ontile( MORSE_desc_t descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlaswpc_ontile(MORSE_desc_t descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlatro(MORSE_enum uplo, MORSE_enum trans,
                 int M, int N,
                 const MORSE_Complex64_t *A, int LDA,
                       MORSE_Complex64_t *B, int LDB);
void CORE_zlauum(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA);
int CORE_zpamm(int op, MORSE_enum side, MORSE_enum storev,
               int M, int N, int K, int L,
               const MORSE_Complex64_t *A1, int LDA1,
                     MORSE_Complex64_t *A2, int LDA2,
               const MORSE_Complex64_t *V, int LDV,
                     MORSE_Complex64_t *W, int LDW);
int  CORE_zparfb(MORSE_enum side, MORSE_enum trans, MORSE_enum direct, MORSE_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       MORSE_Complex64_t *A1, int LDA1,
                       MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                       MORSE_Complex64_t *WORK, int LDWORK);
int CORE_zpemv(MORSE_enum trans, MORSE_enum storev,
               int M, int N, int L,
               MORSE_Complex64_t ALPHA,
               const MORSE_Complex64_t *A, int LDA,
               const MORSE_Complex64_t *X, int INCX,
               MORSE_Complex64_t BETA,
               MORSE_Complex64_t *Y, int INCY,
               MORSE_Complex64_t *WORK);
void CORE_zplghe(double bump, int m, int n, MORSE_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplgsy(MORSE_Complex64_t bump, int m, int n, MORSE_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplrnt(int m, int n, MORSE_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zpotrf(MORSE_enum uplo, int N, MORSE_Complex64_t *A, int LDA, int *INFO);
void CORE_zshift(int s, int m, int n, int L,
                 MORSE_Complex64_t *A);
void CORE_zshiftw(int s, int cl, int m, int n, int L,
                  MORSE_Complex64_t *A, MORSE_Complex64_t *W);
int  CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       MORSE_Complex64_t *A1, int LDA1,
                       MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *L1, int LDL1,
                 const MORSE_Complex64_t *L2, int LDL2,
                 const int *IPIV);
void CORE_zsymm(MORSE_enum side, MORSE_enum uplo,
                int M, int N,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                          const MORSE_Complex64_t *B, int LDB,
                MORSE_Complex64_t beta,        MORSE_Complex64_t *C, int LDC);
void CORE_zsyrk(MORSE_enum uplo, MORSE_enum trans,
                int N, int K,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                MORSE_Complex64_t beta,        MORSE_Complex64_t *C, int LDC);
void CORE_zsyr2k(MORSE_enum uplo, MORSE_enum trans,
                 int N, int K,
                 MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                           const MORSE_Complex64_t *B, int LDB,
                 MORSE_Complex64_t beta,        MORSE_Complex64_t *C, int LDC);
int  CORE_zsyssq(MORSE_enum uplo, int N,
                 const MORSE_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
int CORE_zsytf2_nopiv(MORSE_enum uplo, int n, MORSE_Complex64_t *A, int lda);
void CORE_zswpab(int i, int n1, int n2,
                 MORSE_Complex64_t *A, MORSE_Complex64_t *work);
int  CORE_zswptr_ontile(MORSE_desc_t descA, int i1, int i2, const int *ipiv, int inc,
                        const MORSE_Complex64_t *Akk, int ldak);
int CORE_ztradd(MORSE_enum uplo, MORSE_enum trans, int M, int N,
                      MORSE_Complex64_t alpha,
                const MORSE_Complex64_t *A, int LDA,
                      MORSE_Complex64_t beta,
                      MORSE_Complex64_t *B, int LDB);
void CORE_ztrasm(MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag,
                 int M, int N,
                 const MORSE_Complex64_t *A, int lda, double *work);
void CORE_ztrdalg(MORSE_enum uplo, int N, int NB,
                  const MORSE_desc_t *pA, MORSE_Complex64_t *V, MORSE_Complex64_t *TAU,
                  int i, int j, int m, int grsiz);
void CORE_ztrmm(MORSE_enum side, MORSE_enum uplo,
                MORSE_enum transA, MORSE_enum diag,
                int M, int N,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                                MORSE_Complex64_t *B, int LDB);
void CORE_ztrsm(MORSE_enum side, MORSE_enum uplo,
                MORSE_enum transA, MORSE_enum diag,
                int M, int N,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                                MORSE_Complex64_t *B, int LDB);
int CORE_ztrssq(MORSE_enum uplo, MORSE_enum diag, int M, int N,
                const MORSE_Complex64_t *A, int LDA,
                double *scale, double *sumsq);
void CORE_ztrtri(MORSE_enum uplo, MORSE_enum diag, int N,
                 MORSE_Complex64_t *A, int LDA, int *info);
int  CORE_ztslqt(int M, int N, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU, MORSE_Complex64_t *WORK);
int  CORE_ztsmlq(MORSE_enum side, MORSE_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK, int LDWORK);
int CORE_ztsmlq_hetra1( MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        MORSE_Complex64_t *A1, int lda1,
                        MORSE_Complex64_t *A2, int lda2,
                        const MORSE_Complex64_t *V, int ldv,
                        const MORSE_Complex64_t *T, int ldt,
                        MORSE_Complex64_t *WORK, int ldwork);
int  CORE_ztsmqr(MORSE_enum side, MORSE_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK, int LDWORK);
int CORE_ztsmqr_hetra1( MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        MORSE_Complex64_t *A1, int lda1,
                        MORSE_Complex64_t *A2, int lda2,
                        const MORSE_Complex64_t *V, int ldv,
                        const MORSE_Complex64_t *T, int ldt,
                        MORSE_Complex64_t *WORK, int ldwork);
int  CORE_ztsqrt(int M, int N, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU, MORSE_Complex64_t *WORK);
int  CORE_ztstrf(int M, int N, int IB, int NB,
                 MORSE_Complex64_t *U, int LDU,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *L, int LDL,
                 int *IPIV, MORSE_Complex64_t *WORK,
                 int LDWORK, int *INFO);
int CORE_ztplqt( int M, int N, int L, int IB,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *B, int LDB,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK );
int CORE_ztpmlqt( MORSE_enum side, MORSE_enum trans,
                  int M, int N, int K, int L, int IB,
                  const MORSE_Complex64_t *V, int LDV,
                  const MORSE_Complex64_t *T, int LDT,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB,
                  MORSE_Complex64_t *WORK );
int CORE_ztpmqrt( MORSE_enum side, MORSE_enum trans,
                  int M, int N, int K, int L, int IB,
                  const MORSE_Complex64_t *V, int LDV,
                  const MORSE_Complex64_t *T, int LDT,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB,
                  MORSE_Complex64_t *WORK );
int CORE_ztpqrt( int M, int N, int L, int IB,
                 MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *B, int LDB,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK );
int  CORE_zttmqr(MORSE_enum side, MORSE_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK, int LDWORK);
int  CORE_zttqrt(int M, int N, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU,
                 MORSE_Complex64_t *WORK);
int  CORE_zttmlq(MORSE_enum side, MORSE_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *WORK, int LDWORK);
int  CORE_zttlqt(int M, int N, int IB,
                 MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
                 MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *TAU,
                 MORSE_Complex64_t *WORK);
int  CORE_zunmlq(MORSE_enum side, MORSE_enum trans,
                 int M, int N, int IB, int K,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *C, int LDC,
                 MORSE_Complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(MORSE_enum side, MORSE_enum trans,
                 int M, int N, int K, int IB,
                 const MORSE_Complex64_t *V, int LDV,
                 const MORSE_Complex64_t *T, int LDT,
                 MORSE_Complex64_t *C, int LDC,
                 MORSE_Complex64_t *WORK, int LDWORK);
#endif
