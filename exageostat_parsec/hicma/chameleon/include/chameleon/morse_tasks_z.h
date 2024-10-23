/**
 *
 * @file morse_tasks_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon MORSE_Complex64_t elementary tasks header
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
#ifndef _MORSE_TASKS_Z_H_
#define _MORSE_TASKS_Z_H_


/**
 *  Declarations of QUARK wrappers (called by MORSE) - alphabetical order
 */
void MORSE_TASK_dzasum(const MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, int M, int N,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zgeadd(const MORSE_option_t *options,
                       MORSE_enum trans, int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_Complex64_t beta,  const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlascal(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int m, int n, int nb,
                        MORSE_Complex64_t alpha,
                        const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_zbrdalg(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int N, int NB,
                        const MORSE_desc_t *A,
                        const MORSE_desc_t *C, int Cm, int Cn,
                        const MORSE_desc_t *S, int Sm, int Sn,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void MORSE_TASK_zgelqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zgemm(const MORSE_option_t *options,
                      MORSE_enum transA, MORSE_enum transB,
                      int m, int n, int k, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm2( const MORSE_option_t *options,
                        MORSE_enum transA, MORSE_enum transB,
                        int m, int n, int k, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                        MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm_f2(const MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc,
                         const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                         const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zgemm_p2(const MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_Complex64_t **B, int ldb,
                         MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm_p2f1(const MORSE_option_t *options,
                           MORSE_enum transA, MORSE_enum transB,
                           int m, int n, int k, int nb,
                           MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                           const MORSE_Complex64_t **B, int ldb,
                           MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc,
                           const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1);
void MORSE_TASK_zgemm_p3(const MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         MORSE_Complex64_t beta, MORSE_Complex64_t **C, int ldc);
void MORSE_TASK_zgeqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zgessm(const MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       const MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_zgessq( const MORSE_option_t *options,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn);
void MORSE_TASK_zgetrf(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,

                       MORSE_bool check_info, int iinfo);
void MORSE_TASK_zgetrf_incpiv(const MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo);
void MORSE_TASK_zgetrf_nopiv(const MORSE_option_t *options,
                             int m, int n, int ib, int nb,
                             const MORSE_desc_t *A, int Am, int An, int lda, int iinfo);
void MORSE_TASK_zgetrf_reclap(const MORSE_option_t *options,
                              int m, int n, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              int *IPIV,

                              MORSE_bool check_info, int iinfo,
                              int nbthread);
void MORSE_TASK_zgetrf_rectil(const MORSE_option_t *options,
                              const MORSE_desc_t A, const MORSE_desc_t *Amn, int Amnm, int Amnn, int size,
                              int *IPIV,

                              MORSE_bool check_info, int iinfo,
                              int nbthread);
void MORSE_TASK_zgetrip(const MORSE_option_t *options,
                        int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA);
void MORSE_TASK_zgetrip_f1(const MORSE_option_t *options,
                           int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA,
                           const MORSE_desc_t *fake, int fakem, int faken, int szeF, int paramF);
void MORSE_TASK_zgetrip_f2(const MORSE_option_t *options,
                           int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA,
                           const MORSE_desc_t *fake1, int fake1m, int fake1n, int szeF1, int paramF1,
                           const MORSE_desc_t *fake2, int fake2m, int fake2n, int szeF2, int paramF2);
void MORSE_TASK_zhe2ge(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zhemm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zhegst(const MORSE_option_t *options,
                       int itype, MORSE_enum uplo, int N,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                       int iinfo);
void MORSE_TASK_zherk(const MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      double alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zher2k(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                       double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zherfb(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zlacpy(const MORSE_option_t *options,
                       MORSE_enum uplo, int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlacpyx(const MORSE_option_t *options,
                        MORSE_enum uplo, int m, int n, int mb,
                        int displA, const MORSE_desc_t *A, int Am, int An, int lda,
                        int displB, const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlange(const MORSE_option_t *options,
                       MORSE_enum norm, int M, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlange_max(const MORSE_option_t *options,
                           const MORSE_desc_t *A, int Am, int An,
                           const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zhessq( const MORSE_option_t *options,
                        MORSE_enum uplo, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_zlanhe(const MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlansy(const MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlantr(const MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag,
                       int M, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlaset(const MORSE_option_t *options,
                       MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                       MORSE_Complex64_t beta, const MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea);
void MORSE_TASK_zlaset2(const MORSE_option_t *options,
                        MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                        const MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea);
void MORSE_TASK_zlaswp(const MORSE_option_t *options,
                       int n, const MORSE_desc_t *A, int Am, int An, int lda,
                       int i1,  int i2, int *ipiv, int inc);
void MORSE_TASK_zlaswp_f2(const MORSE_option_t *options,
                          int n, const MORSE_desc_t *A, int Am, int An, int lda,
                          int i1,  int i2, int *ipiv, int inc,
                          const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                          const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zlaswp_ontile(const MORSE_option_t *options,
                              const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                              int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel);
void MORSE_TASK_zlaswp_ontile_f2(const MORSE_option_t *options,
                                 const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                                 int i1,  int i2, int *ipiv, int inc,
                                 const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                                 const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zlaswpc_ontile(const MORSE_option_t *options,
                               const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                               int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel);
void MORSE_TASK_zlatro(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans, int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlauum(const MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_zplghe(const MORSE_option_t *options,
                       double bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zplgsy(const MORSE_option_t *options,
                       MORSE_Complex64_t bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zplrnt(const MORSE_option_t *options,
                       int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zpotrf(const MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,

                       int iinfo);
void MORSE_TASK_zshift( const MORSE_option_t *options,
                        int s, int m, int n, int L,
                        MORSE_Complex64_t *A);
void MORSE_TASK_zshiftw(const MORSE_option_t *options,
                        int s, int cl, int m, int n, int L,
                        const MORSE_desc_t *A, int Am, int An, MORSE_Complex64_t *W);
void MORSE_TASK_zssssm(const MORSE_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                       const MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV);
void MORSE_TASK_zsymm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyrk(const MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyr2k(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyssq( const MORSE_option_t *options,
                        MORSE_enum uplo, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_zsytrf_nopiv(const MORSE_option_t *options,
                             MORSE_enum uplo, int n, int nb,
                             const MORSE_desc_t *A, int Am, int An, int lda,
                             int iinfo);
void MORSE_TASK_zswpab(const MORSE_option_t *options,
                       int i, int n1, int n2,
                       const MORSE_desc_t *A, int Am, int An, int szeA);
void MORSE_TASK_zswptr_ontile(const MORSE_option_t *options,
                              const MORSE_desc_t descA, const MORSE_desc_t *Aij, int Aijm, int Aijn,
                              int i1,  int i2, int *ipiv, int inc,
                              const MORSE_desc_t *Akk, int Akkm, int Akkn, int ldak);
void MORSE_TASK_ztplqt(const MORSE_option_t *options,
                       int m, int n, int l, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void MORSE_TASK_ztpmlqt( const MORSE_option_t *options,
                         MORSE_enum side, MORSE_enum trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void MORSE_TASK_ztpmqrt(const MORSE_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m, int n, int k, int l, int ib, int nb,
                        const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void MORSE_TASK_ztpqrt(const MORSE_option_t *options,
                       int m, int n, int l, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void MORSE_TASK_ztrdalg(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int N, int NB,
                        const MORSE_desc_t *A,
                        const MORSE_desc_t *C, int Cm, int Cn,
                        const MORSE_desc_t *S, int Sm, int Sn,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void MORSE_TASK_ztradd(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans, int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_Complex64_t beta,  const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_ztrasm(const MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_ztrmm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_ztrmm_p2(const MORSE_option_t *options,
                         MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                         int m, int n, int nb,
                         MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                         MORSE_Complex64_t **B, int ldb);
void MORSE_TASK_ztrsm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_ztrssq( const MORSE_option_t *options,
                        MORSE_enum uplo, MORSE_enum diag,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_ztrtri(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,

                       int iinfo);
void MORSE_TASK_ztslqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmlq(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmlq_hetra1(const MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmqr_hetra1(const MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztstrf(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *U, int Um, int Un, int ldu,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       int *IPIV,

                       MORSE_bool check_info, int iinfo);
void MORSE_TASK_zttmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttmlq(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttlqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zpamm(const MORSE_option_t *options,
                       int op, MORSE_enum side, MORSE_enum storev,
                       int m, int n, int k, int l,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *W, int Wm, int Wn, int ldw);
void MORSE_TASK_zplssq( const MORSE_option_t *options,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                        const MORSE_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void MORSE_TASK_zplssq2( const MORSE_option_t *options,
                         const MORSE_desc_t *RESULT, int RESULTm, int RESULTn );
void MORSE_TASK_zunmlq(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zunmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zbuild( const MORSE_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback );

#endif
