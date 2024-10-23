/**
 *
 * @file core_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgemm CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void CORE_zgemm(MORSE_enum transA, int transB,
                int M, int N, int K,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                                         const MORSE_Complex64_t *B, int LDB,
                MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC)
{
#ifdef CBLAS_HAS_ZGEMM3M
void cblas_zgemm3m(const enum CBLAS_ORDER Layout, const CBLAS_TRANSPOSE TransA,
                   const CBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha, const void *A,
                   const int lda, const void *B, const int ldb,
                   const void *beta, void *C, const int ldc);
  if (get_coreblas_gemm3m_enabled())
    cblas_zgemm3m(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
  else
#endif
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}


/**/








