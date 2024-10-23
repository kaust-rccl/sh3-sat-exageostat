/**
 *
 * @file core_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zherk CPU kernel
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
 * @precisions normal z -> c
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void CORE_zherk(MORSE_enum uplo, MORSE_enum trans,
                int N, int K,
                double alpha, const MORSE_Complex64_t *A, int LDA,
                double beta, MORSE_Complex64_t *C, int LDC)
{
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        alpha, A, LDA,
        beta, C, LDC);
}
