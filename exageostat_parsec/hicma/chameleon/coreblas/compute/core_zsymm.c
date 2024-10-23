/**
 *
 * @file core_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zsymm CPU kernel
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
void CORE_zsymm(MORSE_enum side, MORSE_enum uplo,
                int M, int N,
                MORSE_Complex64_t alpha, const MORSE_Complex64_t *A, int LDA,
                const MORSE_Complex64_t *B, int LDB,
                MORSE_Complex64_t beta, MORSE_Complex64_t *C, int LDC)
{
    cblas_zsymm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}


