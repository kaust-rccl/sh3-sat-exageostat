/**
 *
 * @file core_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgetrf CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"


int CORE_zgetrf(int m, int n,
                 MORSE_Complex64_t *A, int lda,
                 int *IPIV, int *info)
{
    *info = LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    return MORSE_SUCCESS;
}


