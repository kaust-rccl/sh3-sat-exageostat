/**
 *
 * @file core_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlag2c CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void CORE_zlag2c(int m, int n,
                 const MORSE_Complex64_t *A, int lda,
                 MORSE_Complex32_t *B, int ldb, int *info)
{
    *info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}
/**/

/**/

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void CORE_clag2z(int m, int n,
                 const MORSE_Complex32_t *A, int lda,
                 MORSE_Complex64_t *B, int ldb)
{
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}



