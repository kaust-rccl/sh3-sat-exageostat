/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @generated d Mon May 11 22:42:39 2020
 *
 */

#include "dplasma.h"
//#include "dplasma/types.h"
//#include "mix_precision_internal.h"
#include "include/dplasmaaux.h"

#include "my_dtrsm_LLN.h"
//#include "my_dtrsm_LLT.h"
/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 *  my_dtrsm_New - Generates parsec taskpool to compute triangular solve
 *     op( A ) * X = B or X * op( A ) = B
 *  WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  op( A ) * X = B
 *          = PlasmaRight: X * op( A ) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          ugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          If uplo = dplasmaUpper, the leading N-by-N upper triangular part of
 *          the array A contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced. If uplo = PlasmaLower,
 *          the leading N-by-N lower triangular part of the array A contains the
 *          lower triangular matrix, and the strictly upper triangular part of A
 *          is not referenced. If diag = PlasmaUnit, the diagonal elements of A
 *          are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          Descriptor of the N-by-NRHS right hand side B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_my_dtrsm_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_my_dtrsm
 * @sa dplasma_my_dtrsm_Destruct
 * @sa dplasma_ctrsm_New
 * @sa my_dtrsm_New
 * @sa dplasma_strsm_New
 *
 ******************************************************************************/
parsec_taskpool_t*
my_dtrsm_New( PLASMA_enum side,  PLASMA_enum uplo,
                   PLASMA_enum trans, PLASMA_enum diag,
                   double alpha,
                   const parsec_tiled_matrix_dc_t *A,
                   parsec_tiled_matrix_dc_t *B )
{
    parsec_taskpool_t *parsec_trsm = NULL;

    if ( side == PlasmaLeft ) {
        if ( uplo == PlasmaLower ) {
            if ( trans == PlasmaNoTrans ) {
                parsec_trsm = (parsec_taskpool_t*)parsec_my_dtrsm_LLN_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            } else { /* trans =! PlasmaNoTrans */
	    }
        }
    }
	int rc;
    rc = parsec_matrix_add2arena(((parsec_my_dtrsm_LLN_taskpool_t*)parsec_trsm)->arenas[PARSEC_my_dtrsm_LLN_FULL_ARENA],
			                parsec_datatype_double_t, matrix_UpperLower,
							1, A->mb, A->nb, A->mb,
							PARSEC_ARENA_ALIGNMENT_SSE, -1 );
    parsec_matrix_add2arena(((parsec_my_dtrsm_LLN_taskpool_t*)parsec_trsm)->arenas[PARSEC_my_dtrsm_LLN_VECTOR_ARENA],
			                parsec_datatype_double_t, matrix_UpperLower,
							1, A->mb, 1, A->mb,
							PARSEC_ARENA_ALIGNMENT_SSE, -1 );
    return parsec_trsm;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 *  dplasma_my_dtrsm_Destruct - Free the data structure associated to an taskpool
 *  created with my_dtrsm_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa my_dtrsm_New
 * @sa dplasma_my_dtrsm
 *
 ******************************************************************************/
void
my_dtrsm_Destruct( parsec_taskpool_t *tp )
{
    parsec_my_dtrsm_LLN_taskpool_t *otrsm = (parsec_my_dtrsm_LLN_taskpool_t *)tp;

    parsec_matrix_del2arena( otrsm->arenas[PARSEC_my_dtrsm_LLN_FULL_ARENA] );
    parsec_matrix_del2arena( otrsm->arenas[PARSEC_my_dtrsm_LLN_VECTOR_ARENA] );
    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 *  dplasma_my_dtrsm - Computes triangular solve
 *     op( A ) * X = B or X * op( A ) = B
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  op( A ) * X = B
 *          = PlasmaRight: X * op( A ) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          ugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          If uplo = dplasmaUpper, the leading N-by-N upper triangular part of
 *          the array A contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced. If uplo = PlasmaLower,
 *          the leading N-by-N lower triangular part of the array A contains the
 *          lower triangular matrix, and the strictly upper triangular part of A
 *          is not referenced. If diag = PlasmaUnit, the diagonal elements of A
 *          are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          Descriptor of the N-by-NRHS right hand side B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa my_dtrsm_New
 * @sa dplasma_my_dtrsm_Destruct
 * @sa dplasma_ctrsm
 * @sa dplasma_my_dtrsm
 * @sa dplasma_strsm
 *
 ******************************************************************************/
int
my_dtrsm( parsec_context_t *parsec,
               PLASMA_enum side,  PLASMA_enum uplo,
               PLASMA_enum trans, PLASMA_enum diag,
               double alpha,
               const parsec_tiled_matrix_dc_t *A,
               parsec_tiled_matrix_dc_t *B)
{
    parsec_taskpool_t *parsec_my_dtrsm = NULL;

    /* Check input arguments */
    if (side != PlasmaLeft ) {
        dplasma_error("dplasma_my_dtrsm", "illegal value of side");
        return -1;
    }
    if (uplo != PlasmaLower) {
        dplasma_error("dplasma_my_dtrsm", "illegal value of uplo");
        return -2;
    }
    if (trans != PlasmaTrans && trans != PlasmaNoTrans && trans != PlasmaTrans ) {
        dplasma_error("dplasma_my_dtrsm", "illegal value of trans");
        return -3;
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        dplasma_error("dplasma_my_dtrsm", "illegal value of diag");
        return -4;
    }

    if ( (A->m != A->n) ||
         (( side == PlasmaLeft )  && (A->n != B->m)) ||
         (( side == PlasmaRight ) && (A->n != B->n)) ) {
        dplasma_error("dplasma_my_dtrsm", "illegal matrix A");
        return -6;
    }

    parsec_my_dtrsm = my_dtrsm_New(side, uplo, trans, diag, alpha, A, B);

    if ( parsec_my_dtrsm != NULL )
    {
        parsec_context_add_taskpool( parsec, parsec_my_dtrsm );
        parsec_context_start( parsec );
        parsec_context_wait( parsec );
        my_dtrsm_Destruct( parsec_my_dtrsm );
        return 0;
    }
    else {
        return -101;
    }
}
