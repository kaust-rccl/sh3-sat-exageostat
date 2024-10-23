/**
 *
 * @file zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zlaset copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = MorseUpperLower: All the matrix A
 *            = MorseUpper: Upper triangular part is set. The lower
 *            triangle is unchanged.
 *            = MorseLower: Lower triangular part is set. The upper
 *            triangle is unchange.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          All the offdiagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa MORSE_zlaset_Tile
 * @sa MORSE_zlaset_Tile_Async
 * @sa MORSE_claset
 * @sa MORSE_dlaset
 * @sa MORSE_slaset
 *
 */
int MORSE_zlaset( MORSE_enum uplo, int M, int N,
                  MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                  MORSE_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlaset", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != MorseUpperLower) &&
         (uplo != MorseUpper) &&
         (uplo != MorseLower) ) {
        morse_error("MORSE_zlaset", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        morse_error("MORSE_zlaset", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_zlaset", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zlaset", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zlaset", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zlaset_Tile_Async( uplo, alpha, beta, &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInout, uplo, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );

    morse_sequence_destroy( morse, sequence );
    return MORSE_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zlaset_Tile - Tile equivalent of MORSE_zlaset().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = MorseUpperLower: All the matrix A
 *            = MorseUpper: Upper triangular part
 *            = MorseLower: Lower triangular part
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zlaset
 * @sa MORSE_zlaset_Tile_Async
 * @sa MORSE_claset_Tile
 * @sa MORSE_dlaset_Tile
 * @sa MORSE_slaset_Tile
 *
 */
int MORSE_zlaset_Tile( MORSE_enum uplo,
                       MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                       MORSE_desc_t *A )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlaset_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zlaset_Tile_Async( uplo, alpha, beta, A, sequence, &request );

    MORSE_Desc_Flush( A, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zlaset_Tile_Async - Non-blocking equivalent of MORSE_zlaset_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MORSE_zlaset
 * @sa MORSE_zlaset_Tile
 * @sa MORSE_claset_Tile_Async
 * @sa MORSE_dlaset_Tile_Async
 * @sa MORSE_slaset_Tile_Async
 *
 */
int MORSE_zlaset_Tile_Async( MORSE_enum uplo,
                             MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                             MORSE_desc_t *A,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlaset_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zlaset_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zlaset_Tile_Async", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS) {
        request->status = MORSE_SUCCESS;
    }
    else {
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zlaset_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zlaset_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != MorseUpperLower) &&
         (uplo != MorseUpper) &&
         (uplo != MorseLower) ) {
        morse_error("MORSE_zlaset_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return MORSE_SUCCESS;
    }

    morse_pzlaset( uplo, alpha, beta, A, sequence, request );

    return MORSE_SUCCESS;
}
