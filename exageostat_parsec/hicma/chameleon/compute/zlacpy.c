/**
 *
 * @file zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy wrappers
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
 *  MORSE_zlacpy copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = MorseUpperLower: All the matrix A
 *            = MorseUpper: Upper triangular part
 *            = MorseLower: Lower triangular part
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = MorseUpper, only the upper trapezium
 *          is accessed; if UPLO = MorseLower, only the lower trapezium is
 *          accessed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa MORSE_zlacpy_Tile
 * @sa MORSE_zlacpy_Tile_Async
 * @sa MORSE_clacpy
 * @sa MORSE_dlacpy
 * @sa MORSE_slacpy
 *
 */
int MORSE_zlacpy( MORSE_enum uplo, int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descBl, descBt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlacpy", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != MorseUpperLower) &&
         (uplo != MorseUpper) &&
         (uplo != MorseLower) ) {
        morse_error("MORSE_zlacpy", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        morse_error("MORSE_zlacpy", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_zlacpy", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zlacpy", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, M)) {
        morse_error("MORSE_zlacpy", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zlacpy", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, uplo,
                     B, NB, NB, LDB, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zlacpy_Tile_Async( uplo, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, uplo, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInout, uplo, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );

    morse_sequence_destroy( morse, sequence );
    return MORSE_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zlacpy_Tile - Tile equivalent of MORSE_zlacpy().
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
 * @param[in] A
 *          The M-by-N matrix A. If uplo = MorseUpper, only the upper trapezium
 *          is accessed; if UPLO = MorseLower, only the lower trapezium is
 *          accessed.
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zlacpy
 * @sa MORSE_zlacpy_Tile_Async
 * @sa MORSE_clacpy_Tile
 * @sa MORSE_dlacpy_Tile
 * @sa MORSE_slacpy_Tile
 *
 */
int MORSE_zlacpy_Tile( MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlacpy_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zlacpy_Tile_Async( uplo, A, B, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( B, sequence );

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
 *  MORSE_zlacpy_Tile_Async - Non-blocking equivalent of MORSE_zlacpy_Tile().
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
 * @sa MORSE_zlacpy
 * @sa MORSE_zlacpy_Tile
 * @sa MORSE_clacpy_Tile_Async
 * @sa MORSE_dlacpy_Tile_Async
 * @sa MORSE_slacpy_Tile_Async
 *
 */
int MORSE_zlacpy_Tile_Async( MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlacpy_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zlacpy_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zlacpy_Tile_Async", "NULL request");
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
        morse_error("MORSE_zlacpy_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zlacpy_Tile_Async", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zlacpy_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != MorseUpperLower) &&
         (uplo != MorseUpper) &&
         (uplo != MorseLower) ) {
        morse_error("MORSE_zlacpy_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return MORSE_SUCCESS;
    }

    morse_pzlacpy( uplo, A, B, sequence, request );

    return MORSE_SUCCESS;
}
