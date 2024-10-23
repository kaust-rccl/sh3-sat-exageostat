/**
 *
 * @file zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon zbuild wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zbuild - Generate a matrix by calling user provided function.
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
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] user_data
 *          The user data used in the matrix generation, it will be passed by chameleon
 *          to the build callback function (see below).
 *
 * @param[in] user_build_callback
 *          The user function to call to generate tiles.
 *          The prototype of the callback is :
 *          void myFcn(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data)
 *          It is expected to build the block of matrix [row_min, row_max] x [col_min, col_max]
 *          (with both min and max values included in the intervals,
 *          index start at 0 like in C, NOT 1 like in Fortran)
 *          and store it at the adresse 'buffer' with leading dimension 'ld'
 *          The argument 'user_data' is an opaque pointer on any user data, it is passed by
 *          the user to Morse_zbuild (see above) and transmitted by chameleon to the callback.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zbuild_Tile
 * @sa MORSE_zbuild_Tile_Async
 * @sa MORSE_cbuild
 * @sa MORSE_dbuild
 * @sa MORSE_sbuild
 *
 */
int MORSE_zbuild( MORSE_enum uplo, int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  void *user_data, void* user_build_callback )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zbuild", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zbuild", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_zbuild", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zbuild", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zbuild", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;
    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescOutput, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zbuild_Tile_Async( uplo, &descAt, user_data, user_build_callback, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescOutput, uplo, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );

    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zbuild_Tile - Generate a matrix by tiles by calling user provided function.
 *  Tile equivalent of MORSE_zbuild().
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
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles (see MORSE_zbuild)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zbuild
 * @sa MORSE_zbuild_Tile_Async
 * @sa MORSE_cbuild_Tile
 * @sa MORSE_dbuild_Tile
 * @sa MORSE_sbuild_Tile
 *
 */
int MORSE_zbuild_Tile( MORSE_enum uplo, MORSE_desc_t *A,
                       void *user_data, void* user_build_callback )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zbuild_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zbuild_Tile_Async( uplo, A, user_data, user_build_callback, sequence, &request );

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
 *  MORSE_zbuild_Tile_Async - Generate a matrix by tiles by calling user provided function.
 *  Non-blocking equivalent of MORSE_zbuild_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
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
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles (see MORSE_zbuild)
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
 * @sa MORSE_zbuild
 * @sa MORSE_zbuild_Tile
 * @sa MORSE_cbuild_Tile_Async
 * @sa MORSE_dbuild_Tile_Async
 * @sa MORSE_sbuild_Tile_Async
 *
 */
int MORSE_zbuild_Tile_Async( MORSE_enum uplo, MORSE_desc_t     *A,
                             void *user_data, void* user_build_callback,
                             MORSE_sequence_t *sequence,
                             MORSE_request_t  *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zbuild_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zbuild_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zbuild_Tile", "NULL request");
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
        morse_error("MORSE_zbuild_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0) {
        return MORSE_SUCCESS;
    }

    morse_pzbuild( uplo, A, user_data, user_build_callback, sequence, request );

    return MORSE_SUCCESS;
}
