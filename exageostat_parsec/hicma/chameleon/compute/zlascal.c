/**
 *
 * @file zlascal.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal wrappers
 *
 * @version 1.0.0
 * @author Dalal Sukkari
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
 *  MORSE_zlascal - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = MorseUpperLower: A is a general matrix.
 *          = MorseUpper: A is an upper trapezoidal matrix.
 *          = MorseLower: A is a lower trapezoidal matrix.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in,out] A
 *          A is a LDA-by-N matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zlascal_Tile
 * @sa MORSE_clascal
 * @sa MORSE_dlascal
 * @sa MORSE_slascal
 *
 */
int MORSE_zlascal( MORSE_enum uplo, int M, int N,
                   MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    MORSE_desc_t descAl, descAt;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlascal", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower) && (uplo != MorseUpperLower)) {
        morse_error("MORSE_zlascal", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        morse_error("MORSE_zlascal", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_zlascal", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zlascal", "illegal value of LDA");
        return -6;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        (alpha == (MORSE_Complex64_t)1.0))
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zlascal", "morse_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zlascal_Tile_Async( uplo, alpha, &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInout, uplo, sequence, &request );

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
 *  MORSE_zlascal_Tile - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = MorseUpperLower: A is a general matrix.
 *          = MorseUpper: A is an upper trapezoidal matrix.
 *          = MorseLower: A is a lower trapezoidal matrix.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-N matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zlascal
 * @sa MORSE_zlascal_Tile_Async
 * @sa MORSE_clascal_Tile
 * @sa MORSE_dlascal_Tile
 * @sa MORSE_slascal_Tile
 *
 */
int MORSE_zlascal_Tile( MORSE_enum uplo,
                        MORSE_Complex64_t alpha, MORSE_desc_t *A )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlascal_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zlascal_Tile_Async( uplo, alpha, A, sequence, &request );

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
 *  MORSE_zlascal_Tile_Async - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *  Non-blocking equivalent of MORSE_zlascal_Tile().
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
 * @sa MORSE_zlascal
 * @sa MORSE_zlascal_Tile
 * @sa MORSE_clascal_Tile_Async
 * @sa MORSE_dlascal_Tile_Async
 * @sa MORSE_slascal_Tile_Async
 *
 */
int MORSE_zlascal_Tile_Async( MORSE_enum uplo,
                              MORSE_Complex64_t alpha, MORSE_desc_t *A,
                              MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t descA;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlascal_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zlascal_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zlascal_Tile_Async", "NULL request");
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
        morse_error("MORSE_zlascal_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower) && (uplo != MorseUpperLower)) {
        morse_error("MORSE_zlascal", "illegal value of uplo");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    if ( (descA.i%descA.mb != 0) || (descA.j%descA.nb != 0) ) {
        morse_error("MORSE_zlascal", "start indexes have to be multiple of tile size");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (descA.m == 0) || (descA.n == 0) ||
         (alpha == (MORSE_Complex64_t)1.0) )
    {
        return MORSE_SUCCESS;
    }

    morse_pzlascal( uplo, alpha, A, sequence, request );

    return MORSE_SUCCESS;
}
