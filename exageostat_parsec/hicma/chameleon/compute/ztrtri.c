/**
 *
 * @file ztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri wrappers
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
 *  MORSE_ztrtri - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = MorseNonUnit: A is non-unit triangular;
 *          = MorseUnit:    A is unit triangular.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_ztrtri_Tile
 * @sa MORSE_ztrtri_Tile_Async
 * @sa MORSE_ctrtri
 * @sa MORSE_dtrtri
 * @sa MORSE_strtri
 * @sa MORSE_zpotri
 *
 */
int MORSE_ztrtri( MORSE_enum uplo, MORSE_enum diag, int N,
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
        morse_fatal_error("MORSE_ztrtri", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_ztrtri", "illegal value of uplo");
        return -1;
    }
    if ((diag != MorseUnit) && (diag != MorseNonUnit)) {
        morse_error("MORSE_ztrtri", "illegal value of diag");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_ztrtri", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_ztrtri", "illegal value of LDA");
        return -5;
    }
    /* Quick return */
    if (chameleon_max(N, 0) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZPOSV, N, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_ztrtri", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    MORSE_ztrtri_Tile_Async( uplo, diag, &descAt, sequence, &request );

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
 *  MORSE_ztrtri_Tile - Computes the inverse of a complex upper or
 *  lower triangular matrix A.
 *  Tile equivalent of MORSE_ztrtri().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = MorseNonUnit: A is non-unit triangular;
 *          = MorseUnit:    A us unit triangular.
 *
 * @param[in] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_ztrtri
 * @sa MORSE_ztrtri_Tile_Async
 * @sa MORSE_ctrtri_Tile
 * @sa MORSE_dtrtri_Tile
 * @sa MORSE_strtri_Tile
 * @sa MORSE_zpotri_Tile
 *
 */
int MORSE_ztrtri_Tile( MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrtri_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_ztrtri_Tile_Async( uplo, diag, A, sequence, &request );

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
 *  MORSE_ztrtri_Tile_Async - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *  Non-blocking equivalent of MORSE_ztrtri_Tile().
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
 * @sa MORSE_ztrtri
 * @sa MORSE_ztrtri_Tile
 * @sa MORSE_ctrtri_Tile_Async
 * @sa MORSE_dtrtri_Tile_Async
 * @sa MORSE_strtri_Tile_Async
 * @sa MORSE_zpotri_Tile_Async
 *
 */
int MORSE_ztrtri_Tile_Async( MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrtri_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztrtri_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztrtri_Tile", "NULL request");
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
        morse_error("MORSE_ztrtri_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_ztrtri_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_ztrtri_Tile", "illegal value of uplo");
        return morse_request_fail(sequence, request, -1);
    }
    if ((diag != MorseUnit) && (diag != MorseNonUnit)) {
        morse_error("MORSE_ztrtri_Tile", "illegal value of diag");
        return morse_request_fail(sequence, request, -2);
    }
    /* Quick return */
    /*
     if (chameleon_max(N, 0) == 0)
     return MORSE_SUCCESS;
     */
    morse_pztrtri( uplo, diag, A, sequence, request );

    return MORSE_SUCCESS;
}
