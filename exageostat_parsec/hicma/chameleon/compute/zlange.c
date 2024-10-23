/**
 *
 * @file zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
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
 *  MORSE_zlange returns the value
 *
 *     zlange = ( max(abs(A(i,j))), NORM = MorseMaxNorm
 *              (
 *              ( norm1(A),         NORM = MorseOneNorm
 *              (
 *              ( normI(A),         NORM = MorseInfNorm
 *              (
 *              ( normF(A),         NORM = MorseFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = MorseMaxNorm: Max norm
 *          = MorseOneNorm: One norm
 *          = MorseInfNorm: Infinity norm
 *          = MorseFrobeniusNorm: Frobenius norm
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0. When M = 0,
 *          the returned value is set to zero.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval the norm described above.
 *
 *******************************************************************************
 *
 * @sa MORSE_zlange_Tile
 * @sa MORSE_zlange_Tile_Async
 * @sa MORSE_clange
 * @sa MORSE_dlange
 * @sa MORSE_slange
 *
 */
double MORSE_zlange(MORSE_enum norm, int M, int N,
                    MORSE_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    double value = -1.;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlange", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (norm != MorseMaxNorm) && (norm != MorseOneNorm)
         && (norm != MorseInfNorm) && (norm != MorseFrobeniusNorm) ) {
        morse_error("MORSE_zlange", "illegal value of norm");
        return -1;
    }
    if (M < 0) {
        morse_error("MORSE_zlange", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_zlange", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zlange", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zlange", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zlange_Tile_Async( norm, &descAt, &value, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );

    morse_sequence_destroy( morse, sequence );
    return value;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zlange_Tile - Tile equivalent of MORSE_zlange().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = MorseMaxNorm: Max norm
 *          = MorseOneNorm: One norm
 *          = MorseInfNorm: Infinity norm
 *          = MorseFrobeniusNorm: Frobenius norm
 *
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zlange
 * @sa MORSE_zlange_Tile_Async
 * @sa MORSE_clange_Tile
 * @sa MORSE_dlange_Tile
 * @sa MORSE_slange_Tile
 *
 */
double MORSE_zlange_Tile(MORSE_enum norm, MORSE_desc_t *A )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;
    double value = -1.;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlange_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zlange_Tile_Async( norm, A, &value, sequence, &request );

    MORSE_Desc_Flush( A, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return ( status == MORSE_SUCCESS ) ? value : (double)status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zlange_Tile_Async - Non-blocking equivalent of MORSE_zlange_Tile().
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
 * @sa MORSE_zlange
 * @sa MORSE_zlange_Tile
 * @sa MORSE_clange_Tile_Async
 * @sa MORSE_dlange_Tile_Async
 * @sa MORSE_slange_Tile_Async
 *
 */
int MORSE_zlange_Tile_Async( MORSE_enum norm, MORSE_desc_t *A, double *value,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zlange_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zlange_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zlange_Tile", "NULL request");
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
        morse_error("MORSE_zlange_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zlange_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ( (norm != MorseMaxNorm) && (norm != MorseOneNorm)
         && (norm != MorseInfNorm) && (norm != MorseFrobeniusNorm) ) {
        morse_error("MORSE_zlange", "illegal value of norm");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        *value = 0.0;
        return MORSE_SUCCESS;
    }

    morse_pzlange( norm, A, value, sequence, request );

    return MORSE_SUCCESS;
}
