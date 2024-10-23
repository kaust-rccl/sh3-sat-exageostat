/**
 *
 * @file zplrnt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt wrappers
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
 *  MORSE_zplrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zplrnt_Tile
 * @sa MORSE_zplrnt_Tile_Async
 * @sa MORSE_cplrnt
 * @sa MORSE_dplrnt
 * @sa MORSE_splrnt
 * @sa MORSE_zplghe
 * @sa MORSE_zplgsy
 *
 */
int MORSE_zplrnt( int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  unsigned long long int seed )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zplrnt", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zplrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_zplrnt", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zplrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zplrnt", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;
    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescOutput, MorseUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zplrnt_Tile_Async( &descAt, seed, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescOutput, MorseUpperLower, sequence, &request );

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
 *  MORSE_zplrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of MORSE_zplrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zplrnt
 * @sa MORSE_zplrnt_Tile_Async
 * @sa MORSE_cplrnt_Tile
 * @sa MORSE_dplrnt_Tile
 * @sa MORSE_splrnt_Tile
 * @sa MORSE_zplghe_Tile
 * @sa MORSE_zplgsy_Tile
 *
 */
int MORSE_zplrnt_Tile( MORSE_desc_t *A,
                       unsigned long long int seed )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zplrnt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zplrnt_Tile_Async( A, seed, sequence, &request );

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
 *  MORSE_zplrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of MORSE_zplrnt_Tile().
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
 * @sa MORSE_zplrnt
 * @sa MORSE_zplrnt_Tile
 * @sa MORSE_cplrnt_Tile_Async
 * @sa MORSE_dplrnt_Tile_Async
 * @sa MORSE_splrnt_Tile_Async
 * @sa MORSE_zplghe_Tile_Async
 * @sa MORSE_zplgsy_Tile_Async
 *
 */
int MORSE_zplrnt_Tile_Async( MORSE_desc_t     *A,
                             unsigned long long int seed,
                             MORSE_sequence_t *sequence,
                             MORSE_request_t  *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zplrnt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zplrnt_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zplrnt_Tile", "NULL request");
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
        morse_error("MORSE_zplrnt_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zplrnt_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return MORSE_SUCCESS;

    morse_pzplrnt( A, seed, sequence, request );

    return MORSE_SUCCESS;
}
