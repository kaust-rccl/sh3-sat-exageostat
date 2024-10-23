/**
 *
 * @file zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe wrappers
 *
 * @version 1.0.0
 * @comment This file is a copy of zplghe.c
 *          wich has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Rade Mathis
 * @author Florent Pruvost
 * @date 2016-08-01
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zplghe - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] uplo
 *          The part of the Matrix wich will be generated.
 *
 * @param[out] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa MORSE_zplghe_Tile
 * @sa MORSE_zplghe_Tile_Async
 * @sa MORSE_cplghe
 * @sa MORSE_dplghe
 * @sa MORSE_splghe
 * @sa MORSE_zplghe
 *
 */
int MORSE_zplghe( double bump, MORSE_enum uplo, int N,
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
        morse_fatal_error("MORSE_zplghe", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        morse_error("MORSE_zplghe", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_zplghe", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_max(0, N) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGEMM, N, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zplghe", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;
    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescOutput, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zplghe_Tile_Async( bump, uplo, &descAt, seed, sequence, &request );

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
 *  MORSE_zplghe_Tile - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 *  Tile equivalent of MORSE_zplghe().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix.
 *
 * @param[in] uplo
 *          The part of the Matrix wich will be generated.
 *
 * @param[in] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa MORSE_zplghe
 * @sa MORSE_zplghe_Tile_Async
 * @sa MORSE_cplghe_Tile
 * @sa MORSE_dplghe_Tile
 * @sa MORSE_splghe_Tile
 * @sa MORSE_zplghe_Tile
 *
 */
int MORSE_zplghe_Tile( double bump, MORSE_enum uplo, MORSE_desc_t *A,
                       unsigned long long int seed )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zplghe_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zplghe_Tile_Async( bump, uplo, A, seed, sequence, &request );

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
 *  MORSE_zplghe_Tile_Async - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 *  Non-blocking equivalent of MORSE_zplghe_Tile().
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
 * @sa MORSE_zplghe
 * @sa MORSE_zplghe_Tile
 * @sa MORSE_cplghe_Tile_Async
 * @sa MORSE_dplghe_Tile_Async
 * @sa MORSE_splghe_Tile_Async
 * @sa MORSE_zplghe_Tile_Async
 *
 */
int MORSE_zplghe_Tile_Async( double                 bump,
                             MORSE_enum             uplo,
                             MORSE_desc_t             *A,
                             unsigned long long int seed,
                             MORSE_sequence_t  *sequence,
                             MORSE_request_t    *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zplghe_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zplghe_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zplghe_Tile", "NULL request");
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
        morse_error("MORSE_zplghe_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zplghe_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return MORSE_SUCCESS;

    morse_pzplghe( bump, uplo, A, seed, sequence, request );

    return MORSE_SUCCESS;
}
