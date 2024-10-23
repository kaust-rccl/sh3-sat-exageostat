/**
 *
 * @file zpotrs.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrs wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
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
 *  MORSE_zpotrs - Solves a system of linear equations A * X = B with a symmetric positive
 *  definite (or Hermitian positive definite in the complex case) matrix A using the Cholesky
 *  factorization A = U**H*U or A = L*L**H computed by MORSE_zpotrf.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by MORSE_zpotrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zpotrs_Tile
 * @sa MORSE_zpotrs_Tile_Async
 * @sa MORSE_cpotrs
 * @sa MORSE_dpotrs
 * @sa MORSE_spotrs
 * @sa MORSE_zpotrf
 *
 */
int MORSE_zpotrs( MORSE_enum uplo, int N, int NRHS,
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
        morse_fatal_error("MORSE_zpotrs", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_zpotrs", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_zpotrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zpotrs", "illegal value of NRHS");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_zpotrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, N)) {
        morse_error("MORSE_zpotrs", "illegal value of LDB");
        return -7;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZPOSV, N, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zpotrs", "morse_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    MORSE_zpotrs_Tile_Async( uplo, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, uplo, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInout, MorseUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zpotrs_Tile - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Tile equivalent of MORSE_zpotrs().
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
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U**H*U or A = L*L**H,
 *          computed by MORSE_zpotrf.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zpotrs
 * @sa MORSE_zpotrs_Tile_Async
 * @sa MORSE_cpotrs_Tile
 * @sa MORSE_dpotrs_Tile
 * @sa MORSE_spotrs_Tile
 * @sa MORSE_zpotrf_Tile
 *
 */
int MORSE_zpotrs_Tile( MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zpotrs_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zpotrs_Tile_Async( uplo, A, B, sequence, &request );

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
 *  MORSE_zpotrs_Tile_Async - Solves a system of linear equations using previously
 *  computed Cholesky factorization.
 *  Non-blocking equivalent of MORSE_zpotrs_Tile().
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
 * @sa MORSE_zpotrs
 * @sa MORSE_zpotrs_Tile
 * @sa MORSE_cpotrs_Tile_Async
 * @sa MORSE_dpotrs_Tile_Async
 * @sa MORSE_spotrs_Tile_Async
 * @sa MORSE_zpotrf_Tile_Async
 *
 */
int MORSE_zpotrs_Tile_Async( MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zpotrs_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zpotrs_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zpotrs_Tile", "NULL request");
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
        morse_error("MORSE_zpotrs_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zpotrs_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zpotrs_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_zpotrs_Tile", "illegal value of uplo");
        return morse_request_fail(sequence, request, -1);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return MORSE_SUCCESS;
     */
    morse_pztrsm( MorseLeft, uplo, uplo == MorseUpper ? MorseConjTrans : MorseNoTrans, MorseNonUnit, 1.0, A, B, sequence, request );

    morse_pztrsm( MorseLeft, uplo, uplo == MorseUpper ? MorseNoTrans : MorseConjTrans, MorseNonUnit, 1.0, A, B, sequence, request );

    return MORSE_SUCCESS;
}
