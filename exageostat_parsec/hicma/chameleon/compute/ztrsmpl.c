/**
 *
 * @file ztrsmpl.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsmpl wrappers
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
 *  MORSE_ztrsmpl - Performs the forward substitution step of solving a system of linear equations
 *  after the tile LU factorization of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] descL
 *          Auxiliary factorization data, related to the tile L factor, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from MORSE_zgetrf_incpiv (not equivalent to LAPACK).
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
 * @sa MORSE_ztrsmpl_Tile
 * @sa MORSE_ztrsmpl_Tile_Async
 * @sa MORSE_ctrsmpl
 * @sa MORSE_dtrsmpl
 * @sa MORSE_strsmpl
 * @sa MORSE_zgetrf_incpiv
 *
 */
int MORSE_ztrsmpl( int N, int NRHS,
                   MORSE_Complex64_t *A, int LDA,
                   MORSE_desc_t *descL, int *IPIV,
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
        morse_fatal_error("MORSE_ztrsmpl", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        morse_error("MORSE_ztrsmpl", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        morse_error("MORSE_ztrsmpl", "illegal value of NRHS");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_ztrsmpl", "illegal value of LDA");
        return -4;
    }
    if (LDB < chameleon_max(1, N)) {
        morse_error("MORSE_ztrsmpl", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGESV, N, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsmpl", "morse_tune() failed");
        return status;
    }

    /* Set Mt, NT & NTRHS */
    NB    = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseLower,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    MORSE_ztrsmpl_Tile_Async( &descAt, descL, IPIV, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseLower, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descL, sequence );

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
 * MORSE_ztrsmpl_Tile - Performs the forward substitution step of solving a system of linear equations
 * after the tile LU factorization of the matrix.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from MORSE_zgetrf_incpiv (not equivalent to LAPACK).
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
 * @sa MORSE_ztrsmpl
 * @sa MORSE_ztrsmpl_Tile_Async
 * @sa MORSE_ctrsmpl_Tile
 * @sa MORSE_dtrsmpl_Tile
 * @sa MORSE_strsmpl_Tile
 * @sa MORSE_zgetrf_incpiv_Tile
 *
 */
int MORSE_ztrsmpl_Tile( MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsmpl_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_ztrsmpl_Tile_Async( A, L, IPIV, B, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( L, sequence );
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
 *  MORSE_ztrsmpl_Tile - Performs the forward substitution step of solving
 *  a system of linear equations after the tile LU factorization of the matrix.
 *  Non-blocking equivalent of MORSE_ztrsmpl_Tile().
 *  Returns control to the user thread before worker threads finish the computation
 *  to allow for pipelined execution of diferent routines.
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
 * @sa MORSE_ztrsmpl
 * @sa MORSE_ztrsmpl_Tile
 * @sa MORSE_ctrsmpl_Tile_Async
 * @sa MORSE_dtrsmpl_Tile_Async
 * @sa MORSE_strsmpl_Tile_Async
 * @sa MORSE_zgetrf_incpiv_Tile_Async
 *
 */
int MORSE_ztrsmpl_Tile_Async( MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B,
                              MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsmpl_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztrsmpl_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztrsmpl_Tile", "NULL request");
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
        morse_error("MORSE_ztrsmpl_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(L) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsmpl_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsmpl_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_ztrsmpl_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return MORSE_SUCCESS;
     */
    morse_pztrsmpl( A, B, L, IPIV, sequence, request );

    return MORSE_SUCCESS;
}
