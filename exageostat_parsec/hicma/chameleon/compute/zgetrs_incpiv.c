/**
 *
 * @file zgetrs_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrs_incpiv wrappers
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
 *  MORSE_zgetrs_incpiv - Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization computed by MORSE_zgetrf_incpiv.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = MorseNoTrans:   A * X = B     (No transpose)
 *          = MorseTrans:     A**T * X = B  (Transpose)
 *          = MorseConjTrans: A**H * X = B  (Conjugate transpose)
 *          Currently only MorseNoTrans is supported.
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by MORSE_zgetrf_incpiv.
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
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zgetrs_incpiv_Tile
 * @sa MORSE_zgetrs_incpiv_Tile_Async
 * @sa MORSE_cgetrs_incpiv
 * @sa MORSE_dgetrs_incpiv
 * @sa MORSE_sgetrs_incpiv
 * @sa MORSE_zgetrf_incpiv
 *
 */
int MORSE_zgetrs_incpiv( MORSE_enum trans, int N, int NRHS,
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
        morse_fatal_error("MORSE_zgetrs_incpiv", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (trans != MorseNoTrans) {
        morse_error("MORSE_zgetrs_incpiv", "only MorseNoTrans supported");
        return MORSE_ERR_NOT_SUPPORTED;
    }
    if (N < 0) {
        morse_error("MORSE_zgetrs_incpiv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zgetrs_incpiv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_zgetrs_incpiv", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, N)) {
        morse_error("MORSE_zgetrs_incpiv", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGESV, N, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgetrs_incpiv", "morse_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    MORSE_zgetrs_incpiv_Tile_Async( &descAt, descL, IPIV, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpperLower, sequence, &request );
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
 *  MORSE_zgetrs_incpiv_Tile - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Tile equivalent of MORSE_zgetrs_incpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by MORSE_zgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from MORSE_zgetrf_incpiv (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgetrs_incpiv
 * @sa MORSE_zgetrs_incpiv_Tile_Async
 * @sa MORSE_cgetrs_incpiv_Tile
 * @sa MORSE_dgetrs_incpiv_Tile
 * @sa MORSE_sgetrs_incpiv_Tile
 * @sa MORSE_zgetrf_incpiv_Tile
 *
 */
int MORSE_zgetrs_incpiv_Tile( MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgetrs_incpiv_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zgetrs_incpiv_Tile_Async( A, L, IPIV, B, sequence, &request );

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
 *  MORSE_zgetrs_incpiv_Tile_Async - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Non-blocking equivalent of MORSE_zgetrs_incpiv_Tile().
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
 * @sa MORSE_zgetrs_incpiv
 * @sa MORSE_zgetrs_incpiv_Tile
 * @sa MORSE_cgetrs_incpiv_Tile_Async
 * @sa MORSE_dgetrs_incpiv_Tile_Async
 * @sa MORSE_sgetrs_incpiv_Tile_Async
 * @sa MORSE_zgetrf_incpiv_Tile_Async
 *
 */
int MORSE_zgetrs_incpiv_Tile_Async( MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B,
                                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgetrs_incpiv_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgetrs_incpiv_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgetrs_incpiv_Tile", "NULL request");
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
        morse_error("MORSE_zgetrs_incpiv_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(L) != MORSE_SUCCESS) {
        morse_error("MORSE_zgetrs_incpiv_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgetrs_incpiv_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zgetrs_incpiv_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return MORSE_SUCCESS;
     */
    morse_pztrsmpl( A, B, L, IPIV, sequence, request );

    morse_pztrsm( MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, 1.0, A, B, sequence, request );

    return MORSE_SUCCESS;
}
