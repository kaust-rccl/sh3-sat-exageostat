/**
 *
 * @file zunglq_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq_param wrappers
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-17
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zunglq_param - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 *  first M rows of a product of the elementary reflectors returned by MORSE_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
 *
 * @param[in] M
 *          The number of rows of the matrix Q. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix Q. N >= M.
 *
 * @param[in] K
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval MORSE_SUCCESS <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zunglq_param_Tile
 * @sa MORSE_zunglq_param_Tile_Async
 * @sa MORSE_cunglq
 * @sa MORSE_dorglq
 * @sa MORSE_sorglq
 * @sa MORSE_zgelqf
 *
 */
int MORSE_zunglq_param( const libhqr_tree_t *qrtree, int M, int N, int K,
                        MORSE_Complex64_t *A, int LDA,
                        MORSE_desc_t *descTS, MORSE_desc_t *descTT,
                        MORSE_Complex64_t *Q, int LDQ )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descQl, descQt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq_param", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zunglq_param", "illegal value of M");
        return -1;
    }
    if (N < M) {
        morse_error("MORSE_zunglq_param", "illegal value of N");
        return -2;
    }
    if (K < 0 || K > M) {
        morse_error("MORSE_zunglq_param", "illegal value of K");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zunglq_param", "illegal value of LDA");
        return -5;
    }
    if (LDQ < chameleon_max(1, M)) {
        morse_error("MORSE_zunglq_param", "illegal value of LDQ");
        return -8;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDQ ) */
    if (chameleon_min(M, chameleon_min(N, K)) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_param", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpper,
                     A, NB, NB, LDA, N, K, N, sequence, &request );
    morse_zlap2tile( morse, &descQl, &descQt, MorseDescInout, MorseUpperLower,
                     Q, NB, NB, LDQ, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zunglq_param_Tile_Async( qrtree, &descAt, descTS, descTT, &descQt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpper, sequence, &request );
    morse_ztile2lap( morse, &descQl, &descQt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descTS, sequence );
    MORSE_Desc_Flush( descTT, sequence );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descQl, &descQt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 * MORSE_zunglq_param_Tile - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 * first M rows of a product of the elementary reflectors returned by MORSE_zgelqf.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zunglq_param
 * @sa MORSE_zunglq_param_Tile_Async
 * @sa MORSE_cunglq_Tile
 * @sa MORSE_dorglq_Tile
 * @sa MORSE_sorglq_Tile
 * @sa MORSE_zgelqf_Tile
 *
 */
int MORSE_zunglq_param_Tile( const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *Q )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zunglq_param_Tile_Async( qrtree, A, TS, TT, Q, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( TS, sequence );
    MORSE_Desc_Flush( TT, sequence );
    MORSE_Desc_Flush( Q, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  Non-blocking equivalent of MORSE_zunglq_param_Tile().
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
 * @sa MORSE_zunglq_param
 * @sa MORSE_zunglq_param_Tile
 * @sa MORSE_cunglq_Tile_Async
 * @sa MORSE_dorglq_Tile_Async
 * @sa MORSE_sorglq_Tile_Async
 * @sa MORSE_zgelqf_Tile_Async
 *
 */
int MORSE_zunglq_param_Tile_Async( const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *Q,
                                   MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zunglq_param_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zunglq_param_Tile", "NULL request");
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
        morse_error("MORSE_zunglq_param_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TS) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_param_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TT) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_param_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(Q) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_param_Tile", "invalid fourth descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || Q->nb != Q->mb) {
        morse_error("MORSE_zunglq_param_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, Q, LDQ ) */
    /*
     if (chameleon_min(M, N) == 0)
     return MORSE_SUCCESS;
     */
#if defined(CHAMELEON_COPY_DIAG)
    {
        int m = chameleon_min(A->mt, A->nt) * A->mb;
        morse_zdesc_alloc(D, A->mb, A->nb, m, A->n, 0, 0, m, A->n, );
        Dptr = &D;
    }
#endif

    morse_pzlaset( MorseUpperLower, 0., 1., Q, sequence, request );
    morse_pzunglq_param( qrtree, A, Q, TS, TT, Dptr, sequence, request );

    if (Dptr != NULL) {
        MORSE_Desc_Flush( A, sequence );
        MORSE_Desc_Flush( Q, sequence );
        MORSE_Desc_Flush( TS, sequence );
        MORSE_Desc_Flush( TT, sequence );
        MORSE_Desc_Flush( Dptr, sequence );
        morse_sequence_wait( morse, sequence );
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return MORSE_SUCCESS;
}
