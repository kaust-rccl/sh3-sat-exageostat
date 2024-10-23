/**
 *
 * @file zgeqrs_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrs_param wrappers
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-17
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zgeqrs_param - Compute a minimum-norm solution min || A*X - B || using the RQ factorization
 *  A = R*Q computed by MORSE_zgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= M >= 0.
 *
 * @param[in] NRHS
 *          The number of columns of B. NRHS >= 0.
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by MORSE_zgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by MORSE_zgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
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
 * @sa MORSE_zgeqrs_param_Tile
 * @sa MORSE_zgeqrs_param_Tile_Async
 * @sa MORSE_cgeqrs
 * @sa MORSE_dgeqrs
 * @sa MORSE_sgeqrs
 * @sa MORSE_zgeqrf
 *
 */
int MORSE_zgeqrs_param( const libhqr_tree_t *qrtree, int M, int N, int NRHS,
                        MORSE_Complex64_t *A, int LDA,
                        MORSE_desc_t *descTS, MORSE_desc_t *descTT,
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
        morse_fatal_error("MORSE_zgeqrs_param", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zgeqrs_param", "illegal value of M");
        return -1;
    }
    if (N < 0 || N > M) {
        morse_error("MORSE_zgeqrs_param", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zgeqrs_param", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zgeqrs_param", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, chameleon_max(1, M))) {
        morse_error("MORSE_zgeqrs_param", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
        return MORSE_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_param", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, NRHS, M, NRHS, sequence, &request );

    /* Call the tile interface */
    MORSE_zgeqrs_param_Tile_Async( qrtree, &descAt, descTS, descTT, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpperLower, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descTS, sequence );
    MORSE_Desc_Flush( descTT, sequence );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zgeqrs_param_Tile - Computes a minimum-norm solution using the tile QR factorization.
 *  Tile equivalent of MORSE_zgeqrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by MORSE_zgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by MORSE_zgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgeqrs_param
 * @sa MORSE_zgeqrs_param_Tile_Async
 * @sa MORSE_cgeqrs_Tile
 * @sa MORSE_dgeqrs_Tile
 * @sa MORSE_sgeqrs_Tile
 * @sa MORSE_zgeqrf_Tile
 *
 */
int MORSE_zgeqrs_param_Tile( const libhqr_tree_t *qrtree, MORSE_desc_t *A,
                             MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgeqrs_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zgeqrs_param_Tile_Async( qrtree, A, TS, TT, B, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( TS, sequence );
    MORSE_Desc_Flush( TT, sequence );
    MORSE_Desc_Flush( B, sequence );

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
 *  MORSE_zgeqrs_param_Tile_Async - Computes a minimum-norm solution using the tile
 *  QR factorization.
 *  Non-blocking equivalent of MORSE_zgeqrs_param_Tile().
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
 * @sa MORSE_zgeqrs_param
 * @sa MORSE_zgeqrs_param_Tile
 * @sa MORSE_cgeqrs_Tile_Async
 * @sa MORSE_dgeqrs_Tile_Async
 * @sa MORSE_sgeqrs_Tile_Async
 * @sa MORSE_zgeqrf_Tile_Async
 *
 */
int MORSE_zgeqrs_param_Tile_Async( const libhqr_tree_t *qrtree,
                                   MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B,
                                   MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_desc_t *subA;
    MORSE_desc_t *subB;
    MORSE_context_t *morse;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgeqrs_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgeqrs_param_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgeqrs_param_Tile", "NULL request");
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
        morse_error("MORSE_zgeqrs_param_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TS) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_param_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TT) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_param_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_param_Tile", "invalid fourth descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zgeqrs_param_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
     return MORSE_SUCCESS;
     }
     */
#if defined(CHAMELEON_COPY_DIAG)
    {
        int n = chameleon_min(A->mt, A->nt) * A->nb;
        morse_zdesc_alloc(D, A->mb, A->nb, A->m, n, 0, 0, A->m, n, );
        Dptr = &D;
    }
#endif

    subB = morse_desc_submatrix(B, 0, 0, A->n, B->n);
    subA = morse_desc_submatrix(A, 0, 0, A->n, A->n);

    morse_pzunmqr_param( qrtree, MorseLeft, MorseConjTrans, A, B, TS, TT, Dptr, sequence, request );
    morse_pztrsm( MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, 1.0, subA, subB, sequence, request );

    free(subA);
    free(subB);

    if (Dptr != NULL) {
        MORSE_Desc_Flush( A, sequence );
        MORSE_Desc_Flush( B, sequence );
        MORSE_Desc_Flush( TS, sequence );
        MORSE_Desc_Flush( TT, sequence );
        MORSE_Desc_Flush( Dptr, sequence );
        morse_sequence_wait( morse, sequence );
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return MORSE_SUCCESS;
}
