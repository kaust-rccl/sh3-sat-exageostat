/**
 *
 * @file zgelqs_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqs_param wrappers
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
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
 *  MORSE_zgelqs_param - Compute a minimum-norm solution min || A*X - B || using the LQ factorization
 *  A = L*Q computed by MORSE_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
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
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] descTS
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in] descTT
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= N.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zgelqs_param_Tile
 * @sa MORSE_zgelqs_param_Tile_Async
 * @sa MORSE_cgelqs
 * @sa MORSE_dgelqs
 * @sa MORSE_sgelqs
 * @sa MORSE_zgelqf
 *
 */
int MORSE_zgelqs_param( const libhqr_tree_t *qrtree, int M, int N, int NRHS,
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
        morse_fatal_error("MORSE_zgelqs_param", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zgelqs_param", "illegal value of M");
        return -1;
    }
    if (N < 0 || M > N) {
        morse_error("MORSE_zgelqs_param", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zgelqs_param", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zgelqs_param", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, chameleon_max(1, N))) {
        morse_error("MORSE_zgelqs_param", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
        return MORSE_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqs_param", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    MORSE_zgelqs_param_Tile_Async( qrtree, &descAt, descTS, descTT, &descBt, sequence, &request );

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
 *  MORSE_zgelqs_param_Tile - Computes a minimum-norm solution using previously computed
 *  LQ factorization.
 *  Tile equivalent of MORSE_zgelqs_param().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] TS
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in] TT
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgelqs_param
 * @sa MORSE_zgelqs_param_Tile_Async
 * @sa MORSE_cgelqs_Tile
 * @sa MORSE_dgelqs_Tile
 * @sa MORSE_sgelqs_Tile
 * @sa MORSE_zgelqf_Tile
 *
 */
int MORSE_zgelqs_param_Tile( const libhqr_tree_t *qrtree, MORSE_desc_t *A,
                             MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgelqs_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zgelqs_param_Tile_Async( qrtree, A, TS, TT, B, sequence, &request );

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
 *  MORSE_zgelqs_param_Tile_Async - Computes a minimum-norm solution using previously
 *  computed LQ factorization.
 *  Non-blocking equivalent of MORSE_zgelqs_param_Tile().
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
 * @sa MORSE_zgelqs_param
 * @sa MORSE_zgelqs_param_Tile
 * @sa MORSE_cgelqs_Tile_Async
 * @sa MORSE_dgelqs_Tile_Async
 * @sa MORSE_sgelqs_Tile_Async
 * @sa MORSE_zgelqf_Tile_Async
 *
 */
int MORSE_zgelqs_param_Tile_Async( const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B,
                                   MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_desc_t *subB;
    MORSE_desc_t *subA;
    MORSE_context_t *morse;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgelqs_param_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgelqs_param_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgelqs_param_Tile", "NULL request");
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
        morse_error("MORSE_zgelqs_param_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TS) != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqs_param_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(TT) != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqs_param_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqs_param_Tile", "invalid fourth descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zgelqs_param_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
     return MORSE_SUCCESS;
     }
     */
    /* subB = morse_desc_submatrix(B, A->m, 0, A->n-A->m, B->n);
     morse_pzlaset( MorseUpperLower, 0., 0., subB, sequence, request );
     free(subB); */

    subB = morse_desc_submatrix(B, 0, 0, A->m, B->n);
    subA = morse_desc_submatrix(A, 0, 0, A->m, A->m);
    morse_pztrsm( MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1.0, subA, subB, sequence, request );
    free(subA);
    free(subB);

#if defined(CHAMELEON_COPY_DIAG)
    {
        int m = chameleon_min(A->mt, A->nt) * A->mb;
        morse_zdesc_alloc(D, A->mb, A->nb, m, A->n, 0, 0, m, A->n, );
        Dptr = &D;
    }
#endif

    morse_pzunmlq_param( qrtree, MorseLeft, MorseConjTrans, A, B, TS, TT, Dptr, sequence, request );
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
