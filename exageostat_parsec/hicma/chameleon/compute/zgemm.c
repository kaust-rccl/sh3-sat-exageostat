/**
 *
 * @file zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm wrappers
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
 *
 * @defgroup MORSE_Complex64_t
 * @brief Linear algebra routines exposed to users. LAPACK matrix data storage
 *
 */

/**
 *
 * @defgroup MORSE_Complex64_t_Tile
 * @brief Linear algebra routines exposed to users. Tile matrix data storage
 *
 */

/**
 *
 * @defgroup MORSE_Complex64_t_Tile_Async
 * @brief Linear algebra routines exposed to users. Tile matrix data storage,
 *  asynchronous interface.
 *
 */

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zgemm - Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   A is not transposed;
 *          = MorseTrans:     A is transposed;
 *          = MorseConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   B is not transposed;
 *          = MorseTrans:     B is transposed;
 *          = MorseConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = MorseNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = MorseNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgemm_Tile
 * @sa MORSE_cgemm
 * @sa MORSE_dgemm
 * @sa MORSE_sgemm
 *
 */
int MORSE_zgemm( MORSE_enum transA, MORSE_enum transB, int M, int N, int K,
                 MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t *B, int LDB,
                 MORSE_Complex64_t beta,  MORSE_Complex64_t *C, int LDC )
{
    int NB;
    int Am, An, Bm, Bn;
    int status;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descBl, descBt;
    MORSE_desc_t descCl, descCt;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgemm", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((transA < MorseNoTrans) || (transA > MorseConjTrans)) {
        morse_error("MORSE_zgemm", "illegal value of transA");
        return -1;
    }
    if ((transB < MorseNoTrans) || (transB > MorseConjTrans)) {
        morse_error("MORSE_zgemm", "illegal value of transB");
        return -2;
    }
    if ( transA == MorseNoTrans ) {
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if ( transB == MorseNoTrans ) {
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }
    if (M < 0) {
        morse_error("MORSE_zgemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        morse_error("MORSE_zgemm", "illegal value of N");
        return -4;
    }
    if (K < 0) {
        morse_error("MORSE_zgemm", "illegal value of N");
        return -5;
    }
    if (LDA < chameleon_max(1, Am)) {
        morse_error("MORSE_zgemm", "illegal value of LDA");
        return -8;
    }
    if (LDB < chameleon_max(1, Bm)) {
        morse_error("MORSE_zgemm", "illegal value of LDB");
        return -10;
    }
    if (LDC < chameleon_max(1, M)) {
        morse_error("MORSE_zgemm", "illegal value of LDC");
        return -13;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (MORSE_Complex64_t)0.0 || K == 0) && beta == (MORSE_Complex64_t)1.0))
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGEMM, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgemm", "morse_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, An, Am, An, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInput, MorseUpperLower,
                     B, NB, NB, LDB, Bn, Bm, Bn, sequence, &request );
    morse_zlap2tile( morse, &descCl, &descCt, MorseDescInout, MorseUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    MORSE_zgemm_Tile_Async( transA, transB, alpha, &descAt, &descBt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpperLower, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInput, MorseUpperLower, sequence, &request );
    morse_ztile2lap( morse, &descCl, &descCt,
                     MorseDescInout, MorseUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );
    morse_ztile2lap_cleanup( morse, &descCl, &descCt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zgemm_Tile - Performs matrix multiplication.
 *  Tile equivalent of MORSE_zgemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   A is not transposed;
 *          = MorseTrans:     A is transposed;
 *          = MorseConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   B is not transposed;
 *          = MorseTrans:     B is transposed;
 *          = MorseConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = MorseNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = MorseNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgemm
 * @sa MORSE_zgemm_Tile_Async
 * @sa MORSE_cgemm_Tile
 * @sa MORSE_dgemm_Tile
 * @sa MORSE_sgemm_Tile
 *
 */
int MORSE_zgemm_Tile( MORSE_enum transA, MORSE_enum transB,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                      MORSE_Complex64_t beta,  MORSE_desc_t *C )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zgemm_Tile_Async( transA, transB, alpha, A, B, beta, C, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( B, sequence );
    MORSE_Desc_Flush( C, sequence );

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
 *  MORSE_zgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of MORSE_zgemm_Tile().
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
 * @sa MORSE_zgemm
 * @sa MORSE_zgemm_Tile
 * @sa MORSE_cgemm_Tile_Async
 * @sa MORSE_dgemm_Tile_Async
 * @sa MORSE_sgemm_Tile_Async
 *
 */
int MORSE_zgemm_Tile_Async( MORSE_enum transA, MORSE_enum transB,
                            MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                            MORSE_Complex64_t beta,  MORSE_desc_t *C,
                            MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "NULL request");
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
        morse_error("MORSE_zgemm_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgemm_Tile_Async", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(C) != MORSE_SUCCESS) {
        morse_error("MORSE_zgemm_Tile_Async", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((transA < MorseNoTrans) || (transA > MorseConjTrans)) {
        morse_error("MORSE_zgemm_Tile_Async", "illegal value of transA");
        return morse_request_fail(sequence, request, -1);
    }
    if ((transB < MorseNoTrans) || (transB > MorseConjTrans)) {
        morse_error("MORSE_zgemm_Tile_Async", "illegal value of transB");
        return morse_request_fail(sequence, request, -2);
    }

    if ( transA == MorseNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }

    if ( transB == MorseNoTrans ) {
        Bm  = B->m;
        Bn  = B->n;
        Bmb = B->mb;
        Bnb = B->nb;
        Bi  = B->i;
        Bj  = B->j;
    } else {
        Bm  = B->n;
        Bn  = B->m;
        Bmb = B->nb;
        Bnb = B->mb;
        Bi  = B->j;
        Bj  = B->i;
    }

    if ( (Amb != C->mb) || (Anb != Bmb) || (Bnb != C->nb) ) {
        morse_error("MORSE_zgemm_Tile_Async", "tile sizes have to match");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != C->m) || (An != Bm) || (Bn != C->n) ) {
        morse_error("MORSE_zgemm_Tile_Async", "sizes of matrices have to match");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != C->i) || (Aj != Bi) || (Bj != C->j) ) {
        morse_error("MORSE_zgemm_Tile_Async", "start indexes have to match");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    M = C->m;
    N = C->n;
    K = An;

    /* Quick return */
    if ( (M == 0) || (N == 0) ||
         (((alpha == (MORSE_Complex64_t)0.0) || (K == 0)) && (beta == (MORSE_Complex64_t)1.0)) )
    {
        return MORSE_SUCCESS;
    }

    morse_pzgemm( transA, transB, alpha, A, B, beta, C, sequence, request );

    return MORSE_SUCCESS;
}
