/**
 *
 * @file zunmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zunmlq - Overwrites the general complex M-by-N matrix C with
 *
 *                  SIDE = 'L'     SIDE = 'R'
 *  TRANS = 'N':      Q * C          C * Q
 *  TRANS = 'C':      Q**H * C       C * Q**H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *        Q = H(1) H(2) . . . H(k)
 *
 *  as returned by MORSE_zgeqrf. Q is of order M if SIDE = MorseLeft
 *  and of order N if SIDE = MorseRight.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = MorseLeft:  apply Q or Q**H from the left;
 *          = MorseRight: apply Q or Q**H from the right.
 *
 * @param[in] trans
 *          Intended usage:
 *          = MorseNoTrans:   no transpose, apply Q;
 *          = MorseConjTrans: conjugate transpose, apply Q**H.
 *
 * @param[in] M
 *          The number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix C. N >= 0.
 *
 * @param[in] K
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          If side == MorseLeft,  M >= K >= 0.
 *          If side == MorseRight, N >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K).
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q**H*C.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zunmlq_Tile
 * @sa MORSE_zunmlq_Tile_Async
 * @sa MORSE_cunmlq
 * @sa MORSE_dormlq
 * @sa MORSE_sormlq
 * @sa MORSE_zgelqf
 *
 */
int MORSE_zunmlq( MORSE_enum side, MORSE_enum trans, int M, int N, int K,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_desc_t *descT,
                  MORSE_Complex64_t *C, int LDC )
{
    int NB, An;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descCl, descCt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunmlq", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    if (side == MorseLeft)
        An = M;
    else
        An = N;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        morse_error("MORSE_zunmlq", "illegal value of side");
        return -1;
    }
    if ((trans != MorseConjTrans) && (trans != MorseNoTrans)){
        morse_error("MORSE_zunmlq", "illegal value of trans");
        return -2;
    }
    if (M < 0) {
        morse_error("MORSE_zunmlq", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        morse_error("MORSE_zunmlq", "illegal value of N");
        return -4;
    }
    if ((K < 0) || (K > An)) {
        morse_error("MORSE_zunmlq", "illegal value of K");
        return -5;
    }
    if (LDA < chameleon_max(1, K)) {
        morse_error("MORSE_zunmlq", "illegal value of LDA");
        return -7;
    }
    if (LDC < chameleon_max(1, M)) {
        morse_error("MORSE_zunmlq", "illegal value of LDC");
        return -10;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
    if (chameleon_min(M, chameleon_min(N, K)) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGELS, M, K, N);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zunmlq", "morse_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
    NB   = MORSE_NB;
    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpper,
                     A, NB, NB, LDA, An, K, An, sequence, &request );
    morse_zlap2tile( morse, &descCl, &descCt, MorseDescInout, MorseUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    MORSE_zunmlq_Tile_Async(  side, trans, &descAt, descT, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpper, sequence, &request );
    morse_ztile2lap( morse, &descCl, &descCt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descT, sequence );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descCl, &descCt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zunmlq_Tile - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by MORSE_zgelqf_Tile Q is of order M.
 *  All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = MorseLeft:  apply Q or Q**H from the left;
 *          = MorseRight: apply Q or Q**H from the right.
 *          Currently only MorseLeft is supported.
 *
 * @param[in] trans
 *          Intended usage:
 *          = MorseNoTrans:   no transpose, apply Q;
 *          = MorseConjTrans: conjugate transpose, apply Q**H.
 *          Currently only MorseConjTrans is supported.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q**H*C.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zunmlq
 * @sa MORSE_zunmlq_Tile_Async
 * @sa MORSE_cunmlq_Tile
 * @sa MORSE_dormlq_Tile
 * @sa MORSE_sormlq_Tile
 * @sa MORSE_zgelqf_Tile
 *
 */
int MORSE_zunmlq_Tile( MORSE_enum side, MORSE_enum trans,
                       MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *C )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunmlq_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zunmlq_Tile_Async(side, trans, A, T, C, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( T, sequence );
    MORSE_Desc_Flush( C, sequence );

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
 *  Non-blocking equivalent of MORSE_zunmlq_Tile().
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
 * @sa MORSE_zunmlq
 * @sa MORSE_zunmlq_Tile
 * @sa MORSE_cunmlq_Tile_Async
 * @sa MORSE_dormlq_Tile_Async
 * @sa MORSE_sormlq_Tile_Async
 * @sa MORSE_zgelqf_Tile_Async
 *
 */
int MORSE_zunmlq_Tile_Async( MORSE_enum side, MORSE_enum trans,
                             MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *C,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunmlq_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zunmlq_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zunmlq_Tile", "NULL request");
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
        morse_error("MORSE_zunmlq_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zunmlq_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(C) != MORSE_SUCCESS) {
        morse_error("MORSE_zunmlq_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || C->nb != C->mb) {
        morse_error("MORSE_zunmlq_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((side != MorseLeft) && (side != MorseRight)) {
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((trans != MorseConjTrans) && (trans != MorseNoTrans)){
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
    /*
     if (chameleon_min(M, chameleon_min(N, K)) == 0)
     return MORSE_SUCCESS;
     */
#if defined(CHAMELEON_COPY_DIAG)
    {
        int m = chameleon_min(A->mt, A->nt) * A->mb;
        morse_zdesc_alloc(D, A->mb, A->nb, m, A->n, 0, 0, m, A->n, );
        Dptr = &D;
    }
#endif

    if (morse->householder == MORSE_FLAT_HOUSEHOLDER) {
        if ( (trans == MorseConjTrans) &&
             (side == MorseLeft) ) {
            morse_pzunmlq( side, trans, A, C, T, Dptr, sequence, request );
        } else {
            morse_pzunmlq( side, trans, A, C, T, Dptr, sequence, request );
        }
    }
    else {
        morse_pzunmlqrh( side, trans, A, C, T, Dptr, MORSE_RHBLK, sequence, request );
    }
    if (Dptr != NULL) {
        MORSE_Desc_Flush( A, sequence );
        MORSE_Desc_Flush( C, sequence );
        MORSE_Desc_Flush( T, sequence );
        MORSE_Desc_Flush( Dptr, sequence );
        morse_sequence_wait( morse, sequence );
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return MORSE_SUCCESS;
}
