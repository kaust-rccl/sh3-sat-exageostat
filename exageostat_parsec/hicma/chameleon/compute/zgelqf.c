/**
 *
 * @file zgelqf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
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
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zgelqf - Computes the tile LQ factorization of a complex M-by-N matrix A: A = L * Q.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data, required by MORSE_zgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zgelqf_Tile
 * @sa MORSE_zgelqf_Tile_Async
 * @sa MORSE_cgelqf
 * @sa MORSE_dgelqf
 * @sa MORSE_sgelqf
 * @sa MORSE_zgelqs
 *
 */
int MORSE_zgelqf( int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_desc_t *descT )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgelqf", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zgelqf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_zgelqf", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("MORSE_zgelqf", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqf", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, MorseUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zgelqf_Tile_Async( &descAt, descT, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descT, sequence );

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
 *  MORSE_zgelqf_Tile - Computes the tile LQ factorization of a matrix.
 *  Tile equivalent of MORSE_zgelqf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by MORSE_zgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgelqf
 * @sa MORSE_zgelqf_Tile_Async
 * @sa MORSE_cgelqf_Tile
 * @sa MORSE_dgelqf_Tile
 * @sa MORSE_sgelqf_Tile
 * @sa MORSE_zgelqs_Tile
 *
 */
int MORSE_zgelqf_Tile( MORSE_desc_t *A, MORSE_desc_t *T )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgelqf_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zgelqf_Tile_Async( A, T, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( T, sequence );

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
 *  MORSE_zgelqf_Tile_Async - Computes the tile LQ factorization of a matrix.
 *  Non-blocking equivalent of MORSE_zgelqf_Tile().
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
 * @sa MORSE_zgelqf
 * @sa MORSE_zgelqf_Tile
 * @sa MORSE_cgelqf_Tile_Async
 * @sa MORSE_dgelqf_Tile_Async
 * @sa MORSE_sgelqf_Tile_Async
 * @sa MORSE_zgelqs_Tile_Async
 *
 */
int MORSE_zgelqf_Tile_Async( MORSE_desc_t *A, MORSE_desc_t *T,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgelqf_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgelqf_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgelqf_Tile", "NULL request");
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
        morse_error("MORSE_zgelqf_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zgelqf_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_zgelqf_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
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

    if (morse->householder == MORSE_FLAT_HOUSEHOLDER) {
        morse_pzgelqf( A, T, Dptr, sequence, request );
    }
    else {
        morse_pzgelqfrh( A, T, Dptr, MORSE_RHBLK, sequence, request );
    }
    if (Dptr != NULL) {
        MORSE_Desc_Flush( A, sequence );
        MORSE_Desc_Flush( T, sequence );
        MORSE_Desc_Flush( Dptr, sequence );
        morse_sequence_wait( morse, sequence );
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return MORSE_SUCCESS;
}
