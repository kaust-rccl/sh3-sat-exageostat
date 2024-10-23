/**
 *
 * @file ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt wrappers
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_ztpqrt - Computes a blocked QR factorization of a
 *  "triangular-pentagonal" matrix C, which is composed of a triangular block A
 *  and a pentagonal block B, using the compact representation for Q.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix B. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix B, and the order of the matrix
 *          A. N >= 0.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of B.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in,out] A
 *          On entry, the upper triangular N-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array
 *          contain the upper triangular matrix R.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
 *          are rectangular, and the last L rows are upper trapezoidal.
 *          On exit, B contains the pentagonal matrix V.  See Further Details.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.  LDB >= max(1,M).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data, required by MORSE_zgeqrs to
 *          solve the system of equations, or by any function to apply the Q.
 *
 * @par Further Details:
 * =====================
 *
 *  The input matrix C is a (N+M)-by-N matrix
 *
 *               C = [ A ]
 *                   [ B ]
 *
 *  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
 *  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
 *  upper trapezoidal matrix B2:
 *
 *               B = [ B1 ]  <- (M-L)-by-N rectangular
 *                   [ B2 ]  <-     L-by-N upper trapezoidal.
 *
 *  The upper trapezoidal matrix B2 consists of the first L rows of a
 *  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
 *  B is rectangular M-by-N; if M=L=N, B is upper triangular.
 *
 *  The matrix W stores the elementary reflectors H(i) in the i-th column
 *  below the diagonal (of A) in the (N+M)-by-N input matrix C
 *
 *               C = [ A ]  <- upper triangular N-by-N
 *                   [ B ]  <- M-by-N pentagonal
 *
 *  so that W can be represented as
 *
 *               W = [ I ]  <- identity, N-by-N
 *                   [ V ]  <- M-by-N, same form as B.
 *
 *  Thus, all of information needed for W is contained on exit in B, which
 *  we call V above.  Note that V has the same form as B; that is,
 *
 *               V = [ V1 ] <- (M-L)-by-N rectangular
 *                   [ V2 ] <-     L-by-N upper trapezoidal.
 *
 *  The columns of V represent the vectors which define the H(i)'s.
 *
 *  The number of blocks is B = ceiling(N/NB), where each
 *  block is of order NB except for the last block, which is of order
 *  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
 *  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
 *  for the last block) T's are stored in the NB-by-N matrix T as
 *
 *               T = [T1 T2 ... TB].
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_ztpqrt_Tile
 * @sa MORSE_ztpqrt_Tile_Async
 * @sa MORSE_ctpqrt
 * @sa MORSE_dtpqrt
 * @sa MORSE_stpqrt
 * @sa MORSE_zgeqrs
 *
 */
int MORSE_ztpqrt( int M, int N, int L,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB,
                  MORSE_desc_t *descT )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descBl, descBt;
    int minMN = chameleon_min( M, N );

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztpqrt", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_ztpqrt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_ztpqrt", "illegal value of N");
        return -2;
    }
    if ((L < 0) || ((L > minMN) && (minMN > 0))) {
        morse_error("MORSE_ztpqrt", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_ztpqrt", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, M)) {
        morse_error("MORSE_ztpqrt", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (minMN == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_ztpqrt", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, MorseUpper,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, MorseDescInout, MorseUpperLower,
                     B, NB, NB, LDB, N, M, N, sequence, &request );

    /* Call the tile interface */
    MORSE_ztpqrt_Tile_Async( L, &descAt, &descBt, descT, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInout, MorseUpper, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     MorseDescInout, MorseUpperLower, sequence, &request );
    MORSE_Desc_Flush( descT, sequence );

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
 *  MORSE_ztpqrt_Tile - Computes the tile QR factorization of a matrix.
 *  Tile equivalent of MORSE_ztpqrt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by MORSE_zgeqrs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_ztpqrt
 * @sa MORSE_ztpqrt_Tile_Async
 * @sa MORSE_ctpqrt_Tile
 * @sa MORSE_dtpqrt_Tile
 * @sa MORSE_stpqrt_Tile
 * @sa MORSE_zgeqrs_Tile
 *
 */
int MORSE_ztpqrt_Tile( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztpqrt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_ztpqrt_Tile_Async( L, A, B, T, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( B, sequence );
    MORSE_Desc_Flush( T, sequence );

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
 *  MORSE_ztpqrt_Tile_Async - Computes the tile QR factorization of a matrix.
 *  Non-blocking equivalent of MORSE_ztpqrt_Tile().
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
 * @sa MORSE_ztpqrt
 * @sa MORSE_ztpqrt_Tile
 * @sa MORSE_ctpqrt_Tile_Async
 * @sa MORSE_dtpqrt_Tile_Async
 * @sa MORSE_stpqrt_Tile_Async
 * @sa MORSE_zgeqrs_Tile_Async
 *
 */
int MORSE_ztpqrt_Tile_Async( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_ztpqrt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztpqrt_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztpqrt_Tile", "NULL request");
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
        morse_error("MORSE_ztpqrt_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpqrt_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpqrt_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_ztpqrt_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ( (L != 0) && (((B->m - L) % B->mb) != 0) ) {
        morse_error("MORSE_ztpqrt_Tile", "Triangular part must be aligned with tiles");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* if (morse->householder == MORSE_FLAT_HOUSEHOLDER) { */
    morse_pztpqrt( L, A, B, T, sequence, request );
    /* } */
    /* else { */
    /*    morse_pztpqrtrh( A, T, MORSE_RHBLK, sequence, request ); */
    /* } */

    return MORSE_SUCCESS;
}
