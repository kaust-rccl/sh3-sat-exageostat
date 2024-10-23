/**
 *
 * @file zheevd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zheevd wrappers
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include <string.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#endif

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zheevd - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex Hermitian matrix A. The matrix A is
 *  preliminary reduced to tridiagonal form using a two-stage
 *  approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = MorseNoVec: computes eigenvalues only;
 *          = MorseVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = MorseLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = MorseLower) or the
 *          upper triangle (if uplo = MorseUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by MORSE_Alloc_Workspace_zheevd
 *          On exit, contains auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa MORSE_zheevd_Tile
 * @sa MORSE_zheevd_Tile_Async
 * @sa MORSE_cheevd
 * @sa MORSE_dsyev
 * @sa MORSE_ssyev
 *
 */
int MORSE_zheevd( MORSE_enum jobz, MORSE_enum uplo, int N,
                  MORSE_Complex64_t *A, int LDA,
                  double *W,
                  MORSE_desc_t *descT )
{
    int NB;
    int status;
    MORSE_context_t  *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t   request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zheevd", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != MorseNoVec && jobz != MorseVec) {
        morse_error("MORSE_zheevd", "illegal value of jobz");
        return -1;
    }
    if ((uplo != MorseLower) && (uplo != MorseUpper)) {
        morse_error("MORSE_zheevd", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        morse_error("MORSE_zheevd", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_zheevd", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (N == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZHEEVD, N, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zheevd", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zheevd_Tile_Async( jobz, uplo, &descAt, W, descT, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInout, uplo, sequence, &request );

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
 *  MORSE_zheevd_Tile - Computes all eigenvalues and, optionally, eigenvectors of a
 *  complex Hermitian matrix A using a two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = MorseNoVec: computes eigenvalues only;
 *          = MorseVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = MorseVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = MorseNoVec, then on exit the lower triangle (if
 *          uplo = MorseLower) or the upper triangle (if uplo =
 *          MorseUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          MORSE_Alloc_Workspace_zheevd
 *          On exit, contains auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa MORSE_zheevd_Tile
 * @sa MORSE_zheevd_Tile_Async
 * @sa MORSE_cheevd
 * @sa MORSE_dsyev
 * @sa MORSE_ssyev
 *
 */
int MORSE_zheevd_Tile( MORSE_enum jobz, MORSE_enum uplo,
                       MORSE_desc_t *A, double *W, MORSE_desc_t *T )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zheevd_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zheevd_Tile_Async( jobz, uplo, A, W, T, sequence, &request );

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
 *  MORSE_zheevd_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = MorseNoVec: computes eigenvalues only;
 *          = MorseVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = MorseVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = MorseNoVec, then on exit the lower triangle (if
 *          uplo = MorseLower) or the upper triangle (if uplo =
 *          MorseUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          MORSE_Alloc_Workspace_zheevd
 *          On exit, contains auxiliary factorization data.
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
 * @sa MORSE_zheevd
 * @sa MORSE_zheevd_Tile
 * @sa MORSE_cheevd_Tile_Async
 * @sa MORSE_dsyev_Tile_Async
 * @sa MORSE_ssyev_Tile_Async
 *
 */
int MORSE_zheevd_Tile_Async( MORSE_enum jobz, MORSE_enum uplo,
                             MORSE_desc_t *A, double *W, MORSE_desc_t *T,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t descA;
    MORSE_desc_t descT;
    MORSE_desc_t D, *Dptr = NULL;
    MORSE_Complex64_t *Q2 = NULL;
    int N, NB, status;
    double *E;
    MORSE_Complex64_t *V;
    MORSE_desc_t descQ2l, descQ2t;
    MORSE_desc_t descVl, descVt;
    MORSE_desc_t *subA, *subQ, *subT;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zheevd_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zheevd_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zheevd_Tile_Async", "NULL request");
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
        morse_error("MORSE_zheevd_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zheevd_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (jobz != MorseNoVec && jobz != MorseVec) {
        morse_error("MORSE_zheevd_Tile_Async", "illegal value of jobz");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != MorseLower) && (uplo != MorseUpper)) {
        morse_error("MORSE_zheevd_Tile_Async", "illegal value of uplo");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        morse_error("MORSE_zheevd_Tile_Async", "matrix need to be square");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        morse_error("MORSE_zheevd_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    N  = descA.m;
    NB = descA.mb;

    /* Allocate data structures for reduction to tridiagonal form */
    E = malloc( (N - 1) * sizeof(double) );
    if (E == NULL) {
        morse_error("MORSE_zheevd_Tile_Async", "malloc(E) failed");
        free(E);
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

    if (jobz == MorseVec){
        /* Have to synchrone right now */
        Q2 = malloc( N * N * sizeof(MORSE_Complex64_t));
        /* For bug in lapacke */
        memset( Q2, 0, N * N * sizeof(MORSE_Complex64_t));
    }

    status = MORSE_zhetrd_Tile_Async( jobz, uplo,
                                      A, W, E, T,
                                      Q2, N,
                                      sequence, request );
    if (status != 0) {
        morse_error("MORSE_zheevd_Tile", "MORSE_zhetrd failed");
    }

    if (jobz == MorseNoVec){
#if !defined(CHAMELEON_SIMULATION)
        /* Tridiagonal eigensolver */
        status = LAPACKE_zstedc( LAPACK_COL_MAJOR,
                                 morse_lapack_const(jobz),
                                 N, W, E,
                                 NULL, N );
        if (status != 0) {
            morse_error("MORSE_zheevd_Tile_Async", "LAPACKE_zstedc failed");
        }
#endif /* !defined(CHAMELEON_SIMULATION) */
        free(E);
        return MORSE_SUCCESS;
    }

    V = malloc( N * N * sizeof(MORSE_Complex64_t) );
    if (V == NULL) {
        morse_error("MORSE_zheevd_Tile_Async", "malloc(V) failed");
        free(V);
        return MORSE_ERR_OUT_OF_RESOURCES;
    }
    /* For bug in lapacke */
    memset(V, 0, N * N * sizeof(MORSE_Complex64_t));

    /*
     * Tridiagonal eigensolver
     * V contains the eigenvectors of the tridiagonal system on exit
     */
#if !defined(CHAMELEON_SIMULATION)
    status = LAPACKE_zstedc( LAPACK_COL_MAJOR,
                             morse_lapack_const(MorseIvec),
                             N, W, E,
                             V, N );
    if (status != 0) {
        morse_error("MORSE_zheevd_Tile_Async", "LAPACKE_zstedc failed");
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    /* Convert matrices in tile format */
    /* A/T from MORSE_zhetrd   refer  to Q1 (tile   layout) */
    /* Q   from MORSE_zhetrd   refers to Q2 (lapack layout) */
    /* V   from LAPACKE_zstedc refers to V  (lapack layout) */
    /* The final eigenvectors are (Q1 Q2 V) or (Q1^h Q2 V)  */
    morse_zlap2tile( morse, &descQ2l, &descQ2t, MorseDescInput, MorseUpperLower,
                     Q2, NB, NB, N, N, N, N, sequence, request );

    morse_zlap2tile( morse, &descVl, &descVt, MorseDescInput, MorseUpperLower,
                     V, NB, NB, N, N, N, N, sequence, request );

    if (uplo == MorseLower)
    {
#if defined(CHAMELEON_COPY_DIAG)
        {
            int n = chameleon_min(A->mt, A->nt) * A->nb;
            morse_zdesc_alloc(D, A->mb, A->nb, A->m, n, 0, 0, A->m, n, );
            Dptr = &D;
        }
#endif
        subA = morse_desc_submatrix(&descA,   descA.mb,   0, descA.m  -descA.mb,   descA.n-descA.nb);
        subQ = morse_desc_submatrix(&descQ2t, descQ2t.mb, 0, descQ2t.m-descQ2t.mb, descQ2t.n );
        subT = morse_desc_submatrix(&descT,   descT.mb,   0, descT.m  -descT.mb,   descT.n-descT.nb);

        /* Compute Q2 = Q1 * Q2 */
        morse_pzunmqr( MorseLeft, MorseNoTrans,
                       subA, subQ, subT, Dptr,
                       sequence, request );

        /* Compute the final eigenvectors A = (Q1 * Q2) * V */
        morse_pzgemm( MorseNoTrans, MorseNoTrans,
                      1.0, &descQ2t, &descVt,
                      0.0, &descA,
                      sequence, request );

    }
    else {
#if defined(CHAMELEON_COPY_DIAG)
        {
            int m = chameleon_min(A->mt, A->nt) * A->mb;
            morse_zdesc_alloc(D, A->mb, A->nb, m, A->n, 0, 0, m, A->n, );
            Dptr = &D;
        }
#endif
        subA = morse_desc_submatrix(&descA,   0,   descA.nb, descA.m  -descA.mb,   descA.n -descA.nb );
        subQ = morse_desc_submatrix(&descQ2t, descQ2t.mb, 0, descQ2t.m-descQ2t.mb, descQ2t.n );
        subT = morse_desc_submatrix(&descT,   0,   descT.nb, descT.m  -descT.mb,   descT.n -descT.nb );

        /* Compute Q2 = Q1^h * Q2 */
        morse_pzunmlq( MorseLeft, MorseConjTrans,
                       subA, subQ, subT, Dptr,
                       sequence, request );

        /* Compute the final eigenvectors A =  (Q1^h * Q2) * V */
        morse_pzgemm( MorseNoTrans, MorseNoTrans,
                      1.0, &descQ2t, &descVt,
                      0.0, &descA,
                      sequence, request );
    }

    morse_ztile2lap( morse, &descQ2l, &descQ2t,
                     MorseDescInput, MorseUpperLower, sequence, request );
    morse_ztile2lap( morse, &descVl, &descVt,
                     MorseDescInput, MorseUpperLower, sequence, request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descQ2l, &descQ2t );
    morse_ztile2lap_cleanup( morse, &descVl, &descVt );

    free(subA); free(subQ); free(subT);
    free(Q2);
    free(V);
    free(E);
    if (Dptr != NULL) {
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return MORSE_SUCCESS;
}
