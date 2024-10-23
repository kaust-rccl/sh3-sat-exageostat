/**
 *
 * @file zhetrd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhetrd wrappers
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#endif

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zhetrd - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary
 *  Q2).  Let Q = Q1 * Q2 be the global unitary transformation; Q**H *
 *  A * Q = S.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = MorseNoVec: computes tridiagonal only;
 *          = MorseVec: computes tridiagonal and generate the orthogonal matrix Q.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = MorseUpper:: Upper triangle of A is stored;
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
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = MorseUpper, E(i) = A(i+1,i) if uplo = MorseLower.
 *
 * @param[out] descT
 *          On entry, descriptor as return by MORSE_Alloc_Workspace_zhetrd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = MorseVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = MorseNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
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
 * @sa MORSE_zhetrd_Tile
 * @sa MORSE_zhetrd_Tile_Async
 * @sa MORSE_chetrd
 * @sa MORSE_dsytrd
 * @sa MORSE_ssytrd
 *
 */
int MORSE_zhetrd( MORSE_enum jobz, MORSE_enum uplo, int N,
                  MORSE_Complex64_t *A, int LDA,
                  double *D,
                  double *E,
                  MORSE_desc_t *descT,
                  MORSE_Complex64_t *Q, int LDQ )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_zhetrd", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != MorseNoVec && jobz != MorseVec) {
        morse_error("MORSE_zhetrd", "illegal value of jobz");
        return -1;
    }
    if ((uplo != MorseLower) && (uplo != MorseUpper)) {
        morse_error("MORSE_zhetrd", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_zhetrd", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_zhetrd", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (N == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZHETRD, N, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zhetrd", "morse_tune() failed");
        return status;
    }
    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    MORSE_zhetrd_Tile_Async( jobz, uplo, &descAt, D, E, descT, Q, LDQ, sequence, &request );

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
 *  MORSE_zhetrd_Tile - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary Q2).
 *  Let Q = Q1 * Q2 be the global unitary transformation;
 *  Q**H * A * Q = S.
 *  Tile equivalent of MORSE_zhetrd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = MorseNoVec: computes tridiagonal only;
 *          = MorseVec: computes tridiagonal and generate the orthogonal matrix Q.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.  If uplo
 *          = MorseUpper, the leading N-by-N upper triangular part of
 *          A contains the upper triangular part of the matrix A, and
 *          the strictly lower triangular part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = MorseVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = MorseNoVec, then on exit the lower triangle (if
 *          uplo = MorseLower) or the upper triangle (if uplo =
 *          MorseUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = MorseUpper,
 *          E(i) = A(i+1,i) if uplo = MorseLower.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = MorseVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = MorseNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
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
 * @sa MORSE_zhetrd
 * @sa MORSE_zhetrd_Tile_Async
 * @sa MORSE_chetrd_Tile
 * @sa MORSE_dsytrd_Tile
 * @sa MORSE_ssytrd_Tile
 * @sa MORSE_zhetrd_Tile
 *
 */
int MORSE_zhetrd_Tile( MORSE_enum jobz, MORSE_enum uplo,
                       MORSE_desc_t *A, double *D, double *E,
                       MORSE_desc_t *T, MORSE_Complex64_t *Q, int LDQ )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zhetrd_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zhetrd_Tile_Async( jobz, uplo, A, D, E, T, Q, LDQ, sequence, &request );

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
 *  MORSE_zhetrd_Tile_Async - Computes all eigenvalues and,
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
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MORSE_zhetrd
 * @sa MORSE_zhetrd_Tile
 * @sa MORSE_chetrd_Tile_Async
 * @sa MORSE_dsytrd_Tile_Async
 * @sa MORSE_ssytrd_Tile_Async
 *
 */
int MORSE_zhetrd_Tile_Async( MORSE_enum jobz,
                             MORSE_enum uplo,
                             MORSE_desc_t *A,
                             double *W,
                             double *E,
                             MORSE_desc_t *T,
                             MORSE_Complex64_t *Q, int LDQ,
                             MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_desc_t descA;
    MORSE_desc_t descAB;
    int N, NB, LDAB;
    int status;
    MORSE_desc_t D, *Dptr = NULL;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zhetrd_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zhetrd_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zhetrd_Tile_Async", "NULL request");
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
        morse_error("MORSE_zhetrd_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zhetrd_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Check input arguments */
    if (jobz != MorseNoVec && jobz != MorseVec) {
        morse_error("MORSE_zhetrd_Tile_Async", "illegal value of jobz");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != MorseLower) && (uplo != MorseUpper)) {
        morse_error("MORSE_zhetrd_Tile_Async", "illegal value of uplo");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        morse_error("MORSE_zhetrd_Tile_Async", "matrix need to be square");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        morse_error("MORSE_zhetrd_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    N  = descA.m;
    NB = descA.mb;
#if defined(CHAMELEON_COPY_DIAG)
    {
        morse_zdesc_alloc_diag(D, A->mb, A->nb, chameleon_min(A->m, A->n), A->nb, 0, 0, chameleon_min(A->m, A->n), A->nb, A->p, A->q);
        Dptr = &D;
    }
#endif
    /* Reduction to band. On exit, T contains reflectors */
    morse_pzhetrd_he2hb( uplo, A, T, Dptr,
                         sequence, request );

    LDAB = NB+1;

    /* Allocate band structure */
    morse_zdesc_alloc_diag( descAB,
                            LDAB, NB, /* mb, nb */
                            LDAB, N,  /* lm, ln */
                            0, 0,     /* i, j */
                            LDAB, N,  /* m, n */
                            1, 1 );

    /* Copy data into band structure */
    morse_pztile2band( uplo, A, &descAB,
                       sequence, request );

    morse_sequence_wait( morse, sequence );

    /* Reduce band matrix to tridiagonal matrix */
#if !defined(CHAMELEON_SIMULATION)
    status = LAPACKE_zhbtrd( LAPACK_COL_MAJOR,
                             morse_lapack_const(jobz),
                             morse_lapack_const(uplo),
                             N, NB,
                             (MORSE_Complex64_t *) descAB.mat, LDAB,
                             W, E, Q, LDQ );
    if (status != 0) {
        morse_error("MORSE_zhetrd_Tile_Async", "LAPACKE_zhbtrd failed");
    }
#endif /* !defined(CHAMELEON_SIMULATION) */
    if (Dptr != NULL) {
        morse_desc_mat_free( Dptr );
    }
    morse_desc_mat_free( &descAB );
    (void)D;
    return MORSE_SUCCESS;
}
