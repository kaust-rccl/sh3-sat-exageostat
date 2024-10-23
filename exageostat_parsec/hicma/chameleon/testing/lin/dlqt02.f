!!!
!
! -- Inria
! -- (C) Copyright 2012
!
! This software is a computer program whose purpose is to process
! Matrices Over Runtime Systems @ Exascale (MORSE). More information
! can be found on the following website: http://www.inria.fr/en/teams/morse.
! 
! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.
!
!!!

      SUBROUTINE DLQT02( M, N, K, A, AF, Q, L, LDA, T, WORK, LWORK,
     $                   RWORK, RESULT )
*
      INCLUDE 'morse_fortran.h'
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M, N
      INTEGER            T( 2 )
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), RESULT( * ), RWORK( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DLQT02 tests DORGLQ, which generates an m-by-n matrix Q with
*  orthonornmal rows that is defined as the product of k elementary
*  reflectors.
*
*  Given the LQ factorization of an m-by-n matrix A, DLQT02 generates
*  the orthogonal matrix Q defined by the factorization of the first k
*  rows of A; it compares L(1:k,1:m) with A(1:k,1:n)*Q(1:m,1:n)', and
*  checks that the rows of Q are orthonormal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q to be generated.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q to be generated.
*          N >= M >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A which was factorized by DLQT01.
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the LQ factorization of A, as returned by DGELQF.
*          See DGELQF for further details.
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (LDA,N)
*
*  L       (workspace) DOUBLE PRECISION array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L. LDA >= N.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (M)
*          The scalar factors of the elementary reflectors corresponding
*          to the LQ factorization in AF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( L - A*Q' ) / ( N * norm(A) * EPS )
*          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ROGUE
      PARAMETER          ( ROGUE = -1.0D+10 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, DORGLQ, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the first k rows of the factorization to the array Q
*
      CALL DLASET( 'Full', M, N, ZERO, ONE, Q, LDA )
*
*     Generate the first n columns of the matrix Q
*
      SRNAMT = 'DORGLQ'
      CALL MORSE_DORGLQ( M, N, K, AF, LDA, T, Q, LDA, INFO )
*
*     Copy L(1:k,1:m)
*
      CALL DLASET( 'Full', K, M, ZERO, ZERO, L, LDA )
      CALL DLACPY( 'Lower', K, M, AF, LDA, L, LDA )
*
*     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
*
      CALL DGEMM( 'No transpose', 'Transpose', K, M, N, -ONE, A, LDA, Q,
     $            LDA, ONE, L, LDA )
*
*     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', K, N, A, LDA, RWORK )
      RESID = DLANGE( '1', K, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL DLASET( 'Full', M, M, ZERO, ONE, L, LDA )
      CALL DSYRK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of DLQT02
*
      END
