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

      SUBROUTINE ZLQT01( M, N, A, AF, Q, L, LDA, T, WORK, LWORK,
     $                   RWORK, RESULT )
*
      INCLUDE 'morse_fortran.h'
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
      INTEGER            T( 2 )
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RESULT( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZLQT01 tests ZGELQF, which computes the LQ factorization of an m-by-n
*  matrix A, and partially tests ZUNGLQ which forms the n-by-n
*  orthogonal matrix Q.
*
*  ZLQT01 compares L with A*Q', and checks that Q is orthogonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) COMPLEX*16 array, dimension (LDA,N)
*          Details of the LQ factorization of A, as returned by ZGELQF.
*          See ZGELQF for further details.
*
*  Q       (output) COMPLEX*16 array, dimension (LDA,N)
*          The n-by-n orthogonal matrix Q.
*
*  L       (workspace) COMPLEX*16 array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L.
*          LDA >= max(M,N).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by ZGELQF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(M,N))
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
      COMPLEX*16         ROGUE
      PARAMETER          ( ROGUE = ( -1.0D+10, -1.0D+10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, MINMN
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGELQF, ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGLQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'ZGELQF'
      CALL MORSE_ZGELQF( M, N, AF, LDA, T, INFO )
*
*     Copy details of Q
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ONE ), Q, LDA )
*
*     Generate the n-by-n matrix Q
*
      SRNAMT = 'ZUNGLQ'
      CALL MORSE_ZUNGLQ( M, N, MINMN, AF, LDA, T, Q, LDA, INFO )
*
*     Copy L
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ZERO ), L,
     $             LDA )
      CALL ZLACPY( 'Lower', M, N, AF, LDA, L, LDA )
*
*     Compute L - A*Q'
*
      CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, M, N,
     $            DCMPLX( -ONE ), A, LDA, Q, LDA, DCMPLX( ONE ), L,
     $            LDA )
*
*     Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      RESID = ZLANGE( '1', M, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA )
      CALL ZHERK( 'Upper', 'No transpose', M, N, ONE, Q, LDA, -ONE, L,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = ZLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of ZLQT01
*
      END
