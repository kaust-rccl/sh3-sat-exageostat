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

      SUBROUTINE CLQT01( M, N, A, AF, Q, L, LDA, T, WORK, LWORK,
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
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CLQT01 tests CGELQF, which computes the LQ factorization of an m-by-n
*  matrix A, and partially tests CUNGLQ which forms the n-by-n
*  orthogonal matrix Q.
*
*  CLQT01 compares L with A*Q', and checks that Q is orthogonal.
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
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) COMPLEX array, dimension (LDA,N)
*          Details of the LQ factorization of A, as returned by CGELQF.
*          See CGELQF for further details.
*
*  Q       (output) COMPLEX array, dimension (LDA,N)
*          The n-by-n orthogonal matrix Q.
*
*  L       (workspace) COMPLEX array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and L.
*          LDA >= max(M,N).
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by CGELQF.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (max(M,N))
*
*  RESULT  (output) REAL array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( L - A*Q' ) / ( N * norm(A) * EPS )
*          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, MINMN
      REAL               ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      EXTERNAL           CLANGE, CLANSY, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGELQF, CGEMM, CHERK, CLACPY, CLASET, CUNGLQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN, REAL
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
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'CGELQF'
      CALL MORSE_CGELQF( M, N, AF, LDA, T, INFO )
*
*     Copy details of Q
*
      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ONE ), Q, LDA )
*
*     Generate the n-by-n matrix Q
*
      SRNAMT = 'CUNGLQ'
      CALL MORSE_CUNGLQ( M, N, MINMN, AF, LDA, T, Q, LDA, INFO )
*
*     Copy L
*
      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA )
      CALL CLACPY( 'Lower', M, N, AF, LDA, L, LDA )
*
*     Compute L - A*Q'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', M, M, N,
     $            CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA )
*
*     Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), L, LDA )
      CALL CHERK( 'Upper', 'No transpose', M, N, ONE, Q, LDA, -ONE, L,
     $            LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of CLQT01
*
      END
