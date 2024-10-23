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

      SUBROUTINE CQRT01( M, N, A, AF, Q, R, LDA, T, WORK, LWORK,
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
      COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
     $                   R( LDA, * ),  WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CQRT01 tests CGEQRF, which computes the QR factorization of an m-by-n
*  matrix A, and partially tests CUNGQR which forms the m-by-m
*  orthogonal matrix Q.
*
*  CQRT01 compares R with Q'*A, and checks that Q is orthogonal.
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
*          Details of the QR factorization of A, as returned by CGEQRF.
*          See CGEQRF for further details.
*
*  Q       (output) COMPLEX array, dimension (LDA,M)
*          The m-by-m orthogonal matrix Q.
*
*  R       (workspace) COMPLEX array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by CGEQRF.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      EXTERNAL           CGEMM, CGEQRF, CHERK, CLACPY, CLASET, CUNGQR
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
      SRNAMT = 'CGEQRF'
      CALL MORSE_CGEQRF( M, N, AF, LDA, T, INFO )
*
*     Copy details of Q
*
      CALL CLASET( 'Full', M, N, CMPLX(ZERO), CMPLX(ONE), Q, LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'CUNGQR'
      CALL MORSE_CUNGQR( M, N, MINMN, AF, LDA, T, Q, LDA, INFO )
*
*     Copy R
*
      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA )
      CALL CLACPY( 'Upper', N, N, AF, LDA, R, LDA )
*
*     Compute R - Q'*A
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', N, N, M,
     $            CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', N, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', N, M, ONE, Q, LDA,
     $            -ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of CQRT01
*
      END
