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

      SUBROUTINE DQRT01( M, N, A, AF, Q, R, LDA, T, WORK, LWORK,
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
      DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
     $                   R( LDA, * ), RESULT( * ), RWORK( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DQRT01 tests DGEQRF, which computes the QR factorization of an m-by-n
*  matrix A, and partially tests DORGQR which forms the m-by-m
*  orthogonal matrix Q.
*
*  DQRT01 compares R with Q'*A, and checks that Q is orthogonal.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the QR factorization of A, as returned by DGEQRF.
*          See DGEQRF for further details.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDA,M)
*          The m-by-m orthogonal matrix Q.
*
*  R       (workspace) DOUBLE PRECISION array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by DGEQRF.
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
*          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      INTEGER            INFO, MINMN
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEQRF, DLACPY, DLASET, DORGQR, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
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
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
*, 
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'DGEQRF'
      CALL MORSE_DGEQRF( M, N, AF, LDA, T, INFO )
*
*     Copy details of Q
*
      CALL DLASET( 'Full', M, N, ZERO, ONE, Q, LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'DORGQR'
      CALL MORSE_DORGQR( M, N, MINMN, AF, LDA, T, Q, LDA, INFO )
*
*     Copy R
*
      CALL DLASET( 'Full', N, N, ZERO, ZERO, R, LDA )
      CALL DLACPY( 'Upper', M, N, AF, LDA, R, LDA )
*
*     Compute R - Q'*A
*
      CALL DGEMM( 'Transpose', 'No transpose', N, N, M, -ONE, Q, LDA,
     $           A, LDA, ONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = DLANGE( '1', M, N, A, LDA, RWORK )
      RESID = DLANGE( '1', N, N, R, LDA, RWORK )
*
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'Transpose', N, M, ONE, Q, LDA, -ONE, R,
     $            LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of DQRT01
*
      END
