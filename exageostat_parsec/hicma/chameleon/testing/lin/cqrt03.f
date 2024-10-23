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

      SUBROUTINE CQRT03( M, N, K, AF, C, CC, Q, LDA, T, WORK, LWORK,
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
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
     $                   Q( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CQRT03 tests CUNMQR, which computes Q*C, Q'*C, C*Q or C*Q'.
*
*  CQRT03 compares the results of a call to CUNMQR with the results of
*  forming Q explicitly by a call to CUNGQR and then performing matrix
*  multiplication by a call to CGEMM.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The order of the orthogonal matrix Q.  M >= 0.
*
*  N       (input) INTEGER
*          The number of rows or columns of the matrix C; C is m-by-n if
*          Q is applied from the left, or n-by-m if Q is applied from
*          the right.  N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          orthogonal matrix Q.  M >= K >= 0.
*
*  AF      (input) COMPLEX array, dimension (LDA,N)
*          Details of the QR factorization of an m-by-n matrix, as
*          returnedby CGEQRF. See CGEQRF for further details.
*
*  C       (workspace) COMPLEX array, dimension (LDA,N)
*
*  CC      (workspace) COMPLEX array, dimension (LDA,N)
*
*  Q       (workspace) COMPLEX array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays AF, C, CC, and Q.
*
*  TAU     (input) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors corresponding
*          to the QR factorization in AF.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK must be at least M, and should be
*          M*NB, where NB is the blocksize for this environment.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (4)
*          The test ratios compare two techniques for multiplying a
*          random matrix C by an m-by-m orthogonal matrix Q.
*          RESULT(1) = norm( Q*C - Q*C )  / ( M * norm(C) * EPS )
*          RESULT(2) = norm( C*Q - C*Q )  / ( M * norm(C) * EPS )
*          RESULT(3) = norm( Q'*C - Q'*C )/ ( M * norm(C) * EPS )
*          RESULT(4) = norm( C*Q' - C*Q' )/ ( M * norm(C) * EPS )
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
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, ISIDE, ITRANS, J, MC, NC
      INTEGER            MORSE_SIDE, MORSE_TRANS
      REAL               CNORM, EPS, RESID
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLARNV, CLASET, CUNGQR, CUNMQR
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, REAL
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEED / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      WORK(1) = ONE
*
*     Copy the first k columns of the factorization to the array Q
*
      IF ( K.EQ.0 ) THEN
          CALL CLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      ELSE
          CALL CLASET( 'Full', M, M, CMPLX(ZERO), CMPLX(ONE), Q, LDA )
      ENDIF
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'CUNGQR'
      CALL MORSE_CUNGQR( M, M, K, AF, LDA, T, Q, LDA, INFO )
*
      DO 30 ISIDE = 1, 2
         IF( ISIDE.EQ.1 ) THEN
            SIDE = 'L'
            MORSE_SIDE = MORSELEFT
            MC = M
            NC = N
         ELSE
            SIDE = 'R'
            MORSE_SIDE = MORSERIGHT
            MC = N
            NC = M
         END IF
*
*        Generate MC by NC matrix C
*
         DO 10 J = 1, NC
            CALL CLARNV( 2, ISEED, MC, C( 1, J ) )
   10    CONTINUE
         CNORM = CLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.ZERO )
     $      CNORM = ONE
*
         DO 20 ITRANS = 1, 2
            IF( ITRANS.EQ.1 ) THEN
               TRANS = 'N'
               MORSE_TRANS = MORSENOTRANS
            ELSE
               TRANS = 'C'
               MORSE_TRANS = MORSECONJTRANS
            END IF
*
*           Copy C
*
            CALL CLACPY( 'Full', MC, NC, C, LDA, CC, LDA )
*
*           Apply Q or Q' to C
*
            SRNAMT = 'CUNMQR'
            CALL MORSE_CUNMQR( MORSE_SIDE, MORSE_TRANS, MC, NC, K,
     $                     AF, LDA, T, CC, LDA, INFO )
*
*           Form explicit product and subtract
*
            IF ( K.EQ.0 ) THEN
               CALL CLASET( 'Full', M, M, CMPLX(ZERO), CMPLX(ONE), 
     $              Q, LDA )
            ENDIF
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL CGEMM( TRANS, 'No transpose', MC, NC, MC,
     $                     CMPLX( -ONE ), Q, LDA, C, LDA, CMPLX( ONE ),
     $                     CC, LDA )
            ELSE
               CALL CGEMM( 'No transpose', TRANS, MC, NC, NC,
     $                     CMPLX( -ONE ), C, LDA, Q, LDA, CMPLX( ONE ),
     $                     CC, LDA )
            END IF
*
*           Compute error in the difference
*
            RESID = CLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID /
     $         ( REAL( MAX( 1, M ) )*CNORM*EPS )
*
   20    CONTINUE
   30 CONTINUE
*
      RETURN
*
*     End of CQRT03
*
      END
