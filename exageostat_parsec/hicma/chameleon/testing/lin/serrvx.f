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

      SUBROUTINE SERRVX( PATH, NUNIT )
*
      INCLUDE 'morse_fortran.h'
*
*  -- LAPACK test routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  SERRVX tests the error exits for the REAL driver routines
*  for solving linear systems of equations.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name for the routines to be tested.
*
*  NUNIT   (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          EQ
      CHARACTER*2        C2
      INTEGER            I, INFO, J
      REAL               RCOND
*     ..
*     .. Local Arrays ..
      INTEGER            HL( 2 ), HPIV( 2 )
      INTEGER            IP( NMAX ), IW( NMAX )
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ),
     $                   W( 2*NMAX ), X( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, SGBSV, SGBSVX, SGESV, SGESVX, SGTSV,
     $                   SGTSVX, SPBSV, SPBSVX, SPOSV, SPOSVX, SPPSV,
     $                   SPPSVX, SPTSV, SPTSVX, SSPSV, SSPSVX, SSYSV,
     $                   SSYSVX
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Disable MORSE warnings/errors
* 
      CALL MORSE_DISABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_DISABLE( MORSE_ERRORS,   INFO )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
            AF( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         B( J ) = 0.
         R1( J ) = 0.
         R2( J ) = 0.
         W( J ) = 0.
         X( J ) = 0.
         C( J ) = 0.
         R( J ) = 0.
         IP( J ) = J
   20 CONTINUE
      EQ = ' '
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'GE' ) ) THEN
*
*        ALLOCATE HL and HPIV
*
         CALL MORSE_ALLOC_WORKSPACE_SGETRF_INCPIV(
     $        2, 1, HL, HPIV, INFO )
*
*        MORSE_SGESV
*
         SRNAMT = 'SGESV '
         INFOT = 1
         CALL MORSE_SGESV_INCPIV( -1, 0, A, 1, HL, HPIV, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 2
         CALL MORSE_SGESV_INCPIV( 0, -1, A, 1, HL, HPIV, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 4
         CALL MORSE_SGESV_INCPIV( 2, 1, A, 1, HL, HPIV, B, 2, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 8
         CALL MORSE_SGESV_INCPIV( 2, 1, A, 2, HL, HPIV, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
*
*        DEALLOCATE HL and HPIV
*
         CALL MORSE_DEALLOC_HANDLE( HL, INFO )
         CALL MORSE_DEALLOC_HANDLE( HPIV, INFO )
*
*
*        SGESV
*
         SRNAMT = 'SGESV '
         INFOT = 1
         CALL MORSE_SGESV( -1, 0, A, 1, IWORK, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 2
         CALL MORSE_SGESV( 0, -1, A, 1, IWORK, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 4
         CALL MORSE_SGESV( 2, 1, A, 1, IWORK, B, 2, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
         INFOT = 7
         CALL MORSE_SGESV( 2, 1, A, 2, IWORK, B, 1, INFO )
         CALL CHKXER( 'SGESV ', INFOT, NOUT, INFO, OK )
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN
*
*        SPOSV
*
         SRNAMT = 'SPOSV '
         INFOT = 1
         CALL MORSE_SPOSV( '/', 0, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'SPOSV ', INFOT, NOUT, INFO, OK )
         INFOT = 2
         CALL MORSE_SPOSV( MORSEUPPER, -1, 0, A, 1, B, 1, INFO )
         CALL CHKXER( 'SPOSV ', INFOT, NOUT, INFO, OK )
         INFOT = 3
         CALL MORSE_SPOSV( MORSEUPPER, 0, -1, A, 1, B, 1, INFO )
         CALL CHKXER( 'SPOSV ', INFOT, NOUT, INFO, OK )
         INFOT = 5
         CALL MORSE_SPOSV( MORSEUPPER, 2, 0, A, 1, B, 2, INFO )
         CALL CHKXER( 'SPOSV ', INFOT, NOUT, INFO, OK )
         INFOT = 7
         CALL MORSE_SPOSV( MORSEUPPER, 2, 0, A, 2, B, 1, INFO )
         CALL CHKXER( 'SPOSV ', INFOT, NOUT, INFO, OK )
*
*        SPOSVX
*
         SRNAMT = 'SPOSVX'
         INFOT = 1
         CALL SPOSVX( '/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 2
         CALL SPOSVX( 'N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 3
         CALL SPOSVX( 'N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 4
         CALL SPOSVX( 'N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 6
         CALL SPOSVX( 'N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 8
         CALL SPOSVX( 'N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 9
         EQ = '/'
         CALL SPOSVX( 'F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 10
         EQ = 'Y'
         CALL SPOSVX( 'F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 12
         CALL SPOSVX( 'N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
         INFOT = 14
         CALL SPOSVX( 'N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1,
     $                RCOND, R1, R2, W, IW, INFO )
         CALL CHKXER( 'SPOSVX', INFOT, NOUT, INFO, OK )
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' drivers passed the tests of the error exits' )
 9998 FORMAT( ' *** ', A3, ' drivers failed the tests of the error ',
     $      'exits ***' )
*
*     Enable MORSE warnings/errors
* 
      CALL MORSE_ENABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_ENABLE( MORSE_ERRORS,   INFO )
*
      RETURN
*
*     End of SERRVX
*
      END
