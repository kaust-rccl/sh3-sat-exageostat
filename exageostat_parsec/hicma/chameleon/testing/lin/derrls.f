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

      SUBROUTINE DERRLS( PATH, NUNIT )
*
      INCLUDE 'morse_fortran.h'
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  DERRLS tests the error exits for the DOUBLE PRECISION least squares
*  driver routines (DGELS, SGELSS, SGELSX, SGELSY, SGELSD).
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
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            INFO, IRNK
      DOUBLE PRECISION   RCOND
      INTEGER            HT( 2 )
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ),
     $                   W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELSX,
     $                   DGELSY
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
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.0D+0
      A( 1, 2 ) = 2.0D+0
      A( 2, 2 ) = 3.0D+0
      A( 2, 1 ) = 4.0D+0
      OK = .TRUE.
*
*     Disable MORSE warnings/errors
* 
      CALL MORSE_DISABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_DISABLE( MORSE_ERRORS,   INFO )
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        Test error exits for the least squares driver routines.
*
*        DGELS
*
         CALL MORSE_ALLOC_WORKSPACE_DGELS( 2, 2, HT, INFO )
*
         SRNAMT = 'DGELS '
         INFOT = 103
         CALL MORSE_DGELS( '/', 0, 0, 0, A, 1, HT, B, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
         INFOT = 2
         CALL MORSE_DGELS( MORSENOTRANS, -1, 0, 0, A, 1, HT,
     $                     B, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
         INFOT = 3
         CALL MORSE_DGELS( MORSENOTRANS, 0, -1, 0, A, 1, HT,
     $                     B, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
         INFOT = 4
         CALL MORSE_DGELS( MORSENOTRANS, 0, 0, -1, A, 1, HT,
     $                     B, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
         INFOT = 6
         CALL MORSE_DGELS( MORSENOTRANS, 2, 0, 0, A, 1, HT,
     $                     B, 2, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
         INFOT = 9
         CALL MORSE_DGELS( MORSENOTRANS, 2, 0, 0, A, 2, HT,
     $                     B, 1, INFO )
         CALL CHKXER( 'DGELS ', INFOT, NOUT, INFO, OK )
*
         CALL MORSE_DEALLOC_HANDLE( HT, INFO )
*
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
*     Enable MORSE warnings/errors
* 
      CALL MORSE_ENABLE( MORSE_WARNINGS, INFO )
      CALL MORSE_ENABLE( MORSE_ERRORS,   INFO )
*
      RETURN
*
*     End of DERRLS
*
      END
