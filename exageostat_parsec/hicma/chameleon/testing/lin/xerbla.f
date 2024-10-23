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

      SUBROUTINE XERBLA( SRNAME, INFO )
*     
*     -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*     
*     Purpose
*     =======
*     
*     This is a special version of XERBLA to be used only as part of
*     the test program for testing error exits from the LAPACK routines.
*     Error messages are printed if INFO.NE.INFOT or if SRNAME.NE.SRMANT,
*     where INFOT and SRNAMT are values stored in COMMON.
*     
*     Arguments
*     =========
*     
*     SRNAME  (input) CHARACTER*(*)
*     The name of the subroutine calling XERBLA.  This name should
*     match the COMMON variable SRNAMT.
*     
*     INFO    (input) INTEGER
*     The error return code from the calling subroutine.  INFO
*     should equal the COMMON variable INFOT.
*     
*     Further Details
*     ======= =======
*     
*     The following variables are passed via the common blocks INFOC and
*     SRNAMC:
*     
*     INFOT   INTEGER      Expected integer return code
*     NOUT    INTEGER      Unit number for printing error messages
*     OK      LOGICAL      Set to .TRUE. if INFO = INFOT and
*     SRNAME = SRNAMT, otherwise set to .FALSE.
*     LERR    LOGICAL      Set to .TRUE., indicating that XERBLA was called
*     SRNAMT  CHARACTER*(*) Expected name of calling subroutine
*     
*     
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*     
      LERR = .TRUE.
      IF( INFO.NE.INFOT ) THEN
         IF( INFOT.NE.0 ) THEN
            WRITE( NOUT, FMT = 9999 )
     $           SRNAMT( 1:LEN_TRIM( SRNAMT ) ), INFO, INFOT
         ELSE
            WRITE( NOUT, FMT = 9997 )
     $           SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
         END IF
         OK = .FALSE.
      END IF
      IF( SRNAME.NE.SRNAMT ) THEN
         WRITE( NOUT, FMT = 9998 )
     $        SRNAME( 1:LEN_TRIM( SRNAME ) ),
     $        SRNAMT( 1:LEN_TRIM( SRNAMT ) )
         OK = .FALSE.
      END IF
      RETURN
*     
 9999 FORMAT( ' *** XERBLA was called from ', A, ' with INFO = ', I6,
     $     ' instead of ', I2, ' ***' )
 9998 FORMAT( ' *** XERBLA was called with SRNAME = ', A,
     $     ' instead of ', A6, ' ***' )
 9997 FORMAT( ' *** On entry to ', A, ' parameter number ', I6,
     $     ' had an illegal value ***' )
*     
*     End of XERBLA
*     
      END
