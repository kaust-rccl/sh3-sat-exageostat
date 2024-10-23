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

      SUBROUTINE ALASUM( TYPE, NOUT, NFAIL, NRUN, NERRS )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3        TYPE
      INTEGER            NFAIL, NOUT, NRUN, NERRS
*     ..
*
*  Purpose
*  =======
*
*  ALASUM prints a summary of results from one of the -CHK- routines.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*3
*          The LAPACK path name.
*
*  NOUT    (input) INTEGER
*          The unit number on which results are to be printed.
*          NOUT >= 0.
*
*  NFAIL   (input) INTEGER
*          The number of tests which did not pass the threshold ratio.
*
*  NRUN    (input) INTEGER
*          The total number of tests.
*
*  NERRS   (input) INTEGER
*          The number of error messages recorded.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      IF( NFAIL.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )TYPE, NFAIL, NRUN
      ELSE
         WRITE( NOUT, FMT = 9998 )TYPE, NRUN
      END IF
      IF( NERRS.GT.0 ) THEN
         WRITE( NOUT, FMT = 9997 )NERRS
      END IF
*
 9999 FORMAT( 1X, A3, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9998 FORMAT( /1X, 'All tests for ', A3,
     $      ' routines passed the threshold (', I6, ' tests run)' )
 9997 FORMAT( 6X, I6, ' error messages recorded' )
      RETURN
*
*     End of ALASUM
*
      END
