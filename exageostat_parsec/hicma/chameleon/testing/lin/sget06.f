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

      REAL             FUNCTION SGET06( RCOND, RCONDC )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL               RCOND, RCONDC
*     ..
*
*  Purpose
*  =======
*
*  SGET06 computes a test ratio to compare two values for RCOND.
*
*  Arguments
*  ==========
*
*  RCOND   (input) REAL
*          The estimate of the reciprocal of the condition number of A,
*          as computed by SGECON.
*
*  RCONDC  (input) REAL
*          The reciprocal of the condition number of A, computed as
*          ( 1/norm(A) ) / norm(inv(A)).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      REAL               EPS, RAT
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      IF( RCOND.GT.ZERO ) THEN
         IF( RCONDC.GT.ZERO ) THEN
            RAT = MAX( RCOND, RCONDC ) / MIN( RCOND, RCONDC ) -
     $            ( ONE-EPS )
         ELSE
            RAT = RCOND / EPS
         END IF
      ELSE
         IF( RCONDC.GT.ZERO ) THEN
            RAT = RCONDC / EPS
         ELSE
            RAT = ZERO
         END IF
      END IF
      SGET06 = RAT
      RETURN
*
*     End of SGET06
*
      END
