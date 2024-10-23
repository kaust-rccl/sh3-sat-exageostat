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

      SUBROUTINE SRSCL( N, SA, SX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SA
*     ..
*     .. Array Arguments ..
      REAL               SX( * )
*     ..
*
*  Purpose
*  =======
*
*  SRSCL multiplies an n-element real vector x by the real scalar 1/a.
*  This is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of components of the vector x.
*
*  SA      (input) REAL
*          The scalar a which is used to divide each component of x.
*          SA must be >= 0, or the subroutine will divide by zero.
*
*  SX      (input/output) REAL array, dimension
*                         (1+(N-1)*abs(INCX))
*          The n-element vector x.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector SX.
*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      REAL               BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLABAD, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply X by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector X by MUL
*
      CALL SSCAL( N, MUL, SX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of SRSCL
*
      END
