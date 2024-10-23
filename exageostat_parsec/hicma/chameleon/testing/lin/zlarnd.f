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

      DOUBLE COMPLEX   FUNCTION ZLARND( IDIST, ISEED )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IDIST
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  ZLARND returns a random complex number from a uniform or normal
*  distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  real and imaginary parts each uniform (0,1)
*          = 2:  real and imaginary parts each uniform (-1,1)
*          = 3:  real and imaginary parts each normal (0,1)
*          = 4:  uniformly distributed on the disc abs(z) <= 1
*          = 5:  uniformly distributed on the circle abs(z) = 1
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine DLARAN to generate a random
*  real number from a uniform (0,1) distribution. The Box-Muller method
*  is used to transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARAN
      EXTERNAL           DLARAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, EXP, LOG, SQRT
*     ..
*     .. Executable Statements ..
*
*     Generate a pair of real random numbers from a uniform (0,1)
*     distribution
*
      T1 = DLARAN( ISEED )
      T2 = DLARAN( ISEED )
*
      IF( IDIST.EQ.1 ) THEN
*
*        real and imaginary parts each uniform (0,1)
*
         ZLARND = DCMPLX( T1, T2 )
      ELSE IF( IDIST.EQ.2 ) THEN
*
*        real and imaginary parts each uniform (-1,1)
*
         ZLARND = DCMPLX( TWO*T1-ONE, TWO*T2-ONE )
      ELSE IF( IDIST.EQ.3 ) THEN
*
*        real and imaginary parts each normal (0,1)
*
         ZLARND = SQRT( -TWO*LOG( T1 ) )*EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.4 ) THEN
*
*        uniform distribution on the unit disc abs(z) <= 1
*
         ZLARND = SQRT( T1 )*EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.5 ) THEN
*
*        uniform distribution on the unit circle abs(z) = 1
*
         ZLARND = EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      END IF
      RETURN
*
*     End of ZLARND
*
      END
