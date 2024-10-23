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

      SUBROUTINE CLAIPD( N, A, INDA, VINDA )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INDA, N, VINDA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  CLAIPD sets the imaginary part of the diagonal elements of a complex
*  matrix A to a large value.  This is used to test LAPACK routines for
*  complex Hermitian matrices, which are not supposed to access or use
*  the imaginary parts of the diagonals.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The number of diagonal elements of A.
*
*  A      (input/output) COMPLEX array, dimension
*                        (1+(N-1)*INDA+(N-2)*VINDA)
*         On entry, the complex (Hermitian) matrix A.
*         On exit, the imaginary parts of the diagonal elements are set
*         to BIGNUM = EPS / SAFMIN, where EPS is the machine epsilon and
*         SAFMIN is the safe minimum.
*
*  INDA   (input) INTEGER
*         The increment between A(1) and the next diagonal element of A.
*         Typical values are
*         = LDA+1:  square matrices with leading dimension LDA
*         = 2:  packed upper triangular matrix, starting at A(1,1)
*         = N:  packed lower triangular matrix, starting at A(1,1)
*
*  VINDA  (input) INTEGER
*         The change in the diagonal increment between columns of A.
*         Typical values are
*         = 0:  no change, the row and column increments in A are fixed
*         = 1:  packed upper triangular matrix
*         = -1:  packed lower triangular matrix
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IA, IXA
      REAL               BIGNUM
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      BIGNUM = SLAMCH( 'Epsilon' ) / SLAMCH( 'Safe minimum' )
      IA = 1
      IXA = INDA
      DO 10 I = 1, N
         A( IA ) = CMPLX( REAL( A( IA ) ), BIGNUM )
         IA = IA + IXA
         IXA = IXA + VINDA
   10 CONTINUE
      RETURN
      END
