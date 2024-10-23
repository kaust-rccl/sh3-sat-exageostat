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

      SUBROUTINE SGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  Purpose
*  =======
*
*  SGET04 computes the difference between a computed solution and the
*  true solution to a system of linear equations.
*
*  RESID =  ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
*  where RCOND is the reciprocal of the condition number and EPS is the
*  machine epsilon.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrices X and XACT.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of the matrices X and XACT.  NRHS >= 0.
*
*  X       (input) REAL array, dimension (LDX,NRHS)
*          The computed solution vectors.  Each vector is stored as a
*          column of the matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  XACT    (input) REAL array, dimension( LDX, NRHS )
*          The exact solution vectors.  Each vector is stored as a
*          column of the matrix XACT.
*
*  LDXACT  (input) INTEGER
*          The leading dimension of the array XACT.  LDXACT >= max(1,N).
*
*  RCOND   (input) REAL
*          The reciprocal of the condition number of the coefficient
*          matrix in the system of equations.
*
*  RESID   (output) REAL
*          The maximum over the NRHS solution vectors of
*          ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      REAL               DIFFNM, EPS, XNORM
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH
      EXTERNAL           ISAMAX, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      EPS = SLAMCH( 'Epsilon' )
      IF( RCOND.LT.ZERO ) THEN
         RESID = 1.0 / EPS
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      RESID = ZERO
      DO 20 J = 1, NRHS
         IX = ISAMAX( N, XACT( 1, J ), 1 )
         XNORM = ABS( XACT( IX, J ) )
         DIFFNM = ZERO
         DO 10 I = 1, N
            DIFFNM = MAX( DIFFNM, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
         IF( XNORM.LE.ZERO ) THEN
            IF( DIFFNM.GT.ZERO )
     $         RESID = 1.0 / EPS
         ELSE
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         END IF
   20 CONTINUE
      IF( RESID*EPS.LT.1.0 )
     $   RESID = RESID / EPS
*
      RETURN
*
*     End of SGET04
*
      END
