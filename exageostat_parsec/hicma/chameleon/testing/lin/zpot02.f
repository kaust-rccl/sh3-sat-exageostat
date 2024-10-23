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

      SUBROUTINE ZPOT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOT02 computes the residual for the solution of a Hermitian system
*  of linear equations  A*x = b:
*
*   RESID = norm( B - A*X ) / ( norm(A) * norm(X) + norm(RHS))* N * EPS )
*
*  where EPS is the machine epsilon.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of B, the matrix of right hand sides.
*          NRHS >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The original Hermitian matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N)
*
*  X       (input) COMPLEX*16 array, dimension (LDX,NRHS)
*          The computed solution vectors for the system of linear
*          equations.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.   LDX >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the system of
*          linear equations.
*          On exit, B is overwritten with the difference B - A*X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          The maximum over the number of right hand sides of
*    norm( B - A*X ) / ( norm(A) * norm(X) + norm(RHS))* N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM, RHSNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANHE, CLANGE
      EXTERNAL           DLAMCH, DZASUM, ZLANHE, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      RHSNORM = CLANGE( '1', N, NRHS, B, LDB, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X
*
      CALL ZHEMM( 'Left', UPLO, N, NRHS, -CONE, A, LDA, X, LDX, CONE, B,
     $            LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DZASUM( N, B( 1, J ), 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( BNORM) / ((ANORM * XNORM + RHSNORM)* 
     $   N *EPS ))
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of ZPOT02
*
      END
