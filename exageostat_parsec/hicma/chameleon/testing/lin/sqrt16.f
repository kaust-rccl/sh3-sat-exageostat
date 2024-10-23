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

      SUBROUTINE SQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB,
     $                   RWORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDA, LDB, LDX, M, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), RWORK( * ),
     $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  SQRT16 computes the residual for a solution of a system of linear
*  equations  A*x = b  or  A'*x = b:
*     RESID = norm(B - A*X) / ( max(m,n) * (norm(A) * norm(X)+ norm(NRHS)) * EPS ),
*  where EPS is the machine epsilon.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A *x = b
*          = 'T':  A'*x = b, where A' is the transpose of A
*          = 'C':  A'*x = b, where A' is the transpose of A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of B, the matrix of right hand sides.
*          NRHS >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The original M x N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  X       (input) REAL array, dimension (LDX,NRHS)
*          The computed solution vectors for the system of linear
*          equations.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  If TRANS = 'N',
*          LDX >= max(1,N); if TRANS = 'T' or 'C', LDX >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the system of
*          linear equations.
*          On exit, B is overwritten with the difference B - A*X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  IF TRANS = 'N',
*          LDB >= max(1,M); if TRANS = 'T' or 'C', LDB >= max(1,N).
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESID   (output) REAL
*          The maximum over the number of right hand sides of
*          norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, N1, N2
      REAL               ANORM, BNORM, EPS, XNORM, RHSNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SASUM, SLAMCH, SLANGE
      EXTERNAL           LSAME, SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if M = 0 or N = 0 or NRHS = 0
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         ANORM = SLANGE( 'I', M, N, A, LDA, RWORK )
         N1 = N
         N2 = M
      ELSE
         ANORM = SLANGE( '1', M, N, A, LDA, RWORK )
         N1 = M
         N2 = N
      END IF
*
      RHSNORM = SLANGE( 'I', N, NRHS, B, LDB, RWORK )
      EPS = SLAMCH( 'Epsilon' )
*
*     Compute  B - A*X  (or  B - A'*X ) and store in B.
*
      CALL SGEMM( TRANS, 'No transpose', N1, NRHS, N2, -ONE, A, LDA, X,
     $            LDX, ONE, B, LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = SASUM( N1, B( 1, J ), 1 )
         XNORM = SASUM( N2, X( 1, J ), 1 )
         IF( ANORM.EQ.ZERO .AND. BNORM.EQ.ZERO ) THEN
            RESID = ZERO
         ELSE IF( ANORM.LE.ZERO .OR. XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, (  BNORM) / (ANORM * XNORM + RHSNORM ) *
     $              ( MAX( M, N )*EPS ) )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of SQRT16
*
      END
