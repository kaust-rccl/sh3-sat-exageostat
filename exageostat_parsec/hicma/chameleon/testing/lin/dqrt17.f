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

      DOUBLE PRECISION FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A,
     $                 LDA, X, LDX, B, LDB, C, WORK, LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDB, * ),
     $                   WORK( LWORK ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DQRT17 computes the ratio
*
*     || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
*
*  where R = op(A)*X - B, op(A) is A or A', and
*
*     alpha = ||B|| if IRESID = 1 (zero-residual problem)
*     alpha = ||R|| if IRESID = 2 (otherwise).
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies whether or not the transpose of A is used.
*          = 'N':  No transpose, op(A) = A.
*          = 'T':  Transpose, op(A) = A'.
*
*  IRESID  (input) INTEGER
*          IRESID = 1 indicates zero-residual problem.
*          IRESID = 2 indicates non-zero residual.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*          If TRANS = 'N', the number of rows of the matrix B.
*          If TRANS = 'T', the number of rows of the matrix X.
*
*  N       (input) INTEGER
*          The number of columns of the matrix  A.
*          If TRANS = 'N', the number of rows of the matrix X.
*          If TRANS = 'T', the number of rows of the matrix B.
*
*  NRHS    (input) INTEGER
*          The number of columns of the matrices X and B.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= M.
*
*  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          If TRANS = 'N', the n-by-nrhs matrix X.
*          If TRANS = 'T', the m-by-nrhs matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.
*          If TRANS = 'N', LDX >= N.
*          If TRANS = 'T', LDX >= M.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          If TRANS = 'N', the m-by-nrhs matrix B.
*          If TRANS = 'T', the n-by-nrhs matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.
*          If TRANS = 'N', LDB >= M.
*          If TRANS = 'T', LDB >= N.
*
*  C       (workspace) DOUBLE PRECISION array, dimension (LDB,NRHS)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= NRHS*(M+N).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, ISCL, NCOLS, NROWS
      DOUBLE PRECISION   BIGNUM, ERR, NORMA, NORMB, NORMRS, NORMX,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
      DQRT17 = ZERO
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         NROWS = N
         NCOLS = M
      ELSE
         CALL XERBLA( 'DQRT17', 1 )
         RETURN
      END IF
*
      IF( LWORK.LT.NCOLS*NRHS ) THEN
         CALL XERBLA( 'DQRT17', 13 )
         RETURN
      END IF
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) THEN
         RETURN
      END IF
*
      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      ISCL = 0
*
*     compute residual and scale it
*
      CALL DLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A,
     $            LDA, X, LDX, ONE, C, LDB )
      NORMRS = DLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB,
     $                INFO )
      END IF
*
*     compute R'*A
*
      CALL DGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB,
     $            A, LDA, ZERO, WORK, NRHS )
*
*     compute and properly scale error
*
      ERR = DLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO )
     $   ERR = ERR / NORMA
*
      IF( ISCL.EQ.1 )
     $   ERR = ERR*NORMRS
*
      IF( IRESID.EQ.1 ) THEN
         NORMB = DLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO )
     $      ERR = ERR / NORMB
      ELSE
         NORMX = DLANGE( 'One-norm', NCOLS, NRHS, X, LDX, RWORK )
         IF( NORMX.NE.ZERO )
     $      ERR = ERR / NORMX
      END IF
*
      DQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
      RETURN
*
*     End of DQRT17
*
      END
