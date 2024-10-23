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

      SUBROUTINE SLAUUM( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLAUUM computes the product U * U' or L' * L, where the triangular
*  factor U or L is stored in the upper or lower triangular part of
*  the array A.
*
*  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*  overwriting the factor U in A.
*  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*  overwriting the factor L in A.
*
*  This is the blocked form of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the triangular factor stored in the array A
*          is upper or lower triangular:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the triangular factor U or L.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the triangular factor U or L.
*          On exit, if UPLO = 'U', the upper triangle of A is
*          overwritten with the upper triangle of the product U * U';
*          if UPLO = 'L', the lower triangle of A is overwritten with
*          the lower triangle of the product L' * L.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLAUU2, SSYRK, STRMM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAUUM', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SLAUUM', UPLO, N, -1, -1, -1 )
*
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL SLAUU2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute the product U * U'.
*
            DO 10 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Non-unit',
     $                     I-1, IB, ONE, A( I, I ), LDA, A( 1, I ),
     $                     LDA )
               CALL SLAUU2( 'Upper', IB, A( I, I ), LDA, INFO )
               IF( I+IB.LE.N ) THEN
                  CALL SGEMM( 'No transpose', 'Transpose', I-1, IB,
     $                        N-I-IB+1, ONE, A( 1, I+IB ), LDA,
     $                        A( I, I+IB ), LDA, ONE, A( 1, I ), LDA )
                  CALL SSYRK( 'Upper', 'No transpose', IB, N-I-IB+1,
     $                        ONE, A( I, I+IB ), LDA, ONE, A( I, I ),
     $                        LDA )
               END IF
   10       CONTINUE
         ELSE
*
*           Compute the product L' * L.
*
            DO 20 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               CALL STRMM( 'Left', 'Lower', 'Transpose', 'Non-unit', IB,
     $                     I-1, ONE, A( I, I ), LDA, A( I, 1 ), LDA )
               CALL SLAUU2( 'Lower', IB, A( I, I ), LDA, INFO )
               IF( I+IB.LE.N ) THEN
                  CALL SGEMM( 'Transpose', 'No transpose', IB, I-1,
     $                        N-I-IB+1, ONE, A( I+IB, I ), LDA,
     $                        A( I+IB, 1 ), LDA, ONE, A( I, 1 ), LDA )
                  CALL SSYRK( 'Lower', 'Transpose', IB, N-I-IB+1, ONE,
     $                        A( I+IB, I ), LDA, ONE, A( I, I ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SLAUUM
*
      END
