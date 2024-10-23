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

      LOGICAL          FUNCTION LSAMEN( N, CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
*     ..
*
*  Purpose
*  =======
*
*  LSAMEN  tests if the first N letters of CA are the same as the
*  first N letters of CB, regardless of case.
*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case
*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
*  or LEN( CB ) is less than N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of characters in CA and CB to be compared.
*
*  CA      (input) CHARACTER*(*)
*  CB      (input) CHARACTER*(*)
*          CA and CB specify two character strings of length at least N.
*          Only the first N characters of each string will be accessed.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN
*     ..
*     .. Executable Statements ..
*
      LSAMEN = .FALSE.
      IF( LEN( CA ).LT.N .OR. LEN( CB ).LT.N )
     $   GO TO 20
*
*     Do for each character in the two strings.
*
      DO 10 I = 1, N
*
*        Test if the characters are equal using LSAME.
*
         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) )
     $      GO TO 20
*
   10 CONTINUE
      LSAMEN = .TRUE.
*
   20 CONTINUE
      RETURN
*
*     End of LSAMEN
*
      END
