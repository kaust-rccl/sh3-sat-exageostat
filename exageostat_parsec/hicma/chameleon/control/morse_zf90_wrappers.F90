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

!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!   
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!     
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
!
!
!  MORSE fortran wrapper for BLAS and LAPACK subroutines.
!  MORSE is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 1.0.0
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @author Emmanuel Agullo
! @author Cedric Castagnede
! @date 2011-09-15
! @precisions normal z -> c d s
!
!
! Wrappers to MORSE functions are provided for the following BLAS
! subroutines since the MORSE and BLAS interfaces match exactly:
! ZGEMM  MORSE_zgemm
! ZHEMM  MORSE_zhemm
! ZHER2K MORSE_zher2k
! ZHERK  MORSE_zherk
! ZSYMM  MORSE_zsymm
! ZSYR2K MORSE_zsyr2k
! ZSYRK  MORSE_zsyrk
! ZTRMM  MORSE_ztrmm
! ZTRSM  MORSE_ztrsm
!
! Wrappers to MORSE functions are provided for the following LAPACK
! subroutines since the MORSE and LAPACK interfaces match exactly:
! ZGESV  MORSE_zgesv
! ZGETRF MORSE_zgetrf
! ZGETRS MORSE_zgetrs
! ZHEGST MORSE_zhegst
! ZLASWP MORSE_zlaswp
! ZLAUUM MORSE_zlauum
! ZPOSV  MORSE_zposv
! ZSYSV  MORSE_zsysv
! ZPOTRF MORSE_zpotrf
! ZSYTRF MORSE_zsytrf
! ZPOTRI MORSE_zpotri
! ZPOTRS MORSE_zpotrs
! ZSYTRS MORSE_zsytrs
! ZTRTRI MORSE_ztrtri
! ZLACPY MORSE_zlacpy
! ZLASET MORSE_zlaset
#define PRECISION_z

      subroutine morse_wrap_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZGESV"
            call MORSE_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine morse_wrap_ZGESV

      subroutine morse_wrap_ZGETRF(M,N,A,LDA,IPIV,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZGETRF"
            call MORSE_ZGETRF(M,N,A,LDA,IPIV,INFO)
      end subroutine morse_wrap_ZGETRF

      subroutine morse_wrap_ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(in), target :: IPIV(*)
            integer, intent(out) :: INFO
            character, intent(in) :: TRANS
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_TRANS
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = MorseNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = MorseTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = MorseConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZGETRS"
            call MORSE_ZGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine morse_wrap_ZGETRS

      subroutine morse_wrap_ZHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: ITYPE
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZHEGST"
            call MORSE_ZHEGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
      end subroutine morse_wrap_ZHEGST

      subroutine morse_wrap_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: INCX
            integer, intent(in) :: K1
            integer, intent(in) :: K2
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(in), target :: IPIV(*)
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_ret
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZLASWP"
            call MORSE_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
      end subroutine morse_wrap_ZLASWP

      subroutine morse_wrap_ZLAUUM(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZLAUUM"
            call MORSE_ZLAUUM(local_UPLO,N,A,LDA,INFO)
      end subroutine morse_wrap_ZLAUUM

      subroutine morse_wrap_ZPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZPOSV"
            call MORSE_ZPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine morse_wrap_ZPOSV

      subroutine morse_wrap_ZSYSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZSYSV"
            call MORSE_ZSYSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine morse_wrap_ZSYSV

      subroutine morse_wrap_ZPOTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZPOTRF"
            call MORSE_ZPOTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine morse_wrap_ZPOTRF

      subroutine morse_wrap_ZSYTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZSYTRF"
            call MORSE_ZSYTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine morse_wrap_ZSYTRF

      subroutine morse_wrap_ZPOTRI(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZPOTRI"
            call MORSE_ZPOTRI(local_UPLO,N,A,LDA,INFO)
      end subroutine morse_wrap_ZPOTRI

      subroutine morse_wrap_ZPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZPOTRS"
            call MORSE_ZPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine morse_wrap_ZPOTRS

      subroutine morse_wrap_ZSYTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZSYTRS"
            call MORSE_ZSYTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine morse_wrap_ZSYTRS

      subroutine morse_wrap_ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: DIAG
            character, intent(in) :: UPLO
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_DIAG
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = MorseUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = MorseNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. morse_initialized) call morse_init(24,INFO)
            ! write(*,*) " Calling MORSE_ZTRTRI"
            call MORSE_ZTRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
      end subroutine morse_wrap_ZTRTRI

      subroutine morse_wrap_ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: TRANSA
            character, intent(in) :: TRANSB
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANSA
            integer :: local_TRANSB
            integer :: local_ret
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = MorseNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = MorseTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = MorseConjTrans
            else
               local_TRANSA = -1
            end if
            if(TRANSB=='N' .or. TRANSB=='n')then
               local_TRANSB = MorseNoTrans
            else if(TRANSB=='T' .or. TRANSB=='t')then
               local_TRANSB = MorseTrans
            else if(TRANSB=='C' .or. TRANSB=='c')then
               local_TRANSB = MorseConjTrans
            else
               local_TRANSB = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZGEMM"
            call MORSE_ZGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine morse_wrap_ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = MorseLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = MorseRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZHEMM"
            call MORSE_ZHEMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZHEMM

      subroutine morse_wrap_ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = MorseNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = MorseTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = MorseConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZHER2K"
            call MORSE_ZHER2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZHER2K

      subroutine morse_wrap_ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            double precision, intent(in) :: ALPHA
            double precision, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = MorseNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = MorseTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = MorseConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZHERK"
            call MORSE_ZHERK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZHERK
#endif

      subroutine morse_wrap_ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = MorseLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = MorseRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZSYMM"
            call MORSE_ZSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZSYMM

      subroutine morse_wrap_ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(in), target :: B(LDB,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = MorseNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = MorseTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = MorseConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZSYR2K"
            call MORSE_ZSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZSYR2K

      subroutine morse_wrap_ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = MorseNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = MorseTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = MorseConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZSYRK"
            call MORSE_ZSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine morse_wrap_ZSYRK

      subroutine morse_wrap_ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = MorseLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = MorseRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = MorseNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = MorseTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = MorseConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = MorseUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = MorseNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZTRMM"
            call MORSE_ZTRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine morse_wrap_ZTRMM

      subroutine morse_wrap_ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = MorseLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = MorseRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = MorseNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = MorseTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = MorseConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = MorseUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = MorseNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZTRSM"
            call MORSE_ZTRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine morse_wrap_ZTRSM

      subroutine morse_wrap_ZLACPY(UPLO,M,N,A,LDA,B,LDB)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in), target :: A(LDA,*)
            complex(kind=wp), intent(out), target :: B(LDB,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZLACPY"
            call MORSE_ZLACPY(local_UPLO,M,N,A,LDA,B,LDB,local_ret)
      end subroutine morse_wrap_ZLACPY

      subroutine morse_wrap_ZLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
            use iso_c_binding
            use morse
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex(kind=wp), intent(in) :: ALPHA
            complex(kind=wp), intent(in) :: BETA
            complex(kind=wp), intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = MorseUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = MorseLower
            else
               local_UPLO = -1
            end if
            if (.not. morse_initialized) call morse_init(24,local_ret)
            ! write(*,*) " Calling MORSE_ZLASET"
            call MORSE_ZLASET(local_UPLO,M,N,ALPHA,BETA,A,LDA,local_ret)
      end subroutine morse_wrap_ZLASET
