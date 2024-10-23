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
!  MORSE Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  MORSE is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 1.0.0
! @author Numerical Algorithms Group
! @author Mathieu Faverge
! @author Emmanuel Agullo
! @author Cedric Castagnede
! @date 2011-12-15
! @precisions normal z -> c d s
!
#define PRECISION_z

module morse_z
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (simple interface)
    !
      interface
         function MORSE_zLapack_to_Tile_c(Af77,LDA,A) &
          & bind(c, name='MORSE_zLapack_to_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zLapack_to_Tile_c
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: A
          end function MORSE_zLapack_to_Tile_c
      end interface

      interface
         function MORSE_zTile_to_Lapack_c(A,Af77,LDA) &
          & bind(c, name='MORSE_zTile_to_Lapack')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zTile_to_Lapack_c
            type(c_ptr), value :: A
            type(c_ptr), value :: Af77
            integer(kind=c_int), value :: LDA
          end function MORSE_zTile_to_Lapack_c
      end interface

      interface
         function MORSE_zgebrd_c(M,N,A,LDA,D,E,descT) &
          & bind(c, name='MORSE_zgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: descT
          end function MORSE_zgebrd_c
      end interface

      interface
         function MORSE_zgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='MORSE_zgecfi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgecfi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: fout
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
          end function MORSE_zgecfi_c
      end interface

      interface
         function MORSE_zgecfi_Async_c(m,n,A,f_in,imb,inb,f_out,omb,onb,sequence,request) &
          & bind(c, name='MORSE_zgecfi_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgecfi_Async_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: f_in
            integer(kind=c_int), value :: imb
            integer(kind=c_int), value :: inb
            integer(kind=c_int), value :: f_out
            integer(kind=c_int), value :: omb
            integer(kind=c_int), value :: onb
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgecfi_Async_c
      end interface

      interface
         function MORSE_zgelqf_c(M,N,A,LDA,descT) &
          & bind(c, name='MORSE_zgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
          end function MORSE_zgelqf_c
      end interface

      interface
         function MORSE_zgelqs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zgelqs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgelqs_c
      end interface

      interface
         function MORSE_zgels_c(trans,M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgels_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgels_c
      end interface

      interface
         function MORSE_zgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='MORSE_zgemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgemm_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zgemm_c
      end interface

      interface
         function MORSE_zgeqrf_c(M,N,A,LDA,descT) &
          & bind(c, name='MORSE_zgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
          end function MORSE_zgeqrf_c
      end interface

      interface
         function MORSE_zgeqrs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zgeqrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrs_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgeqrs_c
      end interface

      interface
         function MORSE_zgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='MORSE_zgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgesv_c
      end interface

      interface
         function MORSE_zgesv_incpiv_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='MORSE_zgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_incpiv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgesv_incpiv_c
      end interface

      interface
         function MORSE_zgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,descT) &
          & bind(c, name='MORSE_zgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesvd_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            integer(kind=c_int), value :: LDU
            type(c_ptr), value :: VT
            integer(kind=c_int), value :: LDVT
            type(c_ptr), value :: descT
          end function MORSE_zgesvd_c
      end interface

      interface
         function MORSE_zgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='MORSE_zgetmi')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetmi_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: fin
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: nb
          end function MORSE_zgetmi_c
      end interface

      interface
         function MORSE_zgetmi_Async_c(m,n,A,f_in,mb,inb,sequence,request) &
          & bind(c, name='MORSE_zgetmi_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetmi_Async_c
            integer(kind=c_int), value :: m
            integer(kind=c_int), value :: n
            type(c_ptr), value :: A
            integer(kind=c_int), value :: f_in
            integer(kind=c_int), value :: mb
            integer(kind=c_int), value :: inb
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetmi_Async_c
      end interface

      interface
         function MORSE_zgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='MORSE_zgetrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
          end function MORSE_zgetrf_c
      end interface

      interface
         function MORSE_zgetrf_incpiv_c(M,N,A,LDA,descL,IPIV) &
          & bind(c, name='MORSE_zgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
          end function MORSE_zgetrf_incpiv_c
      end interface

      interface
         function MORSE_zgetrf_nopiv_c(M,N,A,LDA) &
          & bind(c, name='MORSE_zgetrf_nopiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_nopiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zgetrf_nopiv_c
      end interface

      interface
         function MORSE_zgetri_c(N,A,LDA,IPIV) &
          & bind(c, name='MORSE_zgetri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetri_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
          end function MORSE_zgetri_c
      end interface

      interface
         function MORSE_zgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='MORSE_zgetrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgetrs_c
      end interface

      interface
         function MORSE_zgetrs_incpiv_c(trans,N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='MORSE_zgetrs_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_incpiv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zgetrs_incpiv_c
      end interface

      interface
         function MORSE_zheev_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='MORSE_zheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheev_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zheev_c
      end interface

      interface
         function MORSE_zheevd_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='MORSE_zheevd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheevd_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zheevd_c
      end interface

      interface
         function MORSE_zhegst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zhegst')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegst_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zhegst_c
      end interface

      interface
         function MORSE_zhegv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='MORSE_zhegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegv_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zhegv_c
      end interface

      interface
         function MORSE_zhegvd_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='MORSE_zhegvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegvd_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: W
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zhegvd_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zhemm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='MORSE_zhemm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhemm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zhemm_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zher2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='MORSE_zher2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zher2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zher2k_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zherk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='MORSE_zherk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zherk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zherk_c
      end interface
#endif

      interface
         function MORSE_zhetrd_c(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ) &
          & bind(c, name='MORSE_zhetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhetrd_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: descT
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zhetrd_c
      end interface

      interface
         function MORSE_zlacpy_c(uplo,M,N,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zlacpy')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlacpy_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zlacpy_c
      end interface

      interface
         function MORSE_zlange_c(norm,M,N,A,LDA,work) &
          & bind(c, name='MORSE_zlange')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlange_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
          end function MORSE_zlange_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zlanhe_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='MORSE_zlanhe')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlanhe_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
          end function MORSE_zlanhe_c
      end interface
#endif

      interface
         function MORSE_zlansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='MORSE_zlansy')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlansy_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: work
          end function MORSE_zlansy_c
      end interface

      interface
         function MORSE_zlaset_c(uplo,M,N,alpha,beta,A,LDA) &
          & bind(c, name='MORSE_zlaset')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaset_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zlaset_c
      end interface

      interface
         function MORSE_zlaswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='MORSE_zlaswp')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswp_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function MORSE_zlaswp_c
      end interface

      interface
         function MORSE_zlaswpc_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='MORSE_zlaswpc')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswpc_c
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function MORSE_zlaswpc_c
      end interface

      interface
         function MORSE_zlauum_c(uplo,N,A,LDA) &
          & bind(c, name='MORSE_zlauum')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlauum_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zlauum_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zplghe_c(bump,uplo,N,A,LDA,seed) &
          & bind(c, name='MORSE_zplghe')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplghe_c
            real(kind=c_double), value :: bump
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplghe_c
      end interface
#endif

      interface
         function MORSE_zplgsy_c(bump,uplo,N,A,LDA,seed) &
          & bind(c, name='MORSE_zplgsy')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplgsy_c
            complex(kind=c_double_complex), value :: bump
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplgsy_c
      end interface

      interface
         function MORSE_zplrnt_c(M,N,A,LDA,seed) &
          & bind(c, name='MORSE_zplrnt')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplrnt_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplrnt_c
      end interface

      interface
         function MORSE_zposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zposv_c
      end interface

      interface
         function MORSE_zsysv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zsysv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsysv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zsysv_c
      end interface

      interface
         function MORSE_zpotrf_c(uplo,N,A,LDA) &
          & bind(c, name='MORSE_zpotrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zpotrf_c
      end interface

      interface
         function MORSE_zsytrf_c(uplo,N,A,LDA) &
          & bind(c, name='MORSE_zsytrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrf_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zsytrf_c
      end interface

      interface
         function MORSE_zpotri_c(uplo,N,A,LDA) &
          & bind(c, name='MORSE_zpotri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_zpotri_c
      end interface

      interface
         function MORSE_zpotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zpotrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zpotrs_c
      end interface

      interface
         function MORSE_zsytrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='MORSE_zsytrs')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrs_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zsytrs_c
      end interface

      interface
         function MORSE_zsymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='MORSE_zsymm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsymm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zsymm_c
      end interface

      interface
         function MORSE_zsyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='MORSE_zsyr2k')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyr2k_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zsyr2k_c
      end interface

      interface
         function MORSE_zsyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='MORSE_zsyrk')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyrk_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            integer(kind=c_int), value :: LDC
          end function MORSE_zsyrk_c
      end interface

      interface
         function MORSE_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='MORSE_ztrmm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrmm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_ztrmm_c
      end interface

      interface
         function MORSE_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='MORSE_ztrsm')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsm_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_ztrsm_c
      end interface

      interface
         function MORSE_ztrsmpl_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='MORSE_ztrsmpl')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmpl_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descL
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_ztrsmpl_c
      end interface

      interface
         function MORSE_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='MORSE_ztrsmrv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmrv_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_ztrsmrv_c
      end interface

      interface
         function MORSE_ztrtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='MORSE_ztrtri')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrtri_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            integer(kind=c_int), value :: N
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
          end function MORSE_ztrtri_c
      end interface

      interface
         function MORSE_zunglq_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zunglq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunglq_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zunglq_c
      end interface

      interface
         function MORSE_zungqr_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zungqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zungqr_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zungqr_c
      end interface

      interface
         function MORSE_zunmlq_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zunmlq')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmlq_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zunmlq_c
      end interface

      interface
         function MORSE_zunmqr_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='MORSE_zunmqr')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmqr_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: K
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: descT
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
          end function MORSE_zunmqr_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (native interface)
    !
      interface
         function MORSE_zgebrd_Tile_c(A,D,E,T) &
          & bind(c, name='MORSE_zgebrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgebrd_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
          end function MORSE_zgebrd_Tile_c
      end interface

      interface
         function MORSE_zgelqf_Tile_c(A,T) &
          & bind(c, name='MORSE_zgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
          end function MORSE_zgelqf_Tile_c
      end interface

      interface
         function MORSE_zgelqs_Tile_c(A,T,B) &
          & bind(c, name='MORSE_zgelqs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zgelqs_Tile_c
      end interface

      interface
         function MORSE_zgels_Tile_c(trans,A,T,B) &
          & bind(c, name='MORSE_zgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgels_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zgels_Tile_c
      end interface

      interface
         function MORSE_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='MORSE_zgemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgemm_Tile_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zgemm_Tile_c
      end interface

      interface
         function MORSE_zgeqrf_Tile_c(A,T) &
          & bind(c, name='MORSE_zgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
          end function MORSE_zgeqrf_Tile_c
      end interface

      interface
         function MORSE_zgeqrs_Tile_c(A,T,B) &
          & bind(c, name='MORSE_zgeqrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrs_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zgeqrs_Tile_c
      end interface

      interface
         function MORSE_zgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='MORSE_zgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function MORSE_zgesv_Tile_c
      end interface

      interface
         function MORSE_zgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='MORSE_zgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function MORSE_zgesv_incpiv_Tile_c
      end interface

      interface
         function MORSE_zgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='MORSE_zgesvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesvd_Tile_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
          end function MORSE_zgesvd_Tile_c
      end interface

      interface
         function MORSE_zgetrf_Tile_c(A,IPIV) &
          & bind(c, name='MORSE_zgetrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
          end function MORSE_zgetrf_Tile_c
      end interface

      interface
         function MORSE_zgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='MORSE_zgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
          end function MORSE_zgetrf_incpiv_Tile_c
      end interface

      interface
         function MORSE_zgetrf_nopiv_Tile_c(A) &
          & bind(c, name='MORSE_zgetrf_nopiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_nopiv_Tile_c
            type(c_ptr), value :: A
          end function MORSE_zgetrf_nopiv_Tile_c
      end interface

      interface
         function MORSE_zgetri_Tile_c(A,IPIV) &
          & bind(c, name='MORSE_zgetri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetri_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
          end function MORSE_zgetri_Tile_c
      end interface

      interface
         function MORSE_zgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='MORSE_zgetrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function MORSE_zgetrs_Tile_c
      end interface

      interface
         function MORSE_zgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='MORSE_zgetrs_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_incpiv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function MORSE_zgetrs_incpiv_Tile_c
      end interface

      interface
         function MORSE_zheev_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='MORSE_zheev_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheev_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zheev_Tile_c
      end interface

      interface
         function MORSE_zheevd_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='MORSE_zheevd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheevd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zheevd_Tile_c
      end interface

      interface
         function MORSE_zhegst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='MORSE_zhegst_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegst_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zhegst_Tile_c
      end interface

      interface
         function MORSE_zhegv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='MORSE_zhegv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegv_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
          end function MORSE_zhegv_Tile_c
      end interface

      interface
         function MORSE_zhegvd_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='MORSE_zhegvd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegvd_Tile_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
          end function MORSE_zhegvd_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='MORSE_zhemm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhemm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zhemm_Tile_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='MORSE_zher2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zher2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zher2k_Tile_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zherk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='MORSE_zherk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zherk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zherk_Tile_c
      end interface
#endif

      interface
         function MORSE_zhetrd_Tile_c(jobz,uplo,A,D,E,T,Q,LDQ) &
          & bind(c, name='MORSE_zhetrd_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhetrd_Tile_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
          end function MORSE_zhetrd_Tile_c
      end interface

      interface
         function MORSE_zlacpy_Tile_c(uplo,A,B) &
          & bind(c, name='MORSE_zlacpy_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlacpy_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zlacpy_Tile_c
      end interface

      interface
         function MORSE_zlange_Tile_c(norm,A,work) &
          & bind(c, name='MORSE_zlange_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlange_Tile_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: work
          end function MORSE_zlange_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zlanhe_Tile_c(norm,uplo,A,work) &
          & bind(c, name='MORSE_zlanhe_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlanhe_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
          end function MORSE_zlanhe_Tile_c
      end interface
#endif

      interface
         function MORSE_zlansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='MORSE_zlansy_Tile')
            use iso_c_binding
            implicit none
            real(kind=c_double) :: MORSE_zlansy_Tile_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
          end function MORSE_zlansy_Tile_c
      end interface

      interface
         function MORSE_zlaset_Tile_c(uplo,alpha,beta,A) &
          & bind(c, name='MORSE_zlaset_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaset_Tile_c
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
          end function MORSE_zlaset_Tile_c
      end interface

      interface
         function MORSE_zlaswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='MORSE_zlaswp_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswp_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function MORSE_zlaswp_Tile_c
      end interface

      interface
         function MORSE_zlaswpc_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='MORSE_zlaswpc_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswpc_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
          end function MORSE_zlaswpc_Tile_c
      end interface

      interface
         function MORSE_zlauum_Tile_c(uplo,A) &
          & bind(c, name='MORSE_zlauum_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlauum_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function MORSE_zlauum_Tile_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zplghe_Tile_c(bump,uplo,A,seed) &
          & bind(c, name='MORSE_zplghe_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplghe_Tile_c
            real(kind=c_double), value :: bump
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplghe_Tile_c
      end interface
#endif

      interface
         function MORSE_zplgsy_Tile_c(bump,uplo,A,seed) &
          & bind(c, name='MORSE_zplgsy_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplgsy_Tile_c
            complex(kind=c_double_complex), value :: bump
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplgsy_Tile_c
      end interface

      interface
         function MORSE_zplrnt_Tile_c(A,seed) &
          & bind(c, name='MORSE_zplrnt_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplrnt_Tile_c
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
          end function MORSE_zplrnt_Tile_c
      end interface

      interface
         function MORSE_zposv_Tile_c(uplo,A,B) &
          & bind(c, name='MORSE_zposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zposv_Tile_c
      end interface

      interface
         function MORSE_zsysv_Tile_c(uplo,A,B) &
          & bind(c, name='MORSE_zsysv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsysv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zsysv_Tile_c
      end interface

      interface
         function MORSE_zpotrf_Tile_c(uplo,A) &
          & bind(c, name='MORSE_zpotrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function MORSE_zpotrf_Tile_c
      end interface

      interface
         function MORSE_zsytrf_Tile_c(uplo,A) &
          & bind(c, name='MORSE_zsytrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrf_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function MORSE_zsytrf_Tile_c
      end interface

      interface
         function MORSE_zpotri_Tile_c(uplo,A) &
          & bind(c, name='MORSE_zpotri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotri_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
          end function MORSE_zpotri_Tile_c
      end interface

      interface
         function MORSE_zpotrs_Tile_c(uplo,A,B) &
          & bind(c, name='MORSE_zpotrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zpotrs_Tile_c
      end interface

      interface
         function MORSE_zsytrs_Tile_c(uplo,A,B) &
          & bind(c, name='MORSE_zsytrs_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrs_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_zsytrs_Tile_c
      end interface

      interface
         function MORSE_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='MORSE_zsymm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsymm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zsymm_Tile_c
      end interface

      interface
         function MORSE_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='MORSE_zsyr2k_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyr2k_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zsyr2k_Tile_c
      end interface

      interface
         function MORSE_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='MORSE_zsyrk_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyrk_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
          end function MORSE_zsyrk_Tile_c
      end interface

      interface
         function MORSE_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='MORSE_ztrmm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrmm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_ztrmm_Tile_c
      end interface

      interface
         function MORSE_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='MORSE_ztrsm_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsm_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_ztrsm_Tile_c
      end interface

      interface
         function MORSE_ztrsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='MORSE_ztrsmpl_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmpl_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
          end function MORSE_ztrsmpl_Tile_c
      end interface

      interface
         function MORSE_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='MORSE_ztrsmrv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmrv_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
          end function MORSE_ztrsmrv_Tile_c
      end interface

      interface
         function MORSE_ztrtri_Tile_c(uplo,diag,A) &
          & bind(c, name='MORSE_ztrtri_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrtri_Tile_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
          end function MORSE_ztrtri_Tile_c
      end interface

      interface
         function MORSE_zunglq_Tile_c(A,T,B) &
          & bind(c, name='MORSE_zunglq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunglq_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zunglq_Tile_c
      end interface

      interface
         function MORSE_zungqr_Tile_c(A,T,B) &
          & bind(c, name='MORSE_zungqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zungqr_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zungqr_Tile_c
      end interface

      interface
         function MORSE_zunmlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='MORSE_zunmlq_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmlq_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zunmlq_Tile_c
      end interface

      interface
         function MORSE_zunmqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='MORSE_zunmqr_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmqr_Tile_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
          end function MORSE_zunmqr_Tile_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (asynchronous interface)
    !
      interface
         function MORSE_zgebrd_Tile_Async_c(A,D,E,T,sequence,request) &
          & bind(c, name='MORSE_zgebrd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgebrd_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgebrd_Tile_Async_c
      end interface

      interface
         function MORSE_zgelqf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='MORSE_zgelqf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgelqf_Tile_Async_c
      end interface

      interface
         function MORSE_zgelqs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='MORSE_zgelqs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgelqs_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgelqs_Tile_Async_c
      end interface

      interface
         function MORSE_zgels_Tile_Async_c(trans,A,T,B,sequence,request) &
          & bind(c, name='MORSE_zgels_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgels_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgels_Tile_Async_c
      end interface

      interface
         function MORSE_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='MORSE_zgemm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgemm_Tile_Async_c
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: transB
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgemm_Tile_Async_c
      end interface

      interface
         function MORSE_zgeqrf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='MORSE_zgeqrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgeqrf_Tile_Async_c
      end interface

      interface
         function MORSE_zgeqrs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='MORSE_zgeqrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgeqrs_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgeqrs_Tile_Async_c
      end interface

      interface
         function MORSE_zgesv_Tile_Async_c(A,IPIV,B,sequence,request) &
          & bind(c, name='MORSE_zgesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgesv_Tile_Async_c
      end interface

      interface
         function MORSE_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='MORSE_zgesv_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesv_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgesv_incpiv_Tile_Async_c
      end interface

      interface
         function MORSE_zgesvd_Tile_Async_c(jobu,jobvt,A,S,U,VT,T,sequence,request) &
          & bind(c, name='MORSE_zgesvd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgesvd_Tile_Async_c
            integer(kind=c_int), value :: jobu
            integer(kind=c_int), value :: jobvt
            type(c_ptr), value :: A
            type(c_ptr), value :: S
            type(c_ptr), value :: U
            type(c_ptr), value :: VT
            type(c_ptr), value :: T
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgesvd_Tile_Async_c
      end interface

      interface
         function MORSE_zgetrf_Tile_Async_c(A,IPIV,sequence,request) &
          & bind(c, name='MORSE_zgetrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetrf_Tile_Async_c
      end interface

      interface
         function MORSE_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request) &
          & bind(c, name='MORSE_zgetrf_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetrf_incpiv_Tile_Async_c
      end interface

      interface
         function MORSE_zgetrf_nopiv_Tile_Async_c(A,sequence,request) &
          & bind(c, name='MORSE_zgetrf_nopiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrf_nopiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetrf_nopiv_Tile_Async_c
      end interface

      interface
         function MORSE_zgetri_Tile_Async_c(A,IPIV,W,sequence,request) &
          & bind(c, name='MORSE_zgetri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetri_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: W
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetri_Tile_Async_c
      end interface

      interface
         function MORSE_zgetrs_Tile_Async_c(trans,A,IPIV,B,sequence,request) &
          & bind(c, name='MORSE_zgetrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetrs_Tile_Async_c
      end interface

      interface
         function MORSE_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='MORSE_zgetrs_incpiv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zgetrs_incpiv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zgetrs_incpiv_Tile_Async_c
      end interface

      interface
         function MORSE_zheev_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='MORSE_zheev_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheev_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zheev_Tile_Async_c
      end interface

      interface
         function MORSE_zheevd_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='MORSE_zheevd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zheevd_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zheevd_Tile_Async_c
      end interface

      interface
         function MORSE_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zhegst_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegst_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zhegst_Tile_Async_c
      end interface

      interface
         function MORSE_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='MORSE_zhegv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegv_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zhegv_Tile_Async_c
      end interface

      interface
         function MORSE_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='MORSE_zhegvd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhegvd_Tile_Async_c
            integer(kind=c_int), value :: itype
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: W
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zhegvd_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='MORSE_zhemm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhemm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zhemm_Tile_Async_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='MORSE_zher2k_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zher2k_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zher2k_Tile_Async_c
      end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='MORSE_zherk_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zherk_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            real(kind=c_double), value :: alpha
            type(c_ptr), value :: A
            real(kind=c_double), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zherk_Tile_Async_c
      end interface
#endif

      interface
         function MORSE_zhetrd_Tile_Async_c(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request) &
          & bind(c, name='MORSE_zhetrd_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zhetrd_Tile_Async_c
            integer(kind=c_int), value :: jobz
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: D
            type(c_ptr), value :: E
            type(c_ptr), value :: T
            type(c_ptr), value :: Q
            integer(kind=c_int), value :: LDQ
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zhetrd_Tile_Async_c
      end interface

      interface
         function MORSE_zlacpy_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zlacpy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlacpy_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlacpy_Tile_Async_c
      end interface

      interface
         function MORSE_zlange_Tile_Async_c(norm,A,work,value,sequence,request) &
          & bind(c, name='MORSE_zlange_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlange_Tile_Async_c
            integer(kind=c_int), value :: norm
            type(c_ptr), value :: A
            type(c_ptr), value :: work
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlange_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zlanhe_Tile_Async_c(norm,uplo,A,work,value,sequence,request) &
          & bind(c, name='MORSE_zlanhe_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlanhe_Tile_Async_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlanhe_Tile_Async_c
      end interface
#endif

      interface
         function MORSE_zlansy_Tile_Async_c(norm,uplo,A,work,value,sequence,request) &
          & bind(c, name='MORSE_zlansy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlansy_Tile_Async_c
            integer(kind=c_int), value :: norm
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: work
            type(c_ptr), value :: value
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlansy_Tile_Async_c
      end interface

      interface
         function MORSE_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request) &
          & bind(c, name='MORSE_zlaset_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaset_Tile_Async_c
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlaset_Tile_Async_c
      end interface

      interface
         function MORSE_zlaswp_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='MORSE_zlaswp_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswp_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlaswp_Tile_Async_c
      end interface

      interface
         function MORSE_zlaswpc_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='MORSE_zlaswpc_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlaswpc_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_int), value :: K1
            integer(kind=c_int), value :: K2
            type(c_ptr), value :: IPIV
            integer(kind=c_int), value :: INCX
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlaswpc_Tile_Async_c
      end interface

      interface
         function MORSE_zlauum_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='MORSE_zlauum_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zlauum_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zlauum_Tile_Async_c
      end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
      interface
         function MORSE_zplghe_Tile_Async_c(bump,uplo,A,seed,sequence,request) &
          & bind(c, name='MORSE_zplghe_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplghe_Tile_Async_c
            real(kind=c_double), value :: bump
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zplghe_Tile_Async_c
      end interface
#endif

      interface
         function MORSE_zplgsy_Tile_Async_c(bump,uplo,A,seed,sequence,request) &
          & bind(c, name='MORSE_zplgsy_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplgsy_Tile_Async_c
            complex(kind=c_double_complex), value :: bump
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zplgsy_Tile_Async_c
      end interface

      interface
         function MORSE_zplrnt_Tile_Async_c(A,seed,sequence,request) &
          & bind(c, name='MORSE_zplrnt_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zplrnt_Tile_Async_c
            type(c_ptr), value :: A
            integer(kind=c_long_long), value :: seed
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zplrnt_Tile_Async_c
      end interface

      interface
         function MORSE_zposv_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zposv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zposv_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zposv_Tile_Async_c
      end interface

      interface
         function MORSE_zsysv_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zsysv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsysv_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsysv_Tile_Async_c
      end interface

      interface
         function MORSE_zpotrf_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='MORSE_zpotrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrf_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zpotrf_Tile_Async_c
      end interface

      interface
         function MORSE_zsytrf_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='MORSE_zsytrf_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrf_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsytrf_Tile_Async_c
      end interface

      interface
         function MORSE_zpotri_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='MORSE_zpotri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotri_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zpotri_Tile_Async_c
      end interface

      interface
         function MORSE_zpotrs_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zpotrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zpotrs_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zpotrs_Tile_Async_c
      end interface

      interface
         function MORSE_zsytrs_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='MORSE_zsytrs_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsytrs_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsytrs_Tile_Async_c
      end interface

      interface
         function MORSE_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='MORSE_zsymm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsymm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsymm_Tile_Async_c
      end interface

      interface
         function MORSE_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='MORSE_zsyr2k_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyr2k_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsyr2k_Tile_Async_c
      end interface

      interface
         function MORSE_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='MORSE_zsyrk_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zsyrk_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: trans
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            complex(kind=c_double_complex), value :: beta
            type(c_ptr), value :: C
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zsyrk_Tile_Async_c
      end interface

      interface
         function MORSE_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='MORSE_ztrmm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrmm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_ztrmm_Tile_Async_c
      end interface

      interface
         function MORSE_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='MORSE_ztrsm_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsm_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_ztrsm_Tile_Async_c
      end interface

      interface
         function MORSE_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='MORSE_ztrsmpl_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmpl_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: L
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_ztrsmpl_Tile_Async_c
      end interface

      interface
         function MORSE_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='MORSE_ztrsmrv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrsmrv_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: transA
            integer(kind=c_int), value :: diag
            complex(kind=c_double_complex), value :: alpha
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_ztrsmrv_Tile_Async_c
      end interface

      interface
         function MORSE_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request) &
          & bind(c, name='MORSE_ztrtri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_ztrtri_Tile_Async_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: diag
            type(c_ptr), value :: A
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_ztrtri_Tile_Async_c
      end interface

      interface
         function MORSE_zunglq_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='MORSE_zunglq_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunglq_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zunglq_Tile_Async_c
      end interface

      interface
         function MORSE_zungqr_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='MORSE_zungqr_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zungqr_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zungqr_Tile_Async_c
      end interface

      interface
         function MORSE_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='MORSE_zunmlq_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmlq_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zunmlq_Tile_Async_c
      end interface

      interface
         function MORSE_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='MORSE_zunmqr_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_zunmqr_Tile_Async_c
            integer(kind=c_int), value :: side
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
          end function MORSE_zunmqr_Tile_Async_c
      end interface

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - workspace allocation
    !
      interface
         function MORSE_Alloc_Workspace_zgebrd_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgebrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgebrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgebrd_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgeev_c(N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgeev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgeev_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgeev_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgehrd_c(N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgehrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgehrd_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgehrd_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgelqf_c(M,N,T,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgelqf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgelqf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgelqf_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgelqf_Tile_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgelqf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgelqf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgelqf_Tile_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgels_c(M,N,T,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgels')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgels_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgels_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgels_Tile_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgels_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgels_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgels_Tile_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgeqrf_c(M,N,T,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgeqrf')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgeqrf_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: T ! T is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgeqrf_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgeqrf_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgeqrf_Tile_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgeqrf_Tile_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgesv_incpiv_c(N,descL,IPIV,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgesv_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgesv_incpiv_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgesv_incpiv_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgesv_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgesv_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgesv_incpiv_Tile_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgesvd_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgesvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgesvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgesvd_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgetrf_incpiv_c(M,N,descL,IPIV,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgetrf_incpiv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgetrf_incpiv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgetrf_incpiv_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zgetrf_incpiv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgetrf_incpiv_Tile_c
            integer(kind=c_int), value :: N
            type(c_ptr) :: descL ! descL is **, so pass by reference
            type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zgetrf_incpiv_Tile_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zgetri_Tile_Async_c(A,W) &
          & bind(c, name='MORSE_Alloc_Workspace_zgetri_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zgetri_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: W
          end function MORSE_Alloc_Workspace_zgetri_Tile_Async_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zheev_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zheev')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zheev_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zheev_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zheevd_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zheevd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zheevd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zheevd_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zhegv_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zhegv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zhegv_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zhegv_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zhegvd_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zhegvd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zhegvd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zhegvd_c
      end interface

      interface
         function MORSE_Alloc_Workspace_zhetrd_c(M,N,descT,p,q) &
          & bind(c, name='MORSE_Alloc_Workspace_zhetrd')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Alloc_Workspace_zhetrd_c
            integer(kind=c_int), value :: M
            integer(kind=c_int), value :: N
            type(c_ptr) :: descT ! descT is **, so pass by reference
            integer(kind=c_int), value :: p
            integer(kind=c_int), value :: q
          end function MORSE_Alloc_Workspace_zhetrd_c
      end interface

  contains

       subroutine MORSE_zgebrd(M,N,A,LDA,D,E,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgebrd_c(M,N,c_loc(A),LDA,c_loc(D),c_loc(E),T)
      end subroutine MORSE_zgebrd

      subroutine MORSE_zgelqf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqf_c(M,N,c_loc(A),LDA,T)
      end subroutine MORSE_zgelqf

      subroutine MORSE_zgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zgelqs

      subroutine MORSE_zgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zgels

      subroutine MORSE_zgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = MORSE_zgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine MORSE_zgemm

      subroutine MORSE_zgeqrf(M,N,A,LDA,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrf_c(M,N,c_loc(A),LDA,T)
      end subroutine MORSE_zgeqrf

      subroutine MORSE_zgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zgeqrs

      subroutine MORSE_zgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine MORSE_zgesv

      subroutine MORSE_zgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine MORSE_zgesv_incpiv

      subroutine MORSE_zgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDU
         integer(kind=c_int), intent(in) :: LDVT
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: U(LDU,*)
         complex(kind=c_double_complex), intent(out), target :: VT(LDVT,*)
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
      end subroutine MORSE_zgesvd

      subroutine MORSE_zgetrf(M,N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(out), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine MORSE_zgetrf

      subroutine MORSE_zgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
      end subroutine MORSE_zgetrf_incpiv

      subroutine MORSE_zgetrf_nopiv(M,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zgetrf_nopiv_c(M,N,c_loc(A),LDA)
      end subroutine MORSE_zgetrf_nopiv

      subroutine MORSE_zgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
      end subroutine MORSE_zgetrs

      subroutine MORSE_zgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L    ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine MORSE_zgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zhemm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = MORSE_zhemm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine MORSE_zhemm

      subroutine MORSE_zherk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         info = MORSE_zherk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine MORSE_zherk

      subroutine MORSE_zher2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         real(kind=c_double), intent(in) :: beta
         info = MORSE_zher2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine MORSE_zher2k
#endif

      subroutine MORSE_zheev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zheev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zheev

      subroutine MORSE_zheevd(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zheevd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zheevd

      subroutine MORSE_zhegv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zhegv

      subroutine MORSE_zhegvd(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDQ
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: W(*)
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegvd_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zhegvd

      subroutine MORSE_zhegst(itype,uplo,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zhegst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_zhegst

      subroutine MORSE_zhetrd(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDQ
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: descT ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(inout), target :: Q(LDQ,*)
         info = MORSE_zhetrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),descT,c_loc(Q),LDQ)
      end subroutine MORSE_zhetrd

      function MORSE_zlange(norm,M,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlange
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         MORSE_zlange = MORSE_zlange_c(norm,M,N,c_loc(A),LDA,c_loc(work))
      end function MORSE_zlange

#if defined(PRECISION_z) || defined(PRECISION_c)
      function MORSE_zlanhe(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlanhe
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         MORSE_zlanhe = MORSE_zlanhe_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function MORSE_zlanhe
#endif

      function MORSE_zlansy(norm,uplo,N,A,LDA,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlansy
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         MORSE_zlansy = MORSE_zlansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
      end function MORSE_zlansy

      subroutine MORSE_zlaswp(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zlaswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine MORSE_zlaswp

      subroutine MORSE_zlauum(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zlauum_c(uplo,N,c_loc(A),LDA)
      end subroutine MORSE_zlauum

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zplghe(bump,uplo,N,A,LDA,seed,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = MORSE_zplghe_c(bump,N,c_loc(A),LDA,seed)
       end subroutine MORSE_zplghe
#endif

      subroutine MORSE_zplgsy(bump,uplo,N,A,LDA,seed,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = MORSE_zplgsy_c(bump,N,c_loc(A),LDA,seed)
      end subroutine MORSE_zplgsy
       
      subroutine MORSE_zplrnt(M,N,A,LDA,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_long_long), intent(in) :: seed
         complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
         info = MORSE_zplrnt_c(M,N,c_loc(A),LDA,seed)
      end subroutine MORSE_zplrnt
       
      subroutine MORSE_zposv(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_zposv

      subroutine MORSE_zpotrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zpotrf_c(uplo,N,c_loc(A),LDA)
      end subroutine MORSE_zpotrf

      subroutine MORSE_zsytrf(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zsytrf_c(uplo,N,c_loc(A),LDA)
      end subroutine MORSE_zsytrf

      subroutine MORSE_zpotri(uplo,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zpotri_c(uplo,N,c_loc(A),LDA)
      end subroutine MORSE_zpotri

      subroutine MORSE_zpotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zpotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_zpotrs

      subroutine MORSE_zsytrs(uplo,N,NRHS,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zsytrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_zsytrs

      subroutine MORSE_zsymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = MORSE_zsymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine MORSE_zsymm

      subroutine MORSE_zsyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = MORSE_zsyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
      end subroutine MORSE_zsyrk

      subroutine MORSE_zsyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDC
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
         info = MORSE_zsyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
      end subroutine MORSE_zsyr2k

      subroutine MORSE_ztrmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_ztrmm

      subroutine MORSE_ztrsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_ztrsm

      subroutine MORSE_ztrsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         type(c_ptr), value :: L    ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_ztrsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
      end subroutine MORSE_ztrsmpl

      subroutine MORSE_ztrtri(uplo,diag,N,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_ztrtri_c(uplo,diag,N,c_loc(A),LDA)
      end subroutine MORSE_ztrtri

      subroutine MORSE_zunglq(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = MORSE_zunglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zunglq

      subroutine MORSE_zungqr(M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = MORSE_zungqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zungqr

      subroutine MORSE_zunmlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zunmlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zunmlq

      subroutine MORSE_zunmqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in) :: K
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_zunmqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
      end subroutine MORSE_zunmqr

      subroutine MORSE_zgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = MORSE_zgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
      end subroutine MORSE_zgecfi

      subroutine MORSE_zgetmi(m,n,A,fin,mb,nb,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         info = MORSE_zgetmi_c(m,n,c_loc(A),fin,mb,nb)
      end subroutine MORSE_zgetmi

      subroutine MORSE_zgetri(N,A,LDA,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in), target :: IPIV(*)
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zgetri_c(N,c_loc(A),LDA,c_loc(IPIV))
      end subroutine MORSE_zgetri

      subroutine MORSE_zlacpy(uplo,M,N,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
         info = MORSE_zlacpy_c(uplo,M,N,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_zlacpy

      subroutine MORSE_zlaset(uplo,M,N,alpha,beta,A,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zlaset_c(uplo,M,N,alpha,beta,c_loc(A),LDA)
      end subroutine MORSE_zlaset

      subroutine MORSE_zlaswpc(N,A,LDA,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: N
         complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
         info = MORSE_zlaswpc_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
      end subroutine MORSE_zlaswpc

      subroutine MORSE_ztrsmrv(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
         info = MORSE_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
      end subroutine MORSE_ztrsmrv

      subroutine MORSE_zgebrd_Tile(A,D,E,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgebrd_Tile_c(A,c_loc(D),c_loc(E),T)
      end subroutine MORSE_zgebrd_Tile

      subroutine MORSE_zgelqf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqf_Tile_c(A,T)
      end subroutine MORSE_zgelqf_Tile

      subroutine MORSE_zgelqs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqs_Tile_c(A,T,B)
      end subroutine MORSE_zgelqs_Tile

      subroutine MORSE_zgels_Tile(trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgels_Tile_c(trans,A,T,B)
      end subroutine MORSE_zgels_Tile

      subroutine MORSE_zgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
      end subroutine MORSE_zgemm_Tile

      subroutine MORSE_zgeqrf_Tile(A,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrf_Tile_c(A,T)
      end subroutine MORSE_zgeqrf_Tile

      subroutine MORSE_zgeqrs_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrs_Tile_c(A,T,B)
      end subroutine MORSE_zgeqrs_Tile

      subroutine MORSE_zgesv_Tile(A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesv_Tile_c(A,c_loc(IPIV),B)
      end subroutine MORSE_zgesv_Tile

      subroutine MORSE_zgesv_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesv_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine MORSE_zgesv_incpiv_Tile

      subroutine MORSE_zgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
      end subroutine MORSE_zgesvd_Tile

      subroutine MORSE_zgetrf_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_Tile_c(A,c_loc(IPIV))
      end subroutine MORSE_zgetrf_Tile

      subroutine MORSE_zgetrf_incpiv_Tile(A,L,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_incpiv_Tile_c(A,L,IPIV)
      end subroutine MORSE_zgetrf_incpiv_Tile

      subroutine MORSE_zgetrf_nopiv_Tile(A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_nopiv_Tile_c(A)
      end subroutine MORSE_zgetrf_nopiv_Tile

      subroutine MORSE_zgetrs_Tile(trans,A,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrs_Tile_c(trans,A,c_loc(IPIV),B)
      end subroutine MORSE_zgetrs_Tile

      subroutine MORSE_zgetrs_incpiv_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrs_incpiv_Tile_c(A,L,IPIV,B)
      end subroutine MORSE_zgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zhemm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine MORSE_zhemm_Tile

      subroutine MORSE_zherk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zherk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine MORSE_zherk_Tile

      subroutine MORSE_zher2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine MORSE_zher2k_Tile
#endif

      subroutine MORSE_zheev_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = MORSE_zheev_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zheev_Tile

      subroutine MORSE_zheevd_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = MORSE_zheevd_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
      end subroutine MORSE_zheevd_Tile

      subroutine MORSE_zhegv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine MORSE_zhegv_Tile

      subroutine MORSE_zhegvd_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegvd_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
      end subroutine MORSE_zhegvd_Tile

      subroutine MORSE_zhegst_Tile(itype,uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegst_Tile_c(itype,uplo,A,B)
      end subroutine MORSE_zhegst_Tile

      subroutine MORSE_zhetrd_Tile(jobz,uplo,A,D,E,T,Q,LDQ,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
         info = MORSE_zhetrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
      end subroutine MORSE_zhetrd_Tile

      function MORSE_zlange_Tile(norm,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlange_Tile
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         real(kind=c_double), intent(inout), target :: work(*)
         MORSE_zlange_Tile = MORSE_zlange_Tile_c(norm,A,c_loc(work))
       end function MORSE_zlange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      function MORSE_zlanhe_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlanhe_Tile
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         MORSE_zlanhe_Tile = MORSE_zlanhe_Tile_c(norm,uplo,A,c_loc(work))
      end function MORSE_zlanhe_Tile
#endif

      function MORSE_zlansy_Tile(norm,uplo,A,work)
         use iso_c_binding
         implicit none
         real(kind=c_double) :: MORSE_zlansy_Tile
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         MORSE_zlansy_Tile = MORSE_zlansy_Tile_c(norm,uplo,A,c_loc(work))
      end function MORSE_zlansy_Tile

      subroutine MORSE_zlaswp_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine MORSE_zlaswp_Tile

      subroutine MORSE_zlauum_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlauum_Tile_c(uplo,A)
      end subroutine MORSE_zlauum_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zplghe_Tile(bump,uplo,A,seed,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplghe_Tile_c(bump,A,seed)
      end subroutine MORSE_zplghe_Tile
#endif

      subroutine MORSE_zplgsy_Tile(bump,uplo,A,seed,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplgsy_Tile_c(bump,A,seed)
      end subroutine MORSE_zplgsy_Tile
       
      subroutine MORSE_zplrnt_Tile(A,seed,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplrnt_Tile_c(A,seed)
      end subroutine MORSE_zplrnt_Tile
       
      subroutine MORSE_zposv_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zposv_Tile_c(uplo,A,B)
      end subroutine MORSE_zposv_Tile

      subroutine MORSE_zpotrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotrf_Tile_c(uplo,A)
      end subroutine MORSE_zpotrf_Tile

      subroutine MORSE_zsytrf_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsytrf_Tile_c(uplo,A)
      end subroutine MORSE_zsytrf_Tile

      subroutine MORSE_zpotri_Tile(uplo,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotri_Tile_c(uplo,A)
      end subroutine MORSE_zpotri_Tile

      subroutine MORSE_zpotrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotrs_Tile_c(uplo,A,B)
      end subroutine MORSE_zpotrs_Tile

      subroutine MORSE_zsytrs_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsytrs_Tile_c(uplo,A,B)
      end subroutine MORSE_zsytrs_Tile

      subroutine MORSE_zsymm_Tile(side,uplo,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C)
      end subroutine MORSE_zsymm_Tile

      subroutine MORSE_zsyrk_Tile(uplo,trans,alpha,A,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C)
      end subroutine MORSE_zsyrk_Tile

      subroutine MORSE_zsyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
      end subroutine MORSE_zsyr2k_Tile

      subroutine MORSE_ztrmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine MORSE_ztrmm_Tile

      subroutine MORSE_ztrsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine MORSE_ztrsm_Tile

      subroutine MORSE_ztrsmpl_Tile(A,L,IPIV,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsmpl_Tile_c(A,L,IPIV,B)
      end subroutine MORSE_ztrsmpl_Tile

      subroutine MORSE_ztrtri_Tile(uplo,diag,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrtri_Tile_c(uplo,diag,A)
      end subroutine MORSE_ztrtri_Tile

      subroutine MORSE_zunglq_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunglq_Tile_c(A,T,B)
      end subroutine MORSE_zunglq_Tile

      subroutine MORSE_zungqr_Tile(A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zungqr_Tile_c(A,T,B)
      end subroutine MORSE_zungqr_Tile

      subroutine MORSE_zunmlq_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunmlq_Tile_c(side,trans,A,T,B)
      end subroutine MORSE_zunmlq_Tile

      subroutine MORSE_zunmqr_Tile(side,trans,A,T,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunmqr_Tile_c(side,trans,A,T,B)
      end subroutine MORSE_zunmqr_Tile

      subroutine MORSE_zgetri_Tile(A,IPIV,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetri_Tile_c(A,c_loc(IPIV))
      end subroutine MORSE_zgetri_Tile

      subroutine MORSE_zlacpy_Tile(uplo,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlacpy_Tile_c(uplo,A,B)
      end subroutine MORSE_zlacpy_Tile

      subroutine MORSE_zlaset_Tile(uplo,alpha,beta,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaset_Tile_c(uplo,alpha,beta,A)
      end subroutine MORSE_zlaset_Tile

      subroutine MORSE_zlaswpc_Tile(A,K1,K2,IPIV,INCX,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaswpc_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
      end subroutine MORSE_zlaswpc_Tile

      subroutine MORSE_ztrsmrv_Tile(side,uplo,transA,diag,alpha,A,B,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B)
      end subroutine MORSE_ztrsmrv_Tile

      subroutine MORSE_zgetri_Tile_Async(A,IPIV,W,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: W ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetri_Tile_Async_c(A,c_loc(IPIV),W,sequence,request)
      end subroutine MORSE_zgetri_Tile_Async

      subroutine MORSE_zlange_Tile_Async(norm,A,work,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int) :: info
         real(kind=c_double), intent(out), target :: value
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlange_Tile_Async_c(norm,A,c_loc(work),c_loc(value),sequence,request)
      end subroutine MORSE_zlange_Tile_Async

      subroutine MORSE_zlansy_Tile_Async(norm,uplo,A,work,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int) :: info
         real(kind=c_double), intent(out), target :: value
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlansy_Tile_Async_c(norm,uplo,A,c_loc(work),c_loc(value),sequence,request)
      end subroutine MORSE_zlansy_Tile_Async

      subroutine MORSE_zgebrd_Tile_Async(A,D,E,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgebrd_Tile_Async_c(A,c_loc(D),c_loc(E),T,sequence,request)
      end subroutine MORSE_zgebrd_Tile_Async

      subroutine MORSE_zgecfi_Async(m,n,A,fin,imb,inb,fout,omb,onb,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: fout
         integer(kind=c_int), intent(in) :: imb
         integer(kind=c_int), intent(in) :: inb
         integer(kind=c_int), intent(in) :: omb
         integer(kind=c_int), intent(in) :: onb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgecfi_Async_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb,sequence,request)
      end subroutine MORSE_zgecfi_Async

      subroutine MORSE_zgelqf_Tile_Async(A,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqf_Tile_Async_c(A,T,sequence,request)
      end subroutine MORSE_zgelqf_Tile_Async

      subroutine MORSE_zgelqs_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgelqs_Tile_Async_c(A,T,B,sequence,request)
      end subroutine MORSE_zgelqs_Tile_Async

      subroutine MORSE_zgels_Tile_Async(trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgels_Tile_Async_c(trans,A,T,B,sequence,request)
      end subroutine MORSE_zgels_Tile_Async

      subroutine MORSE_zgemm_Tile_Async(transA,transB,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: transB
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,sequence,request)
      end subroutine MORSE_zgemm_Tile_Async

      subroutine MORSE_zgeqrf_Tile_Async(A,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrf_Tile_Async_c(A,T,sequence,request)
      end subroutine MORSE_zgeqrf_Tile_Async

      subroutine MORSE_zgeqrs_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgeqrs_Tile_Async_c(A,T,B,sequence,request)
      end subroutine MORSE_zgeqrs_Tile_Async

      subroutine MORSE_zgesv_Tile_Async(A,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesv_Tile_Async_c(A,c_loc(IPIV),B,sequence,request)
      end subroutine MORSE_zgesv_Tile_Async

      subroutine MORSE_zgesv_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine MORSE_zgesv_incpiv_Tile_Async

      subroutine MORSE_zgesvd_Tile_Async(jobu,jobvt,A,S,U,VT,T,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobu
         integer(kind=c_int), intent(in) :: jobvt
         real(kind=c_double), intent(out), target :: S(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: U ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: VT ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgesvd_Tile_Async_c(jobu,jobvt,A,c_loc(S),U,VT,T,sequence,request)
      end subroutine MORSE_zgesvd_Tile_Async

      subroutine MORSE_zgetmi_Async(m,n,A,fin,mb,nb,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(inout), target :: A(*)
         integer(kind=c_int), intent(in) :: fin
         integer(kind=c_int), intent(in) :: mb
         integer(kind=c_int), intent(in) :: nb
         integer(kind=c_int), intent(in) :: m
         integer(kind=c_int), intent(in) :: n
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetmi_Async_c(m,n,c_loc(A),fin,mb,nb,sequence,request)
      end subroutine MORSE_zgetmi_Async

      subroutine MORSE_zgetrf_Tile_Async(A,IPIV,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_Tile_Async_c(A,c_loc(IPIV),sequence,request)
      end subroutine MORSE_zgetrf_Tile_Async

      subroutine MORSE_zgetrf_incpiv_Tile_Async(A,L,IPIV,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request)
      end subroutine MORSE_zgetrf_incpiv_Tile_Async

      subroutine MORSE_zgetrf_nopiv_Tile_Async(A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrf_nopiv_Tile_Async_c(A,sequence,request)
      end subroutine MORSE_zgetrf_nopiv_Tile_Async

      subroutine MORSE_zgetrs_Tile_Async(trans,A,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrs_Tile_Async_c(trans,A,c_loc(IPIV),B,sequence,request)
      end subroutine MORSE_zgetrs_Tile_Async

      subroutine MORSE_zgetrs_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine MORSE_zgetrs_incpiv_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zhemm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
      end subroutine MORSE_zhemm_Tile_Async

      subroutine MORSE_zherk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
      end subroutine MORSE_zherk_Tile_Async

      subroutine MORSE_zher2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         real(kind=c_double), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
      end subroutine MORSE_zher2k_Tile_Async
#endif

      subroutine MORSE_zheev_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zheev_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
      end subroutine MORSE_zheev_Tile_Async

      subroutine MORSE_zheevd_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zheevd_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
      end subroutine MORSE_zheevd_Tile_Async

      subroutine MORSE_zhegv_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
      end subroutine MORSE_zhegv_Tile_Async

      subroutine MORSE_zhegvd_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         real(kind=c_double), intent(out), target :: W(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: Q ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
      end subroutine MORSE_zhegvd_Tile_Async

      subroutine MORSE_zhegst_Tile_Async(itype,uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: itype
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request)
      end subroutine MORSE_zhegst_Tile_Async

      subroutine MORSE_zhetrd_Tile_Async(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: jobz
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDQ
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         real(kind=c_double), intent(out), target :: D(*)
         real(kind=c_double), intent(out), target :: E(*)
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zhetrd_Tile_Async_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ,sequence,request)
      end subroutine MORSE_zhetrd_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zlanhe_Tile_Async(norm,uplo,A,work,value,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(inout), target :: work(*)
         integer(kind=c_int), intent(in) :: norm
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(out), target :: value
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlanhe_Tile_Async_c(norm,uplo,A,c_loc(work),c_loc(value),sequence,request)
      end subroutine MORSE_zlanhe_Tile_Async
#endif

      subroutine MORSE_zlaswp_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         integer(kind=c_int), intent(in), target :: IPIV(*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaswp_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
      end subroutine MORSE_zlaswp_Tile_Async

      subroutine MORSE_zlauum_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlauum_Tile_Async_c(uplo,A,sequence,request)
      end subroutine MORSE_zlauum_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine MORSE_zplghe_Tile_Async(bump,uplo,A,seed,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         real(kind=c_double), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplghe_Tile_Async_c(bump,A,seed,sequence,request)
      end subroutine MORSE_zplghe_Tile_Async
#endif

      subroutine MORSE_zplgsy_Tile_Async(bump,uplo,A,seed,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         complex(kind=c_double_complex), intent(in) :: bump
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplgsy_Tile_Async_c(bump,A,seed,sequence,request)
      end subroutine MORSE_zplgsy_Tile_Async
       
      subroutine MORSE_zplrnt_Tile_Async(A,seed,sequence,request,info) 
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_long_long), intent(in) :: seed
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zplrnt_Tile_Async_c(A,seed,sequence,request)
      end subroutine MORSE_zplrnt_Tile_Async

      subroutine MORSE_zposv_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zposv_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine MORSE_zposv_Tile_Async

      subroutine MORSE_zpotrf_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotrf_Tile_Async_c(uplo,A,sequence,request)
      end subroutine MORSE_zpotrf_Tile_Async

      subroutine MORSE_zsytrf_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsytrf_Tile_Async_c(uplo,A,sequence,request)
      end subroutine MORSE_zsytrf_Tile_Async

      subroutine MORSE_zpotri_Tile_Async(uplo,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotri_Tile_Async_c(uplo,A,sequence,request)
      end subroutine MORSE_zpotri_Tile_Async

      subroutine MORSE_zpotrs_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zpotrs_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine MORSE_zpotrs_Tile_Async

      subroutine MORSE_zsytrs_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsytrs_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine MORSE_zsytrs_Tile_Async

      subroutine MORSE_zsymm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
      end subroutine MORSE_zsymm_Tile_Async

      subroutine MORSE_zsyrk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
      end subroutine MORSE_zsyrk_Tile_Async

      subroutine MORSE_zsyr2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: C ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
      end subroutine MORSE_zsyr2k_Tile_Async

      subroutine MORSE_ztrmm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine MORSE_ztrmm_Tile_Async

      subroutine MORSE_ztrsm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine MORSE_ztrsm_Tile_Async

      subroutine MORSE_ztrsmpl_Tile_Async(A,L,IPIV,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: L ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: IPIV ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request)
      end subroutine MORSE_ztrsmpl_Tile_Async

      subroutine MORSE_ztrtri_Tile_Async(uplo,diag,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request)
      end subroutine MORSE_ztrtri_Tile_Async

      subroutine MORSE_zunglq_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunglq_Tile_Async_c(A,T,B,sequence,request)
      end subroutine MORSE_zunglq_Tile_Async

      subroutine MORSE_zungqr_Tile_Async(A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zungqr_Tile_Async_c(A,T,B,sequence,request)
      end subroutine MORSE_zungqr_Tile_Async

      subroutine MORSE_zunmlq_Tile_Async(side,trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request)
      end subroutine MORSE_zunmlq_Tile_Async

      subroutine MORSE_zunmqr_Tile_Async(side,trans,A,T,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request)
      end subroutine MORSE_zunmqr_Tile_Async

      subroutine MORSE_zlacpy_Tile_Async(uplo,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlacpy_Tile_Async_c(uplo,A,B,sequence,request)
      end subroutine MORSE_zlacpy_Tile_Async

      subroutine MORSE_zlaset_Tile_Async(uplo,alpha,beta,A,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         complex(kind=c_double_complex), intent(in) :: beta
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request)
      end subroutine MORSE_zlaset_Tile_Async

      subroutine MORSE_zlaswpc_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: INCX
         integer(kind=c_int), intent(in) :: K1
         integer(kind=c_int), intent(in) :: K2
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zlaswpc_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
      end subroutine MORSE_zlaswpc_Tile_Async

      subroutine MORSE_ztrsmrv_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: diag
         integer(kind=c_int), intent(in) :: side
         integer(kind=c_int), intent(in) :: transA
         integer(kind=c_int), intent(in) :: uplo
         complex(kind=c_double_complex), intent(in) :: alpha
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
      end subroutine MORSE_ztrsmrv_Tile_Async

      subroutine MORSE_Alloc_Workspace_zgelqf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = MORSE_Alloc_Workspace_zgelqf_c(M,N,T)
      end subroutine MORSE_Alloc_Workspace_zgelqf

      subroutine MORSE_Alloc_Workspace_zgels(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = MORSE_Alloc_Workspace_zgels_c(M,N,T)
      end subroutine MORSE_Alloc_Workspace_zgels

      subroutine MORSE_Alloc_Workspace_zgeqrf(M,N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = MORSE_Alloc_Workspace_zgeqrf_c(M,N,T)
      end subroutine MORSE_Alloc_Workspace_zgeqrf

      subroutine MORSE_Alloc_Workspace_zgesv_incpiv(N,L,IPIV,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgesv_incpiv_c(N,L,IPIV,p,q)
      end subroutine MORSE_Alloc_Workspace_zgesv_incpiv

      subroutine MORSE_Alloc_Workspace_zgetrf_incpiv(M,N,L,IPIV,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: L ! L is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgetrf_incpiv_c(M,N,L,IPIV,p,q)
      end subroutine MORSE_Alloc_Workspace_zgetrf_incpiv

      subroutine MORSE_Alloc_Workspace_zgeev(N,T,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         info = MORSE_Alloc_Workspace_zgeev_c(N,T)
      end subroutine MORSE_Alloc_Workspace_zgeev

      subroutine MORSE_Alloc_Workspace_zgebrd(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgebrd_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zgebrd

      subroutine MORSE_Alloc_Workspace_zgehrd(N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgehrd_c(N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zgehrd

      subroutine MORSE_Alloc_Workspace_zgesvd(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgesvd_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zgesvd

      subroutine MORSE_Alloc_Workspace_zheev(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zheev_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zheev

      subroutine MORSE_Alloc_Workspace_zheevd(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zheevd_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zheevd

      subroutine MORSE_Alloc_Workspace_zhegv(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zhegv_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zhegv

      subroutine MORSE_Alloc_Workspace_zhegvd(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zhegvd_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zhegvd

      subroutine MORSE_Alloc_Workspace_zhetrd(M,N,T,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: T ! T is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zhetrd_c(M,N,T,p,q)
      end subroutine MORSE_Alloc_Workspace_zhetrd

      subroutine MORSE_Alloc_Workspace_zgelqf_Tile(M,N,descT,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgelqf_Tile_c(M,N,descT,p,q)
      end subroutine MORSE_Alloc_Workspace_zgelqf_Tile

      subroutine MORSE_Alloc_Workspace_zgels_Tile(M,N,descT,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgels_Tile_c(M,N,descT,p,q)
      end subroutine MORSE_Alloc_Workspace_zgels_Tile

      subroutine MORSE_Alloc_Workspace_zgeqrf_Tile(M,N,descT,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: M
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: descT ! descT is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT,p,q)
      end subroutine MORSE_Alloc_Workspace_zgeqrf_Tile

      subroutine MORSE_Alloc_Workspace_zgesv_incpiv_Tile(N,descL,IPIV,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV,p,q)
      end subroutine MORSE_Alloc_Workspace_zgesv_incpiv_Tile

      subroutine MORSE_Alloc_Workspace_zgetrf_incpiv_Tile(N,descL,IPIV,p,q,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: N
         type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
         type(c_ptr) :: descL ! descL is **, so pass by reference
         integer(kind=c_int), value :: p
         integer(kind=c_int), value :: q
         info = MORSE_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV,p,q)
      end subroutine MORSE_Alloc_Workspace_zgetrf_incpiv_Tile

      subroutine MORSE_Alloc_Workspace_zgetri_Tile_Async(A,W,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         type(c_ptr), value :: W ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_Alloc_Workspace_zgetri_Tile_Async_c(A,W)
      end subroutine MORSE_Alloc_Workspace_zgetri_Tile_Async

      subroutine MORSE_zLapack_to_Tile(Af77,LDA,A,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(in), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zLapack_to_Tile_c(c_loc(Af77),LDA,A)
      end subroutine MORSE_zLapack_to_Tile

      subroutine MORSE_zTile_to_Lapack(A,Af77,LDA,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(in) :: LDA
         complex(kind=c_double_complex), intent(out), target :: Af77(LDA,*)
         type(c_ptr), value :: A ! Arg managed by MORSE: opaque to Fortran
         info = MORSE_zTile_to_Lapack_c(A,c_loc(Af77),LDA)
      end subroutine MORSE_zTile_to_Lapack

end module morse_z
