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
! Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
! - Redistributions of source code must retain the above copyright notice,
!   this list of conditions, and the following disclaimer.
! - Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer listed in
!   this license in the documentation and/or other materials provided with
!   the distribution.
! - Neither the name of the copyright holders nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! This software is provided by the copyright holders and contributors "as
! is" and any express or implied warranties, including, but not limited
! to, the implied warranties of merchantability and fitness for a
! particular purpose are disclaimed. in no event shall the copyright owner
! or contributors be liable for any direct, indirect, incidental, special,
! exemplary, or consequential damages (including, but not limited to,
! procurement of substitute goods or services; loss of use, data, or
! profits; or business interruption) however caused and on any theory of
! liability, whether in contract, strict liability, or tort (including
! negligence or otherwise) arising in any way out of the use of this
! software, even if advised of the possibility of such damage.
!
!
!
!  MORSE fortran 90 interface
!  MORSE is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 1.0.0
! @author Numerical Algorithm Group
! @date 2011-09-15
! @precisions normal z -> c d s
!
module morse

      use morse_s
      use morse_d
      use morse_ds
      use morse_c
      use morse_z
      use morse_zc
      include 'morse_fortran.h'

      logical :: morse_initialized = .false.
      integer, parameter :: sp = kind(0.0)
      integer, parameter :: dp = kind(0.0d0)

      interface
         function MORSE_Init_c(cpus, gpus) &
          & bind(c, name='MORSE_Init')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Init_c
            integer(kind=c_int), value :: cpus, gpus
         end function MORSE_Init_c
      end interface

      interface
         function MORSE_Finalize_c() &
          & bind(c, name='MORSE_Finalize')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Finalize_c
         end function MORSE_Finalize_c
      end interface

      interface
         function MORSE_Set_c(param, pval) &
          & bind(c, name='MORSE_Set')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Set_c
            integer(kind=c_int), value :: param
            integer(kind=c_int), value :: pval
         end function MORSE_Set_c
      end interface

      interface
         function MORSE_Get_c(param, pval) &
          & bind(c, name='MORSE_Get')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Get_c
            integer(kind=c_int), value :: param
            type(c_ptr), value :: pval
         end function MORSE_Get_c
      end interface

      interface
         function MORSE_Enable_c(param) &
          & bind(c, name='MORSE_Enable')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Enable_c
            integer(kind=c_int), value :: param
         end function MORSE_Enable_c
      end interface

      interface
         function MORSE_Disable_c(param) &
          & bind(c, name='MORSE_Disable')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: MORSE_Disable_c
            integer(kind=c_int), value :: param
         end function MORSE_Disable_c
      end interface

      interface
         function MORSE_Lapack_to_Tile_c(a_lpk,lda,a_pma) &
          & bind(c, name='MORSE_Lapack_to_Tile')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Lapack_to_Tile_c
            type(c_ptr), value :: a_lpk, a_pma
            integer(kind=c_int), value :: lda
         end function MORSE_Lapack_to_Tile_c
      end interface

      interface
         function MORSE_Tile_to_Lapack_c(a_pma,a_lpk,lda) &
          & bind(c, name='MORSE_Tile_to_Lapack')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Tile_to_Lapack_c
            type(c_ptr), value :: a_lpk, a_pma
            integer(kind=c_int), value :: lda
         end function MORSE_Tile_to_Lapack_c
      end interface

      interface
         function MORSE_Desc_Create_c(desc, mat, dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q) &
          & bind(c, name='MORSE_Desc_Create')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Desc_Create_c
            type(c_ptr) :: desc
            type(c_ptr), value :: mat
            integer(kind=c_int), value :: dtyp
            integer(kind=c_int), value :: mb, nb, bsiz, lm, ln, i, j, m, n,p, q
         end function MORSE_Desc_Create_c
      end interface

      interface
         function MORSE_Desc_Create_OOC_c(desc, dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q) &
          & bind(c, name='MORSE_Desc_Create_OOC')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Desc_Create_OOC_c
            type(c_ptr) :: desc
            integer(kind=c_int), value :: dtyp
            integer(kind=c_int), value :: mb, nb, bsiz, lm, ln, i, j, m, n,p, q
         end function MORSE_Desc_Create_OOC_c
      end interface

      interface
         function MORSE_Desc_Create_User_c(desc, mat, dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q, get_blkaddr, get_blkldd, get_rankof) &
          & bind(c, name='MORSE_Desc_Create_User')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Desc_Create_User_c
            type(c_ptr) :: desc
            type(c_ptr), value :: mat
            integer(kind=c_int), value :: dtyp
            integer(kind=c_int), value :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
            type(c_ptr) :: get_blkaddr
            type(c_ptr), value :: get_blkldd, get_rankof
         end function MORSE_Desc_Create_User_c
      end interface

      interface
         function MORSE_Desc_Create_OOC_User_c(desc, dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q, get_rankof) &
          & bind(c, name='MORSE_Desc_Create_OOC_User')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Desc_Create_OOC_User_c
            type(c_ptr) :: desc
            integer(kind=c_int), value :: dtyp
            integer(kind=c_int), value :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
            type(c_ptr), value :: get_rankof
         end function MORSE_Desc_Create_OOC_User_c
      end interface

      interface
         function MORSE_Desc_Destroy_c(desc) &
          & bind(c, name='MORSE_Desc_Destroy')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Desc_Destroy_c
            type(c_ptr) :: desc
         end function MORSE_Desc_Destroy_c
      end interface

      interface
         subroutine free_c(ptr) bind(c, name='free')
            use iso_c_binding
            implicit none
            type(c_ptr), value :: ptr
         end subroutine free_c
      end interface

      interface
         function MORSE_Version_c(maj,min,mic) &
          & bind(c, name='MORSE_Version')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Version_c
            type(c_ptr), value ::  maj,min,mic
         end function MORSE_Version_c
      end interface

      interface
         function MORSE_Init_Affinity_c(cores,gpus,bindtab) &
          & bind(c, name='MORSE_Init_Affinity')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Init_Affinity_c
            integer(kind=c_int), value ::  cores, gpus
            type(c_ptr), value :: bindtab
         end function MORSE_Init_Affinity_c
      end interface

      interface
         function MORSE_Dealloc_Handle_c(handle) &
          & bind(c, name='MORSE_Dealloc_Handle')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Dealloc_Handle_c
            type(c_ptr) :: handle
         end function MORSE_Dealloc_Handle_c
      end interface

      interface
         function MORSE_Dealloc_Handle_Tile_c(desc) &
          & bind(c, name='MORSE_Dealloc_Handle_Tile')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Dealloc_Handle_Tile_c
            type(c_ptr) :: desc
         end function MORSE_Dealloc_Handle_Tile_c
      end interface

      interface
         function MORSE_Sequence_Create_c(seq) &
          & bind(c, name='MORSE_Sequence_Create')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Sequence_Create_c
            type(c_ptr) :: seq
         end function MORSE_Sequence_Create_c
      end interface

      interface
         function MORSE_Sequence_Destroy_c(seq) &
          & bind(c, name='MORSE_Sequence_Destroy')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Sequence_Destroy_c
            type(c_ptr), value :: seq
         end function MORSE_Sequence_Destroy_c
      end interface

      interface
         function MORSE_Sequence_Wait_c(seq) &
          & bind(c, name='MORSE_Sequence_Wait')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Sequence_Wait_c
            type(c_ptr), value :: seq
         end function MORSE_Sequence_Wait_c
      end interface

      interface
         function MORSE_Sequence_Flush_c(seq,req) &
          & bind(c, name='MORSE_Sequence_Flush')
            use iso_c_binding
            integer(kind=c_int) :: MORSE_Sequence_Flush_c
            type(c_ptr), value :: seq
            type(c_ptr), value :: req
         end function MORSE_Sequence_Flush_c
      end interface

      interface morse_lapack_to_tile
         module procedure morse_lapack_to_tile_s
         module procedure morse_lapack_to_tile_d
         module procedure morse_lapack_to_tile_cpx
         module procedure morse_lapack_to_tile_z
      end interface morse_lapack_to_tile

      interface morse_tile_to_lapack
         module procedure morse_tile_to_lapack_s
         module procedure morse_tile_to_lapack_d
         module procedure morse_tile_to_lapack_cpx
         module procedure morse_tile_to_lapack_z
      end interface morse_tile_to_lapack

      interface morse_desc_create
         module procedure morse_desc_create_s
         module procedure morse_desc_create_d
         module procedure morse_desc_create_cpx
         module procedure morse_desc_create_z
      end interface morse_desc_create

   contains

   subroutine morse_init(cores,gpus,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: cores, gpus
      integer(kind=c_int), intent(out) :: info
      info = morse_init_c(cores,gpus)
      morse_initialized = .true.
   end subroutine morse_init

   subroutine morse_finalize(info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(out) :: info
      info = morse_finalize_c()
      morse_initialized = .false.
   end subroutine morse_finalize

   subroutine morse_set(param,pval,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(in) :: pval
      integer(kind=c_int), intent(out) :: info
      info = morse_set_c(param,pval)
   end subroutine morse_set

   subroutine morse_get(param,pval,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out), target :: pval
      integer(kind=c_int), intent(out) :: info
      info = morse_get_c(param,c_loc(pval))
   end subroutine morse_get

   subroutine morse_enable(param,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out) :: info
      info = morse_enable_c(param)
   end subroutine morse_enable

   subroutine morse_disable(param,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: param
      integer(kind=c_int), intent(out) :: info
      info = morse_disable_c(param)
   end subroutine morse_disable

! overloaded: single precision
   subroutine morse_lapack_to_tile_s(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine morse_lapack_to_tile_s
! overloaded: double precision
   subroutine morse_lapack_to_tile_d(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine morse_lapack_to_tile_d
! overloaded: single precision complex
   subroutine morse_lapack_to_tile_cpx(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine morse_lapack_to_tile_cpx
! overloaded: double precision complex
   subroutine morse_lapack_to_tile_z(a_lpk,lda,a_pma,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(out) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_lapack_to_tile_c(c_loc(a_lpk),lda,a_pma)
   end subroutine morse_lapack_to_tile_z

! overloaded: single precision
   subroutine morse_tile_to_lapack_s(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine morse_tile_to_lapack_s
! overloaded: double precision
   subroutine morse_tile_to_lapack_d(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      real(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine morse_tile_to_lapack_d
! overloaded: single precision complex
   subroutine morse_tile_to_lapack_cpx(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=sp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine morse_tile_to_lapack_cpx
! overloaded: double precision complex
   subroutine morse_tile_to_lapack_z(a_pma,a_lpk,lda,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: lda
      complex(kind=dp), intent(out), target :: a_lpk(lda,*)
      type(c_ptr), intent(in) ::  a_pma
      integer(kind=c_int), intent(out) :: info
      info = morse_tile_to_lapack_c(a_pma,c_loc(a_lpk),lda)
   end subroutine morse_tile_to_lapack_z

! overloaded: single precision
   subroutine morse_desc_create_s(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
      real(kind=sp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = morse_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q)
   end subroutine morse_desc_create_s
! overloaded: double precision
   subroutine morse_desc_create_d(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
      real(kind=dp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = morse_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q)
   end subroutine morse_desc_create_d
! overloaded: single precision complex
   subroutine morse_desc_create_cpx(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
      complex(kind=sp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = morse_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q)
   end subroutine morse_desc_create_cpx
! overloaded: double precision complex
   subroutine morse_desc_create_z(desc,mat,dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: desc
      integer(kind=c_int), intent(in) :: mb, nb, bsiz, lm, ln, i, j, m, n, p, q
      complex(kind=dp), intent(in), target :: mat(lm,*)
      integer(kind=c_int), intent(in) :: dtyp
      integer(kind=c_int), intent(out) :: info
      info = morse_desc_create_c(desc,c_loc(mat),dtyp,mb,nb,bsiz,lm,ln,i,j,m,n,p,q)
   end subroutine morse_desc_create_z

   subroutine morse_desc_destroy(desc,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: desc
      integer(kind=c_int), intent(out) :: info
      info = morse_desc_destroy_c(desc)
   end subroutine morse_desc_destroy

   subroutine morse_free(ptr)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: ptr
      call free_c(ptr)
   end subroutine morse_free

   subroutine morse_version(ver_major,ver_minor,ver_micro,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(out), target :: ver_major,ver_minor,ver_micro
      integer(kind=c_int), intent(out) :: info
      info = morse_version_c(c_loc(ver_major),c_loc(ver_minor),c_loc(ver_micro))
   end subroutine morse_version

   subroutine morse_init_affinity(cores,bindtab,info)
      use iso_c_binding
      implicit none
      integer(kind=c_int), intent(in) :: cores
      integer(kind=c_int), intent(out), target :: bindtab
      integer(kind=c_int), intent(out) :: info
      info = morse_init_affinity_c(cores,c_loc(bindtab))
   end subroutine morse_init_affinity

   subroutine morse_dealloc_handle(handle,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: handle
      integer(kind=c_int), intent(out) :: info
      info = morse_dealloc_handle_c(handle)
   end subroutine morse_dealloc_handle

   subroutine morse_dealloc_handle_tile(desc,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: desc
      integer(kind=c_int), intent(out) :: info
      info = morse_dealloc_handle_tile_c(desc)
   end subroutine morse_dealloc_handle_tile

   subroutine morse_sequence_create(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = morse_sequence_create_c(sequence)
   end subroutine morse_sequence_create

   subroutine morse_sequence_destroy(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = morse_sequence_destroy_c(sequence)
   end subroutine morse_sequence_destroy

   subroutine morse_sequence_wait(sequence,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence
      integer(kind=c_int), intent(out) :: info
      info = morse_sequence_wait_c(sequence)
   end subroutine morse_sequence_wait

   subroutine morse_sequence_flush(sequence,request,info)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: sequence, request
      integer(kind=c_int), intent(out) :: info
      info = morse_sequence_flush_c(sequence,request)
   end subroutine morse_sequence_flush

end module morse
