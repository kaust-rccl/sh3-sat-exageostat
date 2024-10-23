!**
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
!**
!
! @brief Chameleon Fortran77 header
! @version 1.0.0
! @author Bilel Hadri
! @author Mathieu Faverge
! @author Cedric Castagnede
! @date 2010-11-15
!
!**

!********************************************************************
!   MORSE constants - precisions
!
      integer  MorseByte, MorseInteger, MorseRealFloat
      integer  MorseRealDouble, MorseComplexFloat, MorseComplexDouble
      parameter ( MorseByte          = 0 )
      parameter ( MorseInteger       = 1 )
      parameter ( MorseRealFloat     = 2 )
      parameter ( MorseRealDouble    = 3 )
      parameter ( MorseComplexFloat  = 4 )
      parameter ( MorseComplexDouble = 5 )

!********************************************************************
!   MORSE constants - CBLAS & LAPACK
!
      integer MorseCM, MorseRM, MorseCCRB
      integer MorseCRRB, MorseRCRB, MorseRRRB
      parameter ( MorseCM         = 101 )
      parameter ( MorseRM         = 102 )
      parameter ( MorseCCRB       = 103 )
      parameter ( MorseCRRB       = 104 )
      parameter ( MorseRCRB       = 105 )
      parameter ( MorseRRRB       = 106 )

      integer  MorseNoTrans, MorseTrans, MorseConjTrans
      parameter ( MorseNoTrans    = 111 )
      parameter ( MorseTrans      = 112 )
      parameter ( MorseConjTrans  = 113 )

      integer MorseUpper, MorseLower
      integer MorseUpperLower
      parameter ( MorseUpper      = 121 )
      parameter ( MorseLower      = 122 )
      parameter ( MorseUpperLower = 123 )

      integer MorseNonUnit,MorseUnit
      parameter ( MorseNonUnit    = 131 )
      parameter ( MorseUnit       = 132 )

      integer MorseLeft,MorseRight
      parameter ( MorseLeft       = 141 )
      parameter ( MorseRight      = 142 )

      integer MorseOneNorm, MorseRealOneNorm
      integer MorseTwoNorm, MorseFrobeniusNorm
      integer MorseInfNorm, MorseRealInfNorm
      integer MorseMaxNorm, MorseRealMaxNorm
      parameter ( MorseOneNorm       = 171 )
      parameter ( MorseRealOneNorm   = 172 )
      parameter ( MorseTwoNorm       = 173 )
      parameter ( MorseFrobeniusNorm = 174 )
      parameter ( MorseInfNorm       = 175 )
      parameter ( MorseRealInfNorm   = 176 )
      parameter ( MorseMaxNorm       = 177 )
      parameter ( MorseRealMaxNorm   = 178 )

      integer MorseDistUniform
      integer MorseDistSymmetric
      integer MorseDistNormal
      parameter ( MorseDistUniform   = 201 )
      parameter ( MorseDistSymmetric = 202 )
      parameter ( MorseDistNormal    = 203 )

      integer MorseHermGeev
      integer MorseHermPoev
      integer MorseNonsymPosv
      integer MorseSymPosv
      parameter ( MorseHermGeev    = 241 )
      parameter ( MorseHermPoev    = 242 )
      parameter ( MorseNonsymPosv  = 243 )
      parameter ( MorseSymPosv     = 244 )

      integer MorseNoPacking
      integer MorsePackSubdiag
      integer MorsePackSupdiag
      integer MorsePackColumn
      integer MorsePackLowerBand
      integer MorsePackRow
      integer MorsePackUpeprBand
      integer MorsePackAll
      parameter ( MorseNoPacking     = 291 )
      parameter ( MorsePackSubdiag   = 292 )
      parameter ( MorsePackSupdiag   = 293 )
      parameter ( MorsePackColumn    = 294 )
      parameter ( MorsePackRow       = 295 )
      parameter ( MorsePackLowerBand = 296 )
      parameter ( MorsePackUpeprBand = 297 )
      parameter ( MorsePackAll       = 298 )

      integer MorseNoVec,MorseVec,MorseIvec
      parameter ( MorseNoVec = 301 )
      parameter ( MorseVec   = 302 )
      parameter ( MorseIvec  = 303 )

      integer MorseForward, MorseBackward
      parameter ( MorseForward    = 391 )
      parameter ( MorseBackward   = 392 )

      integer MorseColumnwise,MorseRowwise
      parameter ( MorseColumnwise = 401 )
      parameter ( MorseRowwise    = 402 )

!********************************************************************
!   MORSE constants - boolean
!
      integer MORSE_FALSE, MORSE_TRUE
      parameter ( MORSE_FALSE = 0 )
      parameter ( MORSE_TRUE  = 1 )

!********************************************************************
!   State machine switches
!
      integer MORSE_WARNINGS, MORSE_ERRORS, MORSE_AUTOTUNING
      integer MORSE_DAG, MORSE_PROFILING_MODE, MORSE_PARALLEL_MODE
      integer MORSE_BOUND
      parameter ( MORSE_WARNINGS       = 1 )
      parameter ( MORSE_ERRORS         = 2 )
      parameter ( MORSE_AUTOTUNING     = 3 )
      parameter ( MORSE_DAG            = 4 )
      parameter ( MORSE_PROFILING_MODE = 5 )
      parameter ( MORSE_PARALLEL_MODE  = 6 )
      parameter ( MORSE_BOUND          = 7 )

!********************************************************************
!   MORSE constants - configuration  parameters
!
      integer MORSE_CONCURRENCY, MORSE_TILE_SIZE
      integer MORSE_INNER_BLOCK_SIZE, MORSE_SCHEDULING_MODE
      integer MORSE_HOUSEHOLDER_MODE, MORSE_HOUSEHOLDER_SIZE
      integer MORSE_TRANSLATION_MODE
      parameter ( MORSE_CONCURRENCY      = 1 )
      parameter ( MORSE_TILE_SIZE        = 2 )
      parameter ( MORSE_INNER_BLOCK_SIZE = 3 )
      parameter ( MORSE_SCHEDULING_MODE  = 4 )
      parameter ( MORSE_HOUSEHOLDER_MODE = 5 )
      parameter ( MORSE_HOUSEHOLDER_SIZE = 6 )
      parameter ( MORSE_TRANSLATION_MODE = 7 )

!********************************************************************
!   MORSE constants - householder mode
!
      integer MORSE_FLAT_HOUSEHOLDER, MORSE_TREE_HOUSEHOLDER
      parameter ( MORSE_FLAT_HOUSEHOLDER  = 1 )
      parameter ( MORSE_TREE_HOUSEHOLDER  = 2 )

!*********************************************************************
!   MORSE constants - translation mode
!
      integer MORSE_INPLACE, MORSE_OUTOFPLACE
      parameter ( MORSE_INPLACE     = 1 )
      parameter ( MORSE_OUTOFPLACE  = 2 )

!********************************************************************
!   MORSE constants - success & error codes
!
      integer MORSE_SUCCESS, MORSE_ERR_NOT_INITIALIZED
      integer MORSE_ERR_REINITIALIZED, MORSE_ERR_NOT_SUPPORTED
      integer MORSE_ERR_ILLEGAL_VALUE, MORSE_ERR_NOT_FOUND
      integer MORSE_ERR_OUT_OF_MEMORY, MORSE_ERR_INTERNAL_LIMIT
      integer MORSE_ERR_UNALLOCATED, MORSE_ERR_FILESYSTEM
      integer MORSE_ERR_UNEXPECTED, MORSE_ERR_SEQUENCE_FLUSHED
      parameter ( MORSE_SUCCESS             =    0 )
      parameter ( MORSE_ERR_NOT_INITIALIZED = -101 )
      parameter ( MORSE_ERR_REINITIALIZED   = -102 )
      parameter ( MORSE_ERR_NOT_SUPPORTED   = -103 )
      parameter ( MORSE_ERR_ILLEGAL_VALUE   = -104 )
      parameter ( MORSE_ERR_NOT_FOUND       = -105 )
      parameter ( MORSE_ERR_OUT_OF_MEMORY   = -106 )
      parameter ( MORSE_ERR_INTERNAL_LIMIT  = -107 )
      parameter ( MORSE_ERR_UNALLOCATED     = -108 )
      parameter ( MORSE_ERR_FILESYSTEM      = -109 )
      parameter ( MORSE_ERR_UNEXPECTED      = -110 )
      parameter ( MORSE_ERR_SEQUENCE_FLUSHED= -111 )

!********************************************************************
!   MORSE constants - kernels options
!
      integer MORSE_PRIORITY_MIN, MORSE_PRIORITY_MAX
      parameter ( MORSE_PRIORITY_MIN = 0     )
      parameter ( MORSE_PRIORITY_MAX = 10000 )

!********************************************************************
!   MORSE constants - scheduler properties
!
      integer PRIORITY, CALLBACK, REDUX
      parameter ( PRIORITY = 16 )
      parameter ( CALLBACK = 17 )
      parameter ( REDUX    = 18 )
