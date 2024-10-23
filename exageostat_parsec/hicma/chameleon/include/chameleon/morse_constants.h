/**
 *
 * @file morse_constants.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global constants
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _MORSE_CONSTANTS_H_
#define _MORSE_CONSTANTS_H_

/**
 *
 * @brief Chameleon constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 */
#define MorseByte              0
#define MorseInteger           1
#define MorseRealFloat         2
#define MorseRealDouble        3
#define MorseComplexFloat      4
#define MorseComplexDouble     5

#define MorseCM              101
#define MorseRM              102
#define MorseCCRB            103
#define MorseCRRB            104
#define MorseRCRB            105
#define MorseRRRB            106

#define MorseNoTrans         111
#define MorseTrans           112
#define MorseConjTrans       113

#define MorseUpper           121
#define MorseLower           122
#define MorseUpperLower      123

#define MorseNonUnit         131
#define MorseUnit            132

#define MorseLeft            141
#define MorseRight           142

#define MorseOneNorm         171
#define MorseRealOneNorm     172
#define MorseTwoNorm         173
#define MorseFrobeniusNorm   174
#define MorseInfNorm         175
#define MorseRealInfNorm     176
#define MorseMaxNorm         177
#define MorseRealMaxNorm     178

#define MorseDistUniform     201
#define MorseDistSymmetric   202
#define MorseDistNormal      203

#define MorseHermGeev        241
#define MorseHermPoev        242
#define MorseNonsymPosv      243
#define MorseSymPosv         244

#define MorseNoPacking       291
#define MorsePackSubdiag     292
#define MorsePackSupdiag     293
#define MorsePackColumn      294
#define MorsePackRow         295
#define MorsePackLowerBand   296
#define MorsePackUpeprBand   297
#define MorsePackAll         298

#define MorseNoVec           301
#define MorseVec             302
#define MorseIvec            303

#define MorseForward         391
#define MorseBackward        392

#define MorseColumnwise      401
#define MorseRowwise         402
#define MorseTrd            1001
#define MorseBrd            1002

#define MorseW               501
#define MorseA2              502

#define morse_const_neg(const) (((const-1)^0x01)+1)

/**
 *  MORSE constants - boolean
 */
#define MORSE_FALSE  0
#define MORSE_TRUE   1

#define MORSE_CPU    ((1ULL)<<1)
#define MORSE_CUDA   ((1ULL)<<3)

/**
 *  State machine switches
 */
#define MORSE_WARNINGS        1
#define MORSE_ERRORS          2
#define MORSE_AUTOTUNING      3
#define MORSE_DAG             4
#define MORSE_PROFILING_MODE  5
#define MORSE_PARALLEL_MODE   6
#define MORSE_BOUND           7
#define MORSE_PROGRESS        8
#define MORSE_GEMM3M          9

/**
 *  MORSE constants - configuration parameters
 */
#define MORSE_CONCURRENCY       1
#define MORSE_TILE_SIZE         2
#define MORSE_INNER_BLOCK_SIZE  3
#define MORSE_HOUSEHOLDER_MODE  5
#define MORSE_HOUSEHOLDER_SIZE  6
#define MORSE_TRANSLATION_MODE  7

#define MORSE_FLAT_HOUSEHOLDER  1
#define MORSE_TREE_HOUSEHOLDER  2

#define MORSE_INPLACE           1
#define MORSE_OUTOFPLACE        2

/**
 *  MORSE constants - success & error codes
 */
#define MORSE_SUCCESS                 0
#define MORSE_ERR_NOT_INITIALIZED  -101
#define MORSE_ERR_REINITIALIZED    -102
#define MORSE_ERR_NOT_SUPPORTED    -103
#define MORSE_ERR_ILLEGAL_VALUE    -104
#define MORSE_ERR_NOT_FOUND        -105
#define MORSE_ERR_OUT_OF_RESOURCES -106
#define MORSE_ERR_INTERNAL_LIMIT   -107
#define MORSE_ERR_UNALLOCATED      -108
#define MORSE_ERR_FILESYSTEM       -109
#define MORSE_ERR_UNEXPECTED       -110
#define MORSE_ERR_SEQUENCE_FLUSHED -111

/**
 * Kernels options
 */
#define MORSE_PRIORITY_MIN  0
#define MORSE_PRIORITY_MAX  INT_MAX


/**
 *  Scheduler properties
 */
#define PRIORITY        16
#define CALLBACK        17
#define REDUX           18

/**
 *  MORSE ???
 */
#define MORSE_REQUEST_INITIALIZER {MORSE_SUCCESS}

#endif
