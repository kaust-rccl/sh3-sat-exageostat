/**
 *
 * @file cudaglobal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global cudablas variables and functions
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2017-04-06
 *
 */
#include "cudablas.h"

/**
 *  LAPACK Constants
 */
int morse_cublas_constants[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 100
    0,                      // 101: MorseRowMajor
    0,                      // 102: MorseColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: MorseNoTrans
    CUBLAS_OP_T,            // 112: MorseTrans
    CUBLAS_OP_C,            // 113: MorseConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: MorseUpper
    CUBLAS_FILL_MODE_LOWER, // 122: MorseLower
    0,                      // 123: MorseUpperLower
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: MorseNonUnit
    CUBLAS_DIAG_UNIT,       // 132: MorseUnit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: MorseLeft
    CUBLAS_SIDE_RIGHT,      // 142: MorseRight
    0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 151:
    0,                      // 152:
    0,                      // 153:
    0,                      // 154:
    0,                      // 155:
    0,                      // 156:
    0,                      // 157: MorseEps
    0,                      // 158:
    0,                      // 159:
    0,                      // 160:
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 171: MorseOneNorm
    0,                      // 172: MorseRealOneNorm
    0,                      // 173: MorseTwoNorm
    0,                      // 174: MorseFrobeniusNorm
    0,                      // 175: MorseInfNorm
    0,                      // 176: MorseRealInfNorm
    0,                      // 177: MorseMaxNorm
    0,                      // 178: MorseRealMaxNorm
    0,                      // 179
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 200
    0,                      // 201: MorseDistUniform
    0,                      // 202: MorseDistSymmetric
    0,                      // 203: MorseDistNormal
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 240
    0,                      // 241 MorseHermGeev
    0,                      // 242 MorseHermPoev
    0,                      // 243 MorseNonsymPosv
    0,                      // 244 MorseSymPosv
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 290
    0,                      // 291 MorseNoPacking
    0,                      // 292 MorsePackSubdiag
    0,                      // 293 MorsePackSupdiag
    0,                      // 294 MorsePackColumn
    0,                      // 295 MorsePackRow
    0,                      // 296 MorsePackLowerBand
    0,                      // 297 MorsePackUpeprBand
    0,                      // 298 MorsePackAll
    0,                      // 299
    0,                      // 300
    0,                      // 301 MorseNoVec
    0,                      // 302 MorseVec
    0,                      // 303 MorseIvec
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 390
    0,                      // 391
    0,                      // 392
    0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 401
    0,                      // 402
    0, 0, 0, 0, 0, 0, 0, 0  // Remember to add a coma!
};
