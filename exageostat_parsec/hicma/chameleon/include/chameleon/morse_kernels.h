/**
 *
 * @file morse_kernels.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon elementary kernels enum
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 */
#ifndef _MORSE_KERNELS_H_
#define _MORSE_KERNELS_H_

/**
 * Used to apply operations on specific kernels
 */
typedef enum morse_kernel_e {

  MORSE_GEMM,
  MORSE_HEMM,
  MORSE_HER2K,
  MORSE_HERK,
  MORSE_SYTRF_NOPIV,
  MORSE_SYMM,
  MORSE_SYR2K,
  MORSE_SYRK,
  MORSE_TRSM,
  MORSE_TRMM,

  MORSE_GELQT,
  MORSE_GEQRT,
  MORSE_GESSM,
  MORSE_GETRF,
  MORSE_GETRF_INCPIV,
  MORSE_GETRF_NOPIV,
  MORSE_LAUUM,
  MORSE_ORMLQ,
  MORSE_ORMQR,
  MORSE_POTRF,
  MORSE_SSSSM,
  MORSE_TRTRI,
  MORSE_TSLQT,
  MORSE_TSMLQ,
  MORSE_TSMQR,
  MORSE_TSQRT,
  MORSE_TSTRF,
  MORSE_TTLQT,
  MORSE_TTMLQ,
  MORSE_TTMQR,
  MORSE_TTQRT,
  MORSE_UNMLQ,
  MORSE_UNMQR,

  MORSE_GEADD,
  MORSE_LASCAL,
  MORSE_LACPY,
  MORSE_LAG2C,
  MORSE_LAG2Z,
  MORSE_LANGE,
  MORSE_LANHE,
  MORSE_LANSY,
  MORSE_LASET,
  MORSE_LASET2,
  MORSE_PEMV,
  MORSE_PLGHE,
  MORSE_PLGSY,
  MORSE_PLRNT,
  MORSE_TILE_ZERO,

  MORSE_NBKERNELS
} MORSE_kernel_t;

#endif
