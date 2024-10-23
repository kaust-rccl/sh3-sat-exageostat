/**
 *
 * @file runtime_codelet_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU MORSE_Complex64_t codelets header
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 * @precisions normal z -> c d s
 *
 */
#ifndef _CODELETS_Z_H_
#define _CODELETS_Z_H_

#include <stdio.h>
#include "runtime_codelets.h"

#include "chameleon/morse_tasks_z.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_z.h"
#if defined(CHAMELEON_USE_CUDA)
#include "cudablas.h"
#endif
#endif

/*
 * Management functions
 */
ZCODELETS_HEADER(tile_zero)

/*
 * BLAS 1 functions
 */
ZCODELETS_HEADER(axpy)

/*
 * BLAS 3 functions
 */
ZCODELETS_HEADER(gemm)
ZCODELETS_HEADER(hemm)
ZCODELETS_HEADER(her2k)
ZCODELETS_HEADER(herk)
ZCODELETS_HEADER(symm)
ZCODELETS_HEADER(syr2k)
ZCODELETS_HEADER(syrk)
ZCODELETS_HEADER(trmm)
ZCODELETS_HEADER(trsm)

/*
 * LAPACK functions
 */
ZCODELETS_HEADER(gelqt)
ZCODELETS_HEADER(geqrt)
ZCODELETS_HEADER(gessm)
ZCODELETS_HEADER(gessq)
ZCODELETS_HEADER(getrf)
ZCODELETS_HEADER(getrf_incpiv)
ZCODELETS_HEADER(getrf_nopiv)
ZCODELETS_HEADER(herfb)
ZCODELETS_HEADER(lauum)
ZCODELETS_HEADER(potrf)
ZCODELETS_HEADER(ssssm)
ZCODELETS_HEADER(syssq)
ZCODELETS_HEADER(trasm)
ZCODELETS_HEADER(trssq)
ZCODELETS_HEADER(trtri)
ZCODELETS_HEADER(tplqt)
ZCODELETS_HEADER(tpqrt)
ZCODELETS_HEADER(tpmlqt)
ZCODELETS_HEADER(tpmqrt)
ZCODELETS_HEADER(tslqt)
ZCODELETS_HEADER(tsmlq)
ZCODELETS_HEADER(tsmqr)
ZCODELETS_HEADER(tsmlq_hetra1)
ZCODELETS_HEADER(tsmqr_hetra1)
ZCODELETS_HEADER(tsqrt)
ZCODELETS_HEADER(tstrf)
ZCODELETS_HEADER(ttlqt)
ZCODELETS_HEADER(ttmlq)
ZCODELETS_HEADER(ttmqr)
ZCODELETS_HEADER(ttqrt)
ZCODELETS_HEADER(unmlq)
ZCODELETS_HEADER(unmqr)

/*
 * Auxiliary functions
 */
ZCODELETS_HEADER(geadd)
ZCODELETS_HEADER(he2ge)
ZCODELETS_HEADER(lascal)
ZCODELETS_HEADER(tradd)
ZCODELETS_HEADER(lacpy)
ZCODELETS_HEADER(lange)
ZCODELETS_HEADER(lange_max)
ZCODELETS_HEADER(lansy)
ZCODELETS_HEADER(lantr)
ZCODELETS_HEADER(laset)
ZCODELETS_HEADER(laset2)
ZCODELETS_HEADER(latro)
ZCODELETS_HEADER(plssq)
ZCODELETS_HEADER(plssq2)

/*
 * MIXED PRECISION functions
 */
ZCODELETS_HEADER(lag2c)

/*
 * DZ functions
 */
ZCODELETS_HEADER(asum)

/*
 * CPU only functions
 */
ZCODELETS_HEADER(plrnt)
ZCODELETS_HEADER(build)

#if defined(PRECISION_z) || defined(PRECISION_c)
ZCODELETS_HEADER(hessq)
ZCODELETS_HEADER(lanhe)
ZCODELETS_HEADER(plghe)
ZCODELETS_HEADER(sytrf_nopiv)
#endif
ZCODELETS_HEADER(plgsy)

#endif /* _CODELETS_Z_H_ */
