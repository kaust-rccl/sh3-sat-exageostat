/**
 *
 * @file chameleon_quark.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark runtime main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _MORSE_QUARK_H_
#define _MORSE_QUARK_H_

#include <quark.h>
#include "coreblas.h"
#include "core_blas_dag.h"

#include "control/common.h"

typedef struct quark_option_s {
    Quark_Task_Flags flags;
    Quark *quark;
} quark_option_t;

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (type*)RUNTIME_data_getaddr( desc, m, n ) )

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _MORSE_QUARK_H_ */
