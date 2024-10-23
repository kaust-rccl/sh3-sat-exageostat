/**
 *
 * @file morse_tasks.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon elementary tasks main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 */
#ifndef _MORSE_TASKS_H_
#define _MORSE_TASKS_H_

#include "chameleon/chameleon_config.h"

BEGIN_C_DECLS

/**
 * MORSE Task submission routines
 */
void MORSE_TASK_ztile_zero(const MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_dtile_zero(const MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_ctile_zero(const MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           const MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_stile_zero(const MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           const MORSE_desc_t *A, int Am, int An, int lda);

#include "chameleon/morse_tasks_z.h"
#include "chameleon/morse_tasks_d.h"
#include "chameleon/morse_tasks_c.h"
#include "chameleon/morse_tasks_s.h"
#include "chameleon/morse_tasks_zc.h"
#include "chameleon/morse_tasks_ds.h"

END_C_DECLS

#endif /* _MORSE_TASKS_H_ */
