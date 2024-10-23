/**
 *
 * @file runtime_context.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC context routines
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2017-01-12
 *
 */
#include <stdlib.h>
#include "chameleon_parsec.h"

/**
 *  Create new context
 */
void RUNTIME_context_create( MORSE_context_t *morse )
{
    /* In case of PaRSEC, this is done in init */
    morse->scheduler = RUNTIME_SCHED_PARSEC;
    return;
}

/**
 *  Clean the context
 */
void RUNTIME_context_destroy( MORSE_context_t *morse )
{
    (void)morse;
    return;
}

/**
 *
 */
void RUNTIME_enable(MORSE_enum lever)
{
    switch (lever)
    {
    case MORSE_PROFILING_MODE:
        break;
    default:
        return;
    }
    return;
}

/**
 *
 */
void RUNTIME_disable(MORSE_enum lever)
{
    switch (lever)
    {
    case MORSE_PROFILING_MODE:
        break;
    default:
        return;
    }
    return;
}
