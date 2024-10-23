/**
 *
 * @file context.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#ifndef _MORSE_CONTEXT_H_
#define _MORSE_CONTEXT_H_

#include "chameleon/morse_struct.h"

/**
 *  Routines to handle threads context
 */
#ifdef __cplusplus
extern "C" {
#endif

MORSE_context_t* morse_context_create  ();
MORSE_context_t* morse_context_self    ();
int              morse_context_destroy ();

#ifdef __cplusplus
}
#endif

#endif
