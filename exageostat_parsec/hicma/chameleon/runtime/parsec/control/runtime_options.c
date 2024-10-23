/**
 *
 * @file runtime_options.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC options routines
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2017-01-12
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_parsec.h"

void RUNTIME_options_init( MORSE_option_t *options, MORSE_context_t *morse,
                           MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    options->sequence   = sequence;
    options->request    = request;
    options->profiling  = MORSE_PROFILING == MORSE_TRUE;
    options->parallel   = MORSE_PARALLEL == MORSE_TRUE;
    options->priority   = MORSE_PRIORITY_MIN;
    options->nb         = MORSE_NB;
    options->ws_wsize   = 0;
    options->ws_hsize   = 0;
    options->ws_worker  = NULL;
    options->ws_host    = NULL;
    return;
}

void RUNTIME_options_finalize( MORSE_option_t *options, MORSE_context_t *morse )
{
    (void)options;
    (void)morse;
    return;
}

int RUNTIME_options_ws_alloc( MORSE_option_t *options, size_t worker_size, size_t host_size )
{
    options->ws_wsize = worker_size;
    options->ws_hsize = host_size;
    return MORSE_SUCCESS;
}

int RUNTIME_options_ws_free( MORSE_option_t *options )
{
    options->ws_wsize = 0;
    options->ws_hsize = 0;
    return MORSE_SUCCESS;
}
