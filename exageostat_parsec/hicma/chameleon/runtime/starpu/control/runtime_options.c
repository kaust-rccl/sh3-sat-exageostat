/**
 *
 * @file runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU options routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_starpu.h"

void RUNTIME_options_init( MORSE_option_t *option, MORSE_context_t *morse,
                           MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    option->sequence   = sequence;
    option->request    = request;
    option->profiling  = MORSE_PROFILING == MORSE_TRUE;
    option->parallel   = MORSE_PARALLEL == MORSE_TRUE;
    option->priority   = MORSE_PRIORITY_MIN;
    option->nb         = MORSE_NB;
    option->ws_wsize   = 0;
    option->ws_hsize   = 0;
    option->ws_worker  = NULL;
    option->ws_host    = NULL;
    return;
}

void RUNTIME_options_finalize( MORSE_option_t *option, MORSE_context_t *morse )
{
    (void)option;
    (void)morse;
    return;
}

int RUNTIME_options_ws_alloc( MORSE_option_t *options, size_t worker_size, size_t host_size )
{
    int ret = 0;
    if ( worker_size > 0 ) {
        options->ws_wsize = worker_size;
        starpu_vector_data_register((starpu_data_handle_t*)(&(options->ws_worker)),
                                    -1, (uintptr_t)NULL,
                                    worker_size, sizeof(char));
    }
    if ( host_size > 0 ) {
        options->ws_hsize = host_size;
        ret = RUNTIME_starpu_ws_alloc((MORSE_starpu_ws_t**)&(options->ws_host),
                                      host_size, MORSE_CUDA, MORSE_HOST_MEM);
    }
    return ret;
}

int RUNTIME_options_ws_free( MORSE_option_t *options )
{
    int ret = 0;
    if ( options->ws_worker != NULL ) {
        starpu_data_unregister_submit((starpu_data_handle_t)(options->ws_worker));
        options->ws_worker = NULL;
    }
    if ( options->ws_host != NULL ) {
        starpu_task_wait_for_all();
        ret = RUNTIME_starpu_ws_free( (MORSE_starpu_ws_t*)(options->ws_host) );
        options->ws_host = NULL;
    }
    return ret;
}
