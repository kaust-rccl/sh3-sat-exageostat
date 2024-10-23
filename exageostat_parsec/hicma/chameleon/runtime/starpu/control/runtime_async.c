/**
 *
 * @file runtime_async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU asynchronous routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include "chameleon_starpu.h"

/**
 *  Create a sequence
 */
int RUNTIME_sequence_create( MORSE_context_t  *morse,
                             MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;
    return MORSE_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int RUNTIME_sequence_destroy( MORSE_context_t  *morse,
                              MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;
    return MORSE_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int RUNTIME_sequence_wait( MORSE_context_t  *morse,
                           MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;

    if (morse->progress_enabled) {
        RUNTIME_progress(morse);
    }

    starpu_task_wait_for_all();
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif
    return MORSE_SUCCESS;
}

/**
 *  Terminate a sequence
 */
void RUNTIME_sequence_flush( MORSE_context_t  *morse,
                             MORSE_sequence_t *sequence,
                             MORSE_request_t  *request,
                             int status )
{
    (void)morse;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}
