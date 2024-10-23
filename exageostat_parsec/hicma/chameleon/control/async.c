/**
 *
 * @file async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Sequences
 * @brief Group routines exposed to users to handle asynchronous tasks execution
 *
 */
#include <stdlib.h>
#include "control/common.h"
#include "chameleon/morse_runtime.h"

/**
 *  Register an exception.
 */
int morse_request_fail(MORSE_sequence_t *sequence, MORSE_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return status;
}

/**
 *  Create a sequence
 */
int morse_sequence_create(MORSE_context_t *morse, MORSE_sequence_t **sequence)
{
    if ((*sequence = malloc(sizeof(MORSE_sequence_t))) == NULL) {
        morse_error("MORSE_Sequence_Create", "malloc() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

    RUNTIME_sequence_create( morse, *sequence );

    (*sequence)->status = MORSE_SUCCESS;
    return MORSE_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int morse_sequence_destroy(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    RUNTIME_sequence_destroy( morse, sequence );
    free(sequence);
    return MORSE_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int morse_sequence_wait(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    RUNTIME_sequence_wait( morse, sequence );
    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Sequences
 *
 *  MORSE_Sequence_Create - Create a squence.
 *
 ******************************************************************************
 *
 * @param[out] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Sequence_Create(MORSE_sequence_t **sequence)
{
    MORSE_context_t *morse;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Sequence_Create", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    status = morse_sequence_create(morse, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  MORSE_Sequence_Destroy - Destroy a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Sequence_Destroy(MORSE_sequence_t *sequence)
{
    MORSE_context_t *morse;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Sequence_Destroy", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_Sequence_Destroy", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    status = morse_sequence_destroy(morse, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  MORSE_Sequence_Wait - Wait for the completion of a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Sequence_Wait(MORSE_sequence_t *sequence)
{
    MORSE_context_t *morse;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Sequence_Wait", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_Sequence_Wait", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    status = morse_sequence_wait(morse, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  MORSE_Sequence_Flush - Terminate a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 * @param[in] request
 *          The flush request.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Sequence_Flush(MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Sequence_Flush", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_Sequence_Flush", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }

    RUNTIME_sequence_flush( morse->schedopt, sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    return MORSE_SUCCESS;
}
