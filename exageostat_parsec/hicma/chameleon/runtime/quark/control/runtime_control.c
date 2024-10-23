/**
 *
 * @file runtime_control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark control routines
 *
 * @version 1.0.0
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_quark.h"

/**
 *
 */
int RUNTIME_init( MORSE_context_t *morse,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    int hres = 0;
    if ( ncudas > 0 ) {
        morse_warning( "RUNTIME_init_scheduler(quark)", "GPUs are not supported for now");
    }

    if ( nthreads_per_worker > 0 ) {
        morse_warning( "RUNTIME_init_scheduler(quark)", "Multi-threaded kernels are not supported for now");
    }

    morse->schedopt = (void*)QUARK_New( ncpus );

    return hres;
}

/**
 *
 */
void RUNTIME_finalize( MORSE_context_t *morse )
{
    QUARK_Delete((Quark*)(morse->schedopt));
    return;
}

/**
 *  To suspend the processing of new tasks by workers
 */
void RUNTIME_pause( MORSE_context_t *morse )
{
    (void)morse;
    return;
}

/**
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 */
void RUNTIME_resume( MORSE_context_t *morse )
{
    (void)morse;
    return;
}

/**
 *  Busy-waiting barrier
 */
void RUNTIME_barrier( MORSE_context_t *morse )
{
    QUARK_Barrier((Quark*)(morse->schedopt));
}

/**
 *  Display a progress information when executing the tasks
 */
void RUNTIME_progress( MORSE_context_t *morse )
{
    (void)morse;
    return;
}

/**
 * Thread rank.
 */
int RUNTIME_thread_rank( MORSE_context_t *morse )
{
    return QUARK_Thread_Rank((Quark*)(morse->schedopt));
}

/**
 * Number of threads.
 */
int RUNTIME_thread_size( MORSE_context_t *morse )
{
    (void)morse;
    /*
     * TODO: should add a function to Quark to get the number of thread from the
     * data structure and not from the system function
     */
    return quark_get_numthreads();
}

/**
 *  The process rank
 */
int RUNTIME_comm_rank( MORSE_context_t *morse )
{
    (void)morse;
    return 0;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( MORSE_context_t *morse )
{
    (void)morse;
    return 1;
}
