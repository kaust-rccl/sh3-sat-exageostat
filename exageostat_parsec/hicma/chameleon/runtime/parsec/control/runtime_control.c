/**
 *
 * @file runtime_control.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC control routines
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

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 * Initialize MORSE
 */
int RUNTIME_init( MORSE_context_t *morse,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    int hres = -1, default_ncores = -1;
    int *argc = (int *)malloc(sizeof(int));
    *argc = 0;

    /* Initializing parsec context */
    if( 0 < ncpus ) {
        default_ncores = ncpus;
    }
    morse->parallel_enabled = MORSE_TRUE;
    morse->schedopt = (void *)parsec_init(default_ncores, argc, NULL);

    if(NULL != morse->schedopt) {
        morse->nworkers = ncpus;
        morse->nthreads_per_worker = nthreads_per_worker;
        hres = 0;
    }

    free(argc);

    (void)ncudas;
    return hres;
}

/**
 * Finalize MORSE
 */
void RUNTIME_finalize( MORSE_context_t *morse )
{
    parsec_context_t *parsec = (parsec_context_t*)morse->schedopt;
    parsec_fini(&parsec);
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
 * Barrier MORSE.
 */
void RUNTIME_barrier( MORSE_context_t *morse )
{
    parsec_context_t *parsec = (parsec_context_t*)(morse->schedopt);
    // This will be a problem with the fake tasks inserted to detect end of DTD algorithms
    parsec_context_wait( parsec );
    return;
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
    (void)morse;
    return 0;
}

/**
 * Thread rank.
 */
int RUNTIME_thread_size( MORSE_context_t *morse )
{
    // TODO: fixme
    //return vpmap_get_nb_total_threads();
    (void)morse;
    return 1;
}

/**
 *  This returns the rank of this process
 */
int RUNTIME_comm_rank( MORSE_context_t *morse )
{
    int rank = 0;
#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    (void)morse;
    return rank;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( MORSE_context_t *morse )
{
    int size = 0;
#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    (void)morse;
    return size;
}
