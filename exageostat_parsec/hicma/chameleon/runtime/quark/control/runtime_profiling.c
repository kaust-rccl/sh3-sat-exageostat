/**
 *
 * @file runtime_profiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark profiling routines
 *
 * @version 1.0.0
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2011-10-29
 *
 */
#include "chameleon_quark.h"
#include "chameleon/chameleon_timer.h"

double RUNTIME_get_time(){
    return CHAMELEON_timer();
}

void RUNTIME_start_profiling()
{
    morse_warning("RUNTIME_start_profiling()", "FxT profiling is not available with Quark\n");
}

void RUNTIME_stop_profiling()
{
    morse_warning("RUNTIME_stop_profiling()", "FxT profiling is not available with Quark\n");
}

void RUNTIME_start_stats()
{
    morse_warning("RUNTIME_start_stats()", "pruning stats are not available with Quark\n");
}

void RUNTIME_stop_stats()
{
    morse_warning("RUNTIME_stop_stats()", "pruning stats are not available with Quark\n");
}

void RUNTIME_schedprofile_display(void)
{
    morse_warning("RUNTIME_schedprofile_display(quark)", "Scheduler profiling is not available with Quark\n");
}

void RUNTIME_kernelprofile_display(void)
{
    morse_warning("RUNTIME_kernelprofile_display(quark)", "Kernel profiling is not available with Quark\n");
}

/**
 *  Set iteration numbers for traces
 */
void RUNTIME_iteration_push( MORSE_context_t *morse, unsigned long iteration )
{
    (void)morse; (void)iteration;
    return;
}
void RUNTIME_iteration_pop( MORSE_context_t *morse )
{
    (void)morse;
    return;
}
