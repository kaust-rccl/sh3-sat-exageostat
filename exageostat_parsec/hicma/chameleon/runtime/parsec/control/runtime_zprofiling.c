/**
 *
 * @file runtime_zprofiling.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC MORSE_Complex64_t kernel progiling
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2017-01-12
 *
 */
#include "chameleon_parsec.h"

void RUNTIME_zdisplay_allprofile()
{
    morse_warning("RUNTIME_zdisplay_allprofile(PaRSEC)", "Profiling is not available with PaRSEC");
}

void RUNTIME_zdisplay_oneprofile( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zdisplay_oneprofile(PaRSEC)", "Profiling is not available with PaRSEC\n");
}

