/**
 *
 * @file runtime_zprofiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark MORSE_Complex64_t kernel progiling
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @date 2011-06-01
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_quark.h"

void RUNTIME_zdisplay_allprofile()
{
    morse_warning("RUNTIME_zdisplay_allprofile(quark)", "Profiling is not available with Quark");
}

void RUNTIME_zdisplay_oneprofile( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zdisplay_oneprofile(quark)", "Profiling is not available with Quark\n");
}

