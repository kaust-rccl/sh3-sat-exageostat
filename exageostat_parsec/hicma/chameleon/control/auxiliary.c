/**
 *
 * @file auxiliary.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Auxiliary
 * @brief Group auxiliary routines exposed to users
 *
 */

#include "control/common.h"
#include "control/auxiliary.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while MORSE is being used correctly.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void morse_warning(const char *func_name, const char *msg_text)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL)
        morse_fatal_error("morse_warning", "MORSE not initialized");
    if (morse->warnings_enabled)
        fprintf(stderr, "MORSE WARNING: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of MORSE.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void morse_error(const char *func_name, const char *msg_text)
{
    fprintf(stderr, "MORSE ERROR: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *  Context oblivious.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void morse_fatal_error(const char *func_name, const char *msg_text)
{
    fprintf(stderr, "MORSE FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/**
 *  Returns core id
 */
int morse_rank(MORSE_context_t *morse)
{
    return RUNTIME_thread_rank( morse );
}

/**
 *  Tune block size nb and internal block size ib
 */
int morse_tune(MORSE_enum func, int M, int N, int NRHS)
{
    MORSE_context_t *morse;
    morse = morse_context_self();
    if ( morse && morse->autotuning_enabled == MORSE_TRUE ) {
        morse_warning( "morse_tune", "Autotunning not available for now" );
    }
    (void)func;
    (void)M;
    (void)N;
    (void)NRHS;
    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Version - Reports MORSE version number.
 *
 ******************************************************************************
 *
 * @param[out] ver_major
 *          MORSE major version number.
 *
 * @param[out] ver_minor
 *          MORSE minor version number.
 *
 * @param[out] ver_micro
 *          MORSE micro version number.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Version(int *ver_major, int *ver_minor, int *ver_micro)
{
    if (! ver_major && ! ver_minor && ! ver_micro)
        return  MORSE_ERR_ILLEGAL_VALUE;

    if (ver_major)
        *ver_major = CHAMELEON_VERSION_MAJOR;

    if (ver_minor)
        *ver_minor = CHAMELEON_VERSION_MINOR;

    if (ver_micro)
        *ver_micro = CHAMELEON_VERSION_MICRO;

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Element_Size - Reports the size in bytes of a MORSE precision type
 *  (e.g. MorseInteger, MorseRealFloat, etc).
 *
 ******************************************************************************
 *
 * @param[in] type
 *          MORSE element type, can be one of the following:
 *          - MorseByte
 *          - MorseInteger
 *          - MorseRealFloat
 *          - MorseRealDouble
 *          - MorseComplexFloat
 *          - MorseComplexDouble
 *
 ******************************************************************************
 *
 * @return
 *          \retval Element size in bytes
 *
 */
int MORSE_Element_Size(int type)
{
    switch(type) {
        case MorseByte:          return          1;
        case MorseInteger:       return   sizeof(int);
        case MorseRealFloat:     return   sizeof(float);
        case MorseRealDouble:    return   sizeof(double);
        case MorseComplexFloat:  return 2*sizeof(float);
        case MorseComplexDouble: return 2*sizeof(double);
        default: morse_fatal_error("MORSE_Element_Size", "undefined type");
                 return MORSE_ERR_ILLEGAL_VALUE;

    }
}

/**
 *
 * @ingroup Auxiliary
 *
 *  MORSE_My_Mpi_Rank - Return the MPI rank of the calling process.
 *
 ******************************************************************************
 *
 ******************************************************************************
 *
 * @return
 *          \retval MPI rank
 *
 */
int MORSE_My_Mpi_Rank(void)
{
#if defined(CHAMELEON_USE_MPI)
    MORSE_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Finalize()", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    return MORSE_MPI_RANK;
#else
    return MORSE_SUCCESS;
#endif
}

/**
 *  Display a progress percentage in stderr
 */
void update_progress(int currentValue, int maximumValue) {
    div_t res ;
    static int progress = -1; /* varie de 0 a 100 au cours du calcul concerne */

    if (maximumValue == 0) {
        res.quot = 100;
    }
    else {
        if (currentValue < (INT_MAX / 100) ) {
            res = div(currentValue*100, maximumValue);
        }
        else {
            /* Calcule le quotient de la division */
            res.quot = (int)( (long long)( currentValue * 100 ) / maximumValue );
        }
    }

    // Print the percentage
    if (res.quot > progress) {
        fprintf(stderr, "%3d%%\b\b\b\b", res.quot);
    }
    progress = res.quot;

    if (currentValue >= maximumValue) {
        progress = -1;
    }
}

// A function to display the progress indicator.
// By default it is update_progress()
// The user can change it with MORSE_Set_update_progress_callback()
void (*update_progress_callback)(int, int) = update_progress;

int MORSE_Set_update_progress_callback(void (*p)(int, int)) {
  update_progress_callback = p;
  return MORSE_SUCCESS;
}

