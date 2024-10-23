/**
 *
 * @file context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Options
 * @brief Group routines exposed to users to handle options
 *
 */

#include <stdlib.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include "control/morsewinthread.h"
#else
#include <pthread.h>
#endif

#include "control/common.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "chameleon/morse_runtime.h"

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas.h"
#endif

/**
 *  Global data
 */
/* master threads context lookup table */
static MORSE_context_t *morse_ctxt = NULL;

/**
 *  Create new context
 */
MORSE_context_t *morse_context_create()
{
    MORSE_context_t *morse;

    if ( morse_ctxt != NULL ) {
        morse_error("morse_context_create", "a context is already existing\n");
        return NULL;
    }

    morse = (MORSE_context_t*)malloc(sizeof(MORSE_context_t));
    if (morse == NULL) {
        morse_error("morse_context_create", "malloc() failed");
        return NULL;
    }

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    morse->nb                 = 128;
    morse->ib                 = 32;
    morse->rhblock            = 4;

    morse->nworkers           = 1;
    morse->ncudas             = 0;
    morse->nthreads_per_worker= 1;

    morse->warnings_enabled     = MORSE_TRUE;
    morse->autotuning_enabled   = MORSE_TRUE;
    morse->parallel_enabled     = MORSE_FALSE;
    morse->profiling_enabled    = MORSE_FALSE;
    morse->progress_enabled     = MORSE_FALSE;

    morse->householder        = MORSE_FLAT_HOUSEHOLDER;
    morse->translation        = MORSE_OUTOFPLACE;


    /* Initialize scheduler */
    RUNTIME_context_create(morse);

    morse_ctxt = morse;
    return morse;
}


/**
 *  Return context for a thread
 */
MORSE_context_t *morse_context_self()
{
    return morse_ctxt;
}

/**
 *  Clean the context
 */
int morse_context_destroy(){

    RUNTIME_context_destroy(morse_ctxt);
    free(morse_ctxt);
    morse_ctxt = NULL;

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  MORSE_Enable - Enable MORSE feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be enabled:
 *          @arg MORSE_WARNINGS   printing of warning messages,
 *          @arg MORSE_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg MORSE_PROFILING_MODE  activate profiling of kernels
 *          @arg MORSE_PROGRESS  activate progress indicator
 *          @arg MORSE_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Enable(MORSE_enum option)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Enable", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case MORSE_WARNINGS:
            morse->warnings_enabled = MORSE_TRUE;
            break;
        case MORSE_AUTOTUNING:
            morse->autotuning_enabled = MORSE_TRUE;
            break;
        case MORSE_PROFILING_MODE:
            morse->profiling_enabled = MORSE_TRUE;
            break;
        case MORSE_PROGRESS:
            morse->progress_enabled = MORSE_TRUE;
            break;
        case MORSE_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(1);
#else
            morse_error("MORSE_Enable", "cannot enable GEMM3M (not available in cblas)");
#endif
            break;
        /* case MORSE_PARALLEL: */
        /*     morse->parallel_enabled = MORSE_TRUE; */
        /*     break; */
        default:
            morse_error("MORSE_Enable", "illegal parameter value");
            return MORSE_ERR_ILLEGAL_VALUE;
        case MORSE_BOUND:
            break;
    }

    /* Enable at the lower level if required */
    RUNTIME_enable( option );

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  MORSE_Disable - Disable MORSE feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be disabled:
 *          @arg MORSE_WARNINGS   printing of warning messages,
 *          @arg MORSE_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg MORSE_PROFILING_MODE  deactivate profiling of kernels
 *          @arg MORSE_PROGRESS  deactivate progress indicator
 *          @arg MORSE_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Disable(MORSE_enum option)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Disable", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case MORSE_WARNINGS:
            morse->warnings_enabled = MORSE_FALSE;
            break;
        case MORSE_AUTOTUNING:
            morse->autotuning_enabled = MORSE_FALSE;
            break;
        case MORSE_PROFILING_MODE:
            morse->profiling_enabled = MORSE_FALSE;
            break;
        case MORSE_PROGRESS:
            morse->progress_enabled = MORSE_FALSE;
            break;
        case MORSE_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(0);
#endif
            break;
        case MORSE_PARALLEL_MODE:
            morse->parallel_enabled = MORSE_FALSE;
            break;
        default:
            morse_error("MORSE_Disable", "illegal parameter value");
            return MORSE_ERR_ILLEGAL_VALUE;
    }

    /* Disable at the lower level if required */
    RUNTIME_disable( option );

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  MORSE_Set - Set MORSE parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg MORSE_TILE_SIZE:        size matrix tile,
 *          @arg MORSE_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Set(MORSE_enum param, int value)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Set", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case MORSE_TILE_SIZE:
            if (value <= 0) {
                morse_error("MORSE_Set", "negative tile size");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            morse->nb = value;
            if ( morse->autotuning_enabled ) {
                morse->autotuning_enabled = MORSE_FALSE;
                morse_warning("MORSE_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            morse->ib = chameleon_min( morse->nb, morse->ib );
            break;
        case MORSE_INNER_BLOCK_SIZE:
            if (value <= 0) {
                morse_error("MORSE_Set", "negative inner block size");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            if (value > morse->nb) {
                morse_error("MORSE_Set", "inner block larger than tile");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            /* if (morse->nb % value != 0) { */
            /*     morse_error("MORSE_Set", "inner block does not divide tile"); */
            /*     return MORSE_ERR_ILLEGAL_VALUE; */
            /* } */
            morse->ib = value;

            if ( morse->autotuning_enabled ) {
                morse->autotuning_enabled = MORSE_FALSE;
                morse_warning("MORSE_Set", "autotuning has been automatically disable\n");
            }
            break;
        case MORSE_HOUSEHOLDER_MODE:
            if (value != MORSE_FLAT_HOUSEHOLDER && value != MORSE_TREE_HOUSEHOLDER) {
                morse_error("MORSE_Set", "illegal value of MORSE_HOUSEHOLDER_MODE");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            morse->householder = value;
            break;
        case MORSE_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                morse_error("MORSE_Set", "negative householder size");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            morse->rhblock = value;
            break;
        case MORSE_TRANSLATION_MODE:
            if (value != MORSE_INPLACE && value != MORSE_OUTOFPLACE) {
                morse_error("MORSE_Set", "illegal value of MORSE_TRANSLATION_MODE");
                return MORSE_ERR_ILLEGAL_VALUE;
            }
            morse->translation = value;
            break;
        default:
            morse_error("MORSE_Set", "unknown parameter");
            return MORSE_ERR_ILLEGAL_VALUE;
    }

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  MORSE_Get - Get value of MORSE parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg MORSE_TILE_SIZE:        size matrix tile,
 *          @arg MORSE_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Get(MORSE_enum param, int *value)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Get", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case MORSE_TILE_SIZE:
            *value = morse->nb;
            return MORSE_SUCCESS;
        case MORSE_INNER_BLOCK_SIZE:
            *value = morse->ib;
            return MORSE_SUCCESS;
        case MORSE_HOUSEHOLDER_MODE:
            *value = morse->householder;
            return MORSE_SUCCESS;
        case MORSE_HOUSEHOLDER_SIZE:
            *value = morse->rhblock;
            return MORSE_SUCCESS;
        case MORSE_TRANSLATION_MODE:
            *value = morse->translation;
            return MORSE_SUCCESS;
        default:
            morse_error("MORSE_Get", "unknown parameter");
            return MORSE_ERR_ILLEGAL_VALUE;
    }

    return MORSE_SUCCESS;
}
