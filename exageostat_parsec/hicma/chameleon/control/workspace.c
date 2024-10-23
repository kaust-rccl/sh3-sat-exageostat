/**
 *
 * @file workspace.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon workspace routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Workspace
 * @brief Group routines exposed to users about specific workspaces management
 *
 */

#include <stdlib.h>
#include "control/common.h"
#include "control/auxiliary.h"
#include "control/workspace.h"

/**
 *
 */
int morse_alloc_ibnb_tile(int M, int N, MORSE_enum func, int type, MORSE_desc_t **desc, int p, int q)
{
    int status;
    int IB, NB, MT, NT;
    int64_t lm, ln;
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("morse_alloc_ibnb_tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = morse_tune(func, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("morse_alloc_ibnb_tile", "morse_tune() failed");
        return MORSE_ERR_UNEXPECTED;
    }

    /* Set MT & NT & allocate */
    NB = MORSE_NB;
    IB = MORSE_IB;
    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Size is doubled for RH QR to store the reduction T */
    if ((morse->householder == MORSE_TREE_HOUSEHOLDER) &&
        ((func == MORSE_FUNC_SGELS)  ||
         (func == MORSE_FUNC_DGELS)  ||
         (func == MORSE_FUNC_CGELS)  ||
         (func == MORSE_FUNC_ZGELS)  ||
         (func == MORSE_FUNC_SGESVD) ||
         (func == MORSE_FUNC_DGESVD) ||
         (func == MORSE_FUNC_CGESVD) ||
         (func == MORSE_FUNC_ZGESVD)))
        NT *= 2;

    lm = IB * MT;
    ln = NB * NT;

    /* Allocate and initialize descriptor */
    *desc = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    if (*desc == NULL) {
        morse_error("morse_alloc_ibnb_tile", "malloc() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }
    **desc = morse_desc_init(type, IB, NB, IB*NB, lm, ln, 0, 0, lm, ln, p, q);

    /* Allocate matrix */
    if (morse_desc_mat_alloc(*desc)) {
        morse_error("morse_alloc_ibnb_tile", "malloc() failed");
        free(*desc);
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

    RUNTIME_desc_create( *desc );

    /* Check that everything is ok */
    status = morse_desc_check(*desc);
    if (status != MORSE_SUCCESS) {
        morse_error("morse_alloc_ibnb_tile", "invalid descriptor");
        free(*desc);
        return status;
    }

    return MORSE_SUCCESS;
}

/**
 *
 */
int morse_alloc_ipiv(int M, int N, MORSE_enum func, int type, MORSE_desc_t **desc, void **IPIV, int p, int q)
{
    int status;
    int NB, IB, MT, NT;
    int64_t lm, ln;
    size_t size;
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("morse_alloc_ipiv", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = morse_tune(func, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("morse_alloc_ipiv", "morse_tune() failed");
        return MORSE_ERR_UNEXPECTED;
    }

    /* Set MT & NT & allocate */
    NB = MORSE_NB;
    IB = MORSE_IB;

    NT = (N%NB==0) ? (N/NB) : ((N/NB)+1);
    MT = (M%NB==0) ? (M/NB) : ((M/NB)+1);

    lm = IB * MT;
    ln = NB * NT;

    size = (size_t)(chameleon_min(MT, NT) * NB * NT * sizeof(int));
    if (size == 0) {
        *IPIV = NULL;
        return MORSE_SUCCESS;
    }
    /* TODO: Fix the distribution for IPIV */
    *IPIV = (int*)malloc( size );

    *desc = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    **desc = morse_desc_init(type, IB, NB, IB*NB, lm, ln, 0, 0, lm, ln, p, q );

    if ( morse_desc_mat_alloc(*desc) ) {
        morse_error("morse_alloc_ipiv", "malloc() failed");
        free(*desc);
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

    RUNTIME_desc_create( *desc );

    return MORSE_SUCCESS;
}

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Dealloc_Worksapce - Deallocate workspace descriptor allocated by
 *                            any workspace allocation routine.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Workspace descriptor
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Dealloc_Workspace(MORSE_desc_t **desc)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Dealloc_Workspace", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        morse_error("MORSE_Dealloc_Workspace", "attempting to deallocate a NULL descriptor");
        return MORSE_ERR_UNALLOCATED;
    }
    if ((*desc)->mat == NULL && (*desc)->use_mat == 1) {
        morse_error("MORSE_Dealloc_Worspace", "attempting to deallocate a NULL pointer");
        return MORSE_ERR_UNALLOCATED;
    }
    morse_desc_mat_free( *desc );
    RUNTIME_desc_destroy( *desc );

    free(*desc);
    *desc = NULL;
    return MORSE_SUCCESS;
}
