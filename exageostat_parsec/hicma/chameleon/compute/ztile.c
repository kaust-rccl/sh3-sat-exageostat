/**
 *
 * @file ztile.c
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
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zTile_to_Lapack
 * @sa MORSE_cLapack_to_Tile
 * @sa MORSE_dLapack_to_Tile
 * @sa MORSE_sLapack_to_Tile
 *
 */
int MORSE_zLapack_to_Tile( MORSE_Complex64_t *Af77, int LDA, MORSE_desc_t *A )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request;
    MORSE_desc_t *B;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zLapack_to_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zLapack_to_Tile", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    MORSE_Desc_Create_User( &B, Af77, MorseComplexDouble, A->mb, A->nb, A->bsiz,
                            LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                            morse_getaddr_cm, morse_getblkldd_cm, NULL );

    /* Start the computation */
    morse_sequence_create( morse, &sequence );

    morse_pzlacpy( MorseUpperLower, B, A, sequence, &request );

    MORSE_Desc_Flush( B, sequence );
    MORSE_Desc_Flush( A, sequence );

    morse_sequence_wait( morse, sequence );

    /* Destroy temporary B descriptor */
    MORSE_Desc_Destroy( &B );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          Af77 has to be A->mat, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zLapack_to_Tile
 * @sa MORSE_cTile_to_Lapack
 * @sa MORSE_dTile_to_Lapack
 * @sa MORSE_sTile_to_Lapack
 *
 */
int MORSE_zTile_to_Lapack( MORSE_desc_t *A, MORSE_Complex64_t *Af77, int LDA )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request;
    MORSE_desc_t *B;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zTile_to_Lapack", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zTile_to_Lapack", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    MORSE_Desc_Create_User( &B, Af77, MorseComplexDouble, A->mb, A->nb, A->bsiz,
                            LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                            morse_getaddr_cm, morse_getblkldd_cm, NULL );

    /* Start the computation */
    morse_sequence_create( morse, &sequence );

    morse_pzlacpy( MorseUpperLower, A, B, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( B, sequence );

    morse_sequence_wait( morse, sequence );

    MORSE_Desc_Destroy( &B );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}
