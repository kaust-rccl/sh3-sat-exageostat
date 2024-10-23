/**
 *
 * @file tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon layout conversion wrappers
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 ***
 *
 * @defgroup Tile
 * @brief Group routines exposed to users for matrices conversion LAPACK-Tile
 *
 */
#include "control/common.h"
#include "control/auxiliary.h"

/**
 *
 * @ingroup Tile
 *
 *  MORSE_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 ******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Lapack_to_Tile(void *Af77, int LDA, MORSE_desc_t *A)
{
    switch( A->dtyp ) {
    case MorseComplexDouble:
        return MORSE_zLapack_to_Tile( (MORSE_Complex64_t *)Af77, LDA, A );
        break;
    case MorseComplexFloat:
        return MORSE_cLapack_to_Tile( (MORSE_Complex32_t *)Af77, LDA, A );
        break;
    case MorseRealFloat:
        return MORSE_sLapack_to_Tile( (float *)Af77, LDA, A );
        break;
    case MorseRealDouble:
    default:
        return MORSE_dLapack_to_Tile( (double *)Af77, LDA, A );
    }
    return MORSE_ERR_ILLEGAL_VALUE;
}

/**
 *
 * @ingroup Tile
 *
 *  MORSE_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 ******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 * @param[in] Af77
 *          LAPACK matrix (only needed on proc 0).
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Tile_to_Lapack(MORSE_desc_t *A, void *Af77, int LDA)
{
    switch( A->dtyp ) {
    case MorseComplexDouble:
        return MORSE_zTile_to_Lapack( A, (MORSE_Complex64_t *)Af77, LDA );
        break;
    case MorseComplexFloat:
        return MORSE_cTile_to_Lapack( A, (MORSE_Complex32_t *)Af77, LDA );
        break;
    case MorseRealFloat:
        return MORSE_sTile_to_Lapack( A, (float *)Af77, LDA );
        break;
    case MorseRealDouble:
    default:
        return MORSE_dTile_to_Lapack( A, (double *)Af77, LDA );
    }
    return MORSE_ERR_ILLEGAL_VALUE;
}
