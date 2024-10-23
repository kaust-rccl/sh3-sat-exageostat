/**
 *
 * @file morse_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon mixed precision wrappers header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#ifndef _MORSE_ZC_H_
#define _MORSE_ZC_H_

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
//int MORSE_zcgesv(int N, int NRHS, MORSE_Complex64_t *A, int LDA, int * IPIV, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t *X, int LDX, int *ITER);
//int MORSE_zcposv(MORSE_enum uplo, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t *X, int LDX, int *ITER);
//int MORSE_zcgels(MORSE_enum trans, int M, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t *X, int LDX, int *ITER);
//int MORSE_zcungesv(MORSE_enum trans, int N, int NRHS, MORSE_Complex64_t *A, int LDA, MORSE_Complex64_t *B, int LDB, MORSE_Complex64_t *X, int LDX, int *ITER);

/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
//int MORSE_zcgesv_Tile(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER);
//int MORSE_zcposv_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER);
/* int MORSE_zcgels_Tile(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER); */
//int MORSE_zcungesv_Tile(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER);

/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
//int MORSE_zcgesv_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zcposv_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zcgels_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_zcungesv_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_desc_t *X, int *ITER, MORSE_sequence_t *sequence, MORSE_request_t *request);

END_C_DECLS

#endif
