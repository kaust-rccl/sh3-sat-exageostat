/**
 *
 * @file morse_mf77.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 interface for mixed-precision computational routines
 *
 * @version 1.0.0
 * @author Bilel Hadri
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include "control/common.h"

#define MORSE_ZCGESV           MORSE_FNAME(zcgesv,   ZCGESV)
#define MORSE_DSGESV           MORSE_FNAME(dsgesv,   DSGESV)
#define MORSE_ZCPOSV           MORSE_FNAME(zcposv,   ZCPOSV)
#define MORSE_DSPOSV           MORSE_FNAME(dsposv,   DSPOSV)
#define MORSE_ZCGELS           MORSE_FNAME(zcgels,   ZCGELS)
#define MORSE_DSGELS           MORSE_FNAME(dsgels,   DSGELS)
#define MORSE_ZCUNGESV         MORSE_FNAME(zcungesv, ZCUNGESV)
#define MORSE_DSUNGESV         MORSE_FNAME(dsungesv, DSUNGESV)

#define MORSE_ZCGESV_TILE       MORSE_TILE_FNAME(zcgesv,   ZCGESV)
#define MORSE_DSGESV_TILE       MORSE_TILE_FNAME(dsgesv,   DSGESV)
#define MORSE_ZCPOSV_TILE       MORSE_TILE_FNAME(zcposv,   ZCPOSV)
#define MORSE_DSPOSV_TILE       MORSE_TILE_FNAME(dsposv,   DSPOSV)
#define MORSE_ZCGELS_TILE       MORSE_TILE_FNAME(zcgels,   ZCGELS)
#define MORSE_DSGELS_TILE       MORSE_TILE_FNAME(dsgels,   DSGELS)
#define MORSE_ZCUNGESV_TILE     MORSE_TILE_FNAME(zcungesv, ZCUNGESV)
#define MORSE_DSUNGESV_TILE     MORSE_TILE_FNAME(dsungesv, DSUNGESV)

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  FORTRAN API - math functions (simple interface)
 */
//void MORSE_ZCGESV(int *N, int *NRHS, MORSE_Complex64_t *A, int *LDA, int *IPIV, MORSE_Complex64_t *B, int *LDB, MORSE_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_zcgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_DSGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_dsgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_ZCPOSV(MORSE_enum *uplo, int *N, int *NRHS, MORSE_Complex64_t *A, int *LDA, MORSE_Complex64_t *B, int *LDB, MORSE_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_zcposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_DSPOSV(MORSE_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_dsposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_ZCGELS(MORSE_enum *trans, int *M, int *N, int *NRHS, MORSE_Complex64_t *A, int *LDA, MORSE_Complex64_t **T, MORSE_Complex64_t *B, int *LDB, MORSE_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_zcgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_DSGELS(MORSE_enum *trans, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_dsgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_ZCUNGESV(MORSE_enum *trans, int *N, int *NRHS, MORSE_Complex64_t *A, int *LDA, MORSE_Complex64_t *B, int *LDB, MORSE_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_zcungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void MORSE_DSUNGESV(MORSE_enum *trans, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = MORSE_dsungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

#ifdef __cplusplus
}
#endif

/**
 *  FORTRAN API - math functions (native interface)
 */
//void MORSE_ZCGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_zcgesv_Tile((MORSE_desc_t *)(*A), IPIV, (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_DSGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_zcgesv_Tile((MORSE_desc_t *)(*A), IPIV, (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_ZCPOSV_TILE(MORSE_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_zcposv_Tile(*uplo, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_DSPOSV_TILE(MORSE_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_dsposv_Tile(*uplo, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_ZCGELS_TILE(MORSE_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_zcgels_Tile(*trans, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*T), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_DSGELS_TILE(MORSE_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_dsgels_Tile(*trans, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*T), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_ZCUNGESV_TILE(MORSE_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_zcungesv_Tile(*trans, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*T), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
//
//void MORSE_DSUNGESV_TILE(MORSE_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = MORSE_dsungesv_Tile(*trans, (MORSE_desc_t *)(*A), (MORSE_desc_t *)(*T), (MORSE_desc_t *)(*B), (MORSE_desc_t *)(*X), ITER); }
