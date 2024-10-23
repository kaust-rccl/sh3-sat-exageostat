/**
 *
 * @file workspace_z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon precision dependent workspace routines
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "control/common.h"
#include "control/workspace.h"

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgeev - Allocates workspace for MORSE_zgeev or
 *  MORSE_zgeev_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile Hessenberg.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgeev(int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(N, N, MORSE_FUNC_ZGEEV, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgehrd - Allocates workspace for MORSE_zgehrd or
 *  MORSE_zgehrd_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile Hessenberg.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgehrd(int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(N, N, MORSE_FUNC_ZGEHRD, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgebrd - Allocates workspace for MORSE_zgebrd or MORSE_zgebrd_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile BRD.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgebrd(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGEBRD, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgels - Allocates workspace for MORSE_zgels or
 *  MORSE_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile QR or the tile LQ factorization.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgels(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgels_Tile - Allocates tile workspace for
 *  MORSE_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, MORSE_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgels_Tile(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgeqrf - Allocates workspace for MORSE_zgeqrf or
 *  MORSE_zgeqrf_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile QR factorization.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgeqrf(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgeqrf_Tile - Allocates tile workspace for
 *  MORSE_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, MORSE_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgeqrf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgelqf - Allocates workspace for MORSE_zgelqf or
 *  MORSE_zgelqf_Tile routines.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile LQ
 *          factorization.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgelqf(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgelqf_Tile - Allocates tile workspace for MORSE_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, MORSE_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgelqf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGELS, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgesv - Allocates workspace for MORSE_zgesv or
 *  MORSE_zgesv_Tile routines.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A.
 *          N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors
 *          required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required
 *          by the tile LU factorization (not equivalent to LAPACK).
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgesv_incpiv(int N, MORSE_desc_t **descL, int **IPIV, int p, int q) {
    return morse_alloc_ipiv(N, N, MORSE_FUNC_ZGESV, MorseComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgesv_Tile - Allocates workspace for MORSE_zgesv_Tile
 *  routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A.
 *          N >= 0.
 *
 * @param[out] descL
 *          On exit, MORSE descriptor on workspace handle for storage of the
 *          extra L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by
 *          the tile LU factorization (not equivalent to LAPACK).
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgesv_incpiv_Tile(int N, MORSE_desc_t **descL, int **IPIV, int p, int q)
{
    return morse_alloc_ipiv(N, N, MORSE_FUNC_ZGESV, MorseComplexDouble, descL, (void**)IPIV, p, q);
}
/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgesvd - Allocates workspace for MORSE_zgesvd or
 *  MORSE_zgesvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgesvd(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZGESVD, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgetrf_incpiv - Allocates workspace for
 *  MORSE_zgetrf_incpiv or MORSE_zgetrf_incpiv_Tile or
 *  MORSE_zgetrf_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************
 *
 * @sa MORSE_zgetrf_incpiv
 * @sa MORSE_zgetrf_incpiv_Tile
 * @sa MORSE_zgetrf_incpiv_Tile_Async
 *
 */
int MORSE_Alloc_Workspace_zgetrf_incpiv(int M, int N, MORSE_desc_t **descL, int **IPIV, int p, int q) {
    return morse_alloc_ipiv(M, N, MORSE_FUNC_ZGESV, MorseComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zgetrf_incpiv_Tile - Allocates workspace for
 *  MORSE_zgesv_incpiv_Tile or MORSE_zgesv_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, MORSE descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zgetrf_incpiv_Tile(int N, MORSE_desc_t **descL, int **IPIV, int p, int q) {
    return morse_alloc_ipiv(N, N, MORSE_FUNC_ZGESV, MorseComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zheev - Allocates workspace for MORSE_zheev or MORSE_zheev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zheev(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZHEEV, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zheevd - Allocates workspace for MORSE_zheevd or MORSE_zheevd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zheevd(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZHEEVD, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zhegv - Allocates workspace for MORSE_zhegv or MORSE_zhegv_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zhegv(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZHEGV, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zhegvd - Allocates workspace for MORSE_zhegvd or MORSE_zhegvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zhegvd(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZHEGVD, MorseComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  MORSE_Alloc_Workspace_zhetrd - Allocates workspace for MORSE_zhetrd or MORSE_zhetrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 */
int MORSE_Alloc_Workspace_zhetrd(int M, int N, MORSE_desc_t **descT, int p, int q) {
    return morse_alloc_ibnb_tile(M, N, MORSE_FUNC_ZHETRD, MorseComplexDouble, descT, p, q); }
