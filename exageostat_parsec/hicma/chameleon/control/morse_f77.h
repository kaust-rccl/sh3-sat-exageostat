/**
 *
 * @file morse_f77.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 naming macros
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2017-05-03
 *
 */
#ifndef _MORSE_F77_H_
#define _MORSE_F77_H_

#include "chameleon/morse_mangling.h"

/**
 *  Determine FORTRAN names
 */
#define MORSE_FNAME(lcname, UCNAME) MORSE_GLOBAL(morse_##lcname, MORSE_##UCNAME)
#define MORSE_TILE_FNAME(lcname, UCNAME) MORSE_GLOBAL(morse_##lcname##_tile, MORSE_##UCNAME##_TILE)
#define MORSE_ASYNC_FNAME(lcname, UCNAME) MORSE_GLOBAL(morse_##lcname##_tile_async, MORSE_##UCNAME##_TILE_ASYNC)
#define MORSE_WS_FNAME(lcname, UCNAME) MORSE_GLOBAL(morse_alloc_workspace_##lcname, MORSE_ALLOC_WORKSPACE_##UCNAME)
#define MORSE_WST_FNAME(lcname, UCNAME) MORSE_GLOBAL(morse_alloc_workspace_##lcname##_tile, MORSE_ALLOC_WORKSPACE_##UCNAME##_TILE)

#define MORSE_INIT MORSE_GLOBAL(morse_init, MORSE_INIT)
#define MORSE_FINALIZE MORSE_GLOBAL(morse_finalize, MORSE_FINALIZE)
#define MORSE_ENABLE MORSE_GLOBAL(morse_enable, MORSE_ENABLE)
#define MORSE_DISABLE MORSE_GLOBAL(morse_disable, MORSE_DISABLE)
#define MORSE_SET MORSE_GLOBAL(morse_set, MORSE_SET)
#define MORSE_GET MORSE_GLOBAL(morse_get, MORSE_GET)
#define MORSE_DEALLOC_HANDLE MORSE_GLOBAL(morse_dealloc_handle, MORSE_DEALLOC_HANDLE)
#define MORSE_VERSION MORSE_GLOBAL(morse_version, MORSE_VERSION)
#define MORSE_DESC_CREATE MORSE_GLOBAL(morse_desc_create, MORSE_DESC_CREATE)
#define MORSE_DESC_CREATE_OOC MORSE_GLOBAL(morse_desc_create_ooc, MORSE_DESC_CREATE_OOC)
#define MORSE_DESC_CREATE_USER MORSE_GLOBAL(morse_desc_create_user, MORSE_DESC_CREATE_USER)
#define MORSE_DESC_CREATE_OOC_USER MORSE_GLOBAL(morse_desc_create_ooc_user, MORSE_DESC_CREATE_OOC_USER)
#define MORSE_DESC_DESTROY MORSE_GLOBAL(morse_desc_destroy, MORSE_DESC_DESTROY)
#define MORSE_LAPACK_TO_TILE MORSE_GLOBAL(morse_lapack_to_tile, MORSE_LAPACK_TO_TILE)
#define MORSE_TILE_TO_LAPACK MORSE_GLOBAL(morse_tile_to_lapack, MORSE_TILE_TO_LAPACK)

#endif
