/**
 *
 * @file morse_f77.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 interface
 *
 * @version 1.0.0
 * @author Bilel Hadri
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include "control/common.h"
#include "morse.h"
#include "morse_f77.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**
     *  FORTRAN API - auxiliary function prototypes
     */
    void MORSE_INIT(int *CORES, int *NGPUS, int *INFO)
    {   *INFO = MORSE_Init(*CORES, *NGPUS); }

    void MORSE_FINALIZE(int *INFO)
    {   *INFO = MORSE_Finalize(); }

    void MORSE_ENABLE(MORSE_enum *lever, int *INFO)
    {   *INFO = MORSE_Enable(*lever); }

    void MORSE_DISABLE(MORSE_enum *lever, int *INFO)
    {   *INFO = MORSE_Disable(*lever); }

    void MORSE_SET(MORSE_enum *param, int *value, int *INFO)
    {   *INFO = MORSE_Set(*param, *value); }

    void MORSE_GET(MORSE_enum *param, int *value, int *INFO)
    {   *INFO = MORSE_Get(*param, value); }

    void MORSE_DEALLOC_HANDLE(size_t *sp, int *INFO)
    {   free((void *)(*sp));
        *INFO = MORSE_SUCCESS; }

    void MORSE_VERSION(int *VER_MAJOR, int *VER_MINOR, int *VER_MICRO, int *INFO)
    {
        *VER_MAJOR = CHAMELEON_VERSION_MAJOR;
        *VER_MINOR = CHAMELEON_VERSION_MINOR;
        *VER_MICRO = CHAMELEON_VERSION_MICRO;
        *INFO = MORSE_SUCCESS;
    }

    /**
     *  FORTRAN API - descriptor allocation and deallocation
     */
    void MORSE_DESC_CREATE(MORSE_desc_t **desc, void *mat, MORSE_enum *dtyp,
                           int *mb, int *nb, int *bsiz, int *lm, int *ln,
                           int *i, int *j, int *m, int *n, int *p, int *q,
                           int *INFO)
    {   *INFO = MORSE_Desc_Create(desc, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q); }
    void MORSE_DESC_CREATE_OOC(MORSE_desc_t **desc, MORSE_enum *dtyp,
                               int *mb, int *nb, int *bsiz, int *lm, int *ln,
                               int *i, int *j, int *m, int *n, int *p, int *q,
                               int *INFO)
    {   *INFO = MORSE_Desc_Create_OOC(desc, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q); }
    void MORSE_DESC_CREATE_USER(MORSE_desc_t **descptr, void *mat, MORSE_enum *dtyp,
                                int *mb, int *nb, int *bsiz, int *lm, int *ln,
                                int *i, int *j, int *m, int *n, int *p, int *q,
                                void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                                int   (*get_blkldd) ( const MORSE_desc_t*, int      ),
                                int   (*get_rankof) ( const MORSE_desc_t*, int, int ),
                                int *INFO)
    {   *INFO = MORSE_Desc_Create_User(descptr, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q,
                                       get_blkaddr, get_blkldd, get_rankof); }
    void MORSE_DESC_CREATE_OOC_USER(MORSE_desc_t **descptr, MORSE_enum *dtyp,
                                    int *mb, int *nb, int *bsiz, int *lm, int *ln,
                                    int *i, int *j, int *m, int *n, int *p, int *q,
                                    int (*get_rankof) ( const MORSE_desc_t*, int, int ),
                                    int *INFO)
    {   *INFO = MORSE_Desc_Create_OOC_User(descptr, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q,
                                           get_rankof); }

    void MORSE_DESC_DESTROY(MORSE_desc_t **desc, int *INFO)
    {   *INFO = MORSE_Desc_Destroy(desc); }

    /**
     *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
     */
    void MORSE_LAPACK_TO_TILE(intptr_t *Af77, int *LDA, intptr_t *A, int *INFO)
    {   *INFO = MORSE_Lapack_to_Tile( (void *)Af77, *LDA, (MORSE_desc_t *)(*A)); }

    void MORSE_TILE_TO_LAPACK(intptr_t *A, intptr_t *Af77, int *LDA, int *INFO)
    {   *INFO = MORSE_Tile_to_Lapack((MORSE_desc_t *)(*A), (void *)Af77, *LDA); }

#ifdef __cplusplus
}
#endif
