/**
 *
 * @file morse.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2012-09-15
 *
 */
#ifndef _MORSE_H_
#define _MORSE_H_

/* ****************************************************************************
 * MORSE types and constants
 */
#include "chameleon/chameleon_config.h"
#include "chameleon/morse_types.h"
#include "chameleon/morse_struct.h"
#include "chameleon/morse_constants.h"
#include "chameleon/morse_kernels.h"


/* ****************************************************************************
 * MORSE runtime common API
 */
#include "chameleon/morse_runtime.h"


/* ****************************************************************************
 * MORSE Simulation mode
 */
#include "chameleon/morse_simulate.h"

/* ****************************************************************************
 * Include LibHQR for hierarchical trees QR/LQ factorizations
 */
#include "libhqr.h"

/* ****************************************************************************
 * MORSE Tasks
 */
#include "chameleon/morse_tasks.h"

/* ****************************************************************************
 * MORSE functionnalities
 */
#include "chameleon/morse_z.h"
#include "chameleon/morse_c.h"
#include "chameleon/morse_d.h"
#include "chameleon/morse_s.h"
#include "chameleon/morse_zc.h"
#include "chameleon/morse_ds.h"

/* ****************************************************************************
 * MORSE Functions
 */
BEGIN_C_DECLS

/* Auxiliary */
int MORSE_Version           (int *ver_major, int *ver_minor, int *ver_micro);
int MORSE_My_Mpi_Rank       (void);
int MORSE_Init              (int nworkers, int ncudas);
int MORSE_InitPar           (int nworkers, int ncudas, int nthreads_per_worker);
int MORSE_Finalize          (void);
int MORSE_Pause             (void);
int MORSE_Resume            (void);
int MORSE_Distributed_start (void);
int MORSE_Distributed_stop  (void);
int MORSE_Comm_size         (void);
int MORSE_Comm_rank         (void);
int MORSE_Lapack_to_Tile    (void *Af77, int LDA, MORSE_desc_t *A);
int MORSE_Tile_to_Lapack    (MORSE_desc_t *A, void *Af77, int LDA);
int MORSE_Distributed_start (void);
int MORSE_Distributed_stop  (void);
int MORSE_Distributed_size  (int *size);
int MORSE_Distributed_rank  (int *rank);
int MORSE_GetThreadNbr      (void);

/* Descriptor */
int MORSE_Element_Size(int type);
int MORSE_Desc_Create  (MORSE_desc_t **desc, void *mat, MORSE_enum dtyp,
                        int mb, int nb, int bsiz, int lm, int ln,
                        int i, int j, int m, int n, int p, int q);
int MORSE_Desc_Create_User(MORSE_desc_t **desc, void *mat, MORSE_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                           int (*get_blkldd)( const MORSE_desc_t*, int ),
                           int (*get_rankof)( const MORSE_desc_t*, int, int ));
int MORSE_Desc_Create_OOC(MORSE_desc_t **desc, MORSE_enum dtyp,
                          int mb, int nb, int bsiz, int lm, int ln,
                          int i, int j, int m, int n, int p, int q);
int MORSE_Desc_Create_OOC_User(MORSE_desc_t **desc, MORSE_enum dtyp,
                               int mb, int nb, int bsiz, int lm, int ln,
                               int i, int j, int m, int n, int p, int q,
                               int (*get_rankof)( const MORSE_desc_t*, int, int ));
int MORSE_Desc_Destroy (MORSE_desc_t **desc);
int MORSE_Desc_Acquire (MORSE_desc_t  *desc);
int MORSE_Desc_Release (MORSE_desc_t  *desc);
int MORSE_Desc_Flush   (MORSE_desc_t  *desc, MORSE_sequence_t *sequence);
void MORSE_user_tag_size(int, int) ;

/* Workspaces */
int MORSE_Dealloc_Workspace (MORSE_desc_t **desc);

/* Options */
int MORSE_Enable  (MORSE_enum option);
int MORSE_Disable (MORSE_enum option);
int MORSE_Set     (MORSE_enum param, int  value);
int MORSE_Get     (MORSE_enum param, int *value);
int MORSE_Set_update_progress_callback(void (*p)(int, int)) ;

/* Sequences */
int MORSE_Sequence_Create  (MORSE_sequence_t **sequence);
int MORSE_Sequence_Destroy (MORSE_sequence_t *sequence);
int MORSE_Sequence_Wait    (MORSE_sequence_t *sequence);

END_C_DECLS

#endif
