/**
 *
 * @file morse_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon structures
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _MORSE_STRUCT_H_
#define _MORSE_STRUCT_H_

#include "chameleon/chameleon_config.h"
#include "chameleon/morse_types.h"
#include "chameleon/morse_kernels.h"

BEGIN_C_DECLS

/**
 * RUNTIME headers to include types of :
 *         - QUARK
 *         - PaRSEC
 *         - StarPU
 */
typedef enum morse_sched_e {
  RUNTIME_SCHED_QUARK,
  RUNTIME_SCHED_PARSEC,
  RUNTIME_SCHED_STARPU,
} MORSE_sched_t;


/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
struct morse_desc_s;
typedef struct morse_desc_s MORSE_desc_t;

struct morse_desc_s {
    // function to get matrix tiles address
    void *(*get_blkaddr)( const MORSE_desc_t*, int, int );
    // function to get matrix tiles leading dimension
    int   (*get_blkldd )( const MORSE_desc_t*, int );
    // function to get matrix tiles MPI rank
    int   (*get_rankof) ( const MORSE_desc_t*, int, int );
    void *mat;        // pointer to the beginning of the matrix
    size_t A21;       // pointer to the beginning of the matrix A21
    size_t A12;       // pointer to the beginning of the matrix A12
    size_t A22;       // pointer to the beginning of the matrix A22
    MORSE_enum styp;  // storage layout of the matrix
    MORSE_enum dtyp;  // precision of the matrix
    int mb;           // number of rows in a tile
    int nb;           // number of columns in a tile
    int bsiz;         // size in elements including padding
    int lm;  	      // number of rows of the entire matrix
    int ln;           // number of columns of the entire matrix
    int lmt;          // number of tile rows of the entire matrix - derived parameter
    int lnt;          // number of tile columns of the entire matrix - derived parameter
    int i;            // row index to the beginning of the submatrix
    int j;            // column index to the beginning of the submatrix
    int m;            // number of rows of the submatrix
    int n;            // number of columns of the submatrix
    int mt;           // number of tile rows of the submatrix - derived parameter
    int nt;           // number of tile columns of the submatrix - derived parameter
                      // Data for distributed cases
    int p;            // number of rows of the 2D distribution grid
    int q;            // number of columns of the 2D distribution grid
    int llm;          // number of rows of the 2D distribution grid
    int lln;          // number of columns of the 2D distribution grid
    int llm1;         // number of tile rows of the A11 matrix - derived parameter
    int lln1;         // number of tile columns of the A11 matrix - derived parameter
    int llmt;         // number of tile rows of the local (to a node) matrix
    int llnt;         // number of tile columns of the local (to a node) matrix
    int id;           // identification number of the descriptor
    int occurences;   // identify main matrix desc (occurances=1) or
                      // submatrix desc (occurances>1) to avoid unregistering
                      // GPU data twice
    int use_mat;      // 1 if we have a pointer to the overall data mat - else 0
    int alloc_mat;    // 1 if we handle the allocation of mat - else 0
    int register_mat; // 1 if we have to register mat - else 0 (handled by the application)
    int myrank;       // MPI rank of the descriptor
    int ooc;          // 1 if the matrix is not to fit in memory
    void *schedopt;   // scheduler (QUARK|StarPU) specific structure
};


/**
 *  MORSE request uniquely identifies each asynchronous function call.
 */
typedef struct morse_context_s {
    MORSE_sched_t      scheduler;
    int                nworkers;
    int                ncudas;
    int                nthreads_per_worker;
#if defined(CHAMELEON_USE_MPI)
    int                my_mpi_rank;
    int                mpi_comm_size;
#endif
    int                world_size;
    int                group_size;

    /* Boolean flags */
    MORSE_bool         warnings_enabled;
    MORSE_bool         autotuning_enabled;
    MORSE_bool         parallel_enabled;
    MORSE_bool         profiling_enabled;
    MORSE_bool         progress_enabled;

    MORSE_enum         householder;        // "domino" (flat) or tree-based (reduction) Householder
    MORSE_enum         translation;        // In place or Out of place layout conversion

    int                nb;
    int                ib;
    int                nbnbsize;           // tile size in elements (possibly padded)
    int                ibnbsize;           // T or L tile size in elements (---''---)
    int                rhblock;            // block size for tree-based (reduction) Householder
    void              *schedopt;           // structure for runtimes
    int                mpi_outer_init;     // MPI has been initialized outside our functions
} MORSE_context_t;


/**
 *  MORSE request uniquely identifies each asynchronous function call.
 */
typedef struct morse_request_s {
    MORSE_enum status; // MORSE_SUCCESS or appropriate error code
} MORSE_request_t;


/**
 *  MORSE sequence uniquely identifies a set of asynchronous function calls
 *  sharing common exception handling.
 */
typedef struct morse_sequence_s {
    MORSE_bool       status;    /* MORSE_SUCCESS or appropriate error code */
    MORSE_request_t *request;   /* failed request                          */
    void            *schedopt;
} MORSE_sequence_t;


/**
 *  MORSE options
 */
typedef struct morse_option_s {
    MORSE_sequence_t *sequence;
    MORSE_request_t  *request;
    int               profiling;
    int               parallel;
    int               priority;
    int               nb;
    size_t            ws_wsize;
    size_t            ws_hsize;
    void             *ws_worker;  /*> Workspace located on the worker        */
    void             *ws_host;    /*> Workspace *always* located on the host */
    void             *schedopt;
} MORSE_option_t;

END_C_DECLS

#endif /* __CHAMELEON_H__ */
