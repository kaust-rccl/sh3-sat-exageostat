/**
 *
 * @file libhqr.h
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 *
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#ifndef _LIBHQR_H_
#define _LIBHQR_H_

#ifdef BEGIN_C_DECLS
#undef BEGIN_C_DECLS
#endif

#ifdef END_C_DECLS
#undef END_C_DECLS
#endif

#if defined(c_plusplus) || defined(__cplusplus)
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS          /* empty */
# define END_C_DECLS            /* empty */
#endif

BEGIN_C_DECLS

/**
 * @brief Define which tree level the tile belongs to.
 *
 * Four levels of trees are available:
 *   - The top one (LIBHQR_KILLED_BY_DISTTREE) is meant to be the communication
 *     tree, and should be of size the number of nodes per panel
 *   - The middle one (LIBHQR_KILLED_BY_LOCALTREE) is meant to create parallelism
 *     within a single node between local tiles
 *   - LIBHQR_KILLED_BY_TS is the lowest level tree that is used to keep a high
 *     efficiency kernels within each node.
 *   - LIBHQR_KILLED_BY_DOMINO is the connection tree between the distributed and
 *     the local trees in order to generate a better pipeline for tall and skinny
 *     matrices.
 *
 * @remark LIBHQR_KILLED_BY_TS needs to be 0 for all variant of QR
 * factorizations in order to easily distinguish TT kernels from TS kernels in
 * computations.
 *
 */
typedef enum libhqr_type_ {
    LIBHQR_KILLED_BY_TS        = 0,
    LIBHQR_KILLED_BY_LOCALTREE = 1,
    LIBHQR_KILLED_BY_DOMINO    = 2,
    LIBHQR_KILLED_BY_DISTTREE  = 3,
} libhqr_type_e;

/**
 * @brief Define the type of trees that can be used for the reduction.
 *
 * These are the kind of trees available for low-level reduction in shared
 * memory, or for high-level reduction in distributed memory. Note that all
 * kinds of trees are not available for both levels.
 */
typedef enum libhqr_tree_ {
    LIBHQR_FLAT_TREE      = 0,
    LIBHQR_GREEDY_TREE    = 1,
    LIBHQR_FIBONACCI_TREE = 2,
    LIBHQR_BINARY_TREE    = 3,
    LIBHQR_GREEDY1P_TREE  = 4,
} libhqr_tree_e;

/**
 * @brief Define the type of factorization to apply: QR or LQ.
 */
typedef enum libhqr_facto_ {
    LIBHQR_QR = 0,     /**< QR factorization will be performed, and A shape is considered   */
    LIBHQR_LQ = 1,     /**< LQ factorization will be performed, and A^t shape is considered */
} libhqr_facto_e;

/**
 * @brief Minimal matrix description stucture.
 *
 * This structure describes the dimension of the matrix to factorize in number
 * of tiles: mt-by-nt. When provided, nodes and p provides the total number of
 * nodes involved in the computation, and the P parameter of 2D bloc-cyclic
 * distributions of the tiles. LibHQR doesn't support any other distribtions for
 * now.
 */
typedef struct libhqr_matrix_s {
    int mt;            /**< The number of rows of tiles in the matrix               */
    int nt;            /**< The number of columns of tiles in the matrix            */
    int nodes;         /**< The number of nodes involved in the data distribution   */
    int p;             /**< The number of nodes per column in the data distribution */
} libhqr_matrix_t;

/**
 * @brief Hierarchical tree structure for QR/LQ factorization like kernels.
 *
 * This structure holds the functions that allows one to traverse the tree, or
 * to know form each node which ones are its neigbours in the reduction trees.
 *
 */
struct libhqr_tree_s;
typedef struct libhqr_tree_s libhqr_tree_t;

struct libhqr_tree_s {
    /**
     * @brief Return the number of geqrt/gelqt kernels in the column/row k
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @return The number of geqrt applied to the panel k
     */
    int (*getnbgeqrf)( const libhqr_tree_t *arg, int k );

    /**
     * @brief Compute the row/column index of the i^th geqrt/gelqt kernel in the
     * column/row k
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] i    Index of the geqrt applied on the panel k
     * @return The row index of the i-th geqrt applied on the panel k
     */
    int (*getm)( const libhqr_tree_t *arg, int k, int i );

    /**
     * @brief Compute the index of the row/column m in the list of geqrt/gelqt
     * kernels in the column/row k
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] m    Row index where a geqrt is applied on panel k
     * @return the index in the list of geqrt applied to panel k
     */
    int (*geti)( const libhqr_tree_t *qrtree, int k, int m );
    /**
     * @brief Return the tree level the tile belongs to.
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] m    Row index of the one we request the type
     * @return The type of kernel used to kill the row m at step k:
     *         - 0  if it is a TS kernel
     *         - >0 otherwise. (TT kernel)
     */
    int (*gettype)( const libhqr_tree_t *qrtree, int k, int m );
    /**
     * @brief Return the row/column index of the pivot for m at the step k.
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] m    line you want to eliminate
     * @return The index of the row annihilating the row m at step k
     */
    int (*currpiv)( const libhqr_tree_t *qrtree, int k, int m );
    /**
     * @brief Return the next row/column index for which p is a pivot at the step k.
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] p    line currently used as an annihilator to kill the row m at step k
     * @param[in] m    line actually annihilated by the row p at step k. (k < m <= MT)
     *                 m = MT to find the first time p is used as an annihilator during step k
     *
     * @return the next line that the row p will kill during step k
     *         desc->mt if p will never be used again as an annihilator.
     */
    int (*nextpiv)(const libhqr_tree_t *qrtree, int k, int p, int m);
    /**
     * @brief Return the previous row/column index for which p is a pivot at the step k.
     * @param[in] arg  arguments specific to the reduction tree used
     * @param[in] k    Factorization step
     * @param[in] p    line currently used as an annihilator to kill the row m at step k
     * @param[in] m    line actually annihilated by the row p at step k. (k < m <= MT)
     *                 m = p to find the last time p has been used as an annihilator during step k
     *
     * @return the previous line killed by the row p during step k
     *         desc->mt if p has never been used before as an annihilator.
     */
    int (*prevpiv)(const libhqr_tree_t *qrtree, int k, int p, int m);

    int init;   /**< Internal variable                                               */
    int mt;     /**< Number of rows of tile A if QR factorization, and in A^t if LQ
                     factorization                                                   */
    int nt;     /**< Number of tile reflectors to compute. This is equal to
                     min(A->mt, A->nt)                                               */
    int a;      /**< Size of the local flat TS trees                                 */
    int p;      /**< Size of the highest level tree (distributed one), recommended
                     to be the number of process rows in the 2D-cyclic distribution  */
    int domino; /**< Enable the domino tree to connect high and low level trees
                     in T&S matrices                                                 */
    void *args; /**< Arguments of the qrtree structure (depends on the kind of
                     tree: hqr, systolic, svd, ...)                                  */
};

/**
 * @name Initialization/Finalization functions
 * @{
 */
int  libhqr_init_sys( libhqr_tree_t *qrtree,
                      libhqr_facto_e trans, libhqr_matrix_t *A,
                      int p, int q );

int  libhqr_init_svd( libhqr_tree_t *qrtree,
                      libhqr_facto_e trans, libhqr_matrix_t *A,
                      int type_hlvl, int p, int nbcores_per_node, int ratio );

int  libhqr_init_hqr( libhqr_tree_t *qrtree,
                      libhqr_facto_e trans, libhqr_matrix_t *A,
                      int type_llvl, int type_hlvl,
                      int a, int p, int domino, int tsrr );
void libhqr_finalize( libhqr_tree_t *qrtree );

/**
 * @}
 * @name Function to walk the tree in correct order for STF programing model
 * @{
 */
void libhqr_walk_stepk( const libhqr_tree_t *qrtree, int k, int *tiles );

/**
 * @}
 * @name Drawing functions
 * @{
 */
void libhqr_print_dag( const libhqr_tree_t *qrtree,
                       const char          *filename );
void libhqr_print_svg( const libhqr_tree_t *qrtree,
                       const char          *filename );

/**
 * @}
 * @name Checking/Debugging functions
 * @{
 */
int  libhqr_check        ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree );
void libhqr_print_type   ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree );
void libhqr_print_pivot  ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree );
void libhqr_print_nbgeqrt( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree );
void libhqr_print_perm   ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int *perm );
void libhqr_print_next_k ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k );
void libhqr_print_prev_k ( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k );
void libhqr_print_geqrt_k( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k );

/**
 * @}
 */

END_C_DECLS

#endif /* _LIBHQR_H_ */
