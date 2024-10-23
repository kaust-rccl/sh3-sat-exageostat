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
#ifndef _LIBHQR_INTERNAL_H_
#define _LIBHQR_INTERNAL_H_

#include "libhqr.h"
#include <stdio.h>
#include <assert.h>

#define PRINT_PIVGEN 0
#ifdef PRINT_PIVGEN
#define myassert( test ) {if ( ! (test) ) return -1;}
#else
#define myassert(test) {assert((test)); return -1;}
#endif

typedef enum qrtree_type_ {
    LIBHQR_QRTREE_UNSET = 0,
    LIBHQR_QRTREE_HQR,
    LIBHQR_QRTREE_SVD,
    LIBHQR_QRTREE_SYS,
    LIBHQR_QRTREE_MTX,
} qrtree_type_e;

struct hqr_args_s;
typedef struct hqr_args_s hqr_args_t;

struct hqr_subpiv_s;
typedef struct hqr_subpiv_s hqr_subpiv_t;

/**
 * @brief jhj
 */
struct hqr_args_s {
    int domino;  /**< Switch to enable/disable the domino tree linking high and lw level reduction trees */
    int tsrr;    /**< Switch to enable/disable round-robin on TS to optimise pipelining between TS and local tree */
    hqr_subpiv_t *llvl;
    hqr_subpiv_t *hlvl;
    int *perm;
};

struct hqr_subpiv_s {
    /**
     * currpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] m   line you want to eliminate
     *
     *  @return the annihilator p used with m at step k
     */
    int (*currpiv)(const hqr_subpiv_t *arg, int k, int m);
    /*
     * nextpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] p   line currently used as an annihilator
     *    @param[in] m   line actually annihilated.
     *          m = MT to find the first time p is used as an annihilator during step k
     *
     *  @return the next line m' using the line p as annihilator during step k
     *          mt if p will never be used again as an annihilator.
     */
    int (*nextpiv)(const hqr_subpiv_t *arg, int k, int p, int m);
    /*
     * prevpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] p   line currently used as an annihilator
     *    @param[in] m   line actually annihilated.
     *          m = p to find the last time p has been used as an annihilator during step k
     *
     *  @return the previous line m' using the line p as annihilator during step k
     *          mt if p has never been used before that as an annihilator.
     */
    int (*prevpiv)(const hqr_subpiv_t *arg, int k, int p, int m);
    int *ipiv;
    int minMN;
    int ldd;
    int a;
    int p;
    int domino;
};

/**
 * @brief Minimal structure to store the information related to each tile.
 */
typedef struct libhqr_tile_info_s {
    int type;          /**< The type of the reduction applied to the tile (@sa libhqr_type_e) */
    int currpiv;       /**< The row index of the pivot for this tile                          */
    int nextpiv;       /**< The next tile for which the currpiv is a pivot                    */
    int prevpiv;       /**< The previous tile for which currpiv was a pivot                   */
    int first_nextpiv; /**< If type != 0, the first tile for which this tile is a pivot       */
    int first_prevpiv; /**< If type != 0, the last tile for which this tile is a pivot        */
} libhqr_tile_info_t;

static inline int libhqr_imin(int a, int b){
    return (a > b) ? b : a;
}

static inline int libhqr_imax(int a, int b){
    return (a > b) ? a : b;
}

static inline int libhqr_iceil(int a, int b){
    return (a + b - 1) / b;
}


/* Number of extra tile of type 1 between the tile of type 3 and the first of nb11 */
#define nbextra1_formula (( (k % pa) > (pa - p) ) ? (-k)%pa + pa : 0)

/*
 * Common functions
 */
/* int hqr_getnbgeqrf( const libhqr_tree_t *qrtree, int k ); */
/* int hqr_getm(       const libhqr_tree_t *qrtree, int k, int i   ); */
/* int hqr_geti(       const libhqr_tree_t *qrtree, int k, int m   ); */
/* int hqr_gettype(    const libhqr_tree_t *qrtree, int k, int m   ); */

/*
 * Subtree for low-level
 */
void libhqr_matrix_init(libhqr_tree_t *qrtree, const libhqr_tree_t *qrtree_init);
int  rdmtx_gettype(const libhqr_tree_t *qrtree, int k, int m);
int  rdmtx_currpiv(const libhqr_tree_t *qrtree, int k, int m);
int  rdmtx_nextpiv(const libhqr_tree_t *qrtree, int k, int p, int m);
int  rdmtx_prevpiv(const libhqr_tree_t *qrtree, int k, int p, int m);
void libhqr_matrix_finalize(libhqr_tree_t *qrtree);

/*
 * function for drawing the tree
 */
void draw_rectangle(int k, int p, int m, int step_m, FILE *file);
void draw_lines(const libhqr_tree_t *qrtree, int k, int *tiles, int *step, FILE *file);
void draw_horizontal_line(int k, int p, int m, int step_p, int step_m, FILE *file);
void draw_vertical_line(  int k, int p, int m,             int step_m, FILE *file);

/**
 * @name Low level trees
 * @{
 */
void hqr_low_flat_init     ( hqr_subpiv_t *arg );
void hqr_low_binary_init   ( hqr_subpiv_t *arg );
void hqr_low_fibonacci_init( hqr_subpiv_t *arg, int minMN );
void hqr_low_greedy_init   ( hqr_subpiv_t *arg, int minMN );
void hqr_low_greedy1p_init ( hqr_subpiv_t *arg, int minMN );
void svd_low_adaptiv_init  ( hqr_subpiv_t *arg, int gmt, int gnt, int nbcores, int ratio );

/**
 * @}
 *
 * @name High level trees
 * @{
 */
void hqr_high_flat_init     ( hqr_subpiv_t *arg );
void hqr_high_binary_init   ( hqr_subpiv_t *arg );
void hqr_high_fibonacci_init( hqr_subpiv_t *arg );
void hqr_high_greedy_init   ( hqr_subpiv_t *arg, int minMN );
void hqr_high_greedy1p_init ( hqr_subpiv_t *arg );
/**
 * @}
 */

void libhqr_fct_to_mtx( const libhqr_tree_t *in, libhqr_tree_t *out );

#endif /* _LIBHQR_INTERNAL_H_ */
