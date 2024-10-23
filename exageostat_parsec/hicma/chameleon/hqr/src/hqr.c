/**
 *
 * @file hqr.c
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 * This file contains all the function to describe the dependencies
 * used in the Xgeqrf_param.jdf file.
 * The QR factorization done with this file relies on three levels:
 *     - the first one is using a flat tree with TS kernels. The
 *       height of this tree is defined by the parameter 'a'. If 'a'
 *       is set to A->mt, the factorization is identical to the one
 *       perform by PLASMA_zgeqrf.
 *       For all subdiagonal "macro-tiles", the line reduced is always
 *       the first.  For all diagonal "macro-tiles", the factorization
 *       performed is identical to the one performed by PLASMA_zgeqrf.
 *
 *     - the third level is using a reduction tree of size 'p'. By
 *       default, the parameter 'p' should be equal to the number of
 *       processors used for the computation, but can be set
 *       differently. (see further example). The type of tree used at
 *       this level is defined by the hlvl parameter. It can be flat
 *       or greedy.
 *       CODE DETAILS: This tree and all the function related to it
 *       are performing a QR factorization on a band matrix with 'p'
 *       the size of the band. All the functions take global indices
 *       as input and return global indices as output.
 *
 *     - Finally, a second 'low' level of reduction tree is applied.
 *       The size of this tree is induced by the parameters 'a' and 'p'
 *       from the first and third levels and is A->mt / ( p * a ). This
 *       tree is reproduced p times for each subset of tiles
 *       S_k = {i in [0, A->mt-1] \ i%p*a = k } with k in [0, p-1].
 *       The tree used for the reduction is defined by the llvl
 *       parameter and can be: flat, greedy, fibonacci or binary.
 *       CODE DETAILS: For commodity, the size of this tree is always
 *       ceil(A->mt / (p * a) ) inducing some extra tests in the code.
 *       All the functions related to this level of tree take as input
 *       the local indices in the A->mt / (p*a) matrix and the global
 *       k. They return the local index. The reductions are so
 *       performed on a trapezoidal matrices where the step is defined
 *       by a:
 *                                    <- min( lhlvl_mt, min( mt, nt ) ) ->
 *                                     __a__   a     a
 *                                    |     |_____
 *                                    |           |_____
 *                                    |                 |_____
 *        llvl_mt = ceil(MT/ (a*p))   |                       |_____
 *                                    |                             |_____
 *                                    |___________________________________|
 *
 *
 *
 *   At each step of the factorization, the lines of tiles are divided
 *   in 4 types:
 *     - QRPARAM_TILE_TS: They are the lines annihilated by a TS
 *     kernel, these lines are never used as an annihilator.  They are
 *     the lines i, with 1 < (i/p)%a < a and i > (k+1)*p
 *     - QRPARAM_TILE_LOCALTT: They are the lines used as annhilitor
 *     in the TS kernels annihiling the QRPARAM_TILE_TS lines.  They
 *     are themselves annihilated by the TT kernel of the low level
 *     reduction tree.  The roots of the local trees are the lines i,
 *     with i/p = k.
 *     - QRPARAM_TILE_DOMINO: These are the lines that are
 *     annhilihated with a domino effect in the band defined by (i/p)
 *     <= k and i >= k
 *     - QRPARAM_TILE_DISTTT: These are the lines annihilated by the
 *     high level tree to reduce communications.
 *     These lines are defined by (i-k)/p = 0.
 */
#include "libhqr_internal.h"
#include "libhqr_queue.h"
#include <stdlib.h>

/*
 * Common functions
 */
static int  hqr_getnbgeqrf( const libhqr_tree_t *qrtree, int k );
static int  hqr_getm(       const libhqr_tree_t *qrtree, int k, int i );
static int  hqr_geti(       const libhqr_tree_t *qrtree, int k, int m );
static int  hqr_gettype(    const libhqr_tree_t *qrtree, int k, int m );

/* Permutation */
static void hqr_genperm   (       libhqr_tree_t *qrtree );
static int  hqr_getinvperm( const libhqr_tree_t *qrtree, int k, int m );

/****************************************************
 *             Common ipiv
 ***************************************************
 *
 * Common prameters to the 4 following functions:
 *    a - Parameter a for the tunable QR
 *    p - Parameter p for the tunable QR
 *    k - Step k of the QR factorization
 *
 */

/*
 * Extra parameter:
 *    gmt - Global number of tiles in a column of the complete distributed matrix
 * Return:
 *    The number of geqrt to execute in the panel k
 */
static int
hqr_getnbgeqrf( const libhqr_tree_t *qrtree, int k )
{
    int a = qrtree->a;
    int p = qrtree->p;
    int gmt = qrtree->mt;
    int domino = qrtree->domino;
    int pa = p * a;
    int nb_1, nb_2, nb_3;
    int nb_11, nb_12;

    /* Number of tasks of type 3 */
    nb_3 = p;

    /* Number of tasks of type 2 */
    if ( domino ) {
      nb_2 = k * (p - 1);

      /* First multiple of p*a under the diagonal of step k */
      nb_11 = ( ( (p * (k+1)) + pa - 1 ) / pa ) * pa;
    }
    else {
      /* Number of extra tile of type 1 between the tile of type 3 and the first of nb11 */
      nb_2 = nbextra1_formula;

      /* First multiple of p*a under the diagonal of step 1 */
      nb_11 = ( (k + p + pa - 1 ) / pa ) * pa;
    }

    /* Last multiple of p*a lower than A->mt */
    nb_12 = ( gmt / pa ) * pa;

    /* Number of tasks of type 1 between nb_11 and nb_12 */
    nb_1 = (nb_12 - nb_11) / a;

    /* Add leftover */
    nb_1 += libhqr_imin( p, gmt - nb_12 );

    return libhqr_imin( nb_1 + nb_2 + nb_3, gmt - k);
}

/*
 * Extra parameter:
 *    i - indice of the geqrt in the continuous space
 * Return:
 *    The global indice m of the i th geqrt in the panel k
 */
static int
hqr_getm( const libhqr_tree_t *qrtree, int k, int i )
{
    int a = qrtree->a;
    int p = qrtree->p;
    int domino = qrtree->domino;

    int pos1, j, pa = p * a;
    int nbextra1 = nbextra1_formula;
    int nb23 = p + (domino ? k*(p-1) : nbextra1 );

    /* Tile of type 2 or 3 or the 1 between the diagonal and the multiple after the diagonal */
    if ( i < nb23 )
        return k+i;
    /* Tile of type 1 */
    else {
        j = i - nb23;
        pa = p * a;
        if ( domino )
          pos1 = ( ( (p * (k+1)) + pa - 1 ) / pa ) * pa;
        else
          pos1 = ( ( (p + k    ) + pa - 1 ) / pa ) * pa;
        return pos1 + (j/p) * pa + j%p;
    }
}

/*
 * Extra parameter:
 *    m - The global indice m of a geqrt in the panel k
 * Return:
 *    The index i of the geqrt in the panel k
 */
static int
hqr_geti( const libhqr_tree_t *qrtree, int k, int m )
{
    int a = qrtree->a;
    int p = qrtree->p;
    int domino = qrtree->domino;
    int lm = m;

    int pos1, j, pa = p * a;
    int nbextra1 = nbextra1_formula;
    int nb23 = p + ( domino ? k*(p-1) : nbextra1 );
    int end2 = p + ( domino ? k*p     : k + nbextra1 );

    /* Tile of type 2 or 3 or the 1 between the diagonal and the multiple after the diagonal */
    if ( lm < end2 )
        return lm-k;
    /* Tile of type 1 */
    else {
        if ( domino )
          pos1 = ( ( (p * (k+1)) + pa - 1 ) / pa ) * pa;
        else
          pos1 = ( ( (p + k    ) + pa - 1 ) / pa ) * pa;
        j = lm - pos1;
        return nb23 + (j / pa) * p + j%pa;
    }
}

/*
 * Extra parameter:
 *      m - The global indice m of the row in the panel k
 * Return
 *     -1 - Error
 *      0 - if m is reduced thanks to a TS kernel
 *      1 - if m is reduced thanks to the low level tree
 *      2 - if m is reduced thanks to the bubble tree
 *      3 - if m is reduced in distributed
 */
static int
hqr_gettype( const libhqr_tree_t *qrtree, int k, int m )
{
    int a = qrtree->a;
    int p = qrtree->p;
    int domino = qrtree->domino;

    int lm = m;
    myassert( lm >= k );

    /* Element to be reduce in distributed */
    if (lm < k + p) {
        return 3;
    }
    /* Element on the local diagonal */
    else if ( domino && lm < p * (k+1) )
      return 2;
    /* Lower triangle of the matrix */
    else {
        if( (lm / p) % a == 0 )
            return 1;
        else
            return 0;
    }
}

/****************************************************
 *
 *   Generic functions currpiv,prevpiv,nextpiv
 *
 ***************************************************/
static int
hqr_currpiv(const libhqr_tree_t *qrtree, int k, int m)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, tmpk;
    int lm, rank;
    int a = qrtree->a;
    int p = qrtree->p;
    int domino = qrtree->domino;
    int gmt = qrtree->mt;

    lm   = m / p; /* Local index in the distribution over p domains */
    rank = m % p; /* Staring index in this distribution             */

    /* TS level common to every case */
    if ( domino ) {
        switch( hqr_gettype( qrtree, k, m ) )
        {
        case 0:
            tmp = lm / a;
            if ( tmp == k / a )
                return        k * p + rank ; /* Below to the first bloc including the diagonal */
            else
                return  tmp * a * p + rank ;
            break;
        case 1:
            tmp = arg->llvl->currpiv(arg->llvl, k, m);
            return  ( tmp == k / a ) ? k * p + rank : tmp * a * p + rank;
            break;
        case 2:
            return m - p;
            break;
        case 3:
            if ( arg->hlvl != NULL )
                return arg->hlvl->currpiv(arg->hlvl, k, m);
        default:
            return gmt;
        }
    }
    else {
        switch( hqr_gettype( qrtree, k, m ) )
        {
        case 0:
            tmp = lm / a;
            /* tmpk = (k + p - 1 - m%p) / p / a;  */
            tmpk = k / (p * a);
            return  ( tmp == tmpk ) ? k + (m-k)%p : tmp * a * p + rank ;
            break;
        case 1:
            tmp = arg->llvl->currpiv(arg->llvl, k, m);
            /* tmpk = (k + p - 1 - m%p) / p / a; */
            tmpk = k / (p * a);
            return  ( tmp == tmpk ) ? k + (m-k)%p : tmp * a * p + rank ;
            break;
        case 2:
            return  m - p;
            break;
        case 3:
            if ( arg->hlvl != NULL )
                return arg->hlvl->currpiv(arg->hlvl, k, m);
        default:
            return gmt;
        }
    }
}

/**
 *  hqr_nextpiv - Computes the next row killed by the row p, after
 *  it has kill the row start.
 *
 * @param[in] k
 *         Factorization step
 *
 * @param[in] pivot
 *         Line used as killer
 *
 * @param[in] start
 *         Starting point to search the next line killed by p after start
 *         start must be equal to A.mt to find the first row killed by p.
 *         if start != A.mt, start must be killed by p
 *
 * @return:
 *   - -1 if start doesn't respect the previous conditions
 *   -  m, the following row killed by p if it exists, A->mt otherwise
 */
static int
hqr_nextpiv(const libhqr_tree_t *qrtree, int k, int pivot, int start)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, ls, lp, nextp;
    int lpivot, rpivot, lstart;
    int a = qrtree->a;
    int p = qrtree->p;
    int gmt = qrtree->mt;

    lpivot = pivot / p; /* Local index in the distribution over p domains */
    rpivot = pivot % p; /* Staring index in this distribution             */

    /* Local index in the distribution over p domains */
    lstart = ( start == gmt ) ? arg->llvl->ldd * a : start / p;

    myassert( start > pivot && pivot >= k );
    myassert( start == gmt || pivot == hqr_currpiv( qrtree, k, start ) );

    /* TS level common to every case */
    ls = (start < gmt) ? hqr_gettype( qrtree, k, start ) : -1;
    lp = hqr_gettype( qrtree, k, pivot );

    switch( ls )
        {
        case -1:

            if ( lp == LIBHQR_KILLED_BY_TS ) {
                myassert( start == gmt );
                return gmt;
            }

        case LIBHQR_KILLED_BY_TS:

            /* If the tile is over the diagonal of step k, skip directly to type 2 */
            if ( arg->domino && lpivot < k )
                goto next_2;

            if ( start == gmt )
                nextp = pivot + p;
            else
                nextp = start + p;

            if ( ( nextp < gmt ) &&
                 ( nextp < pivot + a*p ) &&
                 ( (nextp/p)%a != 0 ) )
                return nextp;
            start = gmt;
            lstart = arg->llvl->ldd * a;

        case LIBHQR_KILLED_BY_LOCALTREE:

            /* If the tile is over the diagonal of step k, skip directly to type 2 */
            if ( arg->domino && lpivot < k )
                goto next_2;

            /* Get the next pivot for the low level tree */
            tmp = arg->llvl->nextpiv(arg->llvl, k, pivot, lstart / a );

            if ( (tmp * a * p + rpivot >= gmt)
                 && (tmp == arg->llvl->ldd-1) )
                tmp = arg->llvl->nextpiv(arg->llvl, k, pivot, tmp);

            if ( tmp != arg->llvl->ldd )
                return tmp * a * p + rpivot;

        next_2:
            /* no next of type 1, we reset start to search the next 2 */
            start = gmt;
            lstart = arg->llvl->ldd * a;

        case LIBHQR_KILLED_BY_DOMINO:

            if ( lp < LIBHQR_KILLED_BY_DOMINO ) {
                return gmt;
            }

            /* Type 2 are killed only once if they are strictly in the band */
            if ( arg->domino &&
                 (start == gmt) &&
                 (lpivot < k)             &&
                 (pivot+p < gmt) ) {
                return pivot+p;
            }

            /* no next of type 2, we reset start to search the next 3 */
            start = gmt;
            lstart = arg->llvl->ldd * a;

        case LIBHQR_KILLED_BY_DISTTREE:

            if ( lp < LIBHQR_KILLED_BY_DISTTREE ) {
                return gmt;
            }

            if( arg->hlvl != NULL ) {
                tmp = arg->hlvl->nextpiv( arg->hlvl, k, pivot, start );
                if ( tmp != gmt )
                    return tmp;
            }

        default:
            return gmt;
        }
}

/**
 *  hqr_prevpiv - Computes the previous row killed by the row p, before
 *  to kill the row start.
 *
 * @param[in] k
 *         Factorization step
 *
 * @param[in] pivot
 *         Line used as killer
 *
 * @param[in] start
 *         Starting point to search the previous line killed by p before start
 *         start must be killed by p, and start must be greater or equal to p
 *
 * @return:
 *   - -1 if start doesn't respect the previous conditions
 *   -  m, the previous row killed by p if it exists, A->mt otherwise
 */
static int
hqr_prevpiv(const libhqr_tree_t *qrtree, int k, int pivot, int start)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, ls, lp, nextp;
    int lpivot, rpivot, lstart;
    int a = qrtree->a;
    int p = qrtree->p;
    int gmt = qrtree->mt;

    lpivot = pivot / p; /* Local index in the distribution over p domains */
    rpivot = pivot % p; /* Staring index in this distribution             */
    lstart = start / p; /* Local index in the distribution over p domains */

    myassert( start >= pivot && pivot >= k && start < gmt );
    myassert( start == pivot || pivot == hqr_currpiv( qrtree, k, start ) );

    /* T Slevel common to every case */
    ls = hqr_gettype( qrtree, k, start );
    lp = hqr_gettype( qrtree, k, pivot );

    if ( lp == LIBHQR_KILLED_BY_TS )
      return gmt;

    myassert( lp >= ls );
    switch( ls )
        {
        case LIBHQR_KILLED_BY_DISTTREE:
            if( arg->hlvl != NULL ) {
                tmp = arg->hlvl->prevpiv( arg->hlvl, k, pivot, start );
                if ( tmp != gmt )
                    return tmp;
            }

            start = pivot;
            lstart = pivot / p;

        case LIBHQR_KILLED_BY_DOMINO:
            /* If the tile is over the diagonal of step k, process it as type 2 */
            if ( arg->domino && lpivot < k ) {

                if ( ( start == pivot ) &&
                     (start+p < gmt ) )
                    return start+p;

                if ( lp > LIBHQR_KILLED_BY_LOCALTREE )
                    return gmt;
            }

            start = pivot;
            lstart = pivot / p;

            /* If it is the 'local' diagonal block, we go to 1 */

        case LIBHQR_KILLED_BY_LOCALTREE:
            /* If the tile is over the diagonal of step k and is of type 2,
               it cannot annihilate type 0 or 1 */
            if ( arg->domino && lpivot < k )
                return gmt;

            tmp = arg->llvl->prevpiv(arg->llvl, k, pivot, lstart / a);

            if ( (tmp * a * p + rpivot >= gmt)
                 && (tmp == arg->llvl->ldd-1) )
                tmp = arg->llvl->prevpiv(arg->llvl, k, pivot, tmp);

            if ( tmp != arg->llvl->ldd )
                return tmp * a * p + rpivot;

            start = pivot;

        case LIBHQR_KILLED_BY_TS:
            /* Search for predecessor in TS tree */
            /* if ( ( start+p < gmt ) &&  */
            /*      ( (((start+p) / p) % a) != 0 ) ) */
            /*     return start + p; */

            if ( start == pivot ) {
                tmp = lpivot + a - 1 - lpivot%a;
                nextp = tmp * p + rpivot;

                while( pivot < nextp && nextp >= gmt )
                    nextp -= p;
            } else {
                nextp = start - p; /*(lstart - 1) * p + rpivot;*/
            }
            assert(nextp < gmt);
            if ( pivot < nextp )
                return nextp;

        default:
            return gmt;
        }
}

/****************************************************
 *
 * Generate the permutation required for the round-robin on TS
 *
 ***************************************************/
static void
hqr_genperm( libhqr_tree_t *qrtree )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int m = qrtree->mt;
    int n = qrtree->nt;
    int a = qrtree->a;
    int p = qrtree->p;
    int domino = arg->domino;
    int minMN = libhqr_imin( m, n );
    int pa = p * a;
    int i, j, k;
    int nbextra1;
    int end2;
    int mpa   = m % pa;
    int endpa = m - mpa;
    int *perm;

    arg->perm = (int*)malloc( (m+1) * minMN * sizeof(int) );
    perm = arg->perm;

    if ( arg->tsrr ) {
        for(k=0; k<minMN; k++) {
            for( i=0; i<m+1; i++) {
                perm[i] = -1;
            }
            perm += m+1;
        }
        perm = arg->perm;
        for(k=0; k<minMN; k++) {
            nbextra1 = nbextra1_formula;

            end2 = p + ( domino ? k*p : k + nbextra1 );
            end2 = (( end2 + pa - 1 ) / pa ) * pa;
            end2 = libhqr_imin( end2, m );

            /*
             * All tiles of type 3, 2 and:
             * - 1 when domino is disabled
             * - 0 before the first multiple of pa under the considered diagonal
             */
            for( i=k; i<end2; i++) {
                perm[i] = i;
            }

            /* All permutations in groups of pa tiles */
            assert( i%pa == 0 || i>=endpa);
            for( ; i<endpa; i+=pa ) {
                for(j=0; j<pa; j++) {
                    perm[i+j] = i + ( j + p * (k%a) )%pa;
                }
            }

            /* Last group of tiles */
            for( ; i<m; i++) {
                perm[i] = i;
            }
            perm[m] = m;
            perm += m+1;
        }
    }
    else {
        for(k=0; k<minMN; k++) {
            for( i=0; i<m+1; i++) {
                perm[i] = i;
            }
            perm += m+1;
        }
    }
}

static inline int
hqr_getinvperm( const libhqr_tree_t *qrtree, int k, int m )
{
    int gmt = qrtree->mt + 1;
    int a = qrtree->a;
    int p = qrtree->p;
    int pa = p * a;
    int start = m / pa * pa;
    int stop  = libhqr_imin( start + pa, gmt ) - start;
    int i;

    if (a == 1)
        return m;

    for ( i=m%p; i < stop; i+=p ) {
        //if( perm[i] == m )
        if( i == m )
            return i+start;
    }

    /* We should never arrive here */
    myassert( 0 );
    (void)k;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma
 *
 * libhqr_hqr_init - Creates the tree structure that will describes the
 * operation performed during QR/LQ factorization with parameterized QR/LQ
 * algorithms family.
 *
 * Trees available parameters are described below. It is recommended to:
 *   - set p to the same value than the P-by-Q process grid used to distribute
 *     the data. (P for QR factorization, Q for LQ factorization).
 *   - set the low level tree to LIBHQR_GREEDY_TREE.
 *   - set the high level tree to:
 *         1) LIBHQR_FLAT_TREE when the problem is square, because it divides
 *            by two the volume of communication of any other tree.
 *         2) LIBHQR_FIBONACCI_TREE when the problem is tall and skinny (QR) or
 *            small and fat (LQ), because it reduces the critical path length.
 *   - Disable the domino effect when problem is square, to keep high efficiency
 *     kernel proportion high.
 *   - Enable the domino effect when problem is tall and skinny (QR) or
 *     small and fat (LQ) to increase parallelism and reduce critical path length.
 *   - Round-robin on TS domain (tsrr) option should be disabled. It is
 *     experimental and is not safe.
 *
 * These are the default when a parameter is set to -1;
 *
 * See http://www.netlib.org/lapack/lawnspdf/lawn257.pdf
 *
 *******************************************************************************
 *
 * @param[inout] qrtree
 *          On entry, an allocated structure uninitialized.
 *          On exit, the structure initialized according to the parameter given.
 *
 * @param[in] trans
 *          @arg QR:   Structure is initialized for QR factorization.
 *          @arg LQ:     Structure is initialized for LQ factorization.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A to be factorized, on which
 *          QR/LQ factorization will be performed.
 *          The descriptor is untouched and only mt/nt/P parameters are used.
 *
 * @param[in] type_llvl
 *          Defines the tree used to reduce the main tiles of each local domain
 *          together. The matrix of those tiles has a lower triangular structure
 *          with a diagonal by step a.
 *          @arg LIBHQR_FLAT_TREE: A Flat tree is used to reduce the local
 *          tiles.
 *          @arg LIBHQR_GREEDY_TREE: A Greedy tree is used to reduce the local
 *          tiles.
 *          @arg LIBHQR_FIBONACCI_TREE: A Fibonacci tree is used to reduce the
 *          local tiles.
 *          @arg LIBHQR_BINARY_TREE: A Binary tree is used to reduce the local
 *          tiles.
 *          @arg -1: The default is used (LIBHQR_GREEDY_TREE)
 *
 * @param[in] type_hlvl
 *          Defines the tree used to reduce the main tiles of each domain. This
 *          is a band lower diagonal matrix of width p.
 *          @arg LIBHQR_FLAT_TREE: A Flat tree is used to reduce the tiles.
 *          @arg LIBHQR_GREEDY_TREE: A Greedy tree is used to reduce the tiles.
 *          @arg LIBHQR_FIBONACCI_TREE: A Fibonacci tree is used to reduce the
 *          tiles.
 *          @arg LIBHQR_BINARY_TREE: A Binary tree is used to reduce the tiles.
 *          @arg LIBHQR_GREEDY1P_TREE: A Greedy tree is computed for the first
 *          column and then duplicated on all others.
 *          @arg -1: The default is used (LIBHQR_FIBONACCI_TREE)
 *
 * @param[in] a
 *          Defines the size of the local domains on which a classic flat TS
 *          tree is performed. If a==1, then all local tiles are reduced by the
 *          type_lllvl tree. If a is larger than mt/p, then no local reduction
 *          tree is performed and type_llvl is ignored.
 *          If a == -1, it is set to 4 by default.
 *
 * @param[in] p
 *          Defines the number of distributed domains, ie the width of the high
 *          level reduction tree.  If p == 1, no high level reduction tree is
 *          used. If p == mt, a and type_llvl are ignored since only high level
 *          reduction are performed.
 *          By default, it is recommended to set p to P if trans ==
 *          PlasmaNoTrans, Q otherwise, where P-by-Q is the process grid used to
 *          distributed the data. (p > 0)
 *
 * @param[in] domino
 *          Enable/disable the domino effect that connects the high and low
 *          level reduction trees. Enabling the domino increases the proportion
 *          of TT (Triangle on top of Triangle) kernels that are less efficient,
 *          but increase the pipeline effect between independent factorization
 *          steps, reducin the critical path length.
 *          If disabled, it keeps the proprotion of more efficient TS (Triangle
 *          on top of Square) kernels high, but deteriorates the pipeline
 *          effect, and the critical path length.
 *          If domino == -1, it is enable when ration between M and N is lower
 *          than 1/2, and disabled otherwise.
 *
 * @param[in] tsrr
 *          Enable/Disable a round robin selection of the killer in local
 *          domains reduced by TS kernels. Enabling a round-robin selection of
 *          the killer allows to take benefit of the good pipelining of the flat
 *          trees and the high efficient TS kernels, while having other trees on
 *          top of it to reduce critical path length.
 *          WARNING: This option is under development and should not be enabled
 *          due to problem in corner cases.
 *
 *******************************************************************************
 *
 * @retval -i if the ith parameters is incorrect.
 * @retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa libhqr_hqr_finalize
 * @sa libhqr_systolic_init
 *
 ******************************************************************************/
int
libhqr_initfct_hqr( libhqr_tree_t *qrtree,
                    libhqr_facto_e trans, libhqr_matrix_t *A,
                    int type_llvl, int type_hlvl,
                    int a, int p,
                    int domino, int tsrr )
{
    double ratio = 0.0;
    int low_mt, minMN;
    hqr_args_t *arg;

    if (qrtree == NULL) {
      fprintf(stderr, "libhqr_hqr_init, illegal value of qrtree");
      return -1;
    }
    if ((trans != LIBHQR_QR) &&
        (trans != LIBHQR_LQ)) {
      fprintf(stderr, "libhqr_hqr_init, illegal value of trans");
      return -2;
    }
    if (A == NULL) {
      fprintf(stderr, "libhqr_hqr_init, illegal value of A");
      return -3;
    }

    /* Compute parameters */
    if ( a == -1) {
        a = 4; /* TODO: add automatic computation */
    }
    else {
        a = libhqr_imax( a, 1 );
    }

    if ( p == -1 ) {
        if ( trans == LIBHQR_QR ) {
            p = A->p;
        }
        else {
            p = A->nodes / A->p;
        }
    }
    else {
        p = libhqr_imax( p, 1 );
    }

    /* Domino */
    if ( domino >= 0 ) {
        domino = domino ? 1 : 0;
    }
    else {
        if (trans == LIBHQR_QR) {
            ratio = ((double)(A->nt) / (double)(A->mt));
        } else {
            ratio = ((double)(A->mt) / (double)(A->nt));
        }
        if ( ratio >= 0.5 ) {
            domino = 0;
        } else {
            domino = 1;
        }
    }

    minMN = libhqr_imin(A->mt, A->nt);

    qrtree->domino     = domino;
    qrtree->getnbgeqrf = hqr_getnbgeqrf;
    qrtree->getm       = hqr_getm;
    qrtree->geti       = hqr_geti;
    qrtree->gettype    = hqr_gettype;
    qrtree->currpiv    = hqr_currpiv;
    qrtree->nextpiv    = hqr_nextpiv;
    qrtree->prevpiv    = hqr_prevpiv;

    qrtree->init = LIBHQR_QRTREE_HQR;

    qrtree->mt = (trans == LIBHQR_QR) ? A->mt : A->nt;
    qrtree->nt = minMN;

    a = libhqr_imin( a, qrtree->mt );

    qrtree->a    = a;
    qrtree->p    = p;
    qrtree->args = NULL;

    arg = (hqr_args_t*) malloc( sizeof(hqr_args_t) );
    arg->domino = domino;
    arg->tsrr = tsrr;
    arg->perm = NULL;

    arg->llvl = (hqr_subpiv_t*) malloc( sizeof(hqr_subpiv_t) );
    arg->hlvl = NULL;

    low_mt = (qrtree->mt + p * a - 1) / ( p * a );

    arg->llvl->minMN  = minMN;
    arg->llvl->ldd    = low_mt;
    arg->llvl->a      = a;
    arg->llvl->p      = p;
    arg->llvl->domino = domino;

    switch( type_llvl ) {
    case LIBHQR_FLAT_TREE :
        hqr_low_flat_init(arg->llvl);
        break;
    case LIBHQR_FIBONACCI_TREE :
        hqr_low_fibonacci_init(arg->llvl, minMN);
        break;
    case LIBHQR_BINARY_TREE :
        hqr_low_binary_init(arg->llvl);
        break;
    case LIBHQR_GREEDY1P_TREE :
        hqr_low_greedy1p_init(arg->llvl, minMN);
        break;
    case LIBHQR_GREEDY_TREE :
    default:
        hqr_low_greedy_init(arg->llvl, minMN);
    }

    if ( p > 1 ) {
        arg->hlvl = (hqr_subpiv_t*) malloc( sizeof(hqr_subpiv_t) );

        arg->llvl->minMN  = minMN;
        arg->hlvl->ldd    = qrtree->mt;
        arg->hlvl->a      = a;
        arg->hlvl->p      = p;
        arg->hlvl->domino = domino;

        switch( type_hlvl ) {
        case LIBHQR_FLAT_TREE :
            hqr_high_flat_init(arg->hlvl);
            break;
        case LIBHQR_GREEDY_TREE :
            hqr_high_greedy_init(arg->hlvl, minMN);
            break;
        case LIBHQR_GREEDY1P_TREE :
            hqr_high_greedy1p_init(arg->hlvl);
            break;
        case LIBHQR_BINARY_TREE :
            hqr_high_binary_init(arg->hlvl);
            break;
        case LIBHQR_FIBONACCI_TREE :
            hqr_high_fibonacci_init(arg->hlvl);
            break;
        default:
            if ( ratio >= 0.5 ) {
                hqr_high_flat_init(arg->hlvl);
            } else {
                hqr_high_fibonacci_init(arg->hlvl);
            }
        }
    }

    qrtree->args = (void*)arg;
    hqr_genperm( qrtree );

    return 0;
}

int
libhqr_initmtx_hqr( libhqr_tree_t *qrtree,
                    libhqr_facto_e trans, libhqr_matrix_t *A,
                    int type_llvl, int type_hlvl,
                    int a, int p,
                    int domino, int tsrr )
{
    libhqr_tree_t qrtreefct;
    int rc;

    rc = libhqr_initfct_hqr( &qrtreefct, trans, A,
                             type_llvl, type_hlvl, a, p,
                             domino, tsrr );
    if ( rc != 0 ) {
        return rc;
    }

    libhqr_fct_to_mtx( &qrtreefct, qrtree );

    /* Free the initial qrtree */
    libhqr_finalize( &qrtreefct );

    return 0;
}

int
libhqr_init_hqr( libhqr_tree_t *qrtree,
                 libhqr_facto_e trans, libhqr_matrix_t *A,
                 int type_llvl, int type_hlvl,
                 int a, int p,
                 int domino, int tsrr )
{
    return libhqr_initmtx_hqr( qrtree, trans, A,
                               type_llvl, type_hlvl, a, p,
                               domino, tsrr );

}

/**
 *  libhqr_walk_stepk
 *
 * @param[in] arg
 *         Arguments specific to the reduction tree used
 *
 * @param[in] k
 *         Factorization step
 *
 * @param[in] tiles
 *         Array of size qrtree->mt that stores in the first entries the tiles
 *         to kill in the right order to walk through the tree.
 *
 */
void
libhqr_walk_stepk( const libhqr_tree_t *qrtree,
                   int k, int *tiles )
{
    int tsid = -1, ttid = -1;
    int p = k;
    int pivot = k;
    int m = 0;
    libhqr_queue_tile_t *tt = libhqr_queue_tile_new();
    libhqr_queue_tile_t *ts = libhqr_queue_tile_new();
    libhqr_queue_tile_push(&tt, k);

    while( tt )
    {
        /* Stack all elements killed on the way down */
        p = qrtree->prevpiv(qrtree, k, pivot, p);
        while( p != qrtree->mt ) {
            /* If it's a TS tile, or a TT at the bottom of the tree that kills noone */
            if( (qrtree->gettype(qrtree, k, p) == 0) ||
                (qrtree->prevpiv(qrtree, k, p, p) != qrtree->mt) )
            {
                libhqr_queue_tile_push(&tt,p);
                /* printf("PUSH TT: %d\n", p); */
            }
            libhqr_queue_tile_push(&ts, p);
            /* printf("PUSH TS: %d\n", p); */
            p = qrtree->prevpiv(qrtree, k, pivot, p);
            assert( p != -1 );
        }

        tsid = libhqr_queue_tile_head(&ts);
        ttid = libhqr_queue_tile_pop(&tt);

        /* Unstack all element killed by TS, or by TT and already discovered */
        while( (tsid != -1) &&
               (tsid != ttid) )
        {
            /* printf( "TS=%d, TT=%d\n", tsid, ttid ); */

            tsid = libhqr_queue_tile_pop(&ts);
            /* printf("POP TS: %d\n", tsid); */
            /* printf("Call function on (%d, %d)\n",
             qrtree->currpiv(qrtree, k, tsid), tsid ); */
            assert(m < (qrtree->mt-1));
            tiles[m] = tsid;
            m++;
            tsid = libhqr_queue_tile_head(&ts);
        }
        pivot = p = ttid;
    }

    //assert(m  == (qrtree->mt-1));
    assert(ts == NULL);
    assert(tt == NULL);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma
 *
 * libhqr_hqr_finalize - Cleans the qrtree data structure allocated by call to
 * libhqr_hqr_init().
 *
 *******************************************************************************
 *
 * @param[in,out] qrtree
 *          On entry, an allocated structure to destroy.
 *          On exit, the structure is destroy and cannot be used.
 *
 *******************************************************************************
 *
 * @sa libhqr_hqr_init
 *
 ******************************************************************************/
void
libhqr_finalize( libhqr_tree_t *qrtree )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);

    if ( (qrtree->init != LIBHQR_QRTREE_HQR) &&
         (qrtree->init != LIBHQR_QRTREE_SVD) &&
         (qrtree->init != LIBHQR_QRTREE_SYS) &&
         (qrtree->init != LIBHQR_QRTREE_MTX) )
    {
        fprintf(stderr, "libhqr_finalize: qrtree not initialized\n" );
        return;
    }

    if (arg != NULL) {

        if ( (qrtree->init == LIBHQR_QRTREE_HQR) ||
             (qrtree->init == LIBHQR_QRTREE_HQR) )
        {
            if ( arg->llvl != NULL) {
                if ( arg->llvl->ipiv != NULL )
                    free( arg->llvl->ipiv );
                free( arg->llvl );
            }

            if ( arg->hlvl != NULL) {
                if ( arg->hlvl->ipiv != NULL )
                    free( arg->hlvl->ipiv );
                free( arg->hlvl );
            }

            if ( arg->perm != NULL )
                free(arg->perm);
        }

        free(arg);
    }
}
