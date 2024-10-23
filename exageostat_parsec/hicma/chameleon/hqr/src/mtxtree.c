/**
 *
 * @file mtxtree.c
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
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/****************************************************
 * Reading functions
 **************************************************
 *
 * Common parameters for the following functions
 *    qrtree - tree used for the factorization
 *    k      - Step k of the QR factorization
 *    m      - line anhilated
 */

int  mtxtree_gettype(const libhqr_tree_t *qrtree, int k, int m){
    int id;
    libhqr_tile_info_t *arg = (libhqr_tile_info_t*)(qrtree->args);
    id = (k * qrtree->mt) + m;
    return arg[id].type;
}

int mtxtree_currpiv(const libhqr_tree_t *qrtree, int k, int m){
    int id, perm_m, p, a;
    libhqr_tile_info_t *arg = (libhqr_tile_info_t*)(qrtree->args);
    perm_m = m;
    p = qrtree->p;
    a = qrtree->a;
    myassert( (p==1) || (perm_m / (p*a)) == (m / (p*a)) );
    myassert( (p==1) || (perm_m % p) == (m % p) );
    id = (k * qrtree->mt) + m;
    return arg[id].currpiv;
}

/*
 * Extra parameter:
 *    p - line used as pivot
 */

int mtxtree_nextpiv(const libhqr_tree_t *qrtree, int k, int p, int m){
    int id;
    libhqr_tile_info_t *arg = (libhqr_tile_info_t*)(qrtree->args);
    int gmt = qrtree->mt;
    myassert( m > p && p >= k );
    myassert( m == gmt || p == mtxtree_currpiv( qrtree, k, m ) );
    if(m == qrtree->mt){
        id = (k * qrtree->mt) + p;
        return arg[id].first_nextpiv;
    }
    else{
        id = (k * qrtree->mt) + m;
        return arg[id].nextpiv;
    }
}

int mtxtree_prevpiv(const libhqr_tree_t *qrtree, int k, int p, int m){
    int id;
    libhqr_tile_info_t *arg = (libhqr_tile_info_t*)(qrtree->args);
    int gmt = qrtree->mt;
    myassert( m >= p && p >= k && m < gmt );
    myassert( m == p || p == mtxtree_currpiv( qrtree, k, m ) );

    if(m == p){
        id = (k * qrtree->mt) + p;
        return arg[id].first_prevpiv;
    }
    else{
        id = (k * qrtree->mt) + m;
        return arg[id].prevpiv;
    }
}

void
libhqr_fct_to_mtx( const libhqr_tree_t *in, libhqr_tree_t *out )
{
    libhqr_tile_info_t *tileinfo;
    int i, minMN, p, k;

    minMN = libhqr_imin( in->mt, in->nt );

    /* Copy the input tree to the output one */
    memcpy( out, in, sizeof(libhqr_tree_t) );

    /* Switch to matrix storage format and functions */
    out->init    = LIBHQR_QRTREE_MTX;
    out->gettype = mtxtree_gettype;
    out->currpiv = mtxtree_currpiv;
    out->nextpiv = mtxtree_nextpiv;
    out->prevpiv = mtxtree_prevpiv;
    out->args    = malloc( out->mt * minMN * sizeof(libhqr_tile_info_t) );

    tileinfo = (libhqr_tile_info_t*)(out->args);

    /* Initialize the matrix */
    for (k=0; k<minMN; k++)
    {
        for (i=0; i<in->mt; i++, tileinfo++)
        {
            tileinfo->type          = in->gettype(in, k, i);
            tileinfo->currpiv = p   = in->currpiv(in, k, i);
            tileinfo->first_nextpiv = in->nextpiv(in, k, i, in->mt);
            tileinfo->first_prevpiv = in->prevpiv(in, k, i, i);
            tileinfo->nextpiv       = in->nextpiv(in, k, p, i);
            tileinfo->prevpiv       = in->prevpiv(in, k, p, i);

            assert( ((tileinfo->nextpiv       <= in->mt) && (tileinfo->nextpiv       >= k)) || (tileinfo->nextpiv       == -1) );
            assert( ((tileinfo->prevpiv       <= in->mt) && (tileinfo->prevpiv       >= k)) || (tileinfo->prevpiv       == -1) );
            assert( ((tileinfo->first_nextpiv <= in->mt) && (tileinfo->first_nextpiv >= k)) || (tileinfo->first_nextpiv == -1) );
            assert( ((tileinfo->first_prevpiv <= in->mt) && (tileinfo->first_prevpiv >= k)) || (tileinfo->first_prevpiv == -1) );
        }
    }
}
