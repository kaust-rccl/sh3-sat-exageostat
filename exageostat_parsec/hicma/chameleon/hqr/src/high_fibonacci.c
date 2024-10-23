/**
 *
 * @file high_fibonacci.c
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
 * Functions for high level fibonacci tree, and init for duplicated greedy.
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

/****************************************************
 *          HQR_HIGH_FIBONACCI_TREE
 ***************************************************/
/* Return the pivot to use for the row m at step k */
static inline int
hqr_high_fibonacci_currpiv( const hqr_subpiv_t *qrpiv, int k, int m ) {
    return (qrpiv->ipiv)[ m-k ] + k;
}

/* Return the last row which has used the row m as a pivot in step k before the row start */
static inline int
hqr_high_fibonacci_prevpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start ) {
    int i;
    myassert( p >= k && start >= p && start-k <= qrpiv->p);

    int lp    = p - k;
    int lstart= start - k;
    int end   = libhqr_imin(qrpiv->ldd-k, qrpiv->p);
    for( i=lstart+1; i<end; i++ )
        if ( (qrpiv->ipiv)[i] == lp )
            return i+k;
    return qrpiv->ldd;
}

/* Return the next row which will use the row m as a pivot in step k after it has been used by row start */
static inline int
hqr_high_fibonacci_nextpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start ) {
    int i;
    myassert( p>=k && (start == qrpiv->ldd || start-k <= qrpiv->p) );

    for( i=libhqr_imin(start-k-1, qrpiv->p-1); i>0; i-- )
        if ( (qrpiv->ipiv)[i] == (p-k) )
            return i + k;
    return (qrpiv->ldd);
}

void
hqr_high_fibonacci_init(hqr_subpiv_t *arg) {
    int *ipiv;
    int p;

    arg->currpiv = hqr_high_fibonacci_currpiv;
    arg->nextpiv = hqr_high_fibonacci_nextpiv;
    arg->prevpiv = hqr_high_fibonacci_prevpiv;

    p = arg->p;

    arg->ipiv = (int*)calloc( p, sizeof(int) );
    ipiv = arg->ipiv;

    /*
     * Fibonacci of order 1:
     *    u_(n+1) = u_(n) + 1
     */
    {
        int f1, k, m;

        /* Fill in the first column */
        f1 = 1;
        for (m=1; m < p; ) {
            for (k=0; (k < f1) && (m < p); k++, m++) {
                ipiv[m] = m - f1;
            }
            f1++;
        }
    }
}

/****************************************************
 *                 HQR_HIGH_GREEDY_TREE (1 panel duplicated)
 ***************************************************/
void hqr_high_greedy1p_init(hqr_subpiv_t *arg){
    int *ipiv;
    int mt, p;

    arg->currpiv = hqr_high_fibonacci_currpiv;
    arg->nextpiv = hqr_high_fibonacci_nextpiv;
    arg->prevpiv = hqr_high_fibonacci_prevpiv;

    mt = arg->ldd;
    p = arg->p;

    arg->ipiv = (int*)calloc( p, sizeof(int) );
    ipiv = arg->ipiv;

    {
        int minMN = 1;
        int j, k, height, start, end, firstk = 0;
        int *nT = (int*)calloc(minMN, sizeof(int));
        int *nZ = (int*)calloc(minMN, sizeof(int));

        nT[0] = mt;
        nZ[0] = libhqr_imax( mt - p, 0 );
        for(k=1; k<minMN; k++) {
            height = libhqr_imax(mt-k-p, 0);
            nT[k] = height;
            nZ[k] = height;
        }

        k = 0;
        while ( (!( ( nT[minMN-1] == mt - (minMN - 1) ) &&
                    ( nZ[minMN-1]+1 == nT[minMN-1] ) ) )
                && ( firstk < minMN ) ) {
            height = (nT[k] - nZ[k]) / 2;
            if ( height == 0 ) {
                while ( (firstk < minMN) &&
                        ( nT[firstk] == mt - firstk ) &&
                        ( nZ[firstk]+1 == nT[firstk] ) ) {
                    firstk++;
                }
                k = firstk;
                continue;
            }

            start = mt - nZ[k] - 1;
            end = start - height;
            nZ[k] += height;
            if (k < minMN-1) nT[k+1] = nZ[k];

            for( j=start; j > end; j-- ) {
                ipiv[ k*p + j-k ] = (j - height);
            }

            k++;
            if (k > minMN-1) k = firstk;
        }

        free(nT);
        free(nZ);
    }
}
