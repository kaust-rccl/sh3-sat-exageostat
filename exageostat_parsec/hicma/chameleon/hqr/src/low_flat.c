/**
 *
 * @file low_flat.c
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
 * Functions for low level flat tree
 *
 */
#include "libhqr_internal.h"

static inline int
hqr_low_flat_currpiv(const hqr_subpiv_t *arg, int k, int m)
{
    (void)m;
    if ( arg->domino )
        return k / arg->a;
    else
        return (k + arg->p - 1 - m%(arg->p)) / arg->p / arg->a ;
}

static inline int
hqr_low_flat_nextpiv(const hqr_subpiv_t *arg, int k, int p, int start_pa)
{
    int k_a = arg->domino ? k / arg->a :  (k + arg->p - 1 - p%(arg->p)) / arg->p / arg->a;
    int p_pa = (p / arg->p ) / arg->a;

#ifdef FLAT_UP
    if ( ( p_pa == k_a ) && (start_pa > k_a+1 ) )
        return start_pa-1;
    else
#else /* FLAT_DOWN */
        if ( start_pa <= p_pa )
            return arg->ldd;

    if ( p_pa == k_a && ( arg->ldd - k_a ) > 1 ) {
        if ( start_pa == arg->ldd )
            return p_pa+1;
        else if ( start_pa < arg->ldd )
            return start_pa+1;
    }
#endif
    return arg->ldd;
}

static inline int
hqr_low_flat_prevpiv(const hqr_subpiv_t *arg, int k, int p, int start_pa)
{
    int k_a = arg->domino ? k / arg->a :  (k + arg->p - 1 - p%(arg->p)) / arg->p / arg->a;
    int p_pa = (p / arg->p ) / arg->a;

#ifdef FLAT_UP
    if ( p_pa == k_a && (start_pa+1 < arg->ldd) )
        return start_pa+1;
    else
#else
        if ( p_pa == k_a && ( arg->ldd - k_a ) > 1  ) {
            if ( start_pa == p_pa )
                return arg->ldd - 1;
            else if ( start_pa > p_pa + 1 )
                return start_pa-1;
        }
#endif
    return arg->ldd;
}

void hqr_low_flat_init(hqr_subpiv_t *arg){
    arg->currpiv = hqr_low_flat_currpiv;
    arg->nextpiv = hqr_low_flat_nextpiv;
    arg->prevpiv = hqr_low_flat_prevpiv;
    arg->ipiv = NULL;
}
