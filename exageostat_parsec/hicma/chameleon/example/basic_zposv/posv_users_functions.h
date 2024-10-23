/**
 *
 * @file posv_users_functions.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon posv_users_functions example header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-13
 *
 */
#ifndef POSV_USERS_FUNCTIONS_H
#define POSV_USERS_FUNCTIONS_H

#include "basic_posv.h"

/**
 *  Function to return address of block (m,n)
 */
inline static void* user_getaddr_arrayofpointers(const MORSE_desc_t *A, int m, int n)
{
    MORSE_Complex64_t **matA = (MORSE_Complex64_t **)A->mat;
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    offset = A->mt*nn + mm;
//    printf ("Array address in %d %d: %p\n", mm, nn, *(matA+(offset)));
    return (void*)( *(matA + offset) );
}

/**
 *  Function to return the leading dimension of element A(m,*)
 */
inline static int user_getblkldd_arrayofpointers(const MORSE_desc_t *A, int m)
{
    (void)m;
    return A->mb;
}

/**
 *  Function to return MPI rank of element A(m,n)
 */
inline static int user_getrankof_2d(const MORSE_desc_t *desc, int m, int n)
{
    return (m % desc->p) * desc->q + (n % desc->q);
}

#endif /* POSV_USERS_FUNCTIONS_H */
