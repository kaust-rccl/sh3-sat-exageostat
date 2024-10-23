/**
 *
 * @file runtime_descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark descriptor routines
 *
 * @version 1.0.0
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#include <stdlib.h>
#include "chameleon_quark.h"

void RUNTIME_comm_set_tag_sizes( int user_tag_width,
                                 int user_tag_sep )
{
    (void)user_tag_width;
    (void)user_tag_sep;
}

void *RUNTIME_malloc( size_t size )
{
    return malloc( size );
}

void RUNTIME_free( void  *ptr,
                   size_t size )
{
    (void)size;
    free( ptr );
    return;
}

void RUNTIME_desc_create( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

void RUNTIME_desc_destroy( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

int RUNTIME_desc_acquire( const MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

int RUNTIME_desc_release( const MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

void
RUNTIME_desc_flush( const MORSE_desc_t     *desc,
                    const MORSE_sequence_t *sequence )
{
    (void)desc;
    (void)sequence;
    return;
}


void
RUNTIME_flush( )
{
    return;
}

void
RUNTIME_data_flush( const MORSE_sequence_t *sequence,
                    const MORSE_desc_t *A, int Am, int An )
{
    (void)sequence;
    (void)A;
    (void)Am;
    (void)An;
    return;
}

#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const MORSE_sequence_t *sequence,
                           const MORSE_desc_t *A, int Am, int An, int new_rank )
{
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
}
#endif

void *RUNTIME_data_getaddr( const MORSE_desc_t *desc, int m, int n )
{
    return desc->get_blkaddr( desc, m, n );
}
