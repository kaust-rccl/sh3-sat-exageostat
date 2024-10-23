/**
 *
 * @file runtime_descriptor.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC descriptor routines
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2017-01-12
 *
 */
#include <stdlib.h>
#include "chameleon_parsec.h"
#include <parsec/data.h>
#include <parsec/datatype.h>
#include <parsec/arena.h>

#if defined(CHAMELEON_USE_MPI)

/* Variable parsec_dtd_no_of_arenas is private and cannot be changed */
#define MORSE_PARSEC_DTD_NO_OF_ARENA 16 /**< Number of arenas available per DTD */

typedef struct morse_parsec_arena_s {
    /* int mb; */
    /* int nb; */
    /* MORSE_enum dtype; */
    size_t size;
} morse_parsec_arena_t;

static int morse_parsec_nb_arenas = 0;
static morse_parsec_arena_t morse_parsec_registered_arenas[MORSE_PARSEC_DTD_NO_OF_ARENA] = { { 0 } };

#endif

void RUNTIME_comm_set_tag_sizes( int user_tag_width,
                                 int user_tag_sep )
{
    (void)user_tag_width;
    (void)user_tag_sep;
}

void *RUNTIME_malloc( size_t size )
{
    return malloc(size);
}

void RUNTIME_free( void *ptr, size_t size )
{
    (void)size;
    free(ptr);
    return;
}

static inline void
morse_parsec_key_to_coordinates(parsec_data_collection_t *data_collection, parsec_data_key_t key,
                                int *m, int *n)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    int _m, _n;

    _m = key % mdesc->lmt;
    _n = key / mdesc->lmt;
    *m = _m - mdesc->i / mdesc->mb;
    *n = _n - mdesc->j / mdesc->nb;
}

static inline parsec_data_key_t
morse_parsec_data_key(parsec_data_collection_t *data_collection, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, data_collection);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

    return ((n * mdesc->lmt) + m);
}

static inline uint32_t
morse_parsec_rank_of(parsec_data_collection_t *data_collection, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, data_collection);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

    return mdesc->get_rankof( mdesc, m, n );
}

static inline uint32_t
morse_parsec_rank_of_key(parsec_data_collection_t *data_collection, parsec_data_key_t key)
{
    int m, n;
    morse_parsec_key_to_coordinates(data_collection, key, &m, &n);
    return morse_parsec_rank_of(data_collection, m, n);
}

static inline int32_t
morse_parsec_vpid_of(parsec_data_collection_t *data_collection, ... )
{
    (void)data_collection;
    return 0;
}

static inline int32_t
morse_parsec_vpid_of_key(parsec_data_collection_t *data_collection, parsec_data_key_t key)
{
    int m, n;
    morse_parsec_key_to_coordinates(data_collection, key, &m, &n);
    return morse_parsec_vpid_of(data_collection, m, n);
}

static inline parsec_data_t*
morse_parsec_data_of(parsec_data_collection_t *data_collection, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, data_collection);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

#if defined(CHAMELEON_USE_MPI)
    /* TODO: change displacement in data_map when in distributed */
    //assert( data_collection->nodes == 1 );
#endif
    return parsec_data_create( pdesc->data_map + n * mdesc->lmt + m, data_collection,
                               morse_parsec_data_key( data_collection, m, n ),
                               mdesc->get_blkaddr( mdesc, m, n ),
                               mdesc->bsiz * MORSE_Element_Size(mdesc->dtyp) );
}

static inline parsec_data_t*
morse_parsec_data_of_key(parsec_data_collection_t *data_collection, parsec_data_key_t key)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    int m, n;
    morse_parsec_key_to_coordinates(data_collection, key, &m, &n);

#if defined(CHAMELEON_USE_MPI)
    /* TODO: change displacement in data_map when in distributed */
    //assert( data_collection->nodes == 1 );
#endif
    return parsec_data_create( pdesc->data_map + key, data_collection, key,
                               mdesc->get_blkaddr( mdesc, m, n ),
                               mdesc->bsiz * MORSE_Element_Size(mdesc->dtyp) );
}

#ifdef parsec_PROF_TRACE
static inline int
morse_parsec_key_to_string(parsec_data_collection_t *data_collection, parsec_data_key_t key, char * buffer, uint32_t buffer_size)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)data_collection;
    MORSE_desc_t *mdesc = pdesc->desc;
    int m, n, res;
    morse_parsec_key_to_coordinates(data_collection, key, &m, &n);
    res = snprintf(buffer, buffer_size, "(%d, %d)", m, n);
    if (res < 0)
    {
        printf("error in key_to_string for tile (%u, %u) key: %u\n",
               (unsigned int)m, (unsigned int)n, datakey);
    }
    return res;
}
#endif

/**
 *  Create data descriptor
 */
void RUNTIME_desc_create( MORSE_desc_t *mdesc )
{
    parsec_data_collection_t *data_collection;
    morse_parsec_desc_t *pdesc;
    int comm_size;

    pdesc = malloc( sizeof(morse_parsec_desc_t) );
    data_collection = (parsec_data_collection_t*)pdesc;

    /* Super setup */
    comm_size = RUNTIME_comm_size( NULL );
    data_collection->nodes  = comm_size;
    data_collection->myrank = mdesc->myrank;

    data_collection->data_key    = morse_parsec_data_key;
    data_collection->rank_of     = morse_parsec_rank_of;
    data_collection->rank_of_key = morse_parsec_rank_of_key;
    data_collection->data_of     = morse_parsec_data_of;
    data_collection->data_of_key = morse_parsec_data_of_key;
    data_collection->vpid_of     = morse_parsec_vpid_of;
    data_collection->vpid_of_key = morse_parsec_vpid_of_key;
#if defined(parsec_PROF_TRACE)
    {
        int rc;
        data_collection->key_to_string = morse_parsec_key_to_string;
        data_collection->key           = NULL;
        rc = asprintf(&(data_collection->key_dim), "(%d, %d)", mdesc->lmt, mdesc->lnt);
        (void)rc;
    }
#endif
    data_collection->memory_registration_status = MEMORY_STATUS_UNREGISTERED;

    pdesc->data_map = calloc( mdesc->lmt * mdesc->lnt, sizeof(parsec_data_t*) );

    /* Double linking */
    pdesc->desc     = mdesc;
    mdesc->schedopt = pdesc;

    parsec_dtd_data_collection_init(data_collection);

    /* arena init */
    pdesc->arena_index = 0;

    /* taskpool init to bypass a requirement of PaRSEC  */
#if defined(CHAMELEON_USE_MPI)
    /* Look if an arena already exists for this descriptor */
    {
        morse_parsec_arena_t *arena = morse_parsec_registered_arenas;
        size_t size = mdesc->mb * mdesc->nb * MORSE_Element_Size(mdesc->dtyp);
        int i;

        for(i=0; i<morse_parsec_nb_arenas; i++, arena++) {
            if ( size == arena->size) {
                pdesc->arena_index = i;
                break;
            }
        }

        if (i == morse_parsec_nb_arenas) {
            parsec_datatype_t datatype;

            /* Create a taskpool to make sur the system is initialized */
            if ( i == 0 ) {
                parsec_taskpool_t *tp = parsec_dtd_taskpool_new();
                parsec_taskpool_free( tp );
            }

            /* Internal limitation of PaRSEC */
            assert(morse_parsec_nb_arenas < MORSE_PARSEC_DTD_NO_OF_ARENA);

            switch(mdesc->dtyp) {
            case MorseInteger:       datatype = parsec_datatype_int32_t; break;
            case MorseRealFloat:     datatype = parsec_datatype_float_t; break;
            case MorseRealDouble:    datatype = parsec_datatype_double_t; break;
            case MorseComplexFloat:  datatype = parsec_datatype_complex_t; break;
            case MorseComplexDouble: datatype = parsec_datatype_double_complex_t; break;
            default: morse_fatal_error("MORSE_Element_Size", "undefined type"); break;
            }

            /* Register the new arena */
            parsec_matrix_add2arena( parsec_dtd_arenas[i], datatype, matrix_UpperLower, 1,
                                     mdesc->mb, mdesc->nb, mdesc->mb, PARSEC_ARENA_ALIGNMENT_SSE, -1 );
            arena->size = size;
            pdesc->arena_index = i;
            morse_parsec_nb_arenas++;
        }
    }
#endif
    /* /\* Overwrite the leading dimensions to store the padding *\/ */
    /* mdesc->llm = mdesc->mb * mdesc->lmt; */
    /* mdesc->lln = mdesc->nb * mdesc->lnt; */
    return;
}

/**
 *  Destroy data descriptor
 */
void RUNTIME_desc_destroy( MORSE_desc_t *mdesc )
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)(mdesc->schedopt);
    if ( pdesc == NULL ) {
        return;
    }

    if ( pdesc->data_map != NULL ) {
        parsec_data_t **data = pdesc->data_map;
        int nb_local_tiles = mdesc->lmt * mdesc->lnt;
        int i;

        for(i=0; i<nb_local_tiles; i++, data++)
        {
            if (*data) {
                parsec_data_destroy( *data );
            }
        }

        free( pdesc->data_map );
        pdesc->data_map = NULL;
    }

    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)pdesc );

    free(pdesc);
    mdesc->schedopt = NULL;
    return;
}

/**
 *  Acquire data
 */
int RUNTIME_desc_acquire( const MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

/**
 *  Release data
 */
int RUNTIME_desc_release( const MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

/**
 *  Flush cached data
 */
void RUNTIME_flush()
{
}

void RUNTIME_desc_flush( const MORSE_desc_t     *desc,
                         const MORSE_sequence_t *sequence )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(sequence->schedopt);

    parsec_dtd_data_flush_all( PARSEC_dtd_taskpool, (parsec_data_collection_t*)(desc->schedopt) );
}

void RUNTIME_data_flush( const MORSE_sequence_t *sequence,
                         const MORSE_desc_t *A, int Am, int An )
{
    /*
     * For now, we do nothing in this function as in PaRSEC, once the data is
     * flushed it cannot be reused in the same sequence, when this issue will be
     * fixed, we will uncomment this function
     */
    /* parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(sequence->schedopt); */
    /* parsec_dtd_data_flush( PARSEC_dtd_taskpool, RTBLKADDR( A, MORSE_Complex64_t, Am, An ) ); */

    (void)sequence; (void)A; (void)Am; (void)An;
    return;
}

#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const MORSE_sequence_t *sequence,
                           const MORSE_desc_t *A, int Am, int An, int new_rank )
{
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
}
#endif

/**
 *  Get data addr
 */
void *RUNTIME_desc_getaddr( const MORSE_desc_t *desc, int m, int n )
{
    assert(0); /* This should not be called because we also need the handle to match the address we need. */
    return desc->get_blkaddr( desc, m, n );
}
