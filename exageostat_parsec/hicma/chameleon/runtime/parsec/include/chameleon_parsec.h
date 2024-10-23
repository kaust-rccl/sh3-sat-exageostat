/**
 *
 * @file chameleon_parsec.h
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC runtime header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Reazul Hoque
 * @date 2018-02-06
 *
 */
#ifndef _MORSE_PARSEC_H_
#define _MORSE_PARSEC_H_

#include <parsec.h>
#include <parsec/insert_function.h>
#include <parsec/data_dist/matrix/matrix.h>

/* Undefined PaRSEC definition of BLKLDD */
#undef BLKLDD

#include "control/common.h"

struct morse_parsec_desc_s {
    parsec_data_collection_t super;
    int                      arena_index;
    MORSE_desc_t            *desc;
    parsec_data_t          **data_map;
};

typedef struct morse_parsec_desc_s morse_parsec_desc_t;

static inline int
morse_parsec_get_arena_index(const MORSE_desc_t *desc) {
    return ((morse_parsec_desc_t *)desc->schedopt)->arena_index;
}

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( parsec_dtd_tile_of( (parsec_data_collection_t *) ((desc)->schedopt), \
                                                            ((parsec_data_collection_t *) (desc)->schedopt)->data_key((desc)->schedopt, m, n) ))

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _MORSE_PARSEC_H_ */
