/**
 *
 * @file time_zpotri_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  MORSE_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MORSE_zpotri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRI( N ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRI( N ))

//#define POTRI_SYNC

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    MORSE_enum uplo = MorseLower;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, MORSE_Complex64_t, MorseComplexDouble, LDA, N, N );

    /*
     * Initialize Data
     * It's done in static to avoid having the same sequence than one
     * the function we want to trace
     */
    MORSE_zplghe_Tile( (double)N, MorseUpperLower, descA, 51 );

    /* MORSE ZPOTRF / ZTRTRI / ZLAUUM  */
    /*
     * Example of the different way to combine several asynchonous calls
     */
#if defined(TRACE_BY_SEQUENCE)
    {
        MORSE_sequence_t *sequence;
        MORSE_request_t request[3] = { MORSE_REQUEST_INITIALIZER,
                                       MORSE_REQUEST_INITIALIZER,
                                       MORSE_REQUEST_INITIALIZER };

        MORSE_Sequence_Create(&sequence);

        if ( ! iparam[IPARAM_ASYNC] ) {
            START_TIMING();

            MORSE_zpotrf_Tile_Async(uplo, descA,               sequence, &request[0]);
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);

            MORSE_ztrtri_Tile_Async(uplo, MorseNonUnit, descA, sequence, &request[1]);
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);

            MORSE_zlauum_Tile_Async(uplo, descA,               sequence, &request[2]);
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);

            STOP_TIMING();

        } else {

            START_TIMING();
            MORSE_zpotrf_Tile_Async(uplo, descA,               sequence, &request[0]);
            MORSE_ztrtri_Tile_Async(uplo, MorseNonUnit, descA, sequence, &request[1]);
            MORSE_zlauum_Tile_Async(uplo, descA,               sequence, &request[2]);

            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);
            STOP_TIMING();
        }

        MORSE_Sequence_Destroy(sequence[0]);
        MORSE_Sequence_Destroy(sequence[1]);
        MORSE_Sequence_Destroy(sequence[2]);
    }
#else
    {
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            MORSE_zpotrf_Tile(uplo, descA);
            MORSE_ztrtri_Tile(uplo, MorseNonUnit, descA);
            MORSE_zlauum_Tile(uplo, descA);
            STOP_TIMING();

        } else {

            /* Default: we use Asynchonous call with only one sequence */
            MORSE_sequence_t *sequence;
            MORSE_request_t request[2] = { MORSE_REQUEST_INITIALIZER,
                                           MORSE_REQUEST_INITIALIZER };

            START_TIMING();
            MORSE_Sequence_Create(&sequence);
            MORSE_zpotrf_Tile_Async(uplo, descA, sequence, &request[0]);
            MORSE_zpotri_Tile_Async(uplo, descA, sequence, &request[1]);
            MORSE_Desc_Flush( descA, sequence );
            MORSE_Sequence_Wait(sequence);
            STOP_TIMING();

            MORSE_Sequence_Destroy(sequence);
        }
    }
#endif

    /* Check the solution */
    if ( check )
    {
        dparam[IPARAM_ANORM] = 0.0;
        dparam[IPARAM_XNORM] = 0.0;
        dparam[IPARAM_BNORM] = 0.0;
        dparam[IPARAM_RES]   = 0.0;
    }

    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
