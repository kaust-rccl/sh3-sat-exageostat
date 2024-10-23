/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dtrace.c
 *
 * StarPU codelet to calculate trace of a given matrix (A)
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-09-18
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_dtrace_starpu(void *buffers[],void *cl_arg){
    int m;
    int n;
    double *A;
    int m0;
    int n0;
    double s = 0;
    double *sum = &s;
    double *trace;

    *sum	 = 0;
    A		 = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
    sum      = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
    trace    = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n,  &m0, &n0);

    double local_s = core_dtrace(A, m, n, m0, n0, trace);
    *sum	+= local_s;
}

static struct starpu_codelet cl_dtrace =
{
    .where		= STARPU_CPU,
    .cpu_funcs	= {CORE_dtrace_starpu},
    .nbuffers	= 3,
    .modes		= {STARPU_R, STARPU_RW, STARPU_W},
    .name		= "dtrace"
};


int MORSE_MLE_dtrace_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence,
        MORSE_request_t  *request, MORSE_desc_t *descsum, MORSE_desc_t *desctrace) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0, n0;
    int tempmm;
    MORSE_desc_t A = *descA;
    struct starpu_codelet *cl=&cl_dtrace;


    for(m=0; m < A.mt; m++)
    {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                STARPU_VALUE, &tempmm,  sizeof(int),
                STARPU_VALUE, &tempmm, sizeof(int),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, m),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_VALUE, &n0,   sizeof(int),
                STARPU_RW, EXAGEOSTAT_RTBLKADDR(descsum, MorseRealDouble, 0, 0),
                STARPU_W, EXAGEOSTAT_RTBLKADDR(desctrace, MorseRealDouble, m, 0),
                0);
    }

    //MORSE_TASK_flush_desc( &options, MorseUpperLower, descA);
    //MORSE_TASK_flush_desc( &options, MorseUpperLower, descsum);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    //MORSE_TASK_flush_all();
    //MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}
