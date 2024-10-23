/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dshift.c
 *
 * StarPU codelet to Generate covariance matrix of a set of locations in 2D using Matern kernel or power-exp kernel.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-04-12
 *
 **/
#include "../include/starpu_exageostat.h"

static void cl_dshift_cpu_func(void *buffers[],void *cl_arg){
	int n,  n0;
	double *A;
        void core_dshift(double *,int,int);
	A	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

	starpu_codelet_unpack_args(cl_arg, &n, &n0);

	core_dshift(A, n, n0);
}


#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_dshift_cuda_func(void *buffers[], void *cl_arg)
{
	int n, n0;
	double *A       = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);


	starpu_codelet_unpack_args(cl_arg,  &n,
			&n0);

	cudaStream_t stream = starpu_cuda_get_local_stream();

	dshift_array(A, n, n0, stream);

	cudaStreamSynchronize( stream );
}
#endif


static struct starpu_codelet cl_dshift =
{
	.where		= STARPU_CPU,
	.cpu_func	= cl_dshift_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)

//	.where          = STARPU_CUDA,
	.cuda_func     =  cl_dshift_cuda_func,

#endif
	.nbuffers 	= 1,
	.modes		= STARPU_RW,
	.name		= "dshift"
};

//******************************************************************************
//****************************************************
int MORSE_MLE_dshift_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {

	MORSE_context_t *morse;
	MORSE_option_t options;
	morse = morse_context_self();


	if (sequence->status != MORSE_SUCCESS)
		return -2;
	RUNTIME_options_init(&options, morse, sequence, request);


	int  n, n0;
	int tempmm, tempnn;
	MORSE_desc_t A = *descA;
	struct starpu_codelet *cl=&cl_dshift;


	for (n = 0; n < A.nt; n++) {
		tempnn = n == A.nt -1 ? A.n - n * A.nb : A.nb;
		n0= n * A.nb;
		starpu_insert_task(starpu_mpi_codelet(cl),
				STARPU_VALUE, &tempnn, sizeof(int),
				STARPU_VALUE, &n0, sizeof(int),
				STARPU_RW, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, n, n),
				0);

	}

	//MORSE_TASK_flush_desc( &options, MorseUpperLower, descA );
	RUNTIME_options_ws_free(&options);
	RUNTIME_options_finalize(&options, morse);
	//MORSE_TASK_dataflush_all();
	return MORSE_SUCCESS;
}
