/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dscale.c
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

static void cl_dndscale_cpu_func(void *buffers[],void *cl_arg){
	int m, n, m0, n0;
	double *A;
	A	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
        double *B;
	B       = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
	double *C;
	C       = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);

	starpu_codelet_unpack_args(cl_arg, &m, &n, &m0,
			&n0);

	core_dndscale (A, B, C, m, n,
			m0, n0);

}


static void cl_ddscale_cpu_func(void *buffers[],void *cl_arg){
        int m, n, m0, n0;
        double *A;
        A       = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

        starpu_codelet_unpack_args(cl_arg, &m, &n, &m0,
                        &n0);

        core_ddscale (A, m, n,
                        m0, n0);

}


#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_dndscale_cuda_func(void *buffers[], void *cl_arg)
{
	int m, n, m0, n0;
	double *A;
	A       = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
        double *B;
        B       = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
	double *C;
	C       = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);	

	starpu_codelet_unpack_args(cl_arg, &m, &n, &m0,
			&n0);

	cudaStream_t stream = starpu_cuda_get_local_stream();

	dndscale_array(A, B, C, m, n, m0, n0, stream);

	cudaStreamSynchronize( stream );
}

static void cl_ddscale_cuda_func(void *buffers[], void *cl_arg)
{
        int m, n, m0, n0;
        double *A;
        A       = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

        starpu_codelet_unpack_args(cl_arg, &m, &n, &m0,
                        &n0);

        cudaStream_t stream = starpu_cuda_get_local_stream();

        ddscale_array(A,  m, n, m0, n0, stream);

        cudaStreamSynchronize( stream );
}
#endif


static struct starpu_codelet cl_dndscale =
{
	.where		= STARPU_CPU,
	.cpu_func	= cl_dndscale_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)

	.where          = STARPU_CUDA,
	.cuda_func     =  cl_dndscale_cuda_func,

#endif
	.nbuffers 	= 3,
	.modes		= STARPU_RW, STARPU_R, STARPU_R,
	.name		= "dndscale"
};

static struct starpu_codelet cl_ddscale =
{
        .where          = STARPU_CPU,
        .cpu_func       = cl_ddscale_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)

        .where          = STARPU_CUDA,
        .cuda_func     =  cl_ddscale_cuda_func,

#endif
        .nbuffers       = 1,
        .modes          = STARPU_RW,
        .name           = "ddscale"
};


//******************************************************************************
int MORSE_MLE_dscale_Tile_Async(MORSE_enum uplo, 
		MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {

	MORSE_context_t *morse;
	MORSE_option_t options;
	morse = morse_context_self();

	int tempnn, tempmm;

	if (sequence->status != MORSE_SUCCESS)
		return -2;
	RUNTIME_options_init(&options, morse, sequence, request);


	int m, n, m0, n0;
	MORSE_desc_t A = *descA;
	struct starpu_codelet *cl1=&cl_dndscale;
        struct starpu_codelet *cl2=&cl_ddscale;


	for (n = 0; n < A.nt; n++) {
		tempnn = n == A.nt -1 ? A.n - n * A.nb : A.nb;
		if(uplo == MorseUpperLower)
			m = 0;
		else
			m = A.m == A.n? n : 0;
		for(; m < A.mt; m++)
		{

			tempmm = m == A.mt -1 ? A.m- m* A.mb : A.mb;
			m0= m * A.mb;
			n0= n * A.nb;

			if(m != n)
			{
				starpu_insert_task(starpu_mpi_codelet(cl1),
						STARPU_VALUE, &tempmm, sizeof(int),
						STARPU_VALUE, &tempnn, sizeof(int),
						STARPU_VALUE, &m0, sizeof(int),
						STARPU_VALUE, &n0, sizeof(int),
						STARPU_RW, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, n),
						STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, m),
						STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, n, n),
						0);

			}
			else
			{
				starpu_insert_task(starpu_mpi_codelet(cl2),
						STARPU_VALUE, &tempmm, sizeof(int),
						STARPU_VALUE, &tempnn, sizeof(int),
						STARPU_VALUE, &m0, sizeof(int),
						STARPU_VALUE, &n0, sizeof(int),
						STARPU_RW, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, n),
						0);

			}
		}

	}

	//MORSE_TASK_flush_desc( &options, MorseUpperLower, descA );
	RUNTIME_options_ws_free(&options);
	RUNTIME_options_finalize(&options, morse);
	//MORSE_TASK_dataflush_all();
	return MORSE_SUCCESS;
}

