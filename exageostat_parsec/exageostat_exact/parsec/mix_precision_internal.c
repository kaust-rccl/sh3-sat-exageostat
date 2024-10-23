#include "mix_precision_internal.h"

#include "../include/exageostatcore.h"
#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20

#define NBELEM   1


/* Global to record number of DP10SP90 and DP100 when losing SPD */
int num_iteration_scale_1 = 0;
int num_iteration_scale_2 = 0;
int num_iteration_scale_3 = 0;

static unsigned long long int
my_Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
        unsigned long long int a_k, c_k, ran;
        int i;

        a_k = Rnd64_A;
        c_k = Rnd64_C;

        ran = seed;
        for (i = 0; n; n >>= 1, i++) { 
                if (n & 1)
                        ran = a_k * ran + c_k;
                c_k *= (a_k + 1); 
                a_k *= a_k;
        }

        return ran;
}

void parsec_core_dshift (double *A, int n, int n0, int multiple)
{
        double acc=1e-7;
        unsigned long long int ran, jump; 
        unsigned long long int seed=5956+n0; 
        float val=n0*acc;
    
        double tmp;
        int i, j;
    
        ran = my_Rnd64_jump( NBELEM , seed );
        for (i = 0; i < n; i++) {
                        tmp = (ran * RndF_Mul)/1000000000;
                        tmp=tmp+(val*i/100);
                        A[i + i * n] = A[i + i * n] + tmp * multiple; /*+ 1e-4*/;
        }
}



void convert_s2d_unary_CPU(float *data, int mb, int nb){
    double *data_d = (double *)data;
    for( int j = nb-1; j >= 0; j-- )
        for( int i = mb-1; i >= 0; i-- )
            data_d[j*mb+i] = (double)data[j*mb+i];
}

void convert_d2s_unary_CPU(double *data, int mb, int nb){
    float *data_s = (float *)data;
    for( int j = 0; j < nb; j++ )
        for( int i = 0; i < mb; i++ )
            data_s[j*mb+i] = (float)data[j*mb+i];
}

void convert_s2d_binary_CPU(double *_target, float *_source, int mb, int nb) {
    double *target = (double *)_target;
    float *source = (float *)_source;
    for( int j = 0; j < nb; j++ )
        for( int i = 0; i < mb; i++ )
            target[j*mb+i] = (double)source[j*mb+i];
}

void convert_d2s_binary_CPU(float *_target, double *_source, int mb, int nb) {
    float *target = (float *)_target;
    double *source = (double *)_source;
    for( int j = 0; j < nb; j++ )
        for( int i = 0; i < mb; i++ )
            target[j*mb+i] = (float)source[j*mb+i];
}

void disable_CPU( parsec_taskpool_t * tp ) {
#if defined(EXAGEOSTAT_USE_CUDA)
        for (int i = 0; i < parsec_nb_devices; i++) {
                parsec_device_module_t *device = parsec_mca_device_get(i);
                if( device->type == PARSEC_DEV_CPU )
                        tp->devices_index_mask &= ~(1<<i);
        }
#endif
}

void disable_GPU( parsec_taskpool_t * tp ) {
#if defined(EXAGEOSTAT_USE_CUDA)
	for (int i = 0; i < parsec_nb_devices; i++) {
		parsec_device_module_t *device = parsec_mca_device_get(i);
		if( device->type == PARSEC_DEV_CUDA )
			tp->devices_index_mask &= ~(1<<i);
	}
#endif
}

#if defined(EXAGEOSTAT_USE_CUDA)
/* Lookup GPU workspace */
parsec_potrf_stream_workspace_t lookup_gpu_workspace( parsec_device_cuda_module_t *gpu_device,
                                                      parsec_gpu_exec_stream_t *gpu_stream,
                                                      parsec_potrf_workspace_t *ws ) {
    int i, j;

    /* Look for device */
    for(i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA ) continue;
        parsec_device_cuda_module_t *gpu_device_compare = (parsec_device_cuda_module_t*)device;

        if(gpu_device->cuda_index == gpu_device_compare->cuda_index)
            break;
    }

    /* Look for stream; 0, h2d; 1 d2h*/
    for(j = 2; j < gpu_device->max_exec_streams; j++) {
        if( gpu_stream == &(gpu_device->exec_stream[j]) )
            break;
    }

    return ws->gpu_workspace[i].stream_workspace[j];
}

/* Allocate memory for workspace */
parsec_potrf_workspace_t* workspace_memory_allocate( parsec_potrf_workspace_t *ws ) {
    ws = (parsec_potrf_workspace_t *)malloc( sizeof(parsec_potrf_workspace_t) );
    ws->gpu_workspace = (parsec_potrf_gpu_workspace_t *)malloc( parsec_nb_devices * sizeof(parsec_potrf_gpu_workspace_t) );

    for( int i = 0; i < parsec_nb_devices; i++ ) {
        ws->gpu_workspace[i].stream_workspace = (parsec_potrf_stream_workspace_t *)malloc( PARSEC_MAX_STREAMS * sizeof(parsec_potrf_stream_workspace_t) );
        ws->gpu_workspace[i].gpu_device = (parsec_device_cuda_module_t *)malloc( sizeof(parsec_device_cuda_module_t) );
    }

    return ws;
}

/* Destruct and memory free */ 
void workspace_memory_free( parsec_potrf_workspace_t *ws)
{
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA ) continue;

        for(int j = 0; j < ws->gpu_workspace[i].gpu_device->max_exec_streams; j++) {
            /* j 0, h2d; j 1, d2h */
            if( j <= 1 ) continue;

            /* Free GPU buffer */
            if( NULL != ws->gpu_workspace[i].stream_workspace[j].gpu_buffer ) {
                zone_free( ws->gpu_workspace[i].gpu_device->memory, ws->gpu_workspace[i].stream_workspace[j].gpu_buffer );
                //cudaFree( ws->gpu_workspace[i].stream_workspace[j].gpu_buffer );
            }

            /* Free GPU handle_cusolver */
            if( NULL != ws->gpu_workspace[i].stream_workspace[j].handle_cusolver ) {
                cusolverDnHandle_t handle_cusolver = ws->gpu_workspace[i].stream_workspace[j].handle_cusolver;
                cusolverStatus_t status = cusolverDnDestroy(handle_cusolver);
                assert(status == CUSOLVER_STATUS_SUCCESS);
            }

            /* Free GPU handle_cublas */
            if( NULL != ws->gpu_workspace[i].stream_workspace[j].handle_cublas ) {
                cublasHandle_t handle_cublas = ws->gpu_workspace[i].stream_workspace[j].handle_cublas;
                cublasStatus_t status = cublasDestroy(handle_cublas);
                assert(status == CUBLAS_STATUS_SUCCESS);
            }
        }
    }

    for( int i = 0; i < parsec_nb_devices; i++ ) {
        free( ws->gpu_workspace[i].stream_workspace );
    }

    free( ws->gpu_workspace );
    free( ws );
}

#if GPU_BUFFER_ONCE
#if 0
/* Declare workspace used on GPU */
extern parsec_potrf_workspace_t *ws_handle_cusolver;
extern parsec_potrf_workspace_t *ws_double1;
extern parsec_potrf_workspace_t *ws_double2;
extern parsec_potrf_workspace_t *ws_single1;
extern parsec_potrf_workspace_t *ws_single2;
extern parsec_potrf_workspace_t *ws_handle_cublas_tensor;
extern parsec_potrf_workspace_t *ws_handle_cublas_gpu;
extern parsec_potrf_workspace_t *ws_half1;
extern parsec_potrf_workspace_t *ws_half2;
extern parsec_potrf_workspace_t *ws_half3;
/* For matrix generation */
extern parsec_potrf_workspace_t *ws_l1_x;
extern parsec_potrf_workspace_t *ws_l1_y;
extern parsec_potrf_workspace_t *ws_l1_z;
extern parsec_potrf_workspace_t *ws_l2_x;
extern parsec_potrf_workspace_t *ws_l2_y;
extern parsec_potrf_workspace_t *ws_l2_z;
#endif

/* Declare workspace used on GPU */
parsec_potrf_workspace_t *ws_handle_cusolver, *ws_double1, *ws_double2, *ws_single1, *ws_single2;
parsec_potrf_workspace_t *ws_handle_cublas_tensor, *ws_handle_cublas_gpu, *ws_half1, *ws_half2, *ws_half3;
parsec_potrf_workspace_t *ws_l1_x, *ws_l1_y, *ws_l1_z, *ws_l2_x, *ws_l2_y, *ws_l2_z;

void gpu_temporay_buffer_init( int mb, int nb ) {
    /* Allocate memory */
    ws_handle_cusolver = workspace_memory_allocate( ws_handle_cusolver );
    ws_handle_cublas_tensor = workspace_memory_allocate( ws_handle_cublas_tensor );
    ws_handle_cublas_gpu = workspace_memory_allocate( ws_handle_cublas_gpu);
    ws_double1 = workspace_memory_allocate( ws_double1 );
    ws_double2 = workspace_memory_allocate( ws_double2 );
    ws_single1 = workspace_memory_allocate( ws_single1 );
    ws_single2 = workspace_memory_allocate( ws_single2 );
    ws_half1 = workspace_memory_allocate( ws_half1 );
    ws_half2 = workspace_memory_allocate( ws_half2 );
    ws_half3 = workspace_memory_allocate( ws_half3 );

    ws_l1_x = workspace_memory_allocate( ws_l1_x );
    ws_l1_y = workspace_memory_allocate( ws_l1_y );
    ws_l1_z = workspace_memory_allocate( ws_l1_z );
    ws_l2_x = workspace_memory_allocate( ws_l2_x );
    ws_l2_y = workspace_memory_allocate( ws_l2_y );
    ws_l2_z = workspace_memory_allocate( ws_l2_z );

    int size = mb * nb;

    /* Traverse all gpu device */
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA ) continue;

        parsec_device_cuda_module_t *gpu_device = (parsec_device_cuda_module_t*)device;
        cudaSetDevice(gpu_device->cuda_index);

        ws_handle_cusolver->gpu_workspace[i].gpu_device = gpu_device;
        ws_handle_cublas_tensor->gpu_workspace[i].gpu_device = gpu_device;
        ws_handle_cublas_gpu->gpu_workspace[i].gpu_device = gpu_device;
        ws_double1->gpu_workspace[i].gpu_device = gpu_device;
        ws_double2->gpu_workspace[i].gpu_device = gpu_device;
        ws_single1->gpu_workspace[i].gpu_device = gpu_device;
        ws_single2->gpu_workspace[i].gpu_device = gpu_device;
        ws_half1->gpu_workspace[i].gpu_device = gpu_device;
        ws_half2->gpu_workspace[i].gpu_device = gpu_device;
        ws_half3->gpu_workspace[i].gpu_device = gpu_device;

        ws_l1_x->gpu_workspace[i].gpu_device = gpu_device;
        ws_l1_y->gpu_workspace[i].gpu_device = gpu_device;
        ws_l1_z->gpu_workspace[i].gpu_device = gpu_device;
        ws_l2_x->gpu_workspace[i].gpu_device = gpu_device;
        ws_l2_y->gpu_workspace[i].gpu_device = gpu_device;
        ws_l2_z->gpu_workspace[i].gpu_device = gpu_device;

        /* Traverse all streams */
        for(int j = 0; j < gpu_device->max_exec_streams; j++) {
            /* j 0, h2d; j 1, d2h */
            if( j <= 1 ) continue;

            cublasStatus_t status;
            cusolverStatus_t status_cusolver;
            cudaError_t cudaStatus;
            cusolverDnHandle_t handle_cusolver;
            cublasHandle_t handle_cublas_tensor;
            cublasHandle_t handle_cublas_gpu;

            /* Set unused to NULL */
            {
                ws_handle_cusolver->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_handle_cublas_tensor->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_handle_cublas_tensor->gpu_workspace[i].stream_workspace[j].gpu_buffer = NULL;

                ws_handle_cublas_gpu->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_handle_cublas_gpu->gpu_workspace[i].stream_workspace[j].gpu_buffer = NULL;

                ws_double1->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_double1->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;
                ws_double2->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_double2->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_single1->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_single1->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;
                ws_single2->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_single2->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_half1->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_half1->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;
                ws_half2->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_half2->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;
                ws_half3->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_half3->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l1_x->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l1_x->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l1_y->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l1_y->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l1_z->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l1_z->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l2_x->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l2_x->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l2_y->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l2_y->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;

                ws_l2_z->gpu_workspace[i].stream_workspace[j].handle_cusolver = NULL;
                ws_l2_z->gpu_workspace[i].stream_workspace[j].handle_cublas = NULL;
            }

            /* Create handle_cusolver */
            {
                status_cusolver = cusolverDnCreate(&handle_cusolver);
                assert(CUSOLVER_STATUS_SUCCESS == status_cusolver);
                status_cusolver = cusolverDnSetStream(handle_cusolver, gpu_device->exec_stream[j].cuda_stream);
                assert(CUSOLVER_STATUS_SUCCESS == status_cusolver);
                ws_handle_cusolver->gpu_workspace[i].stream_workspace[j].handle_cusolver = handle_cusolver;
            }

            /* Allocate workspace for potrf handle_cusolver */
            {
                int workspace_size;
                status_cusolver = cusolverDnDpotrf_bufferSize(handle_cusolver, CUBLAS_FILL_MODE_LOWER, nb, NULL, mb, &workspace_size);
                assert(CUSOLVER_STATUS_SUCCESS == status_cusolver);
                ws_handle_cusolver->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) + sizeof(int) );
                assert(NULL != ws_handle_cusolver->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_handle_cusolver->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Create handle_cublas_tensor */
            {
                status = cublasCreate(&handle_cublas_tensor); 
                assert(CUBLAS_STATUS_SUCCESS == status);
                /* Set the math mode to allow cuBLAS to use Tensor Cores */
                status = cublasSetMathMode(handle_cublas_tensor, CUBLAS_TENSOR_OP_MATH);
                assert(CUBLAS_STATUS_SUCCESS == status);
                ws_handle_cublas_tensor->gpu_workspace[i].stream_workspace[j].handle_cublas = handle_cublas_tensor;
            }

            /* Create handle_cublas_gpu */
            {
                status = cublasCreate(&handle_cublas_gpu);
                assert(CUBLAS_STATUS_SUCCESS == status);
                ws_handle_cublas_gpu->gpu_workspace[i].stream_workspace[j].handle_cublas = handle_cublas_gpu;
            }

            /* Temporary buffer for double */
            {
                int workspace_size = size;
                ws_double1->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_double1->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_double1->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_double2->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_double2->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_double2->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for single */
            {
                int workspace_size = size;
                ws_single1->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) );
                assert(NULL != ws_single1->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_single1->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_single2->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) );
                assert(NULL != ws_single2->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_single2->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for half */
            {
                int workspace_size = size;
                ws_half1->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) / 2 );
                assert(NULL != ws_half1->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_half1->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_half2->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) / 2 );
                assert(NULL != ws_half2->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_half2->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_half3->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) / 2 );
                assert(NULL != ws_half3->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_half3->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

            }

            /* Temporary buffer for l1 */
            {
                int workspace_size = nb;
                ws_l1_x->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l1_x->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l1_x->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_l1_y->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l1_y->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l1_y->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_l1_z->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l1_z->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l1_z->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for l2 */
            {
                int workspace_size = mb;
                ws_l2_x->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l2_x->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l2_x->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_l2_y->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l2_y->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l2_y->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_l2_z->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_l2_z->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_l2_z->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

        }
    }

}

void gpu_temporay_buffer_fini( ) {
    /** Find all CUDA devices */
    int nb = 0;
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA == device->type ) {
            nb++;
        }
    }

    if( nb > 0 ) {
        workspace_memory_free( ws_handle_cusolver );
        workspace_memory_free( ws_handle_cublas_tensor );
        workspace_memory_free( ws_handle_cublas_gpu );
        workspace_memory_free( ws_double1 );
        workspace_memory_free( ws_double2 );
        workspace_memory_free( ws_single1 );
        workspace_memory_free( ws_single2 );
        workspace_memory_free( ws_half1 );
        workspace_memory_free( ws_half2 );
        workspace_memory_free( ws_half3 );

        workspace_memory_free( ws_l1_x );
        workspace_memory_free( ws_l1_y );
        workspace_memory_free( ws_l1_z );
        workspace_memory_free( ws_l2_x );
        workspace_memory_free( ws_l2_y );
        workspace_memory_free( ws_l2_z );
    }

}
#endif /* GPU_BUFFER_ONCE */

#endif

#if SWITCH_TO_DP
        /* parameters to swith back to dp for mix-precision */
        int switch_to_dp = 0;
        double loglik_pre = 0.0;
#endif


double parsec_dmle_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data) {
        double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, time_mle = 0.0;
	double logdet_calculate = 0.0, matrix_gen_time=0.0, dzcpy_time=0.0;
        int N, NRHS, success;   
        double flops =0.0;   
	double sum = 0.0;
                                
        MLE_data* data  = ((MLE_data*)MORSE_data);
        data->det       = 0;
        data->dotp      = 0;
        int nodes = data->nodes;
	//double threshold = pow(10, -1.0 * data->opt_tol); 
	double threshold = pow(10, -1.0 * 16); 

#if SWITCH_TO_DP
	/* parameters to swith back to dp for mix-precision */
	double threshold_to_dp = pow(10, -1.0 * data->opt_tol_mix); 
#endif

        parsec_tiled_matrix_dc_t *descC       = (parsec_tiled_matrix_dc_t *) data->descC;
        parsec_tiled_matrix_dc_t *descZ       = (parsec_tiled_matrix_dc_t *) data->descZ;
        parsec_tiled_matrix_dc_t *descZcpy    = (parsec_tiled_matrix_dc_t *) data->descZcpy;
        parsec_tiled_matrix_dc_t *descdet     = (parsec_tiled_matrix_dc_t *) data->descdet;
        parsec_tiled_matrix_dc_t *descproduct = (parsec_tiled_matrix_dc_t *) data->descproduct;
	parsec_context_t *parsec = (parsec_context_t *)data->parsec;
	PLASMA_enum uplo = (PLASMA_enum)data->uplo;
        int rank = data->rank;
	int HNB = data->HNB;

        N       = descC->m;       
        NRHS    = descZ->n;
        int num_params;
        if(strcmp (data->c_fun, "pow-exp-nuggets") == 0)
            num_params=4;
        else
            num_params=3;
	if( 0 == rank ) {
		if(num_params==3) 
			fprintf(stderr, "Optimization Start: %3d- Model Parameters (theta0, theta1, theta2): (%.10lf, %.10lf, %.10lf)\n", data->iter_count,  theta[0], theta[1], theta[2]);
		else
			fprintf(stderr, "Optimization Start: %3d- Model Parameters (theta0, theta1, theta2, theta3): (%.10lf, %.10lf, %.10lf, %.10lf)\n", data->iter_count,  theta[0], theta[1], theta[2], theta[3]);
	}
	START_TIMING(dzcpy_time);
	if(data->iter_count==0) {
		VERBOSE("Save copy of the original Z vector...");
		/* Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration) */
		dplasma_dlacpy(parsec, PlasmaUpperLower, descZ, descZcpy);
		VERBOSE(" Done.\n");
	}
	STOP_TIMING(dzcpy_time);

#ifndef EXAGEOSTAT_ON_SUMMIT
	/* Checking point and recovery */
	if(strcmp(data->recovery_file,"") != 0 && recover(data->recovery_file, data->iter_count, theta, &loglik, 3)) {
		if(0 == rank) fprintf(stderr, "\nIteration %d is recovered\n", data->iter_count);
		data->variance = theta[0];
	} else
#endif
	{
		START_TIMING(dzcpy_time);
		if(data->iter_count > 0) {
			VERBOSE("Re-store the original Z vector...");
			dplasma_dlacpy(parsec, PlasmaUpperLower, descZcpy, descZ);
			VERBOSE(" Done.\n");
		}
		STOP_TIMING(dzcpy_time);

#if DEBUG_INFO
		VERBOSE("Get sum of C...");
		sum = parsec_dmatrix_sum(parsec, descC);
		if( 0 == rank ) fprintf( stderr, "\nC: sum before generation %.17g\n", sum );

		VERBOSE("Get sum of Z...");
		sum = parsec_dZ_sum(parsec, descZ);
		if( 0 == rank ) fprintf( stderr, "\nZ: sum before generation %.17g\n", sum );
#endif

#if defined(EXAGEOSTAT_USE_CUDA) 
                /* Reset cache on GPU before Cholesky */
                parsec_devices_release_memory();
                parsec_devices_reset_load(parsec);
#endif

		VERBOSE("Generate New Covariance Matrix...");
		START_TIMING(matrix_gen_time);  
		parsec_dmatrix_generation( parsec, descC, &data->l1, &data->l1,
				(double *)theta, data->dm, data->c_fun, data->band_size_double );

		STOP_TIMING(matrix_gen_time);
		VERBOSE(" Done.\n");

#if COUNT_VALUE_HALF 
                VERBOSE("Count values exceeding half-precision...");
		// Average of error_ratio per tile
		if( 0 == data->iter_count )
			parsec_dmatrix_ratio_error_half( parsec, descC, theta, data->iter_count );

		double count_value_time = 0.0;
                START_TIMING(count_value_time);
                parsec_dmatrix_ratio_exceed_half( parsec, descC, theta, data->iter_count );
                STOP_TIMING(count_value_time);
                VERBOSE(" Done.\n");

#if defined(EXAGEOSTAT_USE_CUDA) 
                        /* Reset cache on GPU before Cholesky */
                        parsec_devices_release_memory();
                        parsec_devices_reset_load(parsec);
#endif

#endif

		VERBOSE("Cholesky factorization of Sigma...");
		START_TIMING(time_facto);

#if DEBUG_INFO
		sum = parsec_dmatrix_sum(parsec, descC);
		if( 0 == rank ) fprintf( stderr, "\nC: sum before cholesky %.17g\n", sum );
		//dplasma_dprint(parsec, uplo, descC);

                if( 0 && data->iter_count == 5 ) {
                        parsec_dmatrix_avg( parsec, descC );
                        //exit(1);
                }

#endif
		//fprintf(stderr, "After generation:\n");
		//if( data->iter_count >= 28 )
		//parsec_dmatrix_avg( parsec, descC );

#if SINGLE_PO
		if( 0 == rank ) fprintf(stderr, "Single in optinization %d\n", data->iter_count);
                /* Convert d2s */
                parsec_band_convert_d2s(parsec, descC, 0);

                /* Call sdpotrf */
		success = parsec_spotrf( parsec, uplo, descC );

                /* Conver s2d */
                parsec_band_convert_s2d(parsec, descC, 0);
#else
		/* Define band_size */
		if( 0 == data->precision ) {
			int NT = descC->lmt ;
			data->band_size_double = NT; 
			data->band_size_single = NT; 
		}

#if SWITCH_TO_DP
		/* Check whether to switch to DP100% */
		else if( switch_to_dp ) { 
			int NT = descC->lmt ;
			data->band_size_double = NT; 
			data->band_size_single = NT; 
			if( 0 == rank && 0 != data->precision ) fprintf(stderr, "SWITCH_TO_DP at iteration %d\n", data->iter_count);
		}
#endif

		/* Convert d2s */
		parsec_band_convert_d2s(parsec, descC, data->band_size_double);

#if TRSM_CONVERT_HALF
		/* Call sdpotrf */
		success = dplasma_hsdpotrf_trsm_convert( parsec, uplo, descC, data->lookahead, data->sendless,
				data->tensor_gemm, data->band_size_double, data->band_size_single, HNB );
#else
		/* Call sdpotrf */
		success = dplasma_hsdpotrf( parsec, uplo, descC, data->lookahead, data->sendless,
				data->tensor_gemm, data->band_size_double, data->band_size_single, HNB );
#endif

		/* Conver s2d */
		parsec_band_convert_s2d(parsec, descC, data->band_size_double);
#endif /* SINGLE_PO */

                flops = flops + FLOPS_DPOTRF(N);
                VERBOSE(" Done.\n");
                STOP_TIMING(time_facto);

		if( 0 != success ) {
			if( 0 == rank ) {
				fprintf(stderr, "\n\nFactorization cannot be performed: the matrix is not positive definite: %d\n", success);
			}

			if( 0 == data->precision ) {
				if( 0 == rank ) fprintf(stderr, "ERROR: double-precision lose SPD at (%.10lf, %.10lf, %.10lf)\n",
						theta[0], theta[1], theta[2]);
				exit(1);
			}
		}

#if ADAPTIVE_PO
		static int iteration_scale = 0;
		for(iteration_scale = 1; 0 != success; iteration_scale++ ) {
			if( 0 == rank ) {
				fprintf(stderr, "Apply matrix scale to retain s.p.d.\n");
				fprintf(stderr, "**************** iteration_scale %d ****************\n", iteration_scale);
			}

                        if( iteration_scale > 11 ) {
                                if( 0 == rank ) fprintf(stderr, "Not recovery s.p.d after %d iterations, exit\n", iteration_scale);
                                exit(1);
                        }

			/* Only support pow-exp */
			assert( strcmp(data->c_fun, "pow-exp") == 0 );

			/* Re-set success to 0 */
			success = 0;

			VERBOSE("Generate New Covariance Matrix (Matrix_scale) ...");
			parsec_dmatrix_generation( parsec, descC, &data->l1, &data->l1,
					(double *)theta, data->dm, data->c_fun, data->band_size_double );
			VERBOSE(" Done.\n");

#if 0 
			/* First time NPD and first time DP100% loses SPD */
			if( 0 & (1 == iteration_scale || iteration_scale > 3) )
			{
				VERBOSE("Scale Matrix to retain s.p.d (Matrix_scale) ....");
				//parsec_dmatrix_scale( parsec, descC );
				VERBOSE(" Done.\n");

				VERBOSE("Shfit Matrix to retain s.p.d (Matrix_scale)\n");

				int multiple = 1;
				for(int ii = 1; ii <= iteration_scale; ii++) {
					if( 0 == rank ) fprintf(stderr, "************** %.d th shift ************** \n", multiple);
					parsec_dmatrix_shift( parsec, descC, multiple );
					multiple *= 2;
				}
				VERBOSE("Done Shfit Matrix.\n");
			}
#endif

			//fprintf(stderr, "After shifting:\n");
			//parsec_dmatrix_avg( parsec, descC );

#if defined(EXAGEOSTAT_USE_CUDA) 
			/* Reset cache on GPU before Cholesky */
			parsec_devices_release_memory();
			parsec_devices_reset_load(parsec);
#endif

			VERBOSE("Cholesky factorization of Sigma (Matrix_scale) ...");

			/* Define band_size */
			int NT = descC->lmt ;
			int my_band_size_double = data->band_size_double;
			int my_band_size_single = data->band_size_single;

			if( 0 == data->precision ) {
				my_band_size_double = NT;
				my_band_size_single = NT;
			}

			else {
				if( 1 == iteration_scale ) {
					num_iteration_scale_1 += 1;
					double noise = 1.0e-4;
					if( 0 == rank ) fprintf(stderr, " %d iteration add %lf to diagonal : %d %d\n",
							iteration_scale, noise, my_band_size_double, my_band_size_single);
					parsec_dmatrix_set_diagonal(parsec, descC, noise);

				} else if( 2 == iteration_scale ) {
					num_iteration_scale_2 += 1;
					double noise = 1.0e-4;
					if( 0 == rank ) fprintf(stderr, " %d iteration add %lf to diagonal : %d %d\n",
							iteration_scale, noise, my_band_size_double, my_band_size_single);
					parsec_dmatrix_set_diagonal(parsec, descC, noise);

					my_band_size_double = NT/10;
					my_band_size_single = NT;
					if( 0 == rank ) fprintf(stderr, " %d iteration do DP10 SP90 instead : %d %d\n",
							iteration_scale, my_band_size_double, my_band_size_single);

				} else if( 3 == iteration_scale ) {
					num_iteration_scale_3 += 1;
					my_band_size_double = NT;
					my_band_size_single = NT;
					if( 0 == rank ) fprintf(stderr, " %d iteration do DP100 instead : %d %d\n",
							iteration_scale, my_band_size_double, my_band_size_single);

				} else { 
					if( 0 == rank ) fprintf(stderr, "ERROR: double lose SPD at (%.10lf, %.10lf, %.10lf)\n",
					       theta[0], theta[1], theta[2]);
					exit(1);
				}
			}

			/* Convert d2s */
			parsec_band_convert_d2s(parsec, descC, my_band_size_double);

#if TRSM_CONVERT_HALF
			/* Call sdpotrf */
			success = dplasma_hsdpotrf_trsm_convert( parsec, uplo, descC, data->lookahead, data->sendless,
					data->tensor_gemm, my_band_size_double, my_band_size_single, HNB );
#else
			/* Call sdpotrf */
			success = dplasma_hsdpotrf( parsec, uplo, descC, data->lookahead, data->sendless,
					data->tensor_gemm, my_band_size_double, my_band_size_single, HNB );
#endif

			/* Conver s2d */
			parsec_band_convert_s2d(parsec, descC, my_band_size_double);

			VERBOSE(" Done.\n");
		}
#endif /* if 0 */


		//if( 0 == rank ) fprintf(stderr, "\n\n");

#if DEBUG_INFO
		sum = parsec_dmatrix_sum(parsec, descC);
		if( 0 == rank ) fprintf( stderr, "\nC: sum after cholesky %.17g\n", sum );

		//fprintf(stderr, "lookahead %d band_size %d sendless %d\n", data->lookahead, data->band_size_double, data->sendless);
		//dplasma_dprint(parsec, uplo, descC);
#endif

		VERBOSE("Calculating the log determinant ...");
		START_TIMING(logdet_calculate);
		parsec_dmatrix_det( parsec, descC, descdet );
		/* Broadcast determinant */
		int root = descdet->super.rank_of(&descdet->super, 0, 0);
		if( descdet->super.myrank == root )
			data->det = *( (double *)((descdet->super.data_of(&descdet->super, 0, 0))->device_copies[0]->device_private) );
		MPI_Bcast( (void *)&data->det, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 
		logdet = 2 * data->det;
		STOP_TIMING(logdet_calculate);
		VERBOSE(" Done.\n");

#if DEBUG_INFO
		sum = parsec_dZ_sum(parsec, descZ);
		if( 0 == rank ) fprintf( stderr, "\nZ: sum before dtrsm %.17g\n", sum );
		//dplasma_dprint(parsec, uplo, descZ);
#endif

		VERBOSE("Solving the linear system ...");
		START_TIMING(time_solve);
		parsec_taskpool_t *parsec_dtrsm = NULL;
		parsec_dtrsm = my_dtrsm_New( PlasmaLeft,  uplo, PlasmaNoTrans, PlasmaNonUnit, 1.0, descC, descZ );
#if NO_GPU_DPLASMA 
		disable_GPU( parsec_dtrsm );
#endif
		if ( parsec_dtrsm != NULL )
		{
			parsec_context_add_taskpool( parsec, parsec_dtrsm );
			parsec_context_start( parsec );
			parsec_context_wait( parsec );
			my_dtrsm_Destruct( parsec_dtrsm );
		}
		STOP_TIMING(time_solve);
		flops = flops + FLOPS_DTRSM(PlasmaLeft, N, NRHS);
		VERBOSE(" Done.\n"); 

#if DEBUG_INFO
		sum = parsec_dZ_sum(parsec, descZ);
		if( 0 == rank ) fprintf( stderr, "\nZ: sum after dtrsm %.17g\n", sum );
		//dplasma_dprint(parsec, uplo, descZ);
#endif

#if OWN_DOT 
		parsec_dmatrix_dot( parsec, descZ, descZ, descproduct); 
#else
		VERBOSE("Calculating the MLE likelihood function ...");
		START_TIMING(time_mle);
		two_dim_block_cyclic_t descZ0;
		two_dim_block_cyclic_init(&descZ0, matrix_RealDouble, matrix_Tile,
				descZ->super.nodes, descZ->super.myrank, descZ->mb, descZ->nb, descZ->m, 1, 0, 0,
				descZ->m, 1, 1, 1, nodes);
		descZ0.mat = calloc((size_t)descZ0.super.nb_local_tiles *
				(size_t)descZ0.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descZ0.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&descZ0, "descZ0");

		dplasma_dlacpy( parsec, matrix_UpperLower, descZ, (parsec_tiled_matrix_dc_t *)&descZ0 );

		parsec_taskpool_t *parsec_dgemm = NULL;
		parsec_dgemm = dplasma_dgemm_New( PlasmaTrans, PlasmaNoTrans, 1.0, descZ, (parsec_tiled_matrix_dc_t *)&descZ0, 0.0, descproduct );
#if NO_GPU_DPLASMA 
		disable_GPU( parsec_dgemm ); 
#endif
		if ( parsec_dgemm != NULL )
		{
			parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_dgemm);
			parsec_context_start( parsec );
			parsec_context_wait( parsec );
			dplasma_dgemm_Destruct( parsec_dgemm );
		}

		parsec_data_free( descZ0.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&descZ0 );
#endif
		/* Broadcast dot product */
		root = descproduct->super.rank_of(&descproduct->super, 0, 0);
		if( descproduct->super.myrank == root )
			data->dotp = *( (double *)((descproduct->super.data_of(&descproduct->super, 0, 0))->device_copies[0]->device_private) );
		MPI_Bcast( (void *)&data->dotp, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); 

		if(strcmp(data->c_fun, "matern") == 0)
		{
			loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
			data->variance = (1.0/N) * data->dotp;
		}
		else
		{
			loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
			data->variance= theta[0];
		}
		STOP_TIMING(time_mle);
		VERBOSE(" Done.\n");

		/* Distribute the values in the case of MPI */
#if defined(PARSEC_HAVE_MPI)
		MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast((void *)theta, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

#if SWITCH_TO_DP
		/* Check wheter to swith to DP100% from the next iteration */
		if( fabs(loglik_pre-loglik) < threshold_to_dp ) {
			switch_to_dp = 1;
		}

		if( 0 == rank && 0 != data->precision )
			fprintf(stderr, "SWITCH_TO_DP iter %d : %d -- %lf %lf\n", data->iter_count, switch_to_dp, fabs(loglik_pre-loglik), threshold_to_dp);
#endif

#ifndef EXAGEOSTAT_ON_SUMMIT
		if( 0 == rank ) {
			if(strcmp(data->checkpoint_file, "") != 0)
				checkpointing(data->checkpoint_file, data->iter_count, theta, loglik, 3);
		}
#endif

	} /* END if not checking point */

#if defined(PARSEC_HAVE_MPI)
	/* Use one desc */
	if( 0 == rank )
#endif
	{
		fprintf(stderr, "\n------ddotproduct: %.17lf ", data->dotp);
		fprintf(stderr, "\n------logdet: %.17lf ", logdet);
		if(num_params==3){
			fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness): (%.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n", data->iter_count, theta[0], theta[1], theta[2], loglik);
		} else{
			fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness, noise): (%.10lf, %.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n", data->iter_count, theta[0], theta[1], theta[2], theta[3], loglik);
		}
		if(data->log == 1)
			fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n", data->iter_count, theta[0], theta[1], theta[2],loglik);

		fprintf(stderr, " ---- Facto Time: %.2lf\n", time_facto);
		fprintf(stderr, " ---- logdet Time: %.2lf\n", logdet_calculate);
		fprintf(stderr, " ---- dtrsm Time: %.2lf\n", time_solve);
		fprintf(stderr, " ---- MLE Time: %.2lf\n", time_mle);
		fprintf(stderr, " ---- Matrix Generation Time: %.2lf\n", matrix_gen_time);
		fprintf(stderr, " ---- Total Time: %.2lf\n", matrix_gen_time+ time_facto + logdet_calculate + time_solve);
		fprintf(stderr, " ---- Gflop/s: %.2lf\n", (time_facto  + time_solve == 0.0 )? 0.0 : flops / 1e9 / (time_facto  + time_solve));
		fprintf(stderr, "\n\n");
	}

	data->iter_count++;
	/* for experiments */
	data->avg_exec_time_per_iter+=/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
	data->avg_flops_per_iter+=flops / 1e9 / (time_facto +time_solve);
	data->final_loglik=loglik;

#if SWITCH_TO_DP
	/* Force terminate */
	if( fabs(loglik_pre-loglik) < threshold ) {
	//if( data->iter_count > 10 ) {
		if( 0 == rank )
			fprintf(stderr, "\nForce terminate iter %d : %.17lf %.17lf %p\n\n", data->iter_count-1, fabs(loglik_pre-loglik), threshold, *data->opt);
		nlopt_force_stop( *data->opt );
	}

	/* Record the max loglik */ 
	if( 1 || loglik_pre < loglik || 0.0 == loglik_pre ) {
		if( 0 == rank )
			fprintf(stderr, "\nSet loglik_pre %d : %lf %lf\n\n", data->iter_count-1, loglik, loglik_pre);
		loglik_pre = loglik;
	}
#endif

	return loglik;
}

double MLE_alg_parsec(unsigned n, const double * theta, double * grad, void * data) {
	/* Synchronization before optimization */
	MPI_Barrier( MPI_COMM_WORLD );
	if ( 0 == ((MLE_data *)data)->precision || 2 == ((MLE_data *)data)->precision )
		return parsec_dmle_Tile(n, theta,  grad,  data);
	else
		fprintf(stderr, "Wrong precision: it should be 0 (doulbe) or 2 (mix precision)\n");	
}

void parsec_dmle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts,
		int p_grid, int q_grid, int mse_flag) {
	MLE_data     *data              = (MLE_data*) MORSE_data;
	int nodes = data->nodes;
	int rank = data->rank;
	PLASMA_enum uplo = (PLASMA_enum)data->uplo;

	if(nZmiss <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return;
	}

	two_dim_block_cyclic_t *descZobs = (two_dim_block_cyclic_t *)data->descZobs;
	two_dim_block_cyclic_init(descZobs, matrix_RealDouble, matrix_Tile,
			nodes, rank, dts, dts, nZobs, 1, 0, 0,
			nZobs, 1, 1, 1, nodes);
	descZobs->mat = calloc((size_t)descZobs->super.nb_local_tiles *
			(size_t)descZobs->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descZobs->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descZobs, "descZobs");

	two_dim_block_cyclic_t *descZmiss = (two_dim_block_cyclic_t *)data->descZmiss;
	two_dim_block_cyclic_init(descZmiss, matrix_RealDouble, matrix_Tile,
			nodes, rank, dts, dts, nZmiss, 1, 0, 0,
			nZmiss, 1, 1, 1, nodes);
	descZmiss->mat = calloc((size_t)descZmiss->super.nb_local_tiles *
			(size_t)descZmiss->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descZmiss->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descZmiss, "descZmiss");

	two_dim_block_cyclic_t *descC12 = (two_dim_block_cyclic_t *)data->descC12;
	two_dim_block_cyclic_init(descC12, matrix_RealDouble, matrix_Tile,
			nodes, rank, dts, dts, nZmiss, nZobs, 0, 0,
			nZmiss, nZobs, 1, 1, p_grid);
	descC12->mat = calloc((size_t)descC12->super.nb_local_tiles *
			(size_t)descC12->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC12->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC12, "descC12");

	sym_two_dim_block_cyclic_t *descC22 = (sym_two_dim_block_cyclic_t *)data->descC22;
	sym_two_dim_block_cyclic_init(descC22, matrix_RealDouble,
			nodes, rank, dts, dts, nZobs, nZobs, 0, 0,
			nZobs, nZobs, p_grid, uplo);
	descC22->mat = calloc((size_t)descC22->super.nb_local_tiles *
			(size_t)descC22->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC22->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC22, "descC22");

        two_dim_block_cyclic_t *descZtrace = (two_dim_block_cyclic_t *)data->descZtrace;
        two_dim_block_cyclic_init(descZtrace, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZmiss, 1, 0, 0,
                        nZmiss, 1, 1, 1, p_grid);
        descZtrace->mat = calloc((size_t)descZtrace->super.nb_local_tiles *
                        (size_t)descZtrace->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descZtrace->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descZtrace, "descZtrace");

        two_dim_block_cyclic_t *descC11 = (two_dim_block_cyclic_t *)data->descC11;
        two_dim_block_cyclic_init(descC11, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZmiss, nZmiss, 0, 0,
                        nZmiss, nZmiss, 1, 1, p_grid);
        descC11->mat = calloc((size_t)descC11->super.nb_local_tiles *
                        (size_t)descC11->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descC11->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descC11, "descC11");

        two_dim_block_cyclic_t *descC21 = (two_dim_block_cyclic_t *)data->descC21;
        two_dim_block_cyclic_init(descC21, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, nZmiss, 0, 0,
                        nZobs, nZmiss, 1, 1, p_grid);
        descC21->mat = calloc((size_t)descC21->super.nb_local_tiles *
                        (size_t)descC21->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descC21->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descC21, "descC21");

	two_dim_block_cyclic_t *descmse = (two_dim_block_cyclic_t *)data->descmse;
	two_dim_block_cyclic_t *descZactual = (two_dim_block_cyclic_t *)data->descZactual;

	if( mse_flag == 1) {
		two_dim_block_cyclic_init(descZactual, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, nZmiss, 1, 0, 0,
				nZmiss, 1, 1, 1, nodes);
		descZactual->mat = calloc((size_t)descZactual->super.nb_local_tiles *
				(size_t)descZactual->super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descZactual->super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)descZactual, "descZactual");

		two_dim_block_cyclic_init(descmse, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, 1, 1, 0, 0,
				1, 1, 1, 1, nodes);
		descmse->mat = calloc((size_t)descmse->super.nb_local_tiles *
				(size_t)descmse->super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descmse->super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)descmse, "descmse");
	}

}

void parsec_prediction_init(MLE_data *data, int nZmiss, int nZobs, int ts, int p_grid, int q_grid, int mse_flag)
	/* @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] nZmiss: number of missing measurements.
	 * @param[in] nZobs: number of observed measurements.
	 * @param[in] p_grid: p_grid in the case of distributed system.
	 * @param[in] q_grid: q_grid in the case of distributed system.
	 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
	 */                     
{                               
	if ( 0 == data->precision || 2 == data->precision )
		parsec_dmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
	else
		fprintf(stderr, "Wrong precision: it should be 0 (doulbe) or 2 (mix precision)\n");
}

void parsec_pick_random_points(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N) {
	/* Initialization */
	location l;
	location *lmiss;
	location *lobs;
	double *Z;                      
	int i = 0;
	int rank = data->rank;

	/* Memory allocation */
	Z       = (double *) malloc(N * sizeof(double));
	l.x     = (double *) malloc(N * sizeof(double));
	l.y     = (double *) malloc(N * sizeof(double)); 

	if(data->l1.z != NULL)
		l.z     = (double *) malloc(N * sizeof(double));
	else
		l.z = NULL;

	parsec_tiled_matrix_dc_t *descZcpy = (parsec_tiled_matrix_dc_t *)data->descZcpy;
	assert( descZcpy->m == N );

	parsec_context_t *parsec = (parsec_context_t *)data->parsec;

	/* Copy observed measurments */
	parsec_get_zobs( parsec, descZcpy, Z, N );

	/* Broadcast Z */
	int root = descZcpy->super.rank_of(&descZcpy->super, 0, 0);
	MPI_Bcast( Z, N, MPI_DOUBLE, root, MPI_COMM_WORLD );

	for( i = 0; i < N ; i++)        
	{
		l.x[i]=data->l1.x[i];
		l.y[i]=data->l1.y[i];

		if(data->l1.z != NULL)
 			l.z[i]=data->l1.z[i];
	}               

	//print_dmatrix("l.x", N, 1, l.x, N);
	//print_dmatrix("l.y", N, 1, l.y, N);
	//printf("N %d nZmiss %d nZobs %d\n", N, nZmiss, nZobs);
	shuffle(Z, &l, N);

	//print_dmatrix("Z after", N, 1, Z, N);

	for( i = 0; i < nZobs ; i++)
		Zobs[i] = Z[nZmiss+i];

	for ( i = 0; i < nZmiss; i++)
		Zactual[i] = Z[i];

	//print_dmatrix("Zobs", nZobs, 1, Zobs, nZobs);
	//print_dmatrix("Zactual", nZmiss, 1, Zactual, nZmiss);

	lmiss = &(data->lmiss);
	lobs = &(data->lobs);
	lmiss->x = l.x;
        lmiss->y = l.y;
        lobs->x  = &l.x[nZmiss];
        lobs->y  = &l.y[nZmiss];

	if(data->l1.z != NULL)
 	{
 		lmiss->z = l.z;
 		lobs->z  = &l.z[nZmiss];
 	}

        locations_obs_zsort_inplace(nZobs, lobs, Zobs);
	locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

        free(Z);
}


double parsec_dmle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n) {
        double time_solve = 0.0;
        double mat_gen_time = 0.0;
        double time_gemm = 0.0;
        double time_mse = 0.0;
        double flops = 0.0;
	int info = 0;

	parsec_context_t *parsec = NULL;
        parsec_tiled_matrix_dc_t *descZmiss   = NULL;
        parsec_tiled_matrix_dc_t *descC12     = NULL;
        parsec_tiled_matrix_dc_t *descC22     = NULL;
        parsec_tiled_matrix_dc_t *descZtrace  = NULL;
        parsec_tiled_matrix_dc_t *descC11     = NULL;
        parsec_tiled_matrix_dc_t *descC21     = NULL;
        parsec_tiled_matrix_dc_t *descmse     = NULL;
        parsec_tiled_matrix_dc_t *descZactual = NULL;
        parsec_tiled_matrix_dc_t *descZobs    = NULL;
        parsec_tiled_matrix_dc_t *descdet     = NULL;
        MLE_data     *data              = (MLE_data*) MORSE_data;
        data->mserror                   = 0.0;

        if(nZmiss <= 0)
        {
                fprintf(stderr," Number of missing values should be positive value\n");
                return -1;
        }

        descZmiss         = (parsec_tiled_matrix_dc_t *)data->descZmiss;
        descC12           = (parsec_tiled_matrix_dc_t *)data->descC12;
        descC22           = (parsec_tiled_matrix_dc_t *)data->descC22;
        descZtrace        = (parsec_tiled_matrix_dc_t *)data->descZtrace;
        descC11           = (parsec_tiled_matrix_dc_t *)data->descC11;
        descC21           = (parsec_tiled_matrix_dc_t *)data->descC21;
        descmse           = (parsec_tiled_matrix_dc_t *)data->descmse;
        descZactual       = (parsec_tiled_matrix_dc_t *)data->descZactual;
        descZobs          = (parsec_tiled_matrix_dc_t *)data->descZobs;
        descdet           = (parsec_tiled_matrix_dc_t *)data->descdet;
        parsec            = (parsec_context_t *)data->parsec;
	PLASMA_enum uplo = data->uplo;
        int rank = data->rank;
	int HNB = data->HNB;
        int success;

	/* Synchronization before prediction */
	MPI_Barrier( MPI_COMM_WORLD );

        VERBOSE("Copy measurments vector to descZobs descriptor...");
        parsec_Lapack_to_Tile( parsec, descZobs, Zobs, nZobs );
        VERBOSE(" Done.\n");

	//dplasma_dprint(parsec, uplo, descZobs);

        if( Zactual != NULL)
        {
                VERBOSE("Copy actual measurments vector to descZactual descriptor...");
                parsec_Lapack_to_Tile( parsec, descZactual, Zactual, nZmiss );
                VERBOSE(" Done.\n");
        }

#if defined(PARSEC_HAVE_MPI)
	MPI_Bcast( &data->variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
	if(strcmp(data->c_fun, "matern") == 0)
                theta[0] = data->variance;
        if( 0 == rank ) fprintf(stderr, "estimated parameters: %f - %f - %f (%s)\n", theta[0], theta[1], theta[2], data->c_fun);

#if defined(EXAGEOSTAT_USE_CUDA) 
        /* Reset cache on GPU before cholesky */
        parsec_devices_release_memory();
        parsec_devices_reset_load(parsec);
#endif

	/* Generate C22 covariance matrix */
        START_TIMING(mat_gen_time);
        VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	parsec_dmatrix_generation( parsec, descC22, &data->lobs, &data->lobs,
			theta, data->dm, data->c_fun, data->band_size_double );
        VERBOSE(" Done.\n");

	/* Generate C12 covariance matrix */
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");

	parsec_taskpool_t *parsec_dmatrix_generation_C12 = NULL;
	parsec_dmatrix_generation_C12 = parsec_dmatrix_generation_New(
			parsec, descC12, &data->lmiss, &data->lobs,
			theta, data->dm, data->c_fun, data->band_size_double );

	/* TODO disable GPU, it seems issues when on GPU */
        disable_GPU( parsec_dmatrix_generation_C12 );

	if( parsec_dmatrix_generation_C12 != NULL ) {
		parsec_context_add_taskpool(parsec, parsec_dmatrix_generation_C12);
		parsec_context_start(parsec);
		parsec_context_wait(parsec);
		parsec_dmatrix_generation_Destruct(parsec_dmatrix_generation_C12);
	}

        VERBOSE(" Done.\n");

        /* Generate C11 covariance matrix */
        VERBOSE("Generate C11 Covariance Matrix... (Prediction Stage)");
        parsec_dmatrix_generation( parsec, descC11, &data->lmiss, &data->lmiss,
                        theta, data->dm, data->c_fun, data->band_size_double );
        VERBOSE(" Done.\n");

        /* Generate C21 covariance matrix */
        VERBOSE("Generate C21 Covariance Matrix... (Prediction Stage)");

        parsec_taskpool_t *parsec_dmatrix_generation_C21 = NULL;
        parsec_dmatrix_generation_C21 = parsec_dmatrix_generation_New(
                        parsec, descC21, &data->lobs, &data->lmiss,
                        theta, data->dm, data->c_fun, data->band_size_double );
        
        /* TODO disable GPU, it seems issues when on GPU */
        disable_GPU( parsec_dmatrix_generation_C21 );

        if( parsec_dmatrix_generation_C21 != NULL ) {
                parsec_context_add_taskpool(parsec, parsec_dmatrix_generation_C21);
                parsec_context_start(parsec);
                parsec_context_wait(parsec);
                parsec_dmatrix_generation_Destruct(parsec_dmatrix_generation_C21);
        }

        VERBOSE(" Done.\n");

        STOP_TIMING(mat_gen_time);

	/* Start prediction */
        START_TIMING(time_solve);
        VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");

	/* Define band_size */
	if( 0 == data->precision ) {
		int NT = descC22->lmt ;
		data->band_size_double = NT;
		data->band_size_single = NT;
	}

        VERBOSE("Calculate hsdpotrf C22 Covariance Matrix... (Prediction Stage)");
        /* Start prediction */
        START_TIMING(time_solve);
	/* Convert d2s */
	parsec_band_convert_d2s(parsec, descC22, data->band_size_double);

#if TRSM_CONVERT_HALF
	/* hsdpotrf */
	success = dplasma_hsdpotrf_trsm_convert( parsec, uplo, descC22, data->lookahead, data->sendless,
			data->tensor_gemm, data->band_size_double, data->band_size_single, HNB );
#else
	/* hsdpotrf */
	success = dplasma_hsdpotrf( parsec, uplo, descC22, data->lookahead, data->sendless,
			data->tensor_gemm, data->band_size_double, data->band_size_single, HNB );
#endif

	/* Convert s2d */
	parsec_band_convert_s2d(parsec, descC22, data->band_size_double);
        VERBOSE(" Done.\n");

	if( 0 != success && 0 == rank ) {
                fprintf(stderr, "\n\nFactorization cannot be performed..\n The matrix is not positive definite: %d\n", success);
		exit(1);
        }

        VERBOSE("Triangular solve 1... (Prediction Stage)");
        parsec_taskpool_t *parsec_dtrsm = NULL;
        parsec_dtrsm = my_dtrsm_New( PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1.0, descC22, descZobs ); 
#if NO_GPU_DPLASMA 
	disable_GPU( parsec_dtrsm );
#endif
        if ( parsec_dtrsm != NULL )
        {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                my_dtrsm_Destruct( parsec_dtrsm );
        }

        VERBOSE("Triangular solve 2... (Prediction Stage)");
        parsec_dtrsm = dplasma_dtrsm_New( PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit, 1.0, descC22, descZobs );
#if NO_GPU_DPLASMA 
	disable_GPU( parsec_dtrsm );
#endif
        if ( parsec_dtrsm != NULL )
        {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dtrsm_Destruct( parsec_dtrsm );
        }
        VERBOSE(" Done.\n");

        STOP_TIMING(time_solve);

        flops = flops + FLOPS_DPOTRF(nZobs);
        flops = flops + 2 * FLOPS_DTRSM(MorseLeft, nZobs, nZobs);

	//fprintf(stderr, "\nC12: \n");
	//dplasma_dprint(parsec, matrix_UpperLower, descC12);
	//fprintf(stderr, "\nZobs: \n");
	//dplasma_dprint(parsec, uplo, descZobs);
        //fprintf(stderr, "\nZmiss before: \n");
        //dplasma_dprint(parsec, uplo, descZmiss);

        START_TIMING(time_gemm);
        VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");   
        parsec_taskpool_t *parsec_dgemm = NULL;
        parsec_dgemm = dplasma_dgemm_New( PlasmaNoTrans, PlasmaNoTrans, 1.0, descC12, descZobs, 0.0, descZmiss ); 
#if NO_GPU_DPLASMA 
	disable_GPU( parsec_dgemm );
#endif
        if ( parsec_dgemm != NULL )
        {
                parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_dgemm);
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dgemm_Destruct( parsec_dgemm );
        }
        flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
        VERBOSE(" Done.\n");
        STOP_TIMING(time_gemm);

        //fprintf(stderr, "\nZmiss after: \n");
	//dplasma_dprint(parsec, uplo, descZmiss);

        /*****************************************************************************************************/
        /* Schur_complement */
        VERBOSE("Calculate dposv C22 Covariance Matrix and C21...(Prediction Stage - Schur_complement)");
        parsec_dtrsm = dplasma_dtrsm_New( PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1.0, descC22, descC21 );
#if NO_GPU_DPLASMA 
        disable_GPU( parsec_dtrsm );
#endif
        if ( parsec_dtrsm != NULL )
        {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dtrsm_Destruct( parsec_dtrsm );
        }
        VERBOSE(" Done.\n");

        VERBOSE("Calculate dposv C22 Covariance Matrix and C21...(Prediction Stage - Schur_complement)");
        parsec_dtrsm = dplasma_dtrsm_New( PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit, 1.0, descC22, descC21 );
#if NO_GPU_DPLASMA 
        disable_GPU( parsec_dtrsm );
#endif
        if ( parsec_dtrsm != NULL )
        {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dtrsm_Destruct( parsec_dtrsm );
        }

        //dplasma_dtrsm( parsec, PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1.0, descC22, descC21 );
        //dplasma_dtrsm( parsec, PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit, 1.0, descC22, descC21 );
        flops = flops + FLOPS_DTRSM(MorseLeft, nZobs, nZobs);
        VERBOSE(" Done.\n");

        VERBOSE("Calculate dgemm C11 = C12 * C21 Covariance Matrix... (Prediction Stage - Schur_complement)");
        dplasma_dgemm( parsec, PlasmaNoTrans, PlasmaNoTrans, -1.0, descC12, descC21, 1.0, descC11 );
        flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
        VERBOSE(" Done.\n");

        /* Use data->det to store the trace summation */
        data->det=0;
        VERBOSE("Calculate trace estimation... (Prediction Stage - Schur_complement)");
	parsec_MLE_dtrace(parsec, descC11, descdet, descZtrace);
	int root = descdet->super.rank_of(&descdet->super, 0, 0);
	if( descdet->super.myrank == root )
		data->det = *( (double *)((descdet->super.data_of(&descdet->super, 0, 0))->device_copies[0]->device_private) );
#if defined(PARSEC_HAVE_MPI)
	MPI_Bcast( (void *)&data->det, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif

        VERBOSE(" Done.\n");

	double mean = (data->det)/nZmiss;
	double sd = 0;
	double sum_trace = 0.0;
	double* Ztrace = (double *) malloc(nZmiss  * sizeof(double));
	parsec_get_zobs( parsec, descZtrace, Ztrace, nZmiss ); 

#if defined(PARSEC_HAVE_MPI)
         MPI_Bcast( (void *)Ztrace, nZmiss, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif

	for(int i = 0; i < nZmiss; i++) {
		sd += pow((Ztrace[i] - mean), 2);
		sum_trace += Ztrace[i];
        }
	sd = sqrt(sd/nZmiss);
	free( Ztrace );

        data->trace_sum = sum_trace;
        data->trace_mean = mean;
        data->trace_sd = sd;

#if defined(PARSEC_HAVE_MPI)
        if(descC22->super.myrank == 0)
#endif
	{
		fprintf(stderr, "Trace_sum %d : %le %le\n", rank, data->det, sum_trace);
		printf("Trace estimation (trace standard deviation ): %6.4e\n", sd);
		printf("Trace estimation (trace mean ): %6.4e\n", mean);
		printf("Trace estimation (trace sum ): %6.4e\n", data->det);
	}

	/*****************************************************************************************************/

	/* return back descZmiss to zmiss vector */
	parsec_get_zobs( parsec, descZmiss, Zmiss, nZmiss ); 

        root = descZmiss->super.rank_of(&descZmiss->super, 0, 0);
	MPI_Bcast( Zmiss, nZmiss, MPI_DOUBLE, root, MPI_COMM_WORLD );

	/* Estimate Mean Square Error */
        if( Zactual != NULL)
        {
                START_TIMING(time_mse);
                VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
		parsec_mse_calculation(parsec, descZactual, descZmiss, descmse);
		/* Broadcast determinant */
		root = descmse->super.rank_of(&descmse->super, 0, 0);
		if( descmse->super.myrank == root )
			data->mserror = *( (double *)((descmse->super.data_of(&descmse->super, 0, 0))->device_copies[0]->device_private) );

#if defined(PARSEC_HAVE_MPI)
		MPI_Bcast( (void *)&data->mserror, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
		VERBOSE(" Done.\n");    
		STOP_TIMING(time_mse);
		data->mserror /= nZmiss;
	}
        else
                data->mserror = -1;
    
#if defined(PARSEC_HAVE_MPI)
        if(descC22->super.myrank == 0)
#endif
        {
                if(data->log == 1)
                        fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

#ifndef EXAGEOSTAT_ON_SUMMIT
                write_prediction_result("predict_result.dat", n, data->hicma_acc, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));
#endif
        }

        return data->mserror;
}

void parsec_prediction_finalize(MLE_data *data) {
	parsec_data_free( ((two_dim_block_cyclic_t *)data->descZobs)->mat );
	parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descZobs );

        parsec_data_free( ((two_dim_block_cyclic_t *)data->descZmiss)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descZmiss );

        parsec_data_free( ((two_dim_block_cyclic_t *)data->descC12)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descC12 );

        parsec_data_free( ((sym_two_dim_block_cyclic_t *)data->descC22)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descC22 );

        parsec_data_free( ((two_dim_block_cyclic_t *)data->descZtrace)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descZtrace );

        parsec_data_free( ((two_dim_block_cyclic_t *)data->descC11)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descC11 );

        parsec_data_free( ((two_dim_block_cyclic_t *)data->descC21)->mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descC21 );

	if( 1 == data->mse_flag ) {
		parsec_data_free( ((two_dim_block_cyclic_t *)data->descZactual)->mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descZactual );

		parsec_data_free( ((two_dim_block_cyclic_t *)data->descmse)->mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)data->descmse );
	}
}

void parsec_print_matrix( double *A, int m, int n, int lda ) {
#if 0 
    for( int i = 0; i < m; i++) {
        for( int j = 0; j < n; j++)
            fprintf(stderr, "%lf ", A[j*lda+i]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif
}

void parsec_print_matrix_single( float *A, int m, int n, int lda ) {
#if 0 
    for( int i = 0; i < m; i++) {
        for( int j = 0; j < n; j++)
            fprintf(stderr, "%f ", A[j*lda+i]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif
}

int parsec_sdmatrix_init_operator( parsec_execution_stream_t *es,
                         const parsec_tiled_matrix_dc_t *descA,
                         void *_A,
                         PLASMA_enum uplo, int m, int n,
                         void *op_data )
{
    int tempmm, tempnn, ldam;
    matrix_init_args_t     *args = (matrix_init_args_t*)op_data;
    (void)es;
    (void)uplo;

    tempmm = ((m)==((descA->mt)-1)) ? ((descA->m)-(m*(descA->mb))) : (descA->mb);
    tempnn = ((n)==((descA->nt)-1)) ? ((descA->n)-(n*(descA->nb))) : (descA->nb);
    ldam   = descA->mb; 

    if( m-n < args->band_size_double ) {
        double *A    = (double*)_A;
        CORE_dplgsy(
            args->bump, tempmm, tempnn, A, ldam,
            descA->m, m*descA->mb, n*descA->nb, args->seed );
    } else {
        float *A    = (float*)_A;
        CORE_splgsy(
            args->bump, tempmm, tempnn, A, ldam,
            descA->m, m*descA->mb, n*descA->nb, args->seed );
    }

    return 0;
}

void parsec_mloe_mmom_init( MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid ) {
        two_dim_block_cyclic_t *desck_t         = NULL;
        two_dim_block_cyclic_t *desck_a         = NULL;
        two_dim_block_cyclic_t *desck_atmp      = NULL;
        two_dim_block_cyclic_t *desck_ttmp      = NULL;

        two_dim_block_cyclic_t *descK_t         = NULL;
        two_dim_block_cyclic_t *descK_ttmp      = NULL;
        two_dim_block_cyclic_t *descK_a         = NULL;
        
        two_dim_block_cyclic_t *descexpr1       = NULL;
        two_dim_block_cyclic_t *descexpr2       = NULL;
        two_dim_block_cyclic_t *descexpr3       = NULL;
        two_dim_block_cyclic_t *descexpr4       = NULL;
        two_dim_block_cyclic_t *desc_mloe_mmom       = NULL;
        MLE_data     *data              = (MLE_data*) MORSE_data;
        int nodes = data->nodes;
        int rank = data->rank;
        PLASMA_enum uplo = (PLASMA_enum)data->uplo;

        if(nZmiss <= 0)
        {
                if( 0 == rank ) fprintf(stderr," Number of missing values should be positive value\n");
                return;
        }

        two_dim_block_cyclic_init(desc_mloe_mmom, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZmiss, 3, 0, 0,
                        nZmiss, 3, 1, 1, nodes);
        desc_mloe_mmom->mat = calloc((size_t)desc_mloe_mmom->super.nb_local_tiles *
                        (size_t)desc_mloe_mmom->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(desc_mloe_mmom->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)desc_mloe_mmom, "desc_mloe_mmom");

        two_dim_block_cyclic_init(desck_t, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, 1, 0, 0,
                        nZobs, 1, 1, 1, nodes);
        desck_t->mat = calloc((size_t)desck_t->super.nb_local_tiles *
                        (size_t)desck_t->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(desck_t->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)desck_t, "desck_t");

        two_dim_block_cyclic_init(desck_a, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, 1, 0, 0,
                        nZobs, 1, 1, 1, nodes);
        desck_a->mat = calloc((size_t)desck_a->super.nb_local_tiles *
                        (size_t)desck_a->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(desck_a->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)desck_a, "desck_a");

        two_dim_block_cyclic_init(desck_atmp, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, 1, 0, 0,
                        nZobs, 1, 1, 1, nodes);
        desck_atmp->mat = calloc((size_t)desck_atmp->super.nb_local_tiles *
                        (size_t)desck_atmp->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(desck_atmp->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)desck_atmp, "desck_atmp");

        two_dim_block_cyclic_init(desck_ttmp, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, 1, 0, 0,
                        nZobs, 1, 1, 1, nodes);
        desck_ttmp->mat = calloc((size_t)desck_ttmp->super.nb_local_tiles *
                        (size_t)desck_ttmp->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(desck_ttmp->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)desck_ttmp, "desck_ttmp");

        two_dim_block_cyclic_init(descexpr1, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, 1, 1, 0, 0,
                        1, 1, 1, 1, nodes);
        descexpr1->mat = &data->expr1; 
        parsec_data_collection_set_key((parsec_data_collection_t*)descexpr1, "descexpr1");

        two_dim_block_cyclic_init(descexpr2, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, 1, 1, 0, 0,
                        1, 1, 1, 1, nodes);
        descexpr2->mat = &data->expr2;
        parsec_data_collection_set_key((parsec_data_collection_t*)descexpr2, "descexpr2");

        two_dim_block_cyclic_init(descexpr3, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, 1, 1, 0, 0,
                        1, 1, 1, 1, nodes);
        descexpr3->mat = &data->expr3;
        parsec_data_collection_set_key((parsec_data_collection_t*)descexpr3, "descexpr3");

        two_dim_block_cyclic_init(descexpr4, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, 1, 1, 0, 0,
                        1, 1, 1, 1, nodes);
        descexpr4->mat = &data->expr4;
        parsec_data_collection_set_key((parsec_data_collection_t*)descexpr4, "descexpr4");

        two_dim_block_cyclic_init(descK_t, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, nZobs, 0, 0,
                        nZobs, nZobs, 1, 1, nodes);
        descK_t->mat = calloc((size_t)descK_t->super.nb_local_tiles *
                        (size_t)descK_t->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descK_t->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descK_t, "descK_t");

        two_dim_block_cyclic_init(descK_a, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, nZobs, 0, 0,
                        nZobs, nZobs, 1, 1, nodes);
        descK_a->mat = calloc((size_t)descK_a->super.nb_local_tiles *
                        (size_t)descK_a->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descK_a->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descK_a, "descK_a");

        two_dim_block_cyclic_init(descK_ttmp, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, nZobs, nZobs, 0, 0,
                        nZobs, nZobs, 1, 1, nodes);
        descK_ttmp->mat = calloc((size_t)descK_ttmp->super.nb_local_tiles *
                        (size_t)descK_ttmp->super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(descK_ttmp->super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)descK_ttmp, "descK_ttmp");

	/* Initiae data descriptors */
        data->desck_t           = desck_t;
        data->desck_a           = desck_a;
        data->descK_ttmp        = descK_ttmp;
        data->desck_atmp        = desck_atmp;
        data->desck_ttmp        = desck_ttmp;
        data->descK_t           = descK_t;
        data->descK_a           = descK_a;
        data->descexpr1         = descexpr1;
        data->descexpr2         = descexpr2;
        data->descexpr3         = descexpr3;
        data->descexpr4         = descexpr4;
        data->desc_mloe_mmom    = desc_mloe_mmom;
}

location* parsec_readLocsFile(char *locs_file, int n, int dts)
//! Read 2D locations from a given
/*! flatfile
 * Returns location struct.
 * @param[in] locs_file: 2D location file.
 * @param[in]  n : number of spatial locations.
 * */
{

        FILE *fp;
        int i = 0;
        char *line = NULL;
        size_t len  = 0;
        ssize_t read;
        char *pch;
        location *locations;


        fp = fopen(locs_file, "r");
        if (fp == NULL)
        {
                printf("cannot read locations file\n");
                printf("%s: \n",locs_file);
                return NULL;
        }
        else
        {
                //Allocate memory
                locations               = (location *) malloc(sizeof(location));
                locations->x            = (double *) malloc((n+dts) * sizeof(double));
                locations->y            = (double *) malloc((n+dts) * sizeof(double));
        }

        while ((read = getline(&line, &len, fp)) != -1) {
                pch = strtok(line, ",");
                while (pch != NULL)
                {
                        locations->x[i] = atof(pch);
                        pch = strtok (NULL, ",");
                        locations->y[i] = atof(pch);
                        pch = strtok (NULL, ",");
                }
                i++;
        }
       fclose(fp);
        if (line)
                free(line);
    //    zsort_locations(n,locations);
        return locations;
}


location* parsec_readLocsFile3d(char *locs_file, int n, int dts)
//! Read 2D locations from a given
/*! flatfile
 * Returns location struct.
 * @param[in] locs_file: 2D location file.
 * @param[in]  n : number of spatial locations.
 * */
{

        FILE *fp;
        int i = 0;
        char *line = NULL;
        size_t len  = 0;
        ssize_t read;
        char *pch;
        location *locations;


        fp = fopen(locs_file, "r");
        if (fp == NULL)
        {
                printf("cannot read locations file\n");
                printf("%s: \n",locs_file);
                return NULL;
        }
        else
        {
                //Allocate memory
                locations               = (location *) malloc(sizeof(location));
                locations->x            = (double *) malloc((n+dts) * sizeof(double));
                locations->y            = (double *) malloc((n+dts) * sizeof(double));
                locations->z            = (double *) malloc((n+dts) * sizeof(double));
        }

        while ((read = getline(&line, &len, fp)) != -1) {
                pch = strtok(line, ",");
                while (pch != NULL)
                {
                        locations->x[i] = atof(pch);
                        pch = strtok (NULL, ",");
                        locations->y[i] = atof(pch);
                        pch = strtok (NULL, ",");
                        locations->z[i] = atof(pch);
                        pch = strtok (NULL, ",");
                }
                i++;
        }
        fclose(fp);
        if (line)
                free(line);
        //    zsort_locations(n,locations);
        return locations;
}

location* parsec_GenerateXYLoc(int n, int seed, int dts)
        //! Generate XY location for exact computation (MOORSE)        
{       
        //initalization
        int i = 0 ,index = 0, j = 0;
        // unsigned int *seed = &exageostat_seed;
        srand(seed);
        location* locations = (location *) malloc( sizeof(location));
        //Allocate memory
        locations->x            = (double *) malloc((n+dts) * sizeof(double));
        locations->y            = (double *) malloc((n+dts) * sizeof(double));
        locations->z            = NULL;
        // if(strcmp(locs_file, "") == 0)
        // {
        
        int sqrtn = (int)(sqrt(n)+0.1);
        
        //Check if the input is square number or not
        if(pow(sqrtn,2) != n)
        {
                printf("Please use a perfect square number to generate a valid synthetic dataset.....\n\n");
                exit(0);     
        }

        int *grid = (int *) calloc((int)sqrtn, sizeof(int));
        
        for(i = 0; i < sqrtn; i++)
        {
                grid[i] = i+1;
        }
        
        for(i = 0; i < sqrtn; i++)
                for(j = 0; j < sqrtn; j++){
                        locations->x[index] = (grid[i]-0.5+uniform_distribution(-0.4, 0.4))/sqrtn;
                        locations->y[index++] = (grid[j]-0.5+uniform_distribution(-0.4, 0.4))/sqrtn;
                }
        free(grid);
        void zsort_locations(int, location*);
        zsort_locations(n, locations); 
        return locations;
}

location* parsec_GenerateXYZLoc(int n, int seed, int dts)
        //! Generate XY location for exact computation (MOORSE)
{       
        //initalization
        int i = 0 ,index = 0, j = 0, k = 0;
        // unsigned int *seed = &exageostat_seed;
        srand(seed);
        location* locations = (location *) malloc( sizeof(location) );
        //Allocate memory
        locations->x            = (double *) malloc((n+dts) * sizeof(double));
        locations->y            = (double *) malloc((n+dts) * sizeof(double));
        locations->z            = (double *) malloc((n+dts) * sizeof(double));
    
    
        int cbrtn = (int)(cbrt(n)+0.1);

        //Check if the input is square number or not
        if(pow(cbrtn, 3) != n)
        {
                printf("Please use a perfect cubic number to generate a valid synthetic dataset.....\n\n");
                exit(0);
        }

        int *grid = (int *) calloc((int)cbrtn, sizeof(int));
        for(i = 0; i < cbrtn; i++)
        {
                grid[i] = i+1;
        }

        printf("cbrtn:%d\n", cbrtn);

        for(i = 0; i < cbrtn; i++)
                for(j = 0; j < cbrtn; j++)
                        for(k = 0; k < cbrtn; k++)
                        {
                                locations->x[index] = (grid[i]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                                locations->y[index] = (grid[j]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                                locations->z[index] = (grid[k]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                                index++;
                        }
        free(grid);
        void zsort_locations_3d(int, location*);
        zsort_locations_3d(n, locations);
        return locations;
}
