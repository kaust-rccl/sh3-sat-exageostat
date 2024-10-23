/*
 * Copyright (c) 2010-2018 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @generated d Tue Mar 12 15:47:17 2019
 *
 */

#include "mix_precision_internal.h"
#include "hsdpotrf_L_4taskclass.h"
#include "include/dplasmaaux.h"

/* Print more */
static int print_more = 0;

/* Info */
int * potrf_info;

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
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
#endif
#endif

/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 * dplasma_hsdpotrf_New - Generates the taskpool that Computes the Cholesky
 * factorization of a symmetric positive definite (or Hermitian positive
 * definite in the complex case) matrix A, with or without recursive calls.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 * WARNING: The computations are not done by this call.
 *
 * If you want to enable the recursive DAGs, don't forget to set the recursive
 * tile size and to synchonize the taskpool ids after the computations since those
 * are for now local. You can follow the code of dplasma_hsdpotrf_rec() as an
 * example to do this.
 *
 * Hierarchical DAG Scheduling for Hybrid Distributed Systems; Wu, Wei and
 * Bouteiller, Aurelien and Bosilca, George and Faverge, Mathieu and Dongarra,
 * Jack. 29th IEEE International Parallel & Distributed Processing Symposium,
 * May 2015. (https://hal.inria.fr/hal-0107835)
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
 *
 * @param[out] info
 *          Address where to store the output information of the factorization,
 *          this is not synchronized between the nodes, and might not be set
 *          when function exists.
 *          On DAG completion:
 *              - info = 0 on all nodes if successful.
 *              - info > 0 if the leading minor of order i of A is not positive
 *                definite, so the factorization could not be completed, and the
 *                solution has not been computed. Info will be equal to i on the
 *                node that owns the diagonal element (i,i), and 0 on all other
 *                nodes.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_hsdpotrf_Destruct();
 *
 *******************************************************************************
 * @sa dplasma_hsdpotrf
 * @sa dplasma_hsdpotrf_Destruct
 * @sa dplasma_cpotrf_New
 * @sa dplasma_hsdpotrf_New
 * @sa dplasma_spotrf_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_hsdpotrf_4taskclass_New( parsec_context_t *parsec,
                    PLASMA_enum uplo,
                    parsec_tiled_matrix_dc_t *A,
                    int *info,
                    int lookahead, int send_less, int tensor_gemm,
                    int band_size_double, int band_size_single, int HNB )
{
    parsec_hsdpotrf_L_4taskclass_taskpool_t *parsec_dpotrf = NULL;
    parsec_taskpool_t *tp = NULL;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        dplasma_error("dplasma_hsdpotrf_4taskclass_New", "illegal value of uplo");
        return NULL /*-1*/;
    }

#if defined(EXAGEOSTAT_USE_CUDA)
    /** Find all CUDA devices */
    int nb = 0;
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA == device->type ) {
            nb++;
        }
    }

    if(nb == 0) {
        char hostname[256];
        gethostname(hostname, 256);
        if( print_more ) {
            fprintf(stderr, "No CUDA device found on rank %d on %s\n",
                    parsec->my_rank, hostname);
        }

        /* No half precision on CPU */
        assert( band_size_single >= A->lmt );
    }

    int *dev_index = (int*)malloc(nb * sizeof(int));
    nb = 0;
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA == device->type ) {
            dev_index[nb++] = device->device_index;
        }
    }

#if !GPU_BUFFER_ONCE 
    /* Declare workspace used on GPU */
    parsec_potrf_workspace_t *ws_handle_cusolver, *ws_double1, *ws_double2, *ws_single1, *ws_single2;
    parsec_potrf_workspace_t *ws_handle_cublas_tensor, *ws_handle_cublas_gpu, *ws_half1, *ws_half2, *ws_half3;

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
                status_cusolver = cusolverDnDpotrf_bufferSize(handle_cusolver, CUBLAS_FILL_MODE_LOWER, A->nb, NULL, A->mb, &workspace_size);
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
                int workspace_size = A->mb * A->nb;
                ws_double1->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_double1->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_double1->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_double2->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) );
                assert(NULL != ws_double2->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_double2->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for single */
            {
                int workspace_size = A->mb * A->nb;
                ws_single1->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) );
                assert(NULL != ws_single1->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_single1->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;

                ws_single2->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(float) );
                assert(NULL != ws_single2->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_single2->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for half */
            {
                int workspace_size = A->mb * A->nb;
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

        }
    }
#endif /* GPU_BUFFER_ONCE */
#endif

    *info = 0;
    tp = (parsec_taskpool_t*)parsec_hsdpotrf_L_4taskclass_new( uplo, A, info );

    parsec_dpotrf = (parsec_hsdpotrf_L_4taskclass_taskpool_t*)tp;

    parsec_dpotrf->_g_PRI_CHANGE = dplasma_aux_get_priority_limit( "POTRF", A );
    if(0 == parsec_dpotrf->_g_PRI_CHANGE)
          parsec_dpotrf->_g_PRI_CHANGE = A->nt;

    parsec_dpotrf->_g_lookahead = lookahead;
    parsec_dpotrf->_g_send_less = send_less;
    parsec_dpotrf->_g_tensor_gemm = tensor_gemm;
    parsec_dpotrf->_g_band_size_double = band_size_double;
    parsec_dpotrf->_g_band_size_single = band_size_single;
    if( HNB > 0 && HNB < A->nb ) {
        parsec_dpotrf->_g_smallnb = HNB;
    }

    parsec_dpotrf->_g_p_work_double = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_dpotrf->_g_p_work_double, A->mb * A->nb * sizeof(double) );

    parsec_dpotrf->_g_p_work_single = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_dpotrf->_g_p_work_single, A->mb * A->nb * sizeof(float) );

#if defined(EXAGEOSTAT_USE_CUDA)
        parsec_dpotrf->_g_ws_handle_cusolver = (void *)ws_handle_cusolver;
        parsec_dpotrf->_g_ws_handle_cublas_tensor = (void *)ws_handle_cublas_tensor;
        parsec_dpotrf->_g_ws_handle_cublas_gpu = (void *)ws_handle_cublas_gpu;
        parsec_dpotrf->_g_ws_double1 = (void *)ws_double1;
        parsec_dpotrf->_g_ws_double2 = (void *)ws_double2;
        parsec_dpotrf->_g_ws_single1 = (void *)ws_single1;
        parsec_dpotrf->_g_ws_single2 = (void *)ws_single2;
        parsec_dpotrf->_g_ws_half1 = (void *)ws_half1;
        parsec_dpotrf->_g_ws_half2 = (void *)ws_half2;
        parsec_dpotrf->_g_ws_half3 = (void *)ws_half3;
        parsec_dpotrf->_g_nb_cuda_devices = nb;
        parsec_dpotrf->_g_cuda_device_index = dev_index;
#endif

    parsec_matrix_add2arena(parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_DOUBLE_ARENA],
                            parsec_datatype_double_t, matrix_UpperLower,
                            1, A->mb, A->nb, A->mb,
                            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    parsec_matrix_add2arena(parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_SINGLE_ARENA],
                            parsec_datatype_float_t, matrix_UpperLower,
                            1, A->mb, A->nb, A->mb,
                            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    parsec_matrix_add2arena(parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_HALF_ARENA],
                            MPI_BYTE, matrix_UpperLower,
                            1, A->mb*2, A->nb, A->mb*2,
                            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 *  dplasma_hsdpotrf_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_hsdpotrf_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_hsdpotrf_New
 * @sa dplasma_hsdpotrf
 *
 ******************************************************************************/
void
dplasma_hsdpotrf_4taskclass_Destruct( parsec_taskpool_t *tp )
{
    parsec_hsdpotrf_L_4taskclass_taskpool_t *parsec_dpotrf = (parsec_hsdpotrf_L_4taskclass_taskpool_t *)tp;

#if defined(EXAGEOSTAT_USE_CUDA)
    if( parsec_dpotrf->_g_nb_cuda_devices > 0 ) {
#if !GPU_BUFFER_ONCE 
        workspace_memory_free( parsec_dpotrf->_g_ws_handle_cusolver );
        workspace_memory_free( parsec_dpotrf->_g_ws_handle_cublas_tensor );
        workspace_memory_free( parsec_dpotrf->_g_ws_handle_cublas_gpu );
        workspace_memory_free( parsec_dpotrf->_g_ws_double1 );
        workspace_memory_free( parsec_dpotrf->_g_ws_double2 );
        workspace_memory_free( parsec_dpotrf->_g_ws_single1 );
        workspace_memory_free( parsec_dpotrf->_g_ws_single2 );
        workspace_memory_free( parsec_dpotrf->_g_ws_half1 );
        workspace_memory_free( parsec_dpotrf->_g_ws_half2 );
        workspace_memory_free( parsec_dpotrf->_g_ws_half3 );
#endif

        if( NULL != parsec_dpotrf->_g_cuda_device_index )
            free(parsec_dpotrf->_g_cuda_device_index);
    }
#endif

    parsec_matrix_del2arena( parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_DOUBLE_ARENA] );
    parsec_matrix_del2arena( parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_SINGLE_ARENA] );
    parsec_matrix_del2arena( parsec_dpotrf->arenas[PARSEC_hsdpotrf_L_4taskclass_HALF_ARENA] );

    parsec_private_memory_fini( parsec_dpotrf->_g_p_work_double );
    parsec_private_memory_fini( parsec_dpotrf->_g_p_work_single );

    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_double
 *
 * dplasma_hsdpotrf - Computes the Cholesky factorization of a symmetric positive
 * definite (or Hermitian positive definite in the complex case) matrix A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *          \retval > 0 if the leading minor of order i of A is not positive
 *          definite, so the factorization could not be completed, and the
 *          solution has not been computed. Info will be equal to i on the node
 *          that owns the diagonal element (i,i), and 0 on all other nodes.
 *
 *******************************************************************************
 *
 * @sa dplasma_hsdpotrf_New
 * @sa dplasma_hsdpotrf_Destruct
 * @sa dplasma_cpotrf
 * @sa dplasma_hsdpotrf
 * @sa dplasma_spotrf
 *
 ******************************************************************************/
int
dplasma_hsdpotrf_4taskclass( parsec_context_t *parsec,
                PLASMA_enum uplo,
                parsec_tiled_matrix_dc_t *A,
                int lookahead, int send_less, int tensor_gemm,
                int band_size_double, int band_size_single, int HNB )
{
    parsec_taskpool_t *parsec_dpotrf = NULL;
    int info = 0, ginfo = 0 ;
    
    int* gpotrf_info = calloc(A->lmt, sizeof(int));
    potrf_info = calloc(A->lmt, sizeof(int));

    assert( band_size_double <= band_size_single );

    parsec_dpotrf = dplasma_hsdpotrf_4taskclass_New( parsec, uplo, A, &info, lookahead, send_less,
		    tensor_gemm, band_size_double, band_size_single, HNB );

    if ( parsec_dpotrf != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_dpotrf);
        parsec_context_start(parsec);
        parsec_context_wait(parsec);
        dplasma_hsdpotrf_4taskclass_Destruct( parsec_dpotrf );

        if( HNB > 0 && HNB < A->nb ) {
            parsec_taskpool_sync_ids(); /* recursive DAGs are not synchronous on ids */
        }
    }

    /* This covers both cases when we have not compiled with MPI, or we don't need to do the reduce */
    ginfo = info;
    memcpy((void*)gpotrf_info, (void*)potrf_info, A->lmt*sizeof(int));

#if defined(PARSEC_HAVE_MPI)
    /* If we don't need to reduce, don't do it, this way we don't require MPI to be initialized */
    if( A->super.nodes > 1 ) {
        MPI_Allreduce( &info, &ginfo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce( potrf_info, gpotrf_info, A->lmt, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
#endif

    for(int i = 0; i< A->lmt; i++) {
	if(gpotrf_info[i] > 0) {
	    ginfo = gpotrf_info[i] + i*A->mb;
	    break;
	}
    }

    free(potrf_info);
    free(gpotrf_info);

    return ginfo;
}
