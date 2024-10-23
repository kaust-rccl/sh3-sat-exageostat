#ifndef MIX_PRECISION_CHOLESKY_H
#define MIX_PRECISION_CHOLESKY_H

/* Exageostat headers */
#define _MORSE_CORE_DBLAS_H_
#include "examples/examples.h"
#include "misc/include/MLE_misc.h"
#include "src/include/MLE.h"
#include "exageostat_exact/core/include/exageostatcore.h"

/* Single cholesky for optimization */
#define SINGLE_PO 0

/* Count value that exceed half-precion limit */
#define COUNT_VALUE_HALF 0

/* sdpotrf on GPU */
//#define EXAGEOSTAT_USE_CUDA
#undef EXAGEOSTAT_USE_CUDA

/* If swith back to DP after certain threshold (tol -1) */
#define SWITCH_TO_DP 0

/* Not write to file */
#define EXAGEOSTAT_ON_SUMMI

/* Adpative choose kernel when losing spd to go back to DP10% SP90% or DP100% */
#define ADAPTIVE_PO 1

/* Allocate GPU buffer once */
#define GPU_BUFFER_ONCE 1

/* Kernels from dplasm not on GPU */
#define NO_GPU_DPLASMA 1

/* trsm convert half precision */ 
#define TRSM_CONVERT_HALF 0

/* merge potrf_convert to potrf and trsm_convet to trsm */ 
/* Need TRSM_CONVERT_HALF 1 */
/* Need COLLECTIVE OFF in parsec */
#define FOUR_TASKCLASS_COLLECTIVE_OFF 0

/* matrix generation used in performance test */ 
#define MATRIX_GENERATION 0

/* Use own matrix generation */
#define OWN_GENERATION 0

/* Use own my_trmv */ 
#define OWN_TRMV 0

/* Use own converor which used for debug */
#define OWN_CONVERTOR 0

/* NOT use own dot */ 
#define OWN_DOT 0

/* Debug info print */
#define DEBUG_INFO 0

#undef VERBOSE
#define VERBOSE(str)    \
        if (data->verbose == 1 && data->rank == 0){    \
                fprintf(stderr, "%s", str);     \
        }

#undef BLKLDD
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic_band.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic_band.h"
#include "parsec/utils/mca_param.h"
#include "parsec/private_mempool.h"
#include "parsec/runtime.h"
#include "parsec/data_internal.h"
#include "parsec/execution_stream.h"
#include "parsec/data_dist/matrix/matrix.h"
#include "tests/flops.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* dplasma */
#include "dplasma.h"

#if defined(EXAGEOSTAT_USE_CUDA)
#include "parsec/mca/device/cuda/device_cuda.h"
#include "parsec/mca/device/cuda/device_cuda_internal.h"
#include "parsec/utils/zone_malloc.h"
//#define CUBLAS_H_
#include <cublas.h>
#include <cuda_fp16.h>
#include <cuda_runtime.h>
//#define CUBLAS_V2_H_
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Formula for flops */
typedef double DagDouble_t;
#define PASTE_CODE_FLOPS( FORMULA, PARAMS ) \
  double gflops = -1.0, flops = FORMULA PARAMS;

/* Structure to init random matrix */
struct matrix_init_args_s {
    double      bump;
    unsigned long long int seed;
    int band_size_double;
};

typedef struct matrix_init_args_s matrix_init_args_t;


int dplasma_hsdpotrf( parsec_context_t *parsec,
                PLASMA_enum uplo,
                parsec_tiled_matrix_dc_t *A,
                int lookahead, int send_less, int tensor_gemm,
                int band_size_double, int band_size_single, int HNB );

int parsec_spotrf( parsec_context_t *parsec,
               PLASMA_enum uplo,
               parsec_tiled_matrix_dc_t *A);

int dplasma_hsdpotrf_trsm_convert( parsec_context_t *parsec,
                PLASMA_enum uplo,
                parsec_tiled_matrix_dc_t *A,
                int lookahead, int send_less, int tensor_gemm,
                int band_size_double, int band_size_single, int HNB );

int dplasma_hsdpotrf_4taskclass( parsec_context_t *parsec,
                PLASMA_enum uplo,
                parsec_tiled_matrix_dc_t *A,
                int lookahead, int send_less, int tensor_gemm,
                int band_size_double, int band_size_single, int HNB );

int parsec_my_dtrmv(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *A,
                       parsec_tiled_matrix_dc_t *B,
                       parsec_tiled_matrix_dc_t *C);

int parsec_band_convert_d2s(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcY, int band_size);

int parsec_band_convert_s2d(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcY, int band_size);

int parsec_dmatrix_generation(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, location *l1, location *l2,
                         double *theta , char *dm, char *c_fun, int band_size_double);

void parsec_dmatrix_generation_Destruct(parsec_taskpool_t *taskpool);

parsec_taskpool_t* parsec_dmatrix_generation_New(parsec_context_t *parsec,
                parsec_tiled_matrix_dc_t *dcA, location *l1, location *l2,
                double *theta , char *dm, char *c_fun, int band_size_double);

int parsec_smatrix_generation(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, location *l1, location *l2,
                         double *theta , char *dm, char *c_fun, int band_size_double);

int parsec_sdmatrix_generation(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, location *l1, location *l2,
                         double *theta , char *dm, char *c_fun, int band_size_double);

int parsec_dZ_generation(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, double *r);

int parsec_sZ_generation(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, float *r);

int parsec_dmatrix_det(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA,
                       parsec_tiled_matrix_dc_t *dcdet);

int parsec_dmatrix_scale(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA);

int parsec_dmatrix_shift(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, int multiple);

double MLE_alg_parsec(unsigned n, const double * theta, double * grad, void * data);

double parsec_dmle_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data);

double parsec_dmle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n);

int parsec_get_zobs(parsec_context_t *parsec,
                parsec_tiled_matrix_dc_t *dcA, double *Z, int N);

void parsec_prediction_init(MLE_data *data, int nZmiss, int nZobs, int ts, int p_grid, int q_grid, int mse_flag);

int parsec_Lapack_to_Tile(parsec_context_t *parsec,
                parsec_tiled_matrix_dc_t *dcA, double *Zobs, int nZobs);

void parsec_pick_random_points(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N);

int parsec_mse_calculation(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcZpre,
                       parsec_tiled_matrix_dc_t *dcZmiss,
                       parsec_tiled_matrix_dc_t *dcZerror);

void parsec_prediction_finalize(MLE_data *data);

int parsec_dmatrix_init( parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, double alpha, double beta );

double parsec_dmatrix_sum(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA );

double parsec_dZ_sum(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA );

int parsec_dmatrix_dot(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA,
                       parsec_tiled_matrix_dc_t *dcB,
                       parsec_tiled_matrix_dc_t *dcdot);

int parsec_sdmatrix_init_operator( parsec_execution_stream_t *es,
                         const parsec_tiled_matrix_dc_t *descA,
                         void *_A,
                         PLASMA_enum uplo, int m, int n,
                         void *op_data );

int parsec_my_warmup(parsec_context_t *parsec,
                     parsec_tiled_matrix_dc_t *A);

void parsec_mloe_mmom_init( MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid );

void parsec_band_allocate(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA, int band_size_double);

void parsec_dmatrix_avg(parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *dcA);

int parsec_dmatrix_set_diagonal(parsec_context_t *parsec,
                         parsec_tiled_matrix_dc_t *dcA, double noise);

location* parsec_readLocsFile(char *locs_file, int n, int dts);

location* parsec_readLocsFile3d(char *locs_file, int n, int dts);

location* parsec_GenerateXYLoc(int n, int seed, int dts);

location* parsec_GenerateXYZLoc(int n, int seed, int dts);

int parsec_MLE_dtrace(parsec_context_t *parsec,
                      parsec_tiled_matrix_dc_t *dcC11,
                      parsec_tiled_matrix_dc_t *dcdet,
                      parsec_tiled_matrix_dc_t *dcZtrace);
int my_dtrsm( parsec_context_t *parsec,
               PLASMA_enum side,  PLASMA_enum uplo,
               PLASMA_enum trans, PLASMA_enum diag,
               double alpha,
               const parsec_tiled_matrix_dc_t *A,
               parsec_tiled_matrix_dc_t *B);

void my_dtrsm_Destruct( parsec_taskpool_t *tp );

parsec_taskpool_t*
my_dtrsm_New( PLASMA_enum side,  PLASMA_enum uplo,
                   PLASMA_enum trans, PLASMA_enum diag,
                   double alpha,
                   const parsec_tiled_matrix_dc_t *A,
                   parsec_tiled_matrix_dc_t *B );


#if defined(EXAGEOSTAT_USE_CUDA)
typedef struct parsec_potrf_stream_workspace_s {
    cusolverDnHandle_t handle_cusolver;
    cublasHandle_t handle_cublas;
    void *gpu_buffer;
    int buffer_size;
}parsec_potrf_stream_workspace_t;

typedef struct parsec_potrf_gpu_workspace_s {
    parsec_potrf_stream_workspace_t *stream_workspace;
    parsec_device_cuda_module_t *gpu_device;
}parsec_potrf_gpu_workspace_t;

typedef struct parsec_potrf_workspace_s {
    parsec_potrf_gpu_workspace_t *gpu_workspace;
    int info;
}parsec_potrf_workspace_t;

void float2double_GPU(int nrows, int ncols,
                const float *F, int ldf,
                double *H, int ldh ,
                cudaStream_t stream);

void double2float_GPU(int nrows, int ncols,
                const double *H, int ldh,
                float *F, int ldf,
                cudaStream_t stream);

void float2half_GPU( int nrows, int ncols,
                const float *F, int ldf,
                void *H, int ldh,
                cudaStream_t stream );

void half2float_GPU( int nrows, int ncols,
                const void *H, int ldh,
                float *F, int ldf,
                cudaStream_t stream );

void memcpy_float_GPU( int nrows, int ncols, void *_src, void *_dest, cudaStream_t stream );

void memcpy_half_GPU( int nrows, int ncols, void *_src, void *_dest, cudaStream_t stream );

void matrix_print_double_GPU( int nrows, int ncols, double *A, cudaStream_t stream );

void matrix_print_float_GPU( int nrows, int ncols, float *A, cudaStream_t stream );

void dcmg_array_GPU( double *A, int m, int n, int m0,
                int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
                double *localtheta, int distance_metric, int lda, cudaStream_t stream);

parsec_potrf_stream_workspace_t lookup_gpu_workspace( parsec_device_cuda_module_t *gpu_device,
                                                      parsec_gpu_exec_stream_t *gpu_stream,
                                                      parsec_potrf_workspace_t *ws );

parsec_potrf_workspace_t* workspace_memory_allocate( parsec_potrf_workspace_t *ws );

void workspace_memory_free( parsec_potrf_workspace_t *ws);

void disable_GPU( parsec_taskpool_t * tp );

void disable_CPU( parsec_taskpool_t * tp );

#if GPU_BUFFER_ONCE
void gpu_temporay_buffer_init( int mb, int nb );
void gpu_temporay_buffer_fini( );
#endif

static inline int my_gpu_load( int m, int n, int NT, int nb_cuda_devices, int band_size_double, int band_size_single) {
	/* Old formula */
	//return (m * NT + n) % nb_cuda_devices; 

	/* New formula */
	//return (n * NT + m) % nb_cuda_devices; 
	//if( band_size_double < NT && m - n < band_size_single )
	//if( band_size_double < NT )
	//	return m % nb_cuda_devices; 
	//else
		return (n * NT + m) % nb_cuda_devices; 

}

#endif


void convert_s2d_unary_CPU(float *data, int mb, int nb);

void convert_d2s_unary_CPU(double *data, int mb, int nb);

void convert_s2d_binary_CPU(double *target, float *source, int mb, int nb);

void convert_d2s_binary_CPU(float *target, double *source, int mb, int nb);

void parsec_print_matrix( double *A, int m, int n, int lda );

void parsec_print_matrix_single( float *A, int m, int n, int lda );

#define MY_MAX(num1, num2) ((num1 > num2 ) ? num1 : num2)

#ifdef __cplusplus
}
#endif

#undef MAGMA_COMPLEX

#endif
