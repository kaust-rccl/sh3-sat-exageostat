/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_misc.h
 *
 * Header file of  auxiliary functions that are needed by ExaGeoStat.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/
#ifndef _MLE_MISC_H_
#define _MLE_MISC_H_
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "common.h"
#include "descriptor.h"
#include "chameleon/morse_runtime.h"
#include <nlopt.h>
#include <math.h>
#include <morse.h>
#include <lapacke.h>
#include <starpu.h>
#include <cblas.h>
#include "../../include/flops.h"
#include <starpu_profiling.h>
#if defined(CHAMELEON_USE_MPI)
#include <starpu_mpi.h>
#else
#endif
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>
#endif
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include "../../include/morse_starpu.h"
/** ****************************************************************************
 * PI value
 **/
#define PI (3.141592653589793)
/** ****************************************************************************
 * The radius of the  eartch (used by Great Circle Distance (GCD)
 **/
#define earthRadiusKm 6371.0
/** ****************************************************************************
 * Start timing macro
 **/
#define START_TIMING(_t)     _t =- cWtime();
/** ****************************************************************************
 * Stop timing macro
 **/
#define STOP_TIMING(_t)      _t += cWtime();
/*******************************************************************************/
/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */


inline static void *chameleon_getaddr_null(const MORSE_desc_t *A, int m, int n)
{
	(void)A; (void)m; (void)n;
	return NULL;
}


/**
 *  Internal function to return the leading dimension of element A(m,*) with m,n = block indices
 */
inline static int chameleon_getblkldd_ccrb(const MORSE_desc_t *A, int m)
{
	int mm = m + A->i / A->mb;
	return ( ((mm+1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

/**
 *  Internal function to return MPI rank of element A(m,n) with m,n = block indices
 */
inline static int chameleon_getrankof_2d(const MORSE_desc_t *desc, int m, int n)
{

	if(m>=n)
		return (m % desc->p) * desc->q + (n % desc->q);
	else
		return (n % desc->p) * desc->q + (m % desc->q);
}


/** ****************************************************************************
 * Allocate matrix in different modes
 **/
#define EXAGEOSTAT_ALLOCATE_MATRIX_TILE(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_ , _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
	if (data->ooc && _memspace_ == NULL && _mb_ != 1  && _nb_ !=1)                                                         \
MORSE_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
		_p_, _q_);             \
else                                                            \
MORSE_Desc_Create(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
		_p_, _q_); 	\




#define EXAGEOSTAT_ALLOCATE_FULL_MATRIX_TILE(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_ , _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
	if (data->ooc && _memspace_ == NULL && _mb_ != 1  && _nb_ !=1)                                                         \
MORSE_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
		_p_, _q_);             \
else                                                            \
MORSE_Desc_Create_User(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
		_p_, _q_, morse_getaddr_null, morse_getblkldd_ccrb, chameleon_getrankof_2d );	\



#define EXAGEOSTAT_ALLOCATE_DIAG_MATRIX_TILE(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_ , _lda_, _n_, _smb_, _snb_, _m_, _n2_, _p_, _q_) \
	if (data->ooc && _memspace_ == NULL && _mb_ != 1  && _nb_ !=1)                                                         \
MORSE_Desc_Create_OOC(_desc_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_, _smb_, _snb_, _m_, _n2_, \
		_p_, _q_);             \
else                                                            \
MORSE_Desc_Create_User(_desc_, _memspace_, _type2_, _mb_, _nb_, _mbXnb_, _lda_, _n_,_smb_, _snb_ , _m_, _n2_, \
		_p_, _q_, morse_getaddr_null, morse_getblkldd_ccrb, morse_getrankof_2d_diag );   \

/** ****************************************************************************
 *  Structure for  identifies two dimensions struct with two vectors (x and y).
 **/
typedef struct {
	double *x;				///< Values in X dimension.
	double *y;				///< Values in Y dimension.
        double *z;                              ///< Values in Z dimension.
} location;



/** ****************************************************************************
 *  Structure for reading real datasets through STARS-H
 **/
typedef struct {
	double xy;                              ///< Locations (x, y)
	double z;                               ///< Measurements. 
} sdata;

/** ****************************************************************************
 *  Structure for  uniquely identifies three different parameters for accuracy
 *  measurments accuracyDenseAppDiff, normDenseAppDiff, and normA.
 **/
/*typedef struct{
  double accuracyDenseAppDiff; 		///< Accuracy different between dense format and approx format.
  double normDenseAppDiff;		///< Norm difference between dense format and approx format.
  double normA;   			///< normA.
  } acc_struct;
 */

/** ****************************************************************************
 *  Structure for  uniquely identifies different variables that are needed
 *  by  EXAGEOSTAT to  ease arguments pass.
 **/
typedef struct
{
	double variance;		///< Variance parameter. 
	char *computation;		///< Exact or approx computation.
	char *c_fun;	                ///< Matern or pow-exp kernels.	
	int test;			///< Synthetic or real dataset execution.
	int async;			///< Running mode: synchronous or asynchronous.
	int iter_count;			///< Number of iterations to converge.
	location  l1;			///< 2D locations for the first dataset.
	location lmiss;			///< 2D locations for the missing data (prediction stage).
	location lobs;			///< 2D locations for the observed data (prediction stage).
	void *descC;			///< Covariance matrix C descriptor.
	void *descZ;			///< Measurements Z descriptor.
	double *Adense;			///< Dense matrix descriptor in the case of approximation mode - accuracy check.	
	double *Adense2;		///< Dense matrix descriptor2 in the case of approximation mode - accuracy check.
	void *descZcpy;			///< A copy of Measurements Z descriptor.
	void *descdet;			///< Determinant descriptor.
	void *descproduct;		///< Dot product descriptor.
	void *descZmiss;		///< Missing measurements descriptor.
    void *descC11;          ///< Covariance Matrix C11 descriptor.
    void *descC21;          ///< Covariance Matrix C21 descriptor.
	void *descC12;			///< Covariance Matrix C12 descriptor.
	void *descC22;			///< Covariance Matrix C22 descriptor.
	void *descZactual;		///< Actual Measurements Z descriptor.
	void *descZobs;			///< observed Measurements Z descriptor.
	void *descmse;         		///< Mean Square Error (MSE) descriptor.
    void *descZtrace;            ///< Trace descriptor.
	void *sequence;			///< MORSE sequence.
	void *request;			///< MORSE request.
	int verbose;			///< Verbose indicator.
	int check;			///< Check indicator -- approximation mode.
	int log;			///< Log files generation indicator, 0--> no, 1--> yes.
	double avg_exec_time_per_iter;	///< Avergae execution time per iteration (only used in verbose mode).
	double total_exec_time;		///< Total execution time (only used in verbose mode).
	double avg_flops_per_iter;	///< Avergae flops per iteration (only used in verbose mode).
	double final_loglik;		///< Final log likelihood value.
	char *locsFPath;		///< Locations file path -- in the case of real dataset (real mode). 
	char *obsFPath;			///< Observations file path --  in the case of real dataset (real mode).
	char *actualZFPath;		///< Actual observations file path -- in the case of prediction.
	char *actualZLocFPath;		///< Actial locations file path -- in the case of prediction.
	double det;			///< determinant value.
	double  dotp;			///< double dot product value.
	float sdotp;			///< single dot product value.
	double mserror;			///< Mean Square Error (MSE) value.
	char *dm;			///< Distance metric to be used ed->Euclidian Distance -- gcd->Great Circle Distance.
	int diag_thick;		        ///< The thick of used diagonal in the case of diagonal approximation approach.
	char *nFileLog;			///< log file name (only used if log -->1).
	FILE *pFileLog;			///< log file path (only used if log -->1).
	int hicma_maxrank;		///< Max Rank in the case of LR-HiCMA approx
	int hicma_data_type;            ///< To define the type of the problem to HiCMA (HICMA_STARSH_PROB_GEOSTAT (Synthetic) or HICMA_STARSH_PROB_GEOSTAT_POINT (real))
	void *hicma_descC;		///< HiCMA descC descriptor (for accuracy check).
	void *hicma_descZ;		///< HiCMA descZ descriptor.
	void *hicma_descCD;		///< HiCMA descCD descriptor.
	void *hicma_descCUV;		///< HiCMA descCUV descriptor.
	void *hicma_descCrk;		///< HiCMA descCrk descriptor.
	void *hicma_descZcpy;           ///< A copy of Measurements Z descriptor.
	void *hicma_descdet;            ///< Determinant descriptor.
	void *hicma_descproduct;        ///< Dot product descriptor.
	void *hicma_descC12D;		///< HiCMA descCD descriptor.	
	void *hicma_descC12UV;		///< HiCMA descCUV descriptor.
	void *hicma_descC12rk;		///< HiCMA descCrk descriptor.
	void *hicma_descC22D;		///< HiCMA descCD descriptor.	
	void *hicma_descC22UV;		///< HiCMA descCUV descriptor.
	void *hicma_descC22rk;		///< HiCMA descCrk descriptor.
	double hicma_acc;		///< Accuracy in the case of LR-HiCMA approx.
	void *hsequence;                ///< HiCMA sequence.
	void *hrequest;                 ///< HiCMA request.
	int opt_tol; 	                ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
	int opt_tol_mix; 	        ///< (Mix-precision) The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
	int opt_max_iters;	        ///< Maximum number of mle iterations.
	int ooc;                        ///< Support Out-Of-Core execution, 0-->no, 1-->yes.
	int precision;                  ///< (0)Double, (1)Single, and (2)Mixed.

        int band_size_double;           ///< Band size double 
        int band_size_single;           ///< Band size single
        int lookahead;                  ///< Lookahead size 
        int sendless;                   ///< 0 depens on receiver side; 1 depends on sender side 
        int tensor_gemm;                ///< 1 AB16F_C32F_OP32F, 2 AB16F_C16F_OP32F, 3 AB16F_C16F_OP16F
        void *parsec;                   ///< PaRSEC context 
        int uplo;                       ///< Matrix type 
        int nodes;                      ///< Number of nodes 
        int rank;                       ///< process rank 
        int mse_flag;                   ///< Flag of mse 
        int HNB;                        ///< smallnb 
	nlopt_opt *opt;                 ///< nlopt_opt

	//Mixed Precision
	void *desctemp;                 ///< Temporary descriptor for mixed precision Cholesky factorization.


	void *desctemp22;               ///< Temporary descriptor for mixed precision Cholesky factorization.

        //MLOE and MMOM
        void *desck_t;
        void *desck_a;
        void *desck_ttmp;
        void *desck_atmp;
        void *descK_ttmp;
        void *descK_t;
        void *descK_a;
        void *descexpr1;
        void *descexpr2;
        void *descexpr3;
        void *descexpr4;
        void *desc_mloe_mmom;
        double expr1;                 
        double expr2;                 
        double expr3;                 
        double expr4;                 
	int mloe_mmom;

	double trace_sum;
	double trace_sd;
	double trace_mean;

	char* recovery_file;            
	char* checkpoint_file;

} MLE_data;


/** ****************************************************************************
 * Verbose Macro.
 **/
#define VERBOSE(str)	\
	if (data->verbose == 1 && MORSE_My_Mpi_Rank() == 0){	\
		fprintf(stdout, "%s", str);	\
		fflush(stdout);\
	}

/** ****************************************************************************
 *  Success Macro.
 **/
#define SUCCESS(success, str) \
	printf("success = %d\n", success); \
	if (success != MORSE_SUCCESS){ \
		fprintf(stdout, "%s", str);\
		fflush(stdout);\
		exit(EXIT_FAILURE);\
	}



void pick_random_points(MLE_data *data, double *Zobs, double *Zactual,
		int nZmiss, int nZobs, int N);

void generate_interior_points(MLE_data *data, double *Zobs, double *Zactual,
		int nZmiss, int nZobs, int N);

void split_data(MLE_data * data, location *locations, double *Z,   double *Zactual,
		int *N, int nZmiss);

void init_data_values(MLE_data *data);

int locations_obs_zsort_inplace(int n, location *locations, double *z);

void zsort_locations_obs(int n, location *locations, double *z);

double uniform_distribution(double rangeLow, double rangeHigh);

void print_dmatrix(char* desc, int m, int n,
		double* a, int lda);

void print_diagonal(char* desc, int m, double* a,
		int lda);

void print_smatrix(char* desc, int m, int n,
		float* a, int lda);

location* GenerateXYLoc(int n, int seed);

location* GenerateXYZLoc(int n, int seed);

int countlines(char *filename);


void checkpointing(char *path, int iter_count, double* theta, double loglik, int num_params);

bool recover(char *path, int iter_count, double* theta, double* loglik, int num_params);

void write_to_file(char * path, int matrix_size,int ncores,
		int tile_size, int test, char *computation,
		int async, char *obsFPath,double total_exec_time,
		double avg_exec_time_per_iter, double avg_flops_per_iter,
		int p_grid, int q_grid, double final_loglik, int n);

void theta_parser2(double * theta_vec, char * kern, int num_params);

void write_vectors(double * zvec, MLE_data * data, int n);

void write_to_thetafile(char * path, double theta0, double theta1,
		double theta2, double theta3, int n, double time_per_iter,
		int total_no_iters,
		double prediction_error);

//void readObsFile(char *obsfile, int n, double * streamdata);

void shuffle(double *array, location * locations, size_t n);

void theta_parser(double *initial_theta, double *target_theta,
		double *starting_theta, char *ikernel, char *kernel,
		double *lb, double *up, int test, int num_params);

void init_optimizer(nlopt_opt * opt, double *lb, double *up,
		double tol);

void print_summary(int test, int N, int ncores, int gpus,
		int ts, int lts, char *computation,
		int zvecs, int p_grid, int q_grid, int precision);

void print_result(MLE_data *data, double *starting_theta, int N,
		int zvecs, int ncores, int ts, int test,
		double *initial_theta, char *computation,
		int p_grid, int q_grid, double final_loglik,
		double prediction_error);

double cWtime(void);

void readlocfile(char* loc_file, int n,  location* l1);

void write_prediction_result(char *path, int matrix_size, int no_missing,
		double MSE, double solve_time, double flops);

int doesFileExist(const char *filename);

void init_log (MLE_data * data);

void finalize_log (MLE_data * data);
//acc_struct check_acc(MLE_data * HICMA_data, int n, int ts);
double core_matern_vector (double x0, double y0, double x1, double y1, double *localtheta, int distance_metric);
#endif
