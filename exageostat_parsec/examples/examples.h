/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE.h
 *
 * Header file of ExaGeoStat main functions.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/    

#ifndef _EXAGEOSTAT_EXAMPLES_H_
#define _EXAGEOSTAT_EXAMPLES_H_

#include <stdlib.h>
#include <argp.h>
#include "MLE_misc.h"

/** ****************************************************************************
 *  EXAGEOSTAT arguments uniquely identifies a set of input arguments
 **/
typedef struct
{
	int     test;                  ///< The -t flag -- use testing mode (Synthetic Dataset Generator).
	int     check ;                ///< The -c flag -- use check mode (in the case of approximation only).
	char    *zvecs ;               ///< The number of Z vectors to be tested.
	int     verbose;               ///< The -v flag -- verbose mode.
	char    *ncores;               ///< Number of CPU cores.
	char    *gpus;                 ///< Number of GPUs.
	char    *N;                    ///< Problem size -- Only in the case of testing mode.
	char    *p;                    ///< P in distributed grid.
	char    *q;                    ///< q in distributed grid.
	char    *lts;                  ///< HiCMA tile size.
	char    *dts;                  ///< Chameleon tile size.
	char    *kernel;               ///< Target theta vector output (ex., ?:?:?).
	char    *ikernel;              ///< Initial theta vector -- Only in the case of testing mode (ex., 1:0.1:0.5).
	char    *olb;                  ///< Optimizer lower bounds vector (ex., 0.1:0.01:0.01).
	char    *oub;                  ///< Optimizer upper bounds vector (ex., 5:5:5).
	char    *computation;          ///< Approx or exact.
        char    *c_fun;	               ///< matern or power-exp  --- covariance function to be used.
	int     async;                 ///< 0--> tile  1--> tile_async.The --async flag.
	char    *locs_file;            ///< Locations file path -- in the case of real dataset (real mode).
	char    *obs_dir;              ///< Observations file path --  in the case of real dataset (real mode).
	char    *actualZ_file;         ///< Actual observations file path -- in the case of prediction.
	char    *actualZloc_file;      ///< Actial locations file path -- in the case of prediction.
	char    *predict;              ///< Number of missing values  -- in the case of testing mode.
	char    *dm;                   ///< Distance metric to be used ed->Euclidian Distance -- gcd->Greate Circle Distance.
	char    *diag_thick;           ///< The thick of used diagonal in the case of diagonal approximation approch.
	int     log;                   ///< Generate log files -- 0->do not store generated data, 1-->store generated data.
	char    *maxrank;              ///< Max Rank in the case of LR-HiCMA approx.
	char    *acc;                  ///< Accuracy in the case of LR-HiCMA approx.
	int     profile;               ///< profiling the performance of exageostat using FxT.
	char    *opt_tol;              ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
	char    *opt_tol_mix;          ///< (Mix-precision) The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
	char    *opt_max_iters;        ///< Maximum number of mle iterations.
	int     ooc;                   ///< Support Out-Of-Core (OOC) -- 0->do not support, 1-->support.
	int     band_size_double;      ///< Band size double 
 	int     band_size_single;      ///< Band size single
 	int     lookahead;             ///< Lookahead size 
 	int     sendless;              ///< 0 depens on receiver side; 1 depends on sender side 
 	int     tensor_gemm;           ///< 1 AB16F_C32F_OP32F, 2 AB16F_C16F_OP32F, 3 AB16F_C16F_OP16F 
 	int     kernel_type;           ///< 1 mix; 2 double; 3 single 
 	int     HNB;                   ///< Hierarchcal tile size 
 	int     nb_runs;               ///< Number of runs 
 	char    *precision;            ///< Double, Single, and Mixed.
        char    *checkpoint_file;      ///< checkpoint file path.
        char    *recovery_file;        ///< Recovery file path.
        int     mloe_mmom;             ///< Use MLOE and MMOM 
        char    *dim;                  ///< 2D or 3D
} arguments;

void check_args(arguments *arg_values);
void set_args_default( arguments *arg_values);
void init(int *test, int *N,  int *ncores,
		int *gpus, int *p_grid, int *q_grid,
		int *zvecs, int *dts,int *lts,
		int *nZmiss, int *log,  double *initial_theta,
		double *starting_theta, double *target_theta, double *lb,
		double *ub, MLE_data *data, arguments *arguments);
static struct argp_option options[] =
{
	{"test",'t', 0, 0, "Execute in test mode"},
	{"check", 'c', 0, 0, "Produce check output"},
	{"zvecs",   'z', "ZVECS", 0, "number of Z vectors to be tested"},
	{"verbose", 'v', 0, 0, "Produce verbose output"},
	{"ncores",   'n', "NCORES", 0, "Number of cores"},
	{"gpus",   'g', "GPUS", 0, "Number of gpus"},
	{"p",   'd', "P", 0, "p in distributed system"},
	{"q",   'q', "Q", 0, "q in distributed system"},
	{"N",   's', "MATRIX_SIZE", 0, "Synthetic Matrix Size"},
	{"lts",   'L', "HICMA_TILE_SIZE", 0, "Number of tiles in TLR"},
	{"dts",   'T', "Chameleon_TILE_SIZE", 0, "Number of Tiles in dense"},
	{"kernel",   'k', "KERNEL", 0, "Computation model"},
	{"ikernel",   'i', "IKERNEL", 0, "Initial theta(s) used for testing case generation"},
	{"olb",   'y', "LB", 0, "Optimizer Lower Bounds"},
	{"oub",   'j', "UB", 0, "Optimizer Upper Bounds"},
	{"computation",   'b', "COMPUTATION", 0, "Exact or Approx"},
        {"c_fun",   'F', "C_FUN", 0, "matern or power-exp"},
	{"async", 'a', 0, 0, "Asynchronous"},
	{"locs_file",   'l', "LOCATIONS_FILE", 0, "Read Locations from this Location File"},
	{"obs_dir",  'o', "OBSERVATIONS_DIRECTORY", 0, "Read Observations from this directory path"},
	{"actualZ_file",  'f', "ACTUALZ_FILE", 0, "Read actual Z from this observation file"},
	{"actualZloc_file",  'm', "ACTUALZLOC_FILE", 0, "Read actual Z locations from this location file"},
	{"predict",  'p', "PREDICT", 0, "Number of Missing Values"},
	{"dm",  'h', "DISTANCE_METRIC", 0, "Distance Metric"},
	{"diag_thick",  'D', "DIAG_THICK", 0, "Diagonal Thick"},
	{"log",  'r', 0, 0, "Store Generated Data (Test mode only)"},
	{"maxrank",   'x', "MAXRANK", 0, "HiCMA Max RANK"},
	{"acc",   'u', "ACC", 0, "HiCMA Accuracy"},
	{"profile", 'w', 0, 0, "Performance profiling"},
	{"opt_tol",   'O', "OPTIMIZATION_TOLERANCE", 0, "Optimization tolerance"},
	{"opt_tol_mix",   'W', "OPTIMIZATION_TOLERANCE_MIX", 0, "Optimization tolerance mix"},
	{"opt_iters", 'I', "OPTIMIZATION_MAX_ITERS", 0, "Optimization maximum iterations"},
	{"precision",   'Z', "PRECISION", 0, "(0)Double, (1)Single, or (2)Mixed"},
	{"ooc",  'C', 0, 0, "Support Out-Of-Core (OOC) execution"},
	{"band_size_double",  'E', "BAND_SIZE_DOUBLE", 0, "Band size double"},
 	{"band_size_single",  'e', "BAND_SIZE_SINGLE", 0, "Band size single"},
 	{"lookahead",  'G', "LOOKAHEAD", 0, "Lookahead size"},
 	{"sendless",  'S', "SENDLESEE", 0, "send less"},
 	{"tensor_gemm",  'U', "TENSOR_GEMM", 0, "tensor_gemm"},
 	{"kernel_type",  'K', "KERNEL_TYPE", 0, "kernel_type"},
 	{"HNB",  'J', "HNB", 0, "HNB"},
 	{"nb_runs",  'Q', "NB_RUNS", 0, "nb_runs"},
        {"checkpoint_file",  'H', "CHECKPOINT_FILE", 0, "Checkpoint parameters to"},
        {"recovery_file",  'Y', "RECOVERY_FILE", 0, "Recover parameters from"},
        {"mloe_mmom",  'X', 0, 0, "use_mloe_mmom"},
        {"dim",  'M', "DIMENSION", 0, "Dimension"},
	{0}
};

	static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
	arguments *arguments = state->input;

	switch (key)
	{
		case 'Z':
			arguments->precision = arg;
			break;
		case 't':
			arguments->test = 1;
			break;
		case 'c':
			arguments->check = 1;
			break;
		case 'v':
			arguments->verbose = 1;
			break;
		case 'k':
			arguments->kernel = arg;
			break;
		case 'i':
			arguments->ikernel = arg;
			break;
		case 'y':
			arguments->olb = arg;
			break;
		case 'j':
			arguments->oub = arg;
			break;
		case 'z':
			arguments->zvecs = arg;
			break;
		case 'n':
			arguments->ncores = arg;  //non-optional;
			break;
		case 'g':
			arguments->gpus = arg;  //non-optional;
			break;

		case 's':
			arguments->N = arg;
			break;

		case 'd':
			arguments->p = arg;  //non-optional;
			break;

		case 'q':
			arguments->q = arg;
			break;
		case 'T':
			arguments->dts = arg;  //non-optional
			break;
		case 'L':
			arguments->lts = arg;  //optional
			break;
		case 'b':
			arguments->computation = arg;
			break;
                case 'F':
                        arguments->c_fun = arg;
                        break;
		case 'a':
			arguments->async = 1;
			break;

		case 'l':
			arguments->locs_file = arg;
			break;
		case 'o':
			arguments->obs_dir = arg;
			break;
		case 'p':
			arguments->predict = arg;
			break;
		case 'f':
			arguments->actualZ_file = arg;
			break;
		case 'm':
			arguments->actualZloc_file = arg;
			break;
		case 'h':
			arguments->dm = arg;
			break;
		case 'D':
			arguments->diag_thick = arg;
			break;
		case 'r':
			arguments->log= 1;
			break;
		case 'x':
			arguments->maxrank = arg;
			break;
		case 'u':
			arguments->acc = arg;
			break;
		case 'w':
			arguments->profile= 1;
			break;
		case 'O':
			arguments->opt_tol = arg;
			break;
		case 'W':
			arguments->opt_tol_mix = arg;
			break;
		case 'I':
			arguments->opt_max_iters = arg;
			break;
		case 'C':
			arguments->ooc= 1;
			break;
                case 'E':
                        arguments->band_size_double = atoi(arg);
                        break;
                case 'e':
                        arguments->band_size_single = atoi(arg);
                        break;
                case 'G':
                        arguments->lookahead = atoi(arg);
                        break;
                case 'S':
                        arguments->sendless = atoi(arg);
                        break;
                case 'U':
                        arguments->tensor_gemm = atoi(arg);
                        break;
                case 'K':
                        arguments->kernel_type = atoi(arg);
                        break;
                case 'J':
                        arguments->HNB = atoi(arg);
                        break;
                case 'Q':
                        arguments->nb_runs = atoi(arg);
                        break;
                case 'H':
                        arguments->checkpoint_file = arg;
                        break;
                case 'Y':
                         arguments->recovery_file = arg;
                        break;
                case 'X':
                        arguments->mloe_mmom= 1;
                        break;
                case 'M':
                        arguments->dim = arg;
                        break;
		default:
			return ARGP_ERR_UNKNOWN;
	}


	return 0;
}

static char args_doc[] = "";

static char doc[] =
"ExaGeoStat -- A unified geospatial statistic framework to evaluate Maximum Likelihood function using both real and synthetic dataset on CHAMELEON (Dense/DST/mixed-precision) - HiCMA (TLR approximation) - )";

static struct argp argp = {options, parse_opt, args_doc, doc};

#endif
