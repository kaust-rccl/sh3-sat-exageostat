/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file examples.c
 *
 * 
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/    
#include "examples.h"


void set_args_default(arguments *arg_values)
	//! set default values for input
	/*!  arguments
	 * @param[in] arg_values: user arguments
	 * */
{

	arg_values->test                = 0;
	arg_values->check               = 0;
	arg_values->verbose             = 0;
	arg_values->zvecs               = "1";
	arg_values->computation         = "exact";
        arg_values->computation         = "matern";
	arg_values->async               = 0;
	arg_values->kernel              = "";
	arg_values->ikernel             = "";
	arg_values->ncores              = "1";
	arg_values->gpus                = "0";
	arg_values->p                   = "1";
	arg_values->q                   = "1";
	arg_values->N                   = "0";
	arg_values->lts                 = "0";
	arg_values->dts                 = "0";
	arg_values->locs_file           = "";
	arg_values->obs_dir             = "";
	arg_values->actualZ_file        ="";
	arg_values->actualZloc_file     = "";
	arg_values->predict             = "0";
	arg_values->dm                  = "ed";
	arg_values->diag_thick          = "1";
	arg_values->log                 = 0;
	arg_values->maxrank             = "0";
	arg_values->acc                 = "0";
	arg_values->profile             = 0;
	arg_values->opt_tol             = "5";
	arg_values->opt_tol_mix         = "16";
	arg_values->opt_max_iters       = "-1";
	arg_values->ooc                 = 0;
	arg_values->band_size_double    = 1;
	arg_values->band_size_single    = 0;
	arg_values->lookahead           = 0;
	arg_values->sendless            = 0;
	arg_values->tensor_gemm         = 1;
	arg_values->HNB                 = 0;
	arg_values->nb_runs             = 1;
	arg_values->precision           = "0";
        arg_values->recovery_file       = "";
        arg_values->checkpoint_file     = "";
        arg_values->mloe_mmom           = 0;
        arg_values->dim	                = "2d";
	//      arg_values->precision         = "0";
}

void check_args(arguments *arg_values) {
	//! check  values for input
	/*!  arguments
	 * @param[in] arg_values: user arguments
	 * */
	if (arg_values->test == 0) {
		if (strcmp(arg_values->locs_file,"") == 0 || strcmp(arg_values->obs_dir, "") == 0) {
			fprintf(stdout, "(ExaGeoStat Error MSG): Real running mode requires both locs_file and obs_dir path.... \nExaGeoStat terminated \n\n");
			exit(EXIT_FAILURE);
		}
	}    
	else {
		if (strcmp(arg_values->locs_file,"") != 0 || strcmp(arg_values->obs_dir, "") != 0)
			fprintf(stdout, "(ExaGeoStat Warning MSG): Test running mode does not require locs_file and obs_dir paths, any will be ignored-> continue running...\n");
	}

	//if(strcmp(arg_values->computation,"exact") != 0) {
	//        fprintf(stdout, "(ExaGeoStat Error MSG): approximation is not supported yet, please use only exact computation for now.... \nExaGeoStat terminated \n\n");
	//        exit(EXIT_FAILURE);
	//    }    
//	if (atoi(arg_values->predict) == 0) {
//		if (strcmp(arg_values->actualZ_file,"") != 0 || strcmp(arg_values->actualZloc_file,"") != 0)
//			fprintf(stdout, "(ExaGeoStat Warning MSG): Test running mode does not require actualZ_file and actualZloc_file paths, any will be ignored-> continue running...\n");
//	}    
}

void init(int *test, int *N,  int *ncores,
		int *gpus, int *p_grid, int *q_grid, 
		int *zvecs, int *dts,int *lts, 
		int *nZmiss, int *log,  double *initial_theta, 
		double *starting_theta, double *target_theta, double *lb, 
		double *ub, MLE_data *data, arguments *arguments)
	//! initialize exageostat by setting several
	/*! variables from argument
	 * Returns MLE_data struct.
	 * @param[in] arguments: command line arguments.
	 * @param[out] test: if test=1 ->test running mode  if test=0->real running mode.
	 * @param[out] N: number of spatial locations (problem size).
	 * @param[out] ncores: number of CPU computing units.
	 * @param[out] gpus: number of GPU computing units.
	 * @param[out] p_grid: p_grid in the case of distributed system.
	 * @param[out] q_grid: q_grid in the case of distributed system.
	 * @param[out] zvecs: the number of Z vectors that should be generated in the case of test running mode.
	 * @param[out] dts: dense tile size.
	 * @param[out] lts: TLR tile size.
	 * @param[out] nZmiss: number of unknown observation to be predicted in the case of test running mode.
	 * @param[out] log: determine if log files should be generated or not (log files stored on disk)
	 * @param[out] initial_theta: initial_theta Vector with three parameter (Variance, Range, Smoothness)
	 that is used to to generate the Covariance Matrix and initial Z vector.
	 * @param[out] starting_theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 that is used to to generate the Covariance Matrix of the first MLE iteration.
	 * @param[out] target_theta: target theta Vector with three parameter (Variance, Range, Smoothness) unknown theta parameter should be shown as '?'
	 * @param[out] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
	 * @param[out] ub: optimization upper bounds vector ( ub_1, ub_2, ub_3).
	 * @param[out] data: MLE_data struct with different MLE inputs.
	 * @param[out] arguments: command line arguments.
	 * */
{

	init_data_values(data);
	*test                   = arguments->test;
	*ncores                 = atoi(arguments->ncores);
	*gpus                   = atoi(arguments->gpus);
	*p_grid                 = atoi(arguments->p);
	*q_grid                 = atoi(arguments->q);
	*N                      = atoi( arguments->N);
	*zvecs                  = atoi(arguments->zvecs);
	*dts                    = atoi(arguments->dts);
	*lts                    = atoi(arguments->lts);
	*nZmiss                 = atoi(arguments->predict);
	*log                    = arguments->log;
	data->computation       = arguments->computation;// exact or approx.
        data->c_fun 	        = arguments->c_fun;// matern or pow-exp.
	data->async             = arguments->async; // 0-->tile  1-->tile_async.
	data->locsFPath         = arguments->locs_file;
	data->obsFPath          = arguments->obs_dir;
	data->actualZFPath      = arguments->actualZ_file;
	data->actualZLocFPath   = arguments->actualZloc_file;
	data->dm                = arguments->dm;
	data->diag_thick        = atoi(arguments->diag_thick);
	data->log               = arguments->log;
	data->hicma_maxrank     = atoi(arguments->maxrank);
	data->hicma_acc         = atof(arguments->acc);
	data->check             = arguments->check;
	data->verbose           = arguments->verbose;
	data->opt_tol           = atoi(arguments->opt_tol);
	data->opt_tol_mix       = atoi(arguments->opt_tol_mix);
	data->opt_max_iters     = atoi(arguments->opt_max_iters);
	data->ooc               = arguments->ooc;
	data->precision         = atoi(arguments->precision);// double, sinlge, or mixed.
        data->recovery_file     = arguments->recovery_file;
        data->checkpoint_file   = arguments->checkpoint_file;
	data->mloe_mmom		= arguments->mloe_mmom;
	//        data->precision         = atoi(arguments->precision);// double, sinlge, or mixed.
	//data->l2              = data->l1;
	int num_params = 0;
	if(strcmp (arguments->c_fun, "pow-exp-nuggets") == 0)
		num_params=4;
	else
		num_params=3;

	/* Added by LEI */
        data->band_size_double = arguments->band_size_double;
        data->band_size_single = arguments->band_size_single;
        data->lookahead = arguments->lookahead;
        data->sendless = arguments->sendless;
        data->tensor_gemm = arguments->tensor_gemm;
        data->HNB = arguments->HNB;

	theta_parser2(lb, arguments->olb, num_params);
	theta_parser2(ub, arguments->oub, num_params);
	data->ooc               = arguments->ooc;
	//data->l2              = data->l1;
	int i;
	//multiplicative scale
	//for (i=0;i<num_params;i++)
		//starting_theta[i] = lb[i];

	starting_theta[0] = 1; 
	starting_theta[1] = 0.1; 
	starting_theta[2] = 0.5; 

	if( 4 == num_params )
		starting_theta[3] = 0.15; 

/*
	if(strcmp(data->c_fun, "matern") == 0)
	{
		//fprintf(stderr, "Multip[licative scale (variance profiling is used...!\n");
		lb[0] = ub[0] = 1;
	}

*/
	theta_parser(initial_theta, target_theta, starting_theta, arguments->ikernel, arguments->kernel, lb, ub, *test, num_params);

}
