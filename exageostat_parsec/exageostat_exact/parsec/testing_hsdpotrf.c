/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file zgen_mle_test.c
 *
 * A complete example to test ExaGeoStat supported function (i.e., dataset generator, Maximum Likelihood Function (MLE), Prediction)
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2019-08-06
 *
 **/

/* Parsec header file */
#include "mix_precision_internal.h"
#include "common_timing.h"

double time_elapsed = 0.0;
double sync_time_elapsed = 0.0;

extern int num_iteration_scale_1;
extern int num_iteration_scale_2;
extern int num_iteration_scale_3;

int main(int argc, char **argv) {

	//initialization
	/* Copy from parsec */
	parsec_context_t* parsec;
	PLASMA_enum uplo = PlasmaLower;
	int rank, nodes, ch;
	int pargc = 0;
	char **pargv;
	int jj;
	int info = 0;
	int ret = 0;
	double time_dpotrf = 0.0, time_spotrf = 0.0, time_mix_precision = 0.0;

#if defined(PARSEC_HAVE_MPI)
        {
                int provided;
                MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
        }
        MPI_Comm_size(MPI_COMM_WORLD, &nodes);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        nodes = 1;
        rank = 0;
#endif

	int N, dts, lts, log;
	int i = 0;
	int zvecs = 1, nZmiss = 0, test = 0, gpus = 0;
	int p_grid, q_grid, ncores;
	double  opt_f;
	arguments arguments;
	nlopt_opt opt;
	MLE_data data;
	int seed = 0;
	location *locations = NULL;
	double prediction_error = 0.0;
	int iseed[4]={seed, seed, seed, 1};

        pargc = 0; pargv = NULL;
        for(i = 1; i < argc; i++) {
                if( strcmp(argv[i], "--") == 0 ) {
                        pargc = argc - i;
                        pargv = &argv[i];
                        /* New argc */
                        argc = i;
                        break;
                }
        }

	//Arguments default values
	set_args_default(&arguments);
	argp_parse(&argp, argc, argv, 0, 0, &arguments);
	//check_args(&arguments);

        int num_gpus = atoi( arguments.gpus );
#if defined(EXAGEOSTAT_USE_CUDA)
        extern char **environ;
        if(0 == rank ) printf("num_gpus: %d\n", num_gpus);
        if( num_gpus < 1 && 0 == rank ) {
                fprintf(stderr, "Warnning: if run on GPUs, please set --gpus=value bigger than 0\n");
        }
        parsec_setenv_mca_param( "device_cuda_enabled", arguments.gpus, &environ );
#endif

	ncores = atoi(arguments.ncores);
        /* Initialize PaRSEC */
        parsec = parsec_init(ncores, &pargc, &pargv);

        if( NULL == parsec ) {
                /* Failed to correctly initialize. In a correct scenario report
                 * upstream, but in this particular case bail out.
                 */
                exit(-1);
        }

        int num_params = 3;
        if(strcmp (arguments.c_fun, "pow-exp-nuggets") == 0)
                num_params=4;
        else
                num_params=3;

        double *lb = (double *) malloc(num_params * sizeof(double));
        double *up = (double *) malloc(num_params * sizeof(double));

        double *starting_theta  = (double *) malloc(num_params * sizeof(double));
        double *initial_theta   = (double *) malloc(num_params * sizeof(double));
        double *target_theta    = (double *) malloc(num_params * sizeof(double));

	//MLE_data data initialization
	init(&test, &N,  &ncores,
			&gpus, &p_grid, &q_grid,
			&zvecs, &dts, &lts,
			&nZmiss, &log, initial_theta, 
			starting_theta, target_theta, lb,
			up, &data, &arguments);

	if( 0 == rank ) printf("dts: %d, lts: %d\n", dts, lts);

	// Optimizater initialization
	if( 0 == rank ) printf(" %d - %d \n", data.opt_tol, data.opt_max_iters);	

	opt=nlopt_create( NLOPT_LN_BOBYQA, num_params);
	init_optimizer(&opt, lb, up, pow(10, -1.0 * data.opt_tol));
	nlopt_set_maxeval(opt, data.opt_max_iters);

	//nomral random generation of e -- ei~N(0, 1) to generate Z
	double *Nrand	= (double *) malloc (N * zvecs * sizeof(double));
	if( 0 == rank ) printf("Precision: %d, async %d\n", data.precision, data.async);
	// 3 is not related to num_params (3 means sample from standard normal distribution)
	LAPACKE_dlarnv(3, iseed, N*zvecs, Nrand);

        /* uniform random generation for locations / read locations from disk */
        if(strcmp(arguments.dim, "2d") == 0)
        {
                if( 0 == rank ) fprintf(stderr, "2D example......\n");
                locations = parsec_GenerateXYLoc(N, seed, dts);
        }
        else if(strcmp(arguments.dim, "3d") == 0)
        {
                if( 0 == rank ) fprintf(stderr, "3D example......\n");
                locations = parsec_GenerateXYZLoc(N, seed, dts);
        }
        else
        {
                if( 0 == rank ) fprintf(stderr, "Input dimension is not supported\n\n");
                exit(0);
        }

	data.l1   = *locations;	

	int nZobs = strcmp(data.actualZFPath,"") == 0? (N-nZmiss) : N;

        /* Flops */
	PASTE_CODE_FLOPS(FLOPS_DPOTRF, ((DagDouble_t)N));

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
        /* GPU temporary buffer */
        gpu_temporay_buffer_init( dts, dts );
#endif
#endif

        /* stop gsl error handler */
        gsl_set_error_handler_off();

	//for( i = 0; i < zvecs; i++ ) {
	for( i = zvecs-1; i < zvecs; i++ ) {
		if( 0 == rank ) fprintf(stderr, "zvecs : %d of %d\n", i, zvecs);

		/* Data descritor declaration */
		sym_two_dim_block_cyclic_t dcC;
		two_dim_block_cyclic_t dcZ;
		two_dim_block_cyclic_t dcZcpy;
		two_dim_block_cyclic_t dcdet;
		two_dim_block_cyclic_t dcproduct;

		/* Declare for prediction */
		two_dim_block_cyclic_t descZmiss;
		two_dim_block_cyclic_t descC12;
		sym_two_dim_block_cyclic_t descC22;
		two_dim_block_cyclic_t descZtrace;
		two_dim_block_cyclic_t descC11;
		two_dim_block_cyclic_t descC21;
		two_dim_block_cyclic_t descZobs;
		two_dim_block_cyclic_t descmse;
		two_dim_block_cyclic_t descZactual;

		/* Initialize data structure */
		sym_two_dim_block_cyclic_init(&dcC, matrix_RealDouble,
				nodes, rank, dts, dts, N, N, 0, 0,
				N, N, p_grid, uplo);
		two_dim_block_cyclic_init(&dcZ, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, N, 1, 0, 0,
				N, 1, 1, 1, nodes);
		two_dim_block_cyclic_init(&dcZcpy, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, N, 1, 0, 0,
				N, 1, 1, 1, nodes);
		two_dim_block_cyclic_init(&dcdet, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, 1, 1, 0, 0,
				1, 1, 1, 1, nodes);
		two_dim_block_cyclic_init(&dcproduct, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, 1, 1, 0, 0,
				1, 1, 1, 1, nodes);

		/* Allocate memory */
		dcC.mat = calloc((size_t)dcC.super.nb_local_tiles *
				(size_t)dcC.super.bsiz, 
				(size_t)parsec_datadist_getsizeoftype(dcC.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcC, "dcC");

		dcZ.mat = calloc((size_t)dcZ.super.nb_local_tiles *
				(size_t)dcZ.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(dcZ.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcZ, "dcZ");

		dcZcpy.mat = calloc((size_t)dcZcpy.super.nb_local_tiles *
				(size_t)dcZcpy.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(dcZcpy.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcZcpy, "dcZcpy");

		dcdet.mat = calloc((size_t)dcdet.super.nb_local_tiles *
				(size_t)dcdet.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(dcdet.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcdet, "dcdet");

		dcproduct.mat = calloc((size_t)dcproduct.super.nb_local_tiles *
				(size_t)dcproduct.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(dcproduct.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcproduct, "dcproduct");

		/* Pass to data structure */
		data.descC = (void *)&dcC;
		data.descZ = (void *)&dcZ;
		data.descZcpy = (void *)&dcZcpy;
		data.descdet = (void *)&dcdet;
		data.descproduct = (void *)&dcproduct;

		data.descZmiss = (void *)&descZmiss;
		data.descC12 = (void *)&descC12;
		data.descC22 = (void *)&descC22;
		data.descZtrace = (void *)&descZtrace;
		data.descC11 = (void *)&descC11;
		data.descC21 = (void *)&descC21;
		data.descZobs = (void *)&descZobs;
		data.descmse = (void *)&descmse;
		data.descZactual = (void *)&descZactual;

		data.parsec = (void *)parsec;
		data.uplo = (int)uplo;
		data.nodes = nodes;
		data.rank = rank;
		int HNB = arguments.HNB;
		int NT = dcC.super.lmt;

		/* Defalt sed band_size_single to NT */
		if( 0 == data.band_size_single ) {
			data.band_size_single  = MY_MAX(NT, data.band_size_double);
		}

		assert( data.band_size_double > 0 );
		assert( data.band_size_single > 0 );
		assert( data.band_size_single >= data.band_size_double );

		if( 0 == rank ) print_summary(test, N, ncores, gpus, dts, lts,  data.computation, zvecs, p_grid, q_grid, data.precision);

		/* Used to check results */
		int info_potrf = 0;
		int info_trmm = 0;
		double sum = 0.0;

		/* Timer start */
		SYNC_TIME_START();
		parsec_dmatrix_generation( parsec, (parsec_tiled_matrix_dc_t *)&dcC, &data.l1, &data.l1,
				initial_theta , data.dm, data.c_fun, arguments.band_size_double );
		SYNC_TIME_PRINT(rank, ESC("Matrix_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
					p_grid, nodes/p_grid, dts, N, HNB, sync_time_elapsed));

		/* Timer start */
		SYNC_TIME_START();
		parsec_dZ_generation( parsec, (parsec_tiled_matrix_dc_t *)&dcZ, (double *)&Nrand[i*N] );
		SYNC_TIME_PRINT(rank, ESC("Z_vector_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
                                       p_grid, nodes/p_grid, dts, N, HNB, sync_time_elapsed));

#if DEBUG_INFO
		sum = parsec_dmatrix_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcC);
		if( 0 == rank ) fprintf( stderr, "\nin mian C: sum after generation %.17g\n", sum );

		sum = parsec_dZ_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcZ);
		if( 0 == rank ) fprintf( stderr, "\nin main Z: sum after generation %.17g\n", sum );
#endif

		/* Timer start */
                SYNC_TIME_START();
		info_potrf = dplasma_hsdpotrf( parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcC,
				data.lookahead, data.sendless, data.tensor_gemm, NT, NT, HNB);
                SYNC_TIME_PRINT(rank, ESC("Potrf_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d band_size_double= %d "
					"band_size_single= %d sendless= %d : %14f gflops\n",
                                        p_grid, nodes/p_grid, dts, N, HNB, NT, NT, data.sendless, 
                                        gflops=(flops/1e9)/sync_time_elapsed));

#if DEBUG_INFO
                sum = parsec_dmatrix_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcC);
                if( 0 == rank ) fprintf( stderr, "\nin mian C: sum after cholesky %.17g\n", sum );

                sum = parsec_dZ_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcZ);
                if( 0 == rank ) fprintf( stderr, "\nin main Z: sum after cholesky %.17g\n", sum );

                //fprintf(stderr, "Print Z vector before trmv\n");
                //dplasma_dprint(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcZ);
#endif

                /* Timer start */
                SYNC_TIME_START();
#if OWN_TRMV
		/* Data descritor declaration */
		two_dim_block_cyclic_t dcZ0;
		two_dim_block_cyclic_init(&dcZ0, matrix_RealDouble, matrix_Tile,
				nodes, rank, dts, dts, N, 1, 0, 0,
				N, 1, 1, 1, nodes);
		dcZ0.mat = calloc((size_t)dcZ0.super.nb_local_tiles *
				(size_t)dcZ0.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(dcZ0.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcZ0, "dcZ0");

		info_trmm = parsec_my_dtrmv( parsec, (parsec_tiled_matrix_dc_t *)&dcC,
				(parsec_tiled_matrix_dc_t *)&dcZ,
				(parsec_tiled_matrix_dc_t *)&dcZ0 );

		dplasma_dlacpy( parsec, PlasmaUpperLower,
				(parsec_tiled_matrix_dc_t *)&dcZ0,
				(parsec_tiled_matrix_dc_t *)&dcZ );

		parsec_data_free( dcZ0.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcZ0 );
#else
		info_trmm = dplasma_dtrmm( parsec, PlasmaLeft,  uplo, PlasmaNoTrans, PlasmaNonUnit, 1.0,
				(parsec_tiled_matrix_dc_t *)&dcC, (parsec_tiled_matrix_dc_t *)&dcZ);
#endif

		SYNC_TIME_PRINT(rank, ESC("dplasma_dtrmm" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
					p_grid, nodes/p_grid, dts, N, HNB, sync_time_elapsed));

#if DEBUG_INFO
                sum = parsec_dmatrix_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcC);
                if( 0 == rank ) fprintf( stderr, "\nin mian C: sum after trmm %.17g\n", sum );

                sum = parsec_dZ_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcZ);
                if( 0 == rank ) fprintf( stderr, "\nin main Z: sum after trmm %.17g\n", sum );

                //fprintf(stderr, "Print Z vector before trmv\n");
                //dplasma_dprint(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcZ);
#endif

		/* Check results */
		if( 0 == rank && info_potrf != 0 ) {
			fprintf(stderr, "-- Cholesky factorization is suspicious (info_potrf = %d) ! \n", info_potrf);
			ret |= 1;
			//exit(1);
		}

                if( 0 == rank && info_trmm != 0 ) {
                        fprintf(stderr, "-- TRMM is suspicious (info_potrf = %d) ! \n", info_trmm);
                        ret |= 1;
                }

                if( 0 == rank ) {
			fprintf(stderr, "Done Z Vector Generation Phase.\n");
			fprintf(stderr, "************************************************************\n");
		}

		if(log == 1 && test == 1)
			init_log(&data);                 

#if DEBUG_INFO
		/* Check results before optimization */
		//fprintf(stderr, "Print Z vector before optimization\n");
		//dplasma_dprint(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcZ);
#endif

		/* Timer start */
		data.iter_count = 0;
		START_TIMING(data.total_exec_time);
		nlopt_result opt_result;
		nlopt_set_max_objective(opt, MLE_alg_parsec, (void *)&data);

		/* Set opt in data */
		data.opt = &opt;

		opt_result = nlopt_optimize(opt, starting_theta, &opt_f);
		if( 0 == rank ) fprintf( stderr, "Results of nlopt_optimize : %d\n", opt_result );
		STOP_TIMING(data.total_exec_time);

		/* Free memory */
		parsec_data_free( dcC.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC );

		if(nZmiss != 0){

			//initialization
			double *Zobs;
			double *Zactual;
			double *Zmiss;

			//memory allocation
			Zobs 	= (double *) malloc(nZobs * sizeof(double));
			Zactual	= (double *) malloc(nZmiss * sizeof(double));
			Zmiss	= (double *) malloc(nZmiss * sizeof(double));
			data.mse_flag = 1;

			if( 0 == rank ) fprintf(stderr, "nZobs %d nZmiss %d\n", nZobs, nZmiss);

			/* Init */
			parsec_prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, data.mse_flag);

			double avg_pred_value=0.0;
			int pred_samples = 1;
			int j = 0;
			for (j = 0; j < pred_samples; j++)
			{
				parsec_pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);

#if 0 
				fprintf(stderr, "\n Before prediction\n");
				print_dmatrix("Zobs", nZobs, 1, Zobs, nZobs);
				print_dmatrix("Zactual", nZmiss, 1, Zactual, nZmiss);
				print_dmatrix("Zmiss", nZmiss, 1, Zmiss, nZmiss);
#endif

				/* Call predict function */
				prediction_error = parsec_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);

#if 0 
				fprintf(stderr, "\n After prediction\n");
                                print_dmatrix("Zobs", nZobs, 1, Zobs, nZobs);
                                print_dmatrix("Zactual", nZmiss, 1, Zactual, nZmiss);
                                print_dmatrix("Zmiss", nZmiss, 1, Zmiss, nZmiss);

				for (int index=0; index< nZmiss; index++)
					printf ("(%f, %f)\n ", Zactual[index], Zmiss[index]);
#endif

				if( 0 == rank ) {
					fprintf(stderr,"Prediction Error: %lf \n", prediction_error);
				}

				/* Ave predict err */
				avg_pred_value +=prediction_error;
			}

#ifndef EXAGEOSTAT_ON_SUMMIT
			write_to_thetafile("theta-pred.txt", starting_theta[0], starting_theta[1], starting_theta[2], starting_theta[3], 0, 0, 0, (avg_pred_value/=pred_samples));
#endif

			parsec_prediction_finalize(&data);

			//free memory
			free(Zobs);
			free(Zactual);
			free(Zmiss);
		}

                if(rank == 0) {
                        fprintf(stderr, "\n");
                        fprintf(stderr, "Dictionary, zvecs, i, precision, band_size_double, band_size_single, "
                                        "initial_theta0, initial_theta1, initial_theta2, initial_theta3, "
                                        //"lookahead, send_less, tensor_gemm, "
                                        //"M, N, MB, NB, HMB, HNB, nodes, P, Q, ncores, gpus, "
                                        "iteration, opt_time_iter, pred_err, num_iteration_scale_1, num_iteration_scale_2, "
					"num_iteration_scale_3, found_theta0, found_theta1, found_theta2, found_theta3\n\n");
                }

                if(rank == 0) {
                        fprintf(stderr, "OUTPUT ");
                        fprintf(stderr, "%d %d %d %d %d ", zvecs, i, data.precision, data.band_size_double, data.band_size_single); 
                        fprintf(stderr, "%0.3f %0.3f %0.3f %0.3f ", initial_theta[0], initial_theta[1], initial_theta[2], ((num_params == 4)? initial_theta[3]: 0.0)); 
                        //fprintf(stderr, "%d %d %d ", data.lookahead, data.sendless, data.tensor_gemm);
                        //fprintf(stderr, "%d %d %d %d %d %d ", N, N, dts, dts, HNB, HNB);
                        //fprintf(stderr, "%d %d %d %d %d ", nodes, p_grid, nodes/p_grid, ncores, num_gpus);
                        fprintf(stderr, "%d %f %f ", data.iter_count, data.avg_exec_time_per_iter/data.iter_count, prediction_error); 
                        fprintf(stderr, "%d %d %d ", num_iteration_scale_1, num_iteration_scale_2, num_iteration_scale_3); 
                        fprintf(stderr, "%g %g %g %g\n\n", data.variance, starting_theta[1], starting_theta[2], ((num_params == 4)? starting_theta[3]: 0.0));

			if(rank == 0) {
				fprintf(stderr, "\n");
				fprintf(stderr, "Dictionary2, precision, band_size_double, band_size_single, "
						"initial_theta0, initial_theta1, initial_theta2, found_theta0, found_theta1, found_theta2, "
						"iteration, pred_err, sd, sum, mean, time\n\n");
			}

			if(rank == 0) {
				fprintf(stderr, "OUTPUT2 ");
				fprintf(stderr, "%d %d ", data.band_size_double, data.band_size_single);
				fprintf(stderr, "%0.3f %0.3f %0.3f ", initial_theta[0], initial_theta[1], initial_theta[2]);
				fprintf(stderr, "%g %g %g ", data.variance, starting_theta[1], starting_theta[2]);
				fprintf(stderr, "%d %lf %le %le %le ", data.iter_count, prediction_error, data.trace_sd, data.trace_sum, data.trace_mean);
				fprintf(stderr, "%lf \n\n", data.avg_exec_time_per_iter/data.iter_count);
			}


		}

#ifndef EXAGEOSTAT_ON_SUMMIT
		print_result(&data, starting_theta, N, zvecs, ncores, lts, test, initial_theta, data.computation, p_grid, q_grid, data.final_loglik, prediction_error);
#endif

		if(log == 1 && test == 1)
			finalize_log(&data);

		/* Free memory */
		parsec_data_free( dcZ.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcZ );

		parsec_data_free( dcZcpy.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcZcpy );

		parsec_data_free( dcproduct.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcproduct );

		parsec_data_free( dcdet.mat );
		parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcdet );

	}

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
	gpu_temporay_buffer_fini( );
#endif
#endif

	/* Destroy nlopt */
	nlopt_destroy(opt);

	/* Clean up parsec*/
	parsec_fini( &parsec );

#ifdef PARSEC_HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}
