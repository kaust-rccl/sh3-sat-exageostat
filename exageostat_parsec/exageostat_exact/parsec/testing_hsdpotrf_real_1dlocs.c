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

        int N, lts, dts, log;
        int i = 0;
        int zvecs = 1, nZmiss = 0, test = 0, gpus = 0;
        int p_grid, q_grid, ncores;
        double  opt_f;
        arguments arguments;
        nlopt_opt opt;
        double *streamdata_u;
        double *streamdata_v;
        int ncid;
        int ncid_pred;
        double prediction_error = 0.0;
        MLE_data data;
        location *locations;
        location *missing_locations;

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
	check_args(&arguments);

#if defined(EXAGEOSTAT_USE_CUDA)
        extern char **environ;
        int num_gpus = atoi( arguments.gpus );
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

        int num_params = 0;
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

        char *dim1 = (char *) malloc(10 * sizeof(char));
        char *dim2 = (char *) malloc(10 * sizeof(char));

        dim1 = "latitude";
        dim2 = "longitude";
#if defined(EXAGEOSTAT_USE_NETCDF)
        ncid      = openFileNC(&data, data.locsFPath);
        N = countlinesNC(ncid, dim1, dim2);
        readLocsNC_1d(&data, ncid);
#endif
        locations = &data.l1;

	locations->z = NULL;
        assert( strcmp(arguments.dim, "2d") == 0 );

	int nZobs = strcmp(data.actualZFPath,"") == 0? (N-nZmiss) : N;

        /* Flops */
	PASTE_CODE_FLOPS(FLOPS_DPOTRF, ((DagDouble_t)N));

        /* Data descritor declaration */
	sym_two_dim_block_cyclic_t dcC;
        two_dim_block_cyclic_t dcZ;
        two_dim_block_cyclic_t dcZcpy;
	two_dim_block_cyclic_t dcdet;
	two_dim_block_cyclic_t dcproduct;

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

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
        /* GPU temporary buffer */
        gpu_temporay_buffer_init( dts, dts ); 
#endif
#endif

        /* stop gsl error handler */
        gsl_set_error_handler_off();

	/* Read observation file */
	streamdata_u = (double *) malloc(N * sizeof(double));
	streamdata_v = (double *) malloc(N * sizeof(double));
#if defined(EXAGEOSTAT_USE_NETCDF)
	readVarNCs(&data, ncid, "UGRD_10maboveground", streamdata_u, dim1, dim2);
	readVarNCs(&data, ncid, "VGRD_10maboveground", streamdata_v, dim1, dim2);
#endif
	/* Calculate the wind speed */
	for (i = 0; i < N; i++)
		streamdata_u[i] = sqrt(pow(streamdata_u[i], 2) + pow(streamdata_v[i], 2));
	locations_obs_zsort_inplace(N, locations, streamdata_u);

	/* Timer start */
	SYNC_TIME_START();
	parsec_dZ_generation( parsec, (parsec_tiled_matrix_dc_t *)&dcZ, streamdata_u );
	SYNC_TIME_PRINT(rank, ESC("Z_vector_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
				p_grid, nodes/p_grid, dts, N, 1, sync_time_elapsed));

	if(log == 1 && test == 1)
		init_log(&data);       

	START_TIMING(data.total_exec_time);
	nlopt_result opt_result;
	nlopt_set_max_objective(opt, MLE_alg_parsec, (void *)&data);
	opt_result = nlopt_optimize(opt, starting_theta, &opt_f);
	if( 0 == rank ) fprintf( stderr, "Results of nlopt_optimize : %d\n", opt_result );
	STOP_TIMING(data.total_exec_time);

	if (strcmp(data.actualZLocFPath,"")!= 0)
	{               
		printf( "%s ========\n", data.actualZLocFPath);
#if defined(EXAGEOSTAT_USE_NETCDF)		
		ncid_pred       = openFileNC(&data, data.actualZLocFPath);
		nZmiss          = countlinesNC(ncid, dim1, dim2);
		readLocsNC_1d( &data, ncid);
#endif
		missing_locations = &data.l1;
	} 

	/* Destroy nlopt */
	nlopt_destroy(opt);

	/* Free memory */
	parsec_data_free( dcC.mat );
	parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC );

	if(nZmiss != 0) {

		//initialization
		double *Zobs;
		double *Zactual;
		double *Zmiss;

		//memory allocation
		Zobs 	= (double *) malloc(nZobs * sizeof(double));
		Zactual	= (double *) malloc(nZmiss * sizeof(double));
		Zmiss	= (double *) malloc(nZmiss * sizeof(double));
		data.mse_flag = 1;

		if( 0 == rank ) printf("nZobs %d nZmiss %d\n", nZobs, nZmiss);

		/* Init */
		parsec_prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, data.mse_flag);

		double avg_pred_value=0.0;
		int pred_samples = 10;
		int j = 0;
		for (j = 0; j < pred_samples; j++)
		{
			if (strcmp(data.actualZLocFPath,"") == 0)
				parsec_pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);
			else
			{
				streamdata_u = (double *) malloc(nZmiss * sizeof(double));
				streamdata_v = (double *) malloc(nZmiss * sizeof(double));
#ifndef EXAGEOSTAT_ON_NETCDF
				readVarNCs(&data, ncid_pred, "UGRD_10maboveground", streamdata_u, dim1, dim2);
				readVarNCs(&data, ncid_pred, "VGRD_10maboveground", streamdata_v, dim1, dim2);
#endif
				for (int i=0; i < nZmiss; i++)
					Zactual[i] = sqrt(pow(streamdata_u[i], 2) + pow(streamdata_v[i], 2));
				parsec_get_zobs(parsec, (parsec_tiled_matrix_dc_t *)&dcZcpy, Zobs, N);
				data.lmiss = *missing_locations;
				data.lobs  = *locations;
			}
#if DEBUG_INFO
			fprintf(stderr, "\n Before prediction\n");
			print_dmatrix("Zobs", nZobs, 1, Zobs, nZobs);
			print_dmatrix("Zactual", nZmiss, 1, Zactual, nZmiss);
			print_dmatrix("Zmiss", nZmiss, 1, Zmiss, nZmiss);
#endif

			//generate_interior_points(&data, Zobs, NULL, nZmiss, nZobs, N);
			prediction_error = parsec_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);

#if DEBUG_INFO
			fprintf(stderr, "\n After prediction\n");
			print_dmatrix("Zobs", nZobs, 1, Zobs, nZobs);
			print_dmatrix("Zactual", nZmiss, 1, Zactual, nZmiss);
			print_dmatrix("Zmiss", nZmiss, 1, Zmiss, nZmiss);
#endif

			if( 0 == rank ) {
#if DEBUG_INFO
				for (int index=0; index< nZmiss; index++)
					printf ("(%f, %f)\n ", Zactual[index], Zmiss[index]);
#endif

				fprintf(stderr,"Prediction Error: %lf \n", prediction_error);
			}

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
                                "initial_theta0, initial_theta1, initial_theta2, "
                                "iteration, opt_time_iter, pred_err, found_theta0, found_theta1, found_theta2\n\n");
        }       

        if(rank == 0) { 
                fprintf(stderr, "OUTPUT ");
                fprintf(stderr, "%d %d %d %d %d ", zvecs, 0, data.precision, data.band_size_double, data.band_size_single);
                fprintf(stderr, "%0.3f %0.3f %0.3f ", initial_theta[0], initial_theta[1], initial_theta[2]);
                fprintf(stderr, "%d %f %f ", data.iter_count, data.avg_exec_time_per_iter/data.iter_count, prediction_error);
                fprintf(stderr, "%g %g %g\n\n", data.variance, starting_theta[1], starting_theta[2]);
        } 

#ifndef EXAGEOSTAT_ON_SUMMIT
	print_result(&data, starting_theta, N, zvecs, ncores, lts, test, initial_theta, data.computation, p_grid, q_grid, data.final_loglik, prediction_error);
#endif
	if(log == 1 && test == 1)
		finalize_log(&data);
#ifndef EXAGEOSTAT_ON_NETCDF
        closeFileNC(&data, ncid);
#endif
#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
        gpu_temporay_buffer_fini( );
#endif
#endif

	/* Free memory */
        parsec_data_free( dcZ.mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcZ );

        parsec_data_free( dcZcpy.mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcZcpy );

        parsec_data_free( dcproduct.mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcproduct );

        parsec_data_free( dcdet.mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcdet );

	/* Clean up parsec*/
	parsec_fini( &parsec );

#ifdef PARSEC_HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}

