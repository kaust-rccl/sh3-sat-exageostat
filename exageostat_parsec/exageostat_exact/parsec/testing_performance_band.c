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

	double *starting_theta;
	double *target_theta;
	double *initial_theta;  //for testing case
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
	double* lb = (double *) malloc(3 * sizeof(double));
	double* up = (double *) malloc(3 * sizeof(double));
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

	//Memory allocation
	starting_theta	= (double *) malloc(3 * sizeof(double));
	initial_theta	= (double *) malloc(3 * sizeof(double));
	target_theta	= (double *) malloc(3 * sizeof(double));

	//MLE_data data initialization
	init(&test, &N,  &ncores,
			&gpus, &p_grid, &q_grid,
			&zvecs, &dts, &lts,
			&nZmiss, &log, initial_theta, 
			starting_theta, target_theta, lb,
			up, &data, &arguments);

#if MATRIX_GENERATION
	if( 0 == rank ) printf("dts: %d, lts: %d\n", dts, lts);

	// Optimizater initialization
	if( 0 == rank ) printf(" %d - %d \n", data.opt_tol, data.opt_max_iters);	

	//data.precision = 0;
	//nomral random generation of e -- ei~N(0, 1) to generate Z
	double *Nrand	= (double *) malloc (N * zvecs * sizeof(double));
	if( 0 == rank ) printf("Precision: %d, async %d\n", data.precision, data.async);
	LAPACKE_dlarnv(3, iseed, N*zvecs, Nrand);

	//uniform random generation for locations / read locations from disk
	if(strcmp(data.locsFPath, "") == 0)
		locations = GenerateXYLoc(N, seed);
	//	else
	//		locations = readLocsFile(data.locsFPath, N);
	data.l1   = *locations;	
#endif

        data.parsec = (void *)parsec;
        data.uplo = (int)uplo;
        data.nodes = nodes;
        data.rank = rank;
        int band_size_double = arguments.band_size_double;
        int band_size_single = arguments.band_size_single;
        int lookahead = arguments.lookahead;
        int tensor_gemm = arguments.tensor_gemm;
        int kernel_type = arguments.kernel_type;
        int sendless = arguments.sendless;
        int HNB = arguments.HNB;
        int nb_runs = arguments.nb_runs;
        int random_seed = 3872;
        int NT = (N-1)/dts + 1; 

	if( band_size_double != 1 ) {
		if( 0 == rank ) fprintf(stderr, "No support in testing_performance_band\n");
		exit(1);
	}

        /* Flops */
	PASTE_CODE_FLOPS(FLOPS_DPOTRF, ((DagDouble_t)N));

        /* Data descritor declaration, within band_size_double */
        two_dim_block_cyclic_t dcB;

        /* Initialize data structure */
        two_dim_block_cyclic_init(&dcB, matrix_RealDouble, matrix_Tile,
                        nodes, rank, dts, dts, dts*band_size_double, N, 0, 0,
                        dts*band_size_double, N, 1, 1, 1);

        /* Allocate memory */
        dcB.mat = calloc((size_t)dcB.super.nb_local_tiles *
                        (size_t)dcB.super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(dcB.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcB, "dcB");


        /* Data descritor declaration */
        sym_two_dim_block_cyclic_t dcC;

        /* Initialize data structure */
        sym_two_dim_block_cyclic_init(&dcC, matrix_RealFloat,
                        nodes, rank, dts, dts, N, N, 0, 0,
                        N, N, p_grid, uplo);

        /* Allocate memory */
        dcC.mat = calloc((size_t)dcC.super.nb_local_tiles *
                        (size_t)dcC.super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(dcC.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcC, "dcC");


	/* C0 for warmup */
        /* Data descritor declaration */
        sym_two_dim_block_cyclic_t dcC0;

        /* Initialize data structure */
	/* P*10 X Q*10 tiles */
	int NT_C0 = (p_grid >= nodes/p_grid) ? p_grid * 10 : nodes/p_grid * 10; 
        sym_two_dim_block_cyclic_init(&dcC0, matrix_RealDouble,
                        nodes, rank, 10, 10, 10*NT_C0, 10*NT_C0, 0, 0,
                        10*NT_C0, 10*NT_C0, p_grid, uplo);

        /* Allocate memory */
        dcC0.mat = calloc((size_t)dcC0.super.nb_local_tiles *
                        (size_t)dcC0.super.bsiz,
                        (size_t)parsec_datadist_getsizeoftype(dcC0.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcC0, "dcC0");

	if( 0 == band_size_single ) {
		band_size_single = NT;
	}

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
        /* GPU temporary buffer */
        gpu_temporay_buffer_init( dts, dts ); 
#endif
#endif

	/* stop gsl error handler */
	gsl_set_error_handler_off();

	/* Timer start */
	SYNC_TIME_START();

	/* Init C */
	matrix_init_args_t *op_args = (matrix_init_args_t *)malloc( sizeof(matrix_init_args_t) );
	op_args->bump = (double)N;
	op_args->seed = random_seed;
	op_args->band_size_double = 0;

	parsec_apply( parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcC,
			(tiled_matrix_unary_op_t)parsec_sdmatrix_init_operator, op_args );

        /* Init B */   
        dplasma_dplgsy( parsec, (double)(dcB.super.ln*2), matrix_UpperLower,
                        (parsec_tiled_matrix_dc_t *)&dcB, random_seed);

	/* simple performance test, divisible */
	assert( N % dts == 0 );

	/* Init C0 */
	dplasma_dplgsy( parsec, (double)(dcC0.super.ln), uplo,
			(parsec_tiled_matrix_dc_t *)&dcC0, random_seed);

	SYNC_TIME_PRINT(rank, ESC("Matrix_generation_MIX" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d "
				"band_size_double %d band_size_single %d lookahead %d ngpus %d: %lf\n",
				p_grid, nodes/p_grid, dts, N, HNB, band_size_double,
				band_size_single, lookahead, num_gpus, sync_time_elapsed));

	/* Timer start */
	SYNC_TIME_START();
	info = dplasma_hsdpotrf( parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcC0,
			lookahead, sendless, tensor_gemm, NT_C0, NT_C0, HNB);
	SYNC_TIME_PRINT(rank, ESC("Warmup_C0_small" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d "
				"band_size_double %d band_size_single %d lookahead %d ngpus %d: %lf\n",
				p_grid, nodes/p_grid, dcC0.super.nb, dcC0.super.ln, HNB, NT_C0,
				NT_C0, lookahead, num_gpus, sync_time_elapsed));

	for(i = 0; i < nb_runs; i++) {
		info = 0;
		/* Call Kernel */
                int dplasma_hsdpotrf_band( parsec_context_t *, PLASMA_enum , parsec_tiled_matrix_dc_t *, parsec_tiled_matrix_dc_t *,
                int , int , int , int , int , int );

		info = dplasma_hsdpotrf_band( parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcC,
				(parsec_tiled_matrix_dc_t *)&dcB,
				lookahead, sendless, tensor_gemm,
				band_size_double, band_size_single, HNB);
		SYNC_TIME_PRINT(rank, ESC("Potrf_Mix" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d band_size_double %d "
					"band_size_single %d lookahead %d tensor_gemm %d ngpus %d: %lf\n",
					p_grid, nodes/p_grid, dts, N, HNB, band_size_double,
					band_size_single, lookahead, tensor_gemm, num_gpus, (flops/1e9)/sync_time_elapsed));
		/* Record time */
		time_mix_precision = sync_time_elapsed;

		if(rank == 0) {
			fprintf(stderr, "OUTPUT ");
			fprintf(stderr, "%lf %lf %lf ", (flops/1e9)/time_mix_precision,  (flops/1e9)/time_dpotrf,  (flops/1e9)/time_spotrf);
			fprintf(stderr, "%d %d %d %d ", N, N, dts, dts);
			fprintf(stderr, "%d %d ", HNB, HNB);
			fprintf(stderr, "%d %d %d %d %d ", nodes, p_grid, nodes/p_grid, ncores, num_gpus);
			fprintf(stderr, "%d %d %d ", lookahead, band_size_double, band_size_single);
			fprintf(stderr, "%d %d %d ", sendless, tensor_gemm, kernel_type);
			fprintf(stderr, "%lf %lf %lf\n\n", time_mix_precision, time_dpotrf, time_spotrf);
		}

	}

	/* Timer start */
	SYNC_TIME_START();
		//parsec_band_convert_s2d(parsec, (parsec_tiled_matrix_dc_t *)&dcC, band_size_double);
		SYNC_TIME_PRINT(rank, ESC("Convert_s2d" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d "
					"band_size_double %d band_size_single %d lookahead %d ngpus %d: %lf\n",
                                        p_grid, nodes/p_grid, dts, N, HNB, band_size_double,
                                        band_size_single, lookahead, num_gpus, sync_time_elapsed));

#if DEBUG_INFO
                double sum = parsec_dmatrix_sum(parsec, (parsec_tiled_matrix_dc_t *)&dcC);
                if( 0 == rank ) fprintf( stderr, "\nin performance C: sum after Cholesky %.17g\n", sum );
		//dplasma_dprint(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcC);
#endif

		if( 0 == rank && info != 0 ) {
			fprintf(stderr, "-- Factorization is suspicious (info = %d) ! \n", info);
		}

	if(rank == 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Dictionary: FLOPS_mix_precision, FLOPS_double, FLOPS_single, M, N, MB, NB, HMB, HNB, "
				"nodes, P, Q, ncores, gpus, lookahead, band_size_double, band_size_single, "
				"send_less, tensor_gemm, kernel_type, "
				"time_mix_precision, time_double, time_single\n");
	}

#if defined(EXAGEOSTAT_USE_CUDA)
#if GPU_BUFFER_ONCE
        gpu_temporay_buffer_fini( );
#endif
#endif

        parsec_data_free( dcC0.mat );
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC0 );

	parsec_data_free( dcC.mat );
	parsec_data_free( dcB.mat );
	parsec_tiled_matrix_dc_destroy( &dcC.super );
	parsec_tiled_matrix_dc_destroy( &dcB.super );

	/* Clean up parsec*/
	parsec_fini( &parsec );

#ifdef PARSEC_HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}
