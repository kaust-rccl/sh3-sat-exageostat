/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_exact.c
 *
 * ExaGeoStat exact computation main functions (i.e., generate synthetic dataset, evaluate ML function, and predication).
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-17
 *
 **/
#include "../include/MLE_exact.h"
//***************************************************************************************
void MORSE_MLE_dzvg_Tile (MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log)
	//! Generate Observations Vector (Z) for testing Maximum
	/*! Likelihood function -- MORSE-sync 
	 * Returns Z observation vector
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
	 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
	 * 	                     that is used to to generate the Covariance Matrix.
	 * @param[in] n: Problem size (number spatial locations).
	 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
	 * @param[in] log: equals one if the user needs to generate log files for his problem.
	 * */
{
	MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
	MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
	//In the case of testing mode, Z should be generated using Nrand and initial_theta
	//if (test == 1)    
	//{
	//Generate the co-variance matrix C
	VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm, data->c_fun);
	VERBOSE(" Done.\n");
                MORSE_Sequence_Wait(msequence);
//double sum=0;
    //    double *C = (double *) malloc(n * n * sizeof(double));
  //      MORSE_Tile_to_Lapack( data->descC, C, n);

//int i=0;
//for(i=0;i<n*n;i++)
//	sum+=C[i];
//printf("sum= %f\n", sum);
//exit(0);
 //       print_dmatrix("testC", 16, 16, C, 16);
//exit(0);

//double *C = (double *) malloc(n *n *sizeof(double));
//	MORSE_Tile_to_Lapack( data->descC, C, n);
//	print_dmatrix("testC", 16, 16, C, 16);
//	exit(0);

	//Copy Nrand to Z
	VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
	MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
	VERBOSE(" Done.\n");

	//Cholesky factorization for the Co-variance matrix C
	VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
	int success = MORSE_dpotrf_Tile(MorseLower, data->descC);
	//printf(" success=%d \n", success);
	//exit(0);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");
         //double *C = (double *) malloc(n * n * sizeof(double));
        //MORSE_Tile_to_Lapack( data->descC, C, n);
        //print_dmatrix("testC", 16, 16, C, 16);
        //exit(0);

	//Triangular matrix-matrix multiplication    
	VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
	MORSE_dtrmm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ);
	VERBOSE(" Done.\n");

	//if log==1 write vector to disk
	if(log==1)
	{
		double *z;
		MORSE_desc_t *MORSE_descZ = (MORSE_desc_t *)(data->descZ);
#if defined(CHAMELEON_USE_MPI)
		z = (double *) malloc(n * sizeof(double));
		MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
#else
		z = MORSE_descZ->mat;
#endif
		write_vectors(z, data, n);
#if defined(CHAMELEON_USE_MPI)
		free(z);
#endif
	}

	/*	}
		else
		{
		double * streamdata;
		streamdata=(double *) malloc(n * sizeof(double));

	//Reading Observations from disk and copy it to data->descZ
	VERBOSE("Reading Observations from disk .....");
	streamdata = readObsFile(data->obsFPath, n);        
	MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
	MORSE_Sequence_Wait(data->sequence);
	VERBOSE(" Done.\n");
	free(streamdata);
	}
	 */
	MORSE_dlaset_Tile(MorseUpperLower, 0, 0, data->descC);
	VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)\n");
	VERBOSE("************************************************************\n");
}



void MORSE_MLE_dzcpy( MLE_data *data, double *streamdata)
	//! Copy measurements vector from Lapack
	/*! format to Chameleon format.
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] streamdata: measurments vector in lapack format.
	 * */
{
	MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
	MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
	VERBOSE("Copy Z from vector to decriptor.\n");
	MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
	MORSE_Sequence_Wait(msequence);
	VERBOSE("Done Z copying step.\n");
	VERBOSE("************************************************************\n");
}

void MORSE_MLE_dzvg_Tile_Async(MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log)
	//! Generate Observations Vector (Z) for testing Maximum
	/*! Likelihood function -- MORSE-Async
	 * Returns Z observation vector
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
	 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] n: Problem size (number spatial locations).
	 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
	 * @param[in] log: equals one if the user needs to generate log files for his problem.
	 * */
{
	MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
	MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
	//In the case of testing mode, Z should be generated using Nrand and initial_theta
	//       if (test ==1)
	//      {
	//Generate the co-variance matrix C
	VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm, data->c_fun);
	VERBOSE(" Done.\n");

	//Copy Nrand to Z
	VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
	MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
	VERBOSE(" Done.\n");

	//Cholesky factorization for the Co-variance matrix C
	VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
	int success = MORSE_dpotrf_Tile_Async(MorseLower, data->descC, msequence, mrequest);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");

	//Triangular matrix-matrix multiplication
	VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
	MORSE_dtrmm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ, msequence, mrequest);
	VERBOSE(" Done.\n");

	//if log == 1 write vector to disk
	if(log == 1)
	{
		double *z;
		MORSE_desc_t *MORSE_descZ = (MORSE_desc_t *)(data->descZ);
#if defined(CHAMELEON_USE_MPI)
		z = (double *) malloc(n * sizeof(double));
		MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
#else
		z = MORSE_descZ->mat;
#endif
		write_vectors(z, data, n);

#if defined(CHAMELEON_USE_MPI)
		free(z);
#endif
	}

	/*       }
		 else
		 {
		 double * streamdata;
		 streamdata=(double *) malloc(n * sizeof(double));

	//Reading Observations from disk and copy it to data->descZ
	VERBOSE("Reading Observations from disk .....");
	streamdata = readObsFile(data->obsFPath, n);
	MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
	MORSE_Sequence_Wait(data->sequence);
	VERBOSE(" Done.\n");
	free(streamdata);
	}
	 */
	VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous)\n");
	VERBOSE("************************************************************\n");
}



double MORSE_dmle_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data) {
	//! Maximum Likelihood Evaluation (MLE)
	/*!  -- MORSE-sync
	 * Returns the loglikelihhod value for the given theta.
	 * @param[in] n: unsigned variable used by NLOPT library.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] grad: double variable used by NLOPT library. 
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * */
	//Initialization
	double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, dzcpy_time=0.0;
	int N, NRHS, success;
	double flops =0.0;	

	MLE_data* data	= ((MLE_data*)MORSE_data);
	data->det	= 0;
	data->dotp	= 0;

	MORSE_desc_t *MORSE_descC	= (MORSE_desc_t *) data->descC;
	MORSE_desc_t *MORSE_descZ	= (MORSE_desc_t *) data->descZ;
	MORSE_desc_t *MORSE_descZcpy	= (MORSE_desc_t *) data->descZcpy; 
	MORSE_desc_t *MORSE_descdet	= (MORSE_desc_t *) data->descdet;
	MORSE_desc_t *MORSE_descproduct	= (MORSE_desc_t *) data->descproduct;
	MORSE_sequence_t *msequence	= (MORSE_sequence_t *) data->sequence;
	MORSE_request_t  *mrequest	= (MORSE_request_t *) data->request;
	int num_params = 0;
	N	= MORSE_descC->m;
	NRHS	= MORSE_descZ->n;

        if(strcmp (data->c_fun, "pow-exp-nuggets") == 0)
                num_params=4;
        else
                num_params=3;


	START_TIMING(dzcpy_time);
	if(data->iter_count==0)
		//Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
		MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy);
	if(strcmp(data->recovery_file,"") != 0 && recover(data->recovery_file, data->iter_count, theta, &loglik, 3));
	else
	{
		START_TIMING(dzcpy_time);
		if(data->iter_count==0)
			//Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
			MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy); 
		else
		{	
			VERBOSE("Re-store the original Z vector...");
			MORSE_dlacpy_Tile(MorseUpperLower ,MORSE_descZcpy,MORSE_descZ);
			VERBOSE(" Done.\n");
		}
		STOP_TIMING(dzcpy_time);	

		// double *C = (double *) malloc(N * N * sizeof(double));

		//Generate new co-variance matrix C based on new theta	
		VERBOSE("Generate New Covariance Matrix...");
		START_TIMING(matrix_gen_time);	
		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm, data->c_fun);    
		MORSE_Sequence_Wait(msequence);
		STOP_TIMING(matrix_gen_time);
		VERBOSE(" Done.\n");

		//VERBOSE("Scale the Covariance Matrix...");
		//START_TIMING(matrix_gen_time);
		//MORSE_MLE_dscale_Tile_Async(MorseUpperLower, MORSE_descC, msequence, &mrequest[0]);
		//MORSE_Sequence_Wait(msequence);
		//STOP_TIMING(matrix_gen_time);
		//VERBOSE(" Done.\n");

		int i;
		//VERBOSE("Shift the Covariance Matrix...");
		//START_TIMING(matrix_gen_time);
		//MORSE_MLE_dshift_Tile_Async(MORSE_descC, msequence, &mrequest[0]);
		//MORSE_Sequence_Wait(msequence);
		//STOP_TIMING(matrix_gen_time);
		//VERBOSE(" Done.\n");
		//MORSE_MLE_dprint_Tile_Async(MORSE_descC, msequence, &mrequest[0]);
		// MORSE_Sequence_Wait(msequence);
		//exit(0);

		//double *C = (double *) malloc(N * N * sizeof(double));
		//MORSE_Tile_to_Lapack( MORSE_descC, C, N);
		//print_dmatrix("testC", 16, 16, C, 16);

		//Calculate Cholesky Factorization (C=LL-1)
		VERBOSE("Cholesky factorization of Sigma...");
		START_TIMING(time_facto);
		success = MORSE_dpotrf_Tile(MorseLower, MORSE_descC);
		STOP_TIMING(time_facto);
		SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
		flops = flops + FLOPS_DPOTRF(N);
		VERBOSE(" Done.\n");

		//*********************************************
		//you need to generate the full matrix
		/*MORSE_desc_t *MORSE_descC2       = NULL;
		  MORSE_desc_t *MORSE_descC3       = NULL;
		  MORSE_desc_t *MORSE_descC4       = NULL;
		  MORSE_Sequence_Wait(msequence);
		  EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC2,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
		  EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC4,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
		  MORSE_dlaset_Tile(MorseUpperLower, 0, 0, MORSE_descC4);
		  MORSE_dlacpy_Tile(MorseLower ,MORSE_descC,MORSE_descC4);
		  EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC3,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
		  MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC2, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm);
		  MORSE_Sequence_Wait(msequence);
		  MORSE_dgemm_Tile (MorseNoTrans, MorseTrans, 1, MORSE_descC4, MORSE_descC4, 0, MORSE_descC3);

		  double error=0;
		  double norm_c=0;
		  MORSE_dgeadd_Tile_Async( MorseTrans, 1.0, MORSE_descC2
		  ,-1.0, MORSE_descC3, msequence, &mrequest[0] );
		  MORSE_Sequence_Wait(msequence);
		  MORSE_dlange_Tile_Async( MorseFrobeniusNorm,
		  MORSE_descC3, &error, msequence, &mrequest[0]
		  );
		  MORSE_Sequence_Wait(msequence);
		  MORSE_dlange_Tile_Async( MorseFrobeniusNorm,
		  MORSE_descC2, &norm_c, msequence, &mrequest[0]
		  );


		  MORSE_Sequence_Wait(msequence);
		  printf("error: %e\n", (error/norm_c));
		  exit(0);
		 */
		//***************************************
		//MORSE_Tile_to_Lapack( MORSE_descC, C, N);
		//print_dmatrix("testC", 16, 16, C, 16);
		//exit(0);
		//MORSE_MLE_dprint_Tile_Async(MORSE_descC, msequence, &mrequest[0]);
		//MORSE_Sequence_Wait(msequence);
		//exit(0);

		//Calculate log(|C|) --> log(square(|L|))
		VERBOSE("Calculating the log determinant ...");
		START_TIMING(logdet_calculate);
		MORSE_MLE_dmdet_Tile_Async(MORSE_descC, msequence, &mrequest[0], MORSE_descdet);
		MORSE_Sequence_Wait(msequence);
		// printf("det: %f\n", data->det);
		logdet= 2*data->det;
		STOP_TIMING(logdet_calculate);
		VERBOSE(" Done.\n");

		//        printf("logdet: %f\n",logdet);        


		//Solving Linear System (L*X=Z)--->inv(L)*Z
		VERBOSE("Solving the linear system ...\n");
		START_TIMING(time_solve);
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);
		STOP_TIMING(time_solve);
		flops = flops + FLOPS_DTRSM(MorseLeft,N, NRHS);
		VERBOSE(" Done.\n");    



		//Calculate MLE likelihood
		VERBOSE("Calculating the MLE likelihood function ...");
		//dotp=0;
		//MORSE_MLE_core_ddotp_Async(MORSE_descZ,MORSE_descproduct,msequence, &mrequest[0]);
		//MORSE_Sequence_Wait(msequence);        

		MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct); 
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
		//loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
		VERBOSE(" Done.\n");
	}
	//Distribute the values in the case of MPI

#if defined(CHAMELEON_USE_MPI)
	MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast(theta, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		if(strcmp(data->checkpoint_file,"") != 0)	
			checkpointing(data->checkpoint_file, data->iter_count, theta, loglik, 3);
		//Print Iteration Summary
		//fprintf(stderr,"***************************************************\n");
		fprintf(stderr,"\n------ddotproduct: %.17g ", data->dotp);
		fprintf(stderr,"\n------logdet: %.17g ", logdet);
		//fprintf(stderr,"------det: %.*e ", det);
		//fprintf(stderr,"\n------expr2: %2.6f \n",((double) (N / 2) * log(2 * PI)));
		//fprintf(stderr," ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----LogLi: %2.6f\n", theta[0], theta[1], theta[2],loglik);
		//reformat
		if(num_params == 4)
			printf(" %3d- Model Parameters (variance, range, smoothness): (%2.8f, %2.8f, %2.8f, %2.8f) ----> LogLi: %.17g\n", data->iter_count+1,  data->variance, theta[1], theta[2], theta[3], loglik);

		else
			printf(" %3d- Model Parameters (variance, range, smoothness): (%2.8f, %2.8f, %2.8f) ----> LogLi: %.17g\n", data->iter_count+1,  data->variance, theta[1], theta[2],loglik);

		if(data->log == 1)
			fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%2.8f, %2.8, %2.8f) ----> LogLi: %.17g\n", data->iter_count+1, data->variance, theta[1], theta[2],loglik);

		printf(" ---- Facto Time: %6.2f\n", time_facto);
		printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
		printf(" ---- dtrsm Time: %6.2f\n", time_solve);
		printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
		//fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", zcpy_time);
		printf(" ---- Total Time: %6.2f\n", matrix_gen_time+ time_facto + logdet_calculate + time_solve);
		//fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );
		printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
		//fprintf(stderr," ---- Peak Performance: %6.2f Gflops/s\n",  (ncores*p_grid*q_grid*16*2.3) );
		//fprintf(stderr,"***************************************************\n");

#if defined(CHAMELEON_USE_MPI)
	}
#endif

	data->iter_count++;
	// for experiments
	data->avg_exec_time_per_iter+=/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
	data->avg_flops_per_iter+=flops / 1e9 / (time_facto +time_solve);
	data->final_loglik=loglik;

	return loglik;
}

double MORSE_dmle_Tile_Async(unsigned n, const double * theta, double * grad, void * MORSE_data) {
	//! Maximum Likelihood Evaluation (MLE)
	/*!  -- MORSE-Async
	 * Returns the loglikelihhod value for the given theta.
	 * @param[in] n: unsigned variable used by NLOPT library.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] grad: double variable used by NLOPT library.
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * */
	//Initialization
	double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, dzcpy_time=0.0, flops = 0.0;
	int N, NRHS, success;

	MLE_data* data	= ((MLE_data*)MORSE_data);
	data->det	= 0;
	data->dotp	= 0;

	MORSE_desc_t *MORSE_descC	= (MORSE_desc_t *) data->descC;
	MORSE_desc_t *MORSE_descZ	= (MORSE_desc_t *) data->descZ;
	MORSE_desc_t *MORSE_descZcpy	= (MORSE_desc_t *) data->descZcpy; 
	MORSE_desc_t *MORSE_descdet	= (MORSE_desc_t *) data->descdet;
	MORSE_desc_t *MORSE_descproduct	= (MORSE_desc_t *) data->descproduct;
	MORSE_sequence_t *msequence	= (MORSE_sequence_t *) data->sequence;
	MORSE_request_t  *mrequest	= (MORSE_request_t *) data->request;

	N	= MORSE_descC->m;
	NRHS	= MORSE_descZ->n;
	START_TIMING(dzcpy_time);
	if(data->iter_count == 0)
		//Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZ, MORSE_descZcpy, msequence, mrequest); 
	else
	{	
		VERBOSE("re-store the original Z vector...");
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZcpy, MORSE_descZ, msequence, mrequest);
		VERBOSE(" Done.\n");
	}
	STOP_TIMING(dzcpy_time);	


	//Generate new co-variance matrix C based on new theta	
	VERBOSE("Generate New Covariance Matrix...");
	START_TIMING(matrix_gen_time);	
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest, &data->l1, &data->l1,(double*) theta,  data->dm, data->c_fun);    
	MORSE_Sequence_Wait(msequence);
	STOP_TIMING(matrix_gen_time);
	VERBOSE(" Done.\n");

	//Calculate Cholesky Factorization (C=LL-1)
	VERBOSE("Cholesky factorization of Sigma...");
	START_TIMING(time_facto);
	success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest);
	STOP_TIMING(time_facto);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	flops = flops + FLOPS_DPOTRF(N);
	VERBOSE(" Done.\n");

	//Calculate log(|C|) --> log(square(|L|))
	VERBOSE("Calculating the log determinant ...");
	START_TIMING(logdet_calculate);
	MORSE_MLE_dmdet_Tile_Async(MORSE_descC, msequence, &mrequest[0],MORSE_descdet);
	MORSE_Sequence_Wait(msequence);
	logdet= 2*data->det;
	STOP_TIMING(logdet_calculate);
	VERBOSE(" Done.\n");

	//Solving Linear System (L*X=Z)--->inv(L)*Z
	VERBOSE("Solving the linear system ...\n");
	START_TIMING(time_solve);
	MORSE_dtrsm_Tile_Async(MorseLeft,MorseLower,MorseNoTrans,MorseNonUnit,1,MORSE_descC,MORSE_descZ, msequence, mrequest);
	STOP_TIMING(time_solve);
	flops = flops + FLOPS_DTRSM(MorseLeft,N, NRHS);
	VERBOSE(" Done.\n");    

	//Claculate MLE likelihood
	VERBOSE("Calculating the MLE likelihood function ...");
	MORSE_dgemm_Tile_Async (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct, msequence, mrequest); 
	loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
	VERBOSE(" Done.\n");

	//Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
	MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		//Print Iteration Summary
		//fprintf(stderr,"***************************************************\n");
		fprintf(stderr,"------ddotproduct: %2.6f ", data->dotp);
		fprintf(stderr,"------logdet: %2.6f ", logdet);
		//fprintf(stderr,"------det: %.*e ", det);
		fprintf(stderr,"------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
		//fprintf(stderr," ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----LogLi: %2.6f\n", theta[0], theta[1], theta[2],loglik);
		//reformat
		fprintf(stderr," %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

		if(data->log == 1)
			fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

		fprintf(stderr," ---- Facto Time: %6.2f\n", time_facto);
		fprintf(stderr," ---- logdet Time: %6.2f\n", logdet_calculate);
		fprintf(stderr," ---- dtrsm Time: %6.2f\n", time_solve);
		fprintf(stderr," ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
		//fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", dzcpy_time);
		fprintf(stderr," ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
		//fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );    
		//fprintf(stderr," ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
		//fprintf(stderr," ---- Peak Performance: %6.2f Gflops/s\n",  (ncores*p_grid*q_grid*16*2.3) );
		//fprintf(stderr,"***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
	}
#endif

	data->iter_count++;
	// for experiments
	data->avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
	data->avg_flops_per_iter+=flops / 1e9 / (time_facto +time_solve);
	data->final_loglik=loglik;

	return loglik;
}


void MORSE_dmle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid, int mse_flag)
	//! Allocate prediction operation descriptors.
	/*!  
	 * Returns MLE_data data with initial values and new descriptors locations.
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] dts: tile size (MB).
	 * @param[in] p_grid: p_grid in the case of distributed system.
	 * @param[in] q_grid: q_grid in the case of distributed system.
	 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
	 * */
{

	MORSE_desc_t *MORSE_descZmiss   = NULL;
    MORSE_desc_t *MORSE_descZtrace  = NULL;
    MORSE_desc_t *MORSE_descC11     = NULL;
    MORSE_desc_t *MORSE_descC21     = NULL;
	MORSE_desc_t *MORSE_descC12     = NULL;
	MORSE_desc_t *MORSE_descC22     = NULL;
	MORSE_desc_t *MORSE_descmse     = NULL;
	MORSE_desc_t *MORSE_descZactual = NULL;
	MORSE_desc_t *MORSE_descZobs    = NULL;
	MLE_data     *data              = (MLE_data*) MORSE_data;

	if(nZmiss <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return;
	}
	//Descriptors Creation
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
	if( mse_flag == 1)
	{
		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);

		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	}
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZtrace, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC11, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZmiss, 0, 0, nZmiss, nZmiss, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC21, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZmiss, 0, 0, nZobs, nZmiss, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

	//Initiate data descriptors
	data->descZmiss         = MORSE_descZmiss;
	data->descC12           = MORSE_descC12;
	data->descC22           = MORSE_descC22;
	data->descmse           = MORSE_descmse;
	data->descZactual       = MORSE_descZactual;
	data->descZobs          = MORSE_descZobs;
    data->descZtrace        = MORSE_descZtrace;
    data->descC11           = MORSE_descC11;
    data->descC21           = MORSE_descC21;
}


double MORSE_dmle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-sync
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] Zobs: observed values vector (known observations).
	 * @param[in] Zmiss missing values vector (unknown observations).
	 * @param[in] Zactual: actual missing values vector (in the case of testing MSE).
	 * @param[in] n: number of spatial locations.
	 * */
{

	//initialization	
	//double *z = NULL, *streamdata = NULL;
	double time_solve = 0.0;
	double mat_gen_time = 0.0;
	double time_gemm = 0.0;
	double time_mse = 0.0;
	double flops = 0.0;

	MORSE_desc_t *MORSE_descZmiss   = NULL;
    MORSE_desc_t *MORSE_descZtrace  = NULL;
    MORSE_desc_t *MORSE_descC11     = NULL;
    MORSE_desc_t *MORSE_descC21     = NULL;
	MORSE_desc_t *MORSE_descC12     = NULL;
	MORSE_desc_t *MORSE_descC22     = NULL;
	MORSE_desc_t *MORSE_descmse     = NULL;
	MORSE_desc_t *MORSE_descZactual = NULL;
	MORSE_desc_t *MORSE_descZobs    = NULL;
    MORSE_desc_t *MORSE_descdet     = NULL;
	MLE_data     *data              = (MLE_data*) MORSE_data;
	MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
	MORSE_request_t *mrequest       = (MORSE_request_t  *) data->request;
	data->mserror                   = 0;

	if(nZmiss <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return -1;
	}

	//Descriptors Creation
	//EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
	//if( Zactual != NULL)
	//{
	//        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
	//        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	//}
	//EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
	//EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
	//EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

	//Initiate data descriptors
	MORSE_descZmiss		= data->descZmiss;
	MORSE_descC12		= data->descC12;
	MORSE_descC22		= data->descC22;
	MORSE_descmse		= data->descmse;
	MORSE_descZactual	= data->descZactual;
	MORSE_descZobs		= data->descZobs;
    MORSE_descZtrace    = data->descZtrace;
    MORSE_descC11       = data->descC11;
    MORSE_descC21       = data->descC21;
	MORSE_descdet       = data->descdet;

	//Copy data to vectors 
	VERBOSE("Copy measurments vector to descZobs descriptor...");
	//MORSE_MLE_dzcpy_Tile_Async(MORSE_descZobs, Zobs, msequence, mrequest);
	MORSE_Lapack_to_Tile( Zobs, nZobs, MORSE_descZobs);
	VERBOSE(" Done.\n");

	if( Zactual != NULL)	
	{
		VERBOSE("Copy actual measurments vector to descZactual descriptor...");
		//MORSE_MLE_dzcpy_Tile_Async(MORSE_descZactual, Zactual, msequence, mrequest);
		MORSE_Lapack_to_Tile( Zactual, nZmiss, MORSE_descZactual);
		VERBOSE(" Done.\n");
	}

	MORSE_Sequence_Wait(msequence);


	//        int i=0;
	//      for (i=0;i<nZmiss;i++)
	//    printf("%f, %f, %f\n", data->lmiss.x[i], data->lmiss.y[i], Zactual[i]);

	//printf("\n\n");

	//	for (i=0;i<100;i++)
	//	printf("%f, %f, %f\n", data->lobs.x[i], data->lobs.y[i], Zobs[i]);

#if defined(CHAMELEON_USE_MPI)
	MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
	if(strcmp(data->c_fun, "matern") == 0)
		theta[0] = data->variance;
	printf("estimated parameters: %f - %f - %f (%s)\n", theta[0], theta[1], theta[2], data->c_fun);

    //Generate C22 covariance matrix
    VERBOSE("Generate C11 Covariance Matrix... (Prediction Stage)");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC11, msequence, mrequest, &data->lmiss, &data->lmiss, theta, data->dm, data->c_fun);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");


    //Generate C22 covariance matrix
    VERBOSE("Generate C21 Covariance Matrix... (Prediction Stage)");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC21, msequence, mrequest, &data->lobs, &data->lmiss, theta, data->dm, data->c_fun);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");

	START_TIMING(mat_gen_time);
	//Generate C22 covariance matrix
	VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm, data->c_fun);
	MORSE_Sequence_Wait(msequence);
	//flops = flops + FLOPS_DPOTRF(nZobs);
	VERBOSE(" Done.\n");


	//Generate C12 covariance matrix, TODO: no need for this, we can use C21 transpose
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm, data->c_fun);
	MORSE_Sequence_Wait(msequence);
	//flops = flops + FLOPS_DPOTRF(nZmiss);
	VERBOSE(" Done.\n");
	STOP_TIMING(mat_gen_time);



	START_TIMING(time_solve);
	//Start prediction
	VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
	MORSE_dposv_Tile(MorseLower, MORSE_descC22, MORSE_descZobs);
	flops = flops + FLOPS_DPOTRF(nZobs);
	flops = flops + FLOPS_DTRSM(MorseLeft, nZobs, nZobs);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_solve);


	START_TIMING(time_gemm);
	VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
	MORSE_dgemm_Tile (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss);
	flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_gemm);

    //************************************************************************
    //Schur_complement
    VERBOSE("Calculate dposv C22 Covariance Matrix and C21... (Prediction Stage - Schur_complement)");
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC22, MORSE_descC21);
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descC22, MORSE_descC21);
    VERBOSE(" Done.\n");

    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage - Schur_complement)");
    MORSE_dgemm_Tile (MorseNoTrans, MorseNoTrans, -1, MORSE_descC12, MORSE_descC21, 1, MORSE_descC11);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");

	//Use data->det to store the trace summation
    data->det=0;
    VERBOSE("Calculate trace estimation... (Prediction Stage - Schur_complement)");
    MORSE_MLE_dtrace_Tile_Async(MORSE_descC11, msequence, mrequest, MORSE_descdet, MORSE_descZtrace);
    MORSE_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    double mean = (data->det)/nZmiss;
    double sd = 0;
    double* Ztrace = (double *) malloc(nZmiss  * sizeof(double));
    MORSE_Tile_to_Lapack( MORSE_descZtrace, Ztrace, nZmiss);
    int i=0;
    for(i=0;i<nZmiss;i++)
        sd +=pow((Ztrace[i] - mean), 2);
    sd = sqrt(sd/nZmiss);
    ////**********************************************************************

	//return back descZmiss to zmiss vector
	MORSE_Tile_to_Lapack( MORSE_descZmiss, Zmiss, nZmiss);

	//Estimate Mean Square Error
	if( Zactual != NULL)
	{
		START_TIMING(time_mse);
		VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
		MORSE_MLE_dmse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
		MORSE_Sequence_Wait(msequence);
		VERBOSE(" Done.\n");	
		STOP_TIMING(time_mse);
		data->mserror /= nZmiss;
	}
	else
		data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		if(data->log == 1)
			fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

		write_prediction_result("predict_result.dat", n, data->hicma_acc, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));

        printf("Trace estimation (trace standard deviation ): %6.4e\n", sd);
        printf("Trace estimation (trace mean ): %6.4e\n", mean);
        printf("Trace estimation (trace sum ): %6.4e\n", data->det);
#if defined(CHAMELEON_USE_MPI)
	}
#endif

	return data->mserror;

}



double MORSE_dmle_Predict_Tile_Async(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-Async
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] n: number of spatial locations.
	 * */
{

	//initialization	
	location *l1 = NULL, *l2 = NULL;
	location temp_loc;
	double mat_gen_time = 0.0;
	double time_solve = 0.0;
	double time_mse = 0.0;
	double flops = 0.0; 
	MLE_data * data			=  (MLE_data*)MORSE_data;
	MORSE_desc_t *MORSE_descZ    	= (MORSE_desc_t *)(data->descZcpy);
	MORSE_desc_t *MORSE_descZobs	= (MORSE_desc_t *)(data->descZobs);	
	MORSE_desc_t *MORSE_descZactual	= (MORSE_desc_t *)(data->descZactual);
	MORSE_desc_t *MORSE_descZmiss 	= (MORSE_desc_t *)(data->descZmiss);	
	MORSE_desc_t *MORSE_descC12 	= (MORSE_desc_t *)(data->descC12);
	MORSE_desc_t *MORSE_descC22 	= (MORSE_desc_t *)(data->descC22);
	MORSE_desc_t *MORSE_descmse 	= (MORSE_desc_t *)(data->descmse);
	MORSE_sequence_t *msequence 	= (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest 	= (MORSE_request_t *)data->request;

	if(strcmp(data->actualZFPath,"")==0)
	{
		double *z = NULL;
#if defined(CHAMELEON_USE_MPI)
		z = (double *) malloc(n * sizeof(double));
		MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
#else
		z = MORSE_descZ->mat;
#endif

		//random  shuffle
		//	shuffle(z, &data->l1, n);

#if defined(CHAMELEON_USE_MPI)
		MORSE_Lapack_to_Tile( z, n, MORSE_descZ);
#endif

		l1 = &data->l1;
		temp_loc.x=&l1->x[nZmiss];
		temp_loc.y=&l1->y[nZmiss];
		l2 = &temp_loc;
	}
	else
	{
		double *streamdata = NULL;
		l1 = &data->l1;
		temp_loc.x=&l1->x[nZmiss];
		temp_loc.y=&l1->y[nZmiss];
		l2 = &temp_loc;

		//l1 = (location *) malloc(sizeof(location));

		//l1->x=(double *) malloc(nZmiss * sizeof(double));
		//l1->y=(double *) malloc(nZmiss * sizeof(double));

		VERBOSE("Reading ActualZ locations for prediction from disk .....");
		l1 = readLocsFile(data->actualZLocFPath, nZmiss);	
		VERBOSE(" Done.\n");

		//streamdata=(double *) malloc(nZmiss * sizeof(double));
		VERBOSE("Reading ActualZ for prediction from disk .....");
		streamdata = readObsFile(data->actualZFPath, nZmiss);
		MORSE_MLE_dzcpy_Tile_Async(MORSE_descZactual, streamdata, msequence, mrequest);
		MORSE_Sequence_Wait(data->sequence);
		VERBOSE(" Done.\n");
	}


	//MORSE_dposv_Tile_Async(MorseLower, MORSE_descC22, MORSE_descZobs, data->sequence, &data->request[0]);
	START_TIMING(mat_gen_time);

	//Generate C22 covariance matrix
	VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  l2, l2, theta, data->dm, data->c_fun);
	//flops = flops + FLOPS_DPOTRF(nZobs);
	VERBOSE(" Done.\n");

	//Generate C12 covariance matrix
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  l1, l2, theta, data->dm, data->c_fun);
	//flops = flops + FLOPS_DPOTRF(nZmiss);
	VERBOSE(" Done.\n");
	STOP_TIMING(mat_gen_time);

	START_TIMING(time_solve);
	//Start prediction
	VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
	MORSE_dposv_Tile_Async(MorseLower, MORSE_descC22, MORSE_descZobs, msequence, mrequest);
	flops = flops + FLOPS_DPOTRF(nZobs);
	flops = flops + FLOPS_DTRSM(MorseLeft, nZobs, nZobs);
	VERBOSE(" Done.\n");

	VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
	MORSE_dgemm_Tile_Async (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss, msequence, mrequest);
	flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_solve);


	//Estimate Mean Square Error
	START_TIMING(time_mse);
	VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
	MORSE_MLE_dmse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_mse);


	//if you do not have actual value to compare with
	if(data->descZactual==NULL)
		return -1;

	data->mserror /= nZmiss;

#if defined(CHAMELEON_USE_MPI)
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		if(data->log == 1)
			fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

		// write_prediction_result("predict_result.dat", n, nZmiss, data->mserror, (mat_gen_time+time_solve+ time_mse), (flops / 1e9 / (time_solve )));

#if defined(CHAMELEON_USE_MPI)
	}
#endif

	return data->mserror;
}



void MORSE_dmle_mloe_mmom_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid)


	//! Allocate prediction operation descriptors.
	/*!
	 * Returns MLE_data data with initial values and new descriptors locations.
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] dts: tile size (MB).
	 * @param[in] p_grid: p_grid in the case of distributed system.
	 * @param[in] q_grid: q_grid in the case of distributed system.
	 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
	 * */
{

	MORSE_desc_t *MORSE_desck_t         = NULL;
	MORSE_desc_t *MORSE_desck_a         = NULL;
	MORSE_desc_t *MORSE_desck_atmp      = NULL;
	MORSE_desc_t *MORSE_desck_ttmp      = NULL;
	//MORSE_desc_t *MORSE_desck_atmp3     = NULL;

	MORSE_desc_t *MORSE_descK_t         = NULL;
	MORSE_desc_t *MORSE_descK_ttmp      = NULL;
	MORSE_desc_t *MORSE_descK_a         = NULL;

	MORSE_desc_t *MORSE_descexpr1       = NULL;
	MORSE_desc_t *MORSE_descexpr2       = NULL;
	MORSE_desc_t *MORSE_descexpr3       = NULL;
	MORSE_desc_t *MORSE_descexpr4       = NULL;
	MORSE_desc_t *MORSE_desc_mloe_mmom       = NULL;
	MLE_data     *data              = (MLE_data*) MORSE_data;

	if(nZmiss <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return;
	}

	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desc_mloe_mmom, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 3,  0, 0, nZmiss, 3, p_grid, q_grid);

	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_t,    NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1,  0, 0, nZobs, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_a,    NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1,  0, 0, nZobs, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_atmp, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1,  0, 0, nZobs, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_ttmp, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1,  0, 0, nZobs, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr1, &data->expr1, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr2, &data->expr2, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr3, &data->expr3, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr4, &data->expr4, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_t, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_a, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_ttmp, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

	//Initiae data descriptors
	data->desck_t		= MORSE_desck_t;
	data->desck_a		= MORSE_desck_a;
	data->descK_ttmp	= MORSE_descK_ttmp;
	data->desck_atmp        = MORSE_desck_atmp;
	data->desck_ttmp        = MORSE_desck_ttmp;
	data->descK_t		= MORSE_descK_t;
	data->descK_a		= MORSE_descK_a;
	data->descexpr1		= MORSE_descexpr1;
	data->descexpr2		= MORSE_descexpr2;
	data->descexpr3		= MORSE_descexpr3;
	data->descexpr4		= MORSE_descexpr4;
	data->desc_mloe_mmom    = MORSE_desc_mloe_mmom;
}


void MORSE_dmle_mloe_mmom_Tile_Async(MLE_data *MORSE_data, double * truth_theta, double* estimated_theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-sync
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] Zobs: observed values vector (known observations).
	 * @param[in] Zmiss missing values vector (unknown observations).
	 * @param[in] Zactual: actual missing values vector (in the case of testing MSE).
	 * @param[in] n: number of spatial locations.
	 * */
{
	// truth_theta[0]=1; truth_theta[1]=0.1; truth_theta[2]=0.5;
	// estimated_theta[0]=1.01; estimated_theta[1]=0.09; estimated_theta[2]=0.49;

	printf("%f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2]);
	printf("%f, %f, %f\n", estimated_theta[0], estimated_theta[1], estimated_theta[2]);
	double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
	double loe_sum  = 0.0;
	double mom_sum  = 0.0;
	int i           = 0;
	int p           = 0;
	double all_time = 0.0;
	//************************************************************************
	double *loe=(double *) malloc(nZmiss * sizeof(double));
	double *mom=(double *) malloc(nZmiss * sizeof(double));

	MLE_data     *data                  = (MLE_data*) MORSE_data;
	MORSE_desc_t *MORSE_desck_t         = data->desck_t;
	MORSE_desc_t *MORSE_desck_a         = data->desck_a;
	MORSE_desc_t *MORSE_descK_t         = data->descK_t;
	MORSE_desc_t *MORSE_descK_a         = data->descK_a;
	MORSE_desc_t *MORSE_descK_ttmp      = data->descK_ttmp;
	MORSE_desc_t *MORSE_desck_atmp      = data->desck_atmp;
	MORSE_desc_t *MORSE_desck_ttmp      = data->desck_ttmp;
	MORSE_desc_t *MORSE_descexpr1       = data->descexpr1;
	MORSE_desc_t *MORSE_descexpr2       = data->descexpr2;
	MORSE_desc_t *MORSE_descexpr3       = data->descexpr3;
	MORSE_desc_t *MORSE_descexpr4       = data->descexpr4;

	MORSE_desc_t *MORSE_desc_mloe_mmom   = data->desc_mloe_mmom;
	MORSE_sequence_t *msequence         = (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest           = (MORSE_request_t *)data->request;
	location lmiss;
	lmiss.x            = (double *) malloc(sizeof(double));
	lmiss.y            = (double *) malloc(sizeof(double));


	START_TIMING(all_time);
	//(1)Generate the co-variance matrix descK_a
	VERBOSE("Create MORSE_descK_a Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descK_a, msequence, mrequest, &data->lobs, &data->lobs, estimated_theta, data->dm, data->c_fun);
	VERBOSE(" Done.\n");

	//(2)Generate the co-variance matrix descK_t
	VERBOSE("Create MORSE_descK_t Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descK_t, msequence, mrequest, &data->lobs, &data->lobs, truth_theta, data->dm, data->c_fun);
	VERBOSE(" Done.\n");

	//(3)Cholesky factorization for the Co-variance matrix MORSE_descK_a
	VERBOSE("Cholesky factorization of MORSE_descK_a (MLOE-MMOM) .....");
	int success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descK_a, msequence, mrequest);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");

	//(4)Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM)
	VERBOSE("Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM).....");
	MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descK_t, MORSE_descK_ttmp, msequence, mrequest);
	VERBOSE(" Done.\n");
	//(5)Cholesky factorization for the Co-variance matrix MORSE_descK_t
	VERBOSE("Cholesky factorization of MORSE_descK_t (MLOE-MMOM) .....");
	success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descK_t, msequence, mrequest);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");


	for(p=0; p<nZmiss; p++)
	{
		//      for (i=0; i<nZobs; i++)
		//      {
		//              k_t[i]= core_matern_vector (data->lobs.x[i],  data->lobs.y[i],
		//                              data->lmiss.x[p], data->lmiss.y[p], truth_theta, data->dm);

		//              k_a[i]= core_matern_vector (data->lobs.x[i],  data->lobs.y[i],
		//                              data->lmiss.x[p], data->lmiss.y[p], estimated_theta, data->dm);

		//      }
		//      printf("%f- %f\n", data->lmiss.x[p], data->lmiss.y[p]);

		//      MORSE_Lapack_to_Tile( k_t, nZobs, MORSE_desck_t);
		//      MORSE_Lapack_to_Tile( k_a, nZobs, MORSE_desck_a);


		lmiss.x[0]            = data->lmiss.x[p];
		lmiss.y[0]            = data->lmiss.y[p];

		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_desck_t, msequence, mrequest, &data->lobs, &lmiss, truth_theta, data->dm, data->c_fun);
		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_desck_a, msequence, mrequest, &data->lobs, &lmiss, estimated_theta, data->dm, data->c_fun);


		//(5)Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_desck_t, MORSE_desck_ttmp, msequence, mrequest );
		VERBOSE(" Done.\n");

		//(6)Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_desck_a, MORSE_desck_atmp, msequence, mrequest);
		VERBOSE(" Done.\n");

		//(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-1, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a, msequence, mrequest);
		VERBOSE(" Done.\n");

		//(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
		VERBOSE("Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
		VERBOSE("Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(12) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("Calculate dgemm MORSE_descexpr4 = MORSE_desck_a^T * MORSE_desck_t... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_t, 0, MORSE_descexpr3, msequence, mrequest);
		VERBOSE(" Done.\n");

		//(11) Calculate dgemm value= MORSE_desck_t^T * MORSE_desck_a
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_t^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_a, 0, MORSE_descexpr1, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(8) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_atmp
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_a, 0, MORSE_descexpr4, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(14) Calculate dgemm MORSE_desck_a= MORSE_descK_t * MORSE_desck_a (use k_t as k_a)
		VERBOSE("Calculate dgemm MORSE_desck_a = MORSE_descK_ttmp * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseNoTrans, MorseNoTrans, 1, MORSE_descK_ttmp, MORSE_desck_a, 0, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");

		//(15) Triangular Solve (TRSM) k_atmp = TRSM(K_a^-1, k_atmp)
		VERBOSE("Solving the linear system k_atmp = TRSM(K_t^-1, k_atmp) ...\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");

		//(16) Triangular Solve (TRSM) k_a = TRSM(K_a^-T, k_a)
		VERBOSE("Solving the linear system k_atmp = TRSM(K_t^-T, k_atmp) ...\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");


		//(13) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_t, 0, MORSE_descexpr2, msequence, mrequest);
		VERBOSE(" Done.\n");


		//MORSE_dstemp_Tile_Async(ORSE_descexpr1, MORSE_descexpr2, MORSE_descexpr2, MORSE_descexpr2, MORSE_desc_mloe_mmom, p);
		MORSE_Sequence_Wait(data->sequence);
		temp1 = truth_theta[0]-  2* data->expr1 + data->expr2;
		temp2 = truth_theta[0]-  data->expr3;
		temp3 = estimated_theta[0]- data->expr4;

		//printf ("%f, %f, %f, %f, %f, %f, %f\n", data->expr1, data->expr2, data->expr3, data->expr4, temp1, temp2, temp3);
		//*******************************
		loe[p]  = temp1/temp2-1.0;
		mom[p]  = temp3/temp1-1.0;
		loe_sum += loe[p];
		mom_sum += mom[p];
	}

	printf("\nMLOE = %f", loe_sum/nZmiss);
	printf("\nMMOM = %f\n\n\n", mom_sum/nZmiss);
	free(loe);
	free(mom);
	free(lmiss.x);
	free(lmiss.y);
	STOP_TIMING(all_time);

	fprintf(stderr," ---- mloe_mmom Time: %6.2f\n\n", all_time);
}
void MORSE_dmle_mloe_mmom_Tile(MLE_data *MORSE_data, double * truth_theta, double* estimated_theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-sync
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] Zobs: observed values vector (known observations).
	 * @param[in] Zmiss missing values vector (unknown observations).
	 * @param[in] Zactual: actual missing values vector (in the case of testing MSE).
	 * @param[in] n: number of spatial locations.
	 * */
{
	truth_theta[0]=1; truth_theta[1]=0.1; truth_theta[2]=0.5;
	estimated_theta[0]=1.01; estimated_theta[1]=0.09; estimated_theta[2]=0.49;
	printf("%f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2]);
	printf("%f, %f, %f\n", estimated_theta[0], estimated_theta[1], estimated_theta[2]);
	//exit(0);     
	double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
	double loe_sum	= 0.0;
	double mom_sum	= 0.0;
	int i		= 0;
	int p		= 0;
	double all_time = 0.0;
	//************************************************************************
	double *loe=(double *) malloc(nZmiss * sizeof(double));
	double *mom=(double *) malloc(nZmiss * sizeof(double));

	MLE_data     *data                  = (MLE_data*) MORSE_data;
	MORSE_desc_t *MORSE_desck_t         = data->desck_t;
	MORSE_desc_t *MORSE_desck_a         = data->desck_a;
	MORSE_desc_t *MORSE_descK_t         = data->descK_t;
	MORSE_desc_t *MORSE_descK_a         = data->descK_a;
	MORSE_desc_t *MORSE_descK_ttmp      = data->descK_ttmp;
	MORSE_desc_t *MORSE_desck_atmp      = data->desck_atmp;
	MORSE_desc_t *MORSE_desck_ttmp      = data->desck_ttmp;
	MORSE_desc_t *MORSE_descexpr1       = data->descexpr1;
	MORSE_desc_t *MORSE_descexpr2       = data->descexpr2;
	MORSE_desc_t *MORSE_descexpr3       = data->descexpr3;
	MORSE_desc_t *MORSE_descexpr4       = data->descexpr4;
	MORSE_sequence_t *msequence  	    = (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest           = (MORSE_request_t *)data->request;
	location lmiss; 
	lmiss.x            = (double *) malloc(sizeof(double));
	lmiss.y            = (double *) malloc(sizeof(double));

	START_TIMING(all_time);
	//(1)Generate the co-variance matrix descK_a
	VERBOSE("Create MORSE_descK_a Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descK_a, msequence, mrequest, &data->lobs, &data->lobs, estimated_theta, data->dm, data->c_fun);
	VERBOSE(" Done.\n");

	//(2)Generate the co-variance matrix descK_t
	VERBOSE("Create MORSE_descK_t Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descK_t, msequence, mrequest, &data->lobs, &data->lobs, truth_theta, data->dm, data->c_fun);
	MORSE_Sequence_Wait(msequence);
	VERBOSE(" Done.\n");

	//(3)Cholesky factorization for the Co-variance matrix MORSE_descK_a
	VERBOSE("Cholesky factorization of MORSE_descK_a (MLOE-MMOM) .....");
	int success = MORSE_dpotrf_Tile(MorseLower, MORSE_descK_a);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");

	//(4)Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM)
	VERBOSE("Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM).....");
	MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descK_t, MORSE_descK_ttmp);
	VERBOSE(" Done.\n");

	//(5)Cholesky factorization for the Co-variance matrix MORSE_descK_t
	VERBOSE("Cholesky factorization of MORSE_descK_t (MLOE-MMOM) .....");
	success = MORSE_dpotrf_Tile(MorseLower, MORSE_descK_t);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");


	for(p=0; p<nZmiss; p++)
	{
		//	for (i=0; i<nZobs; i++)
		//	{
		//		k_t[i]= core_matern_vector (data->lobs.x[i],  data->lobs.y[i],
		//				data->lmiss.x[p], data->lmiss.y[p], truth_theta, data->dm);

		//		k_a[i]= core_matern_vector (data->lobs.x[i],  data->lobs.y[i],
		//				data->lmiss.x[p], data->lmiss.y[p], estimated_theta, data->dm);

		//	}
		//      printf("%f- %f\n", data->lmiss.x[p], data->lmiss.y[p]);

		//	MORSE_Lapack_to_Tile( k_t, nZobs, MORSE_desck_t);
		//	MORSE_Lapack_to_Tile( k_a, nZobs, MORSE_desck_a);


		lmiss.x[0]            = data->lmiss.x[p];
		lmiss.y[0]            = data->lmiss.y[p];

		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_desck_t, msequence, mrequest, &data->lobs, &lmiss, truth_theta, data->dm, data->c_fun);
		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_desck_a, msequence, mrequest, &data->lobs, &lmiss, estimated_theta, data->dm, data->c_fun);

		MORSE_Sequence_Wait(msequence);

		//(5)Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile(MorseUpperLower, MORSE_desck_t, MORSE_desck_ttmp);
		VERBOSE(" Done.\n");

		//(6)Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile(MorseUpperLower, MORSE_desck_a, MORSE_desck_atmp);
		VERBOSE(" Done.\n");

		//(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-1, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a);
		VERBOSE(" Done.\n");


		//(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a);
		VERBOSE(" Done.\n");

		//(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
		VERBOSE("Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t);
		VERBOSE(" Done.\n");


		//(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
		VERBOSE("Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t);
		VERBOSE(" Done.\n");


		//(12) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("Calculate dgemm MORSE_descexpr4 = MORSE_desck_a^T * MORSE_desck_t... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_t, 0, MORSE_descexpr3);
		VERBOSE(" Done.\n");

		//(11) Calculate dgemm value= MORSE_desck_t^T * MORSE_desck_a
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_t^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_a, 0, MORSE_descexpr1);
		VERBOSE(" Done.\n");


		//(8) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_atmp
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_a, 0, MORSE_descexpr4);
		VERBOSE(" Done.\n");


		//(14) Calculate dgemm MORSE_desck_a= MORSE_descK_t * MORSE_desck_a (use k_t as k_a)
		VERBOSE("Calculate dgemm MORSE_desck_a = MORSE_descK_ttmp * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseNoTrans, MorseNoTrans, 1, MORSE_descK_ttmp, MORSE_desck_a, 0, MORSE_desck_t);
		VERBOSE(" Done.\n");

		//(15) Triangular Solve (TRSM) k_atmp = TRSM(K_a^-1, k_atmp)
		VERBOSE("Solving the linear system k_atmp = TRSM(K_t^-1, k_atmp) ...\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t);
		VERBOSE(" Done.\n");

		//(16) Triangular Solve (TRSM) k_a = TRSM(K_a^-T, k_a)
		VERBOSE("Solving the linear system k_atmp = TRSM(K_t^-T, k_atmp) ...\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t);
		VERBOSE(" Done.\n");


		//(13) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_t, 0, MORSE_descexpr2);
		VERBOSE(" Done.\n");


		temp1 = truth_theta[0]-  2* data->expr1 + data->expr2;
		temp2 = truth_theta[0]-  data->expr3;
		temp3 = estimated_theta[0]- data->expr4;


		//		printf ("%f, %f, %f, %f, %f, %f, %f\n", data->expr1, data->expr2, data->expr3, data->expr4, temp1, temp2, temp3);

		//*******************************
		loe[p]  = temp1/temp2-1.0;
		mom[p]  = temp3/temp1-1.0;
		loe_sum += loe[p];
		mom_sum += mom[p];
	}

	printf("\nMLOE = %f", loe_sum/nZmiss);
	printf("\nMMOM = %f\n\n\n", mom_sum/nZmiss);
	free(loe);
	free(mom);
	free(lmiss.x);
	free(lmiss.y);
	STOP_TIMING(all_time);

	fprintf(stderr," ---- mloe_mmom Time: %6.2f\n\n", all_time);
}




//init Chameleon descriptors
void MORSE_dmle_Call(MLE_data  *data, int ncores,int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
	//! //Initiate MORSE and allocate different descriptors for
	/*!  CHAMELEON
	 * Returns MLE_data data with initial values and new descriptors locations.
	 * @param[in] data: MLE_data struct with different MLE inputs.
	 * @param[in] ncores: number of CPU workers.
	 * @param[in] gpus: number of GPU workers.
	 * @param[in] dts: tile size (MB).
	 * @param[in] p_grid: p_grid in the case of distributed system.
	 * @param[in] q_grid: q_grid in the case of distributed system.
	 * @param[in] N: number of spatial locations.
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * */
{

	MORSE_sequence_t *msequence;
	MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };
	MORSE_desc_t *MORSE_descC	= NULL;
	MORSE_desc_t *MORSE_descZ	= NULL;
	MORSE_desc_t *MORSE_descZcpy	= NULL;
	MORSE_desc_t *MORSE_descproduct	= NULL;
	MORSE_desc_t *MORSE_descdet	= NULL;
	//MORSE_desc_t *MORSE_descZmiss	= NULL;
	//MORSE_desc_t *MORSE_descC12	= NULL;
	//MORSE_desc_t *MORSE_descC22	= NULL;
	//MORSE_desc_t *MORSE_descmse	= NULL;
	//MORSE_desc_t *MORSE_descZactual	= NULL;
	//MORSE_desc_t *MORSE_descZobs	= NULL;



	// For ditributed system and should be removed
	double *Zcpy = (double *) malloc(N * sizeof(double));

	//Identifies a set of routines sharing common exception handling.
	MORSE_Sequence_Create(&msequence);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC,NULL , MorseRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ, NULL, MorseRealDouble, dts, dts, dts * dts, N, 1,  0, 0, N, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZcpy, Zcpy, MorseRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct, &data->dotp, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descdet, &data->det, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);


	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	//	}
	//Fill data struct
	data->descC		= MORSE_descC;
	data->descZ		= MORSE_descZ;
	data->descZcpy		= MORSE_descZcpy;
	data->descdet		= MORSE_descdet;
	data->descproduct	= MORSE_descproduct;
	//data->descZmiss		= MORSE_descZmiss;
	//data->descC12		= MORSE_descC12;
	//data->descC22		= MORSE_descC22;
	//data->descmse		= MORSE_descmse;
	//data->descZactual	= MORSE_descZactual;
	//data->descZobs		= MORSE_descZobs;
	data->sequence		= msequence;
	data->request		= mrequest;
	//stop gsl error handler
	gsl_set_error_handler_off () ;

	}


