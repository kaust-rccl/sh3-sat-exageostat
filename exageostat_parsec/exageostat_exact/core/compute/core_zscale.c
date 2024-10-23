/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dscale.c
 *
 * Scaling to avoid NPD matrix problem..
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-04-12
 *
 **/
#include "../include/exageostatcore.h"
/***************************************************************************//**
 *
 *  core_dscale - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.	
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to compute the covariance matrix.
 *
 * @param[in] m
 *          The number of rows in the tile A. 
 *
 * @param[in] n
 *          The number of cols in the tile A. 
 *
 * @param[in] m0
 *          Global row index of the tile A.
 *
 * @param[in] n0
 *          Global col index of the tile A.
 *
 * @param[in] l1
 *          Location struct of the first input.
 *
 * @param[in] l2
 *          Location struct of the second input.
 *
 * @param[in] localtheta
 *          Parameter vector that is used to generate the output covariance matrix.
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dndscale (double *A, double *B, double *C,
		int m, int n,
		int m0, int n0)
{

	//TODO export accuracy

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0, tmp;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
				A[i + j * m] *= B[i + i * m] * C[j + j * n]; /*+ 1e-4*/;
		
	}
}

}

/******************************************************************************/
void core_ddscale (double *A,
                int m, int n,
                int m0, int n0)
{

	//TODO export accuracy

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0, tmp;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			if(i != j) // avoid changing the diagonal.
				A[i + j * m] *= A[i + i * m] * A[j + j * m]; /*+ 1e-4*/;
		}

	}
}
