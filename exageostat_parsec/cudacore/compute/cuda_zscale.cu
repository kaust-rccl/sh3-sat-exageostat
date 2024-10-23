/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file cudaconv.c
 *
 * Cuda datatypes conversion.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-03-17
 *
 **/

#define CHUNKSIZE 32

#include <cublas.h>
#include <stdio.h>
#include "../include/exageostatcudacore.h"


__global__ void dndscale_array_kernel(double *A, double *B, double *C, 
		int m, int n, int m0, int n0)
{

	const int tx  = threadIdx.x;
	const int ty  = threadIdx.y;
	const int idx = blockIdx.x * blockDim.x + tx;
	const int idy = blockIdx.y * blockDim.y + ty;

	if(idx>=m || idy >=n){return;}

	A[idx + idy * m] *= B[idx + idx *m] * C[idy +idy * m]; // power-exp kernel

}


__global__ void ddscale_array_kernel(double *A, 
                int m, int n, int m0, int n0)
{

        const int tx  = threadIdx.x;
        const int ty  = threadIdx.y;
        const int idx = blockIdx.x * blockDim.x + tx;
        const int idy = blockIdx.y * blockDim.y + ty;

        if(idx>=m || idy >=n){return;}

	if( idx != idy )
		A[idx + idy * m] *= A[idx + idx *m] * A[idy +idy * m]; // power-exp kernel

}


void dndscale_array( double *A, double *B, double *C, int m, int n, int m0,
		int n0,  cudaStream_t stream){

	int nBlockx= (m+CHUNKSIZE-1)/CHUNKSIZE;
	int nBlocky= (n+CHUNKSIZE-1)/CHUNKSIZE;
	dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
	dim3 dimGrid(nBlockx,nBlocky);

	dndscale_array_kernel<<<dimGrid,dimBlock,0,stream>>>(A, B, C, m, n, m0, n0);
}

void ddscale_array( double *A, int m, int n, int m0,
                int n0,  cudaStream_t stream){

        int nBlockx= (m+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (n+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
        dim3 dimGrid(nBlockx,nBlocky);

        ddscale_array_kernel<<<dimGrid,dimBlock,0,stream>>>(A, m, n, m0, n0);
}
