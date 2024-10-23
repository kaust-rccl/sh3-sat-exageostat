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

#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20
#define NBELEM   1
static unsigned long long int
Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, i++) {
    if (n & 1)
      ran = a_k * ran + c_k;
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

__global__ void dshift_array_kernel(double *A, int n, int n0, unsigned long long int ran, float val)
{

	const int tx  = threadIdx.x;
	const int ty  = threadIdx.y;
	const int idx = blockIdx.x * blockDim.x + tx;
	const int idy = blockIdx.y * blockDim.y + ty;

	if(idy >=n){return;}

	double tmp;

	tmp = (ran * RndF_Mul)/1000000000;
	tmp=tmp+(val*idy/100);
	A[idy + idy * n] = A[idy + idy * n] +tmp; /*+ 1e-4*/;

}


void dshift_array( double *A, int n,
		int n0, cudaStream_t stream){

	int nBlockx= (n+CHUNKSIZE-1)/CHUNKSIZE;
	int nBlocky= (n+CHUNKSIZE-1)/CHUNKSIZE;
	dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
	dim3 dimGrid(nBlockx,nBlocky);

	double acc=1e-7;
	unsigned long long int ran, jump;
	unsigned long long int seed=5956+n0;
	float val=n0*acc;

	ran = Rnd64_jump( NBELEM , seed );

	dshift_array_kernel<<<dimGrid,dimBlock,0,stream>>>(A, n, n0,  ran, val);



}
