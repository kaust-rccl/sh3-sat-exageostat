/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dcmg.c
 *
 * Generate covariance matrix of a set of locations in 2D using Matern kernel or power exponential kernel.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-03-29
 *
 **/
#include "../include/exageostatcore.h"
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

void core_dshift (double *A, int n, int n0)
{
	//TODO export accuracy
	double acc=1e-7;
	unsigned long long int ran, jump;
	unsigned long long int seed=5956+n0; 
	float val=n0*acc;

	double tmp;
	int i, j;

	ran = Rnd64_jump( NBELEM , seed );
	for (i = 0; i < n; i++) {
			tmp = (ran * RndF_Mul)/1000000000;
			tmp=tmp+(val*i/100);
			A[i + i * n] = A[i + i * n] +tmp; /*+ 1e-4*/;
	}
}


