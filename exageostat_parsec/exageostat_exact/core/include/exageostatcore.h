/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file exageostatcore.h
 *
 * Core functions header file.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/
#ifndef _EXAGEOSTATCORE_H_
#define _EXAGEOSTATCORE_H_
#include "../../../misc/include/MLE_misc.h"
#include "../../../misc/include/flat_file.h"


void core_dndscale (double *A, double *B, double *c,
                int m, int n,
                int m0, int n0);

void core_ddscale (double *A,
                int m, int n,
                int m0, int n0);

//Generate the covariance matrix.
void   core_scmg(float *A, int m, int n,
		int m0, int n0,
		location *l1, location *l2,
                double *localtheta, int distance_metric);

void   core_dcmg(double *A, int m, int n,
                int m0, int n0,
		location *l1, location *l2,
		double *localtheta, int distance_metric);


void   core_sdcmg(double *A, int m, int n,
		int m0, int n0,
		location  *l1, location *l2,
		double *localtheta, int distance_metric);



void   core_scmg_pow_exp(float *A, int m, int n,
		int m0, int n0,
		location *l1, location *l2,
		double *localtheta, int distance_metric);

void   core_dcmg_pow_exp(double *A, int m, int n,
		int m0, int n0,
		location *l1, location *l2,
		double *localtheta, int distance_metric);

void   core_dcmg_pow_exp_nuggets(double *A, int m, int n,
                int m0, int n0,
                location *l1, location *l2,
                double *localtheta, int distance_metric);


void   core_sdcmg_pow_exp(double *A, int m, int n,
		int m0, int n0,
		location  *l1, location *l2,
		double *localtheta, int distance_metric);


float  core_smdet(float * A, int m, int n,
		int m0, int n0);

double core_dmdet(double *A, int m, int n,
		int m0, int n0);

void   core_szcpy(float  *Z,  int m,
		int m0, float  *r);

void   core_dzcpy(double *Z,  int m,
		int m0, double *r);

float  core_sdotp(float  *Z, float *dotproduct,
		int n);

double core_ddotp(double *Z, double *dotproduct,
		int n);

/*void core_sdconv(float *A, double *B,
  int m, int n);

  void core_dsconv(double *A, float *B,
  int m, int n);
 */

void core_dlag2s(int m, int n,
		const double *A, int lda,
		float *B, int ldb);

void core_slag2d(int m, int n,
		const float *A, int lda,
		double *B, int ldb);

void core_sprint(float *A,
		int m, int n,
		int m0, int n0);

void core_dprint(double *A,
		int m, int n,
		int m0, int n0);

double core_dtrace (double *A, int m, int n,
        int m0, int n0, double* trace);
#endif
