/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_ztrace.c
 *
 * Calculate trace of a given matrix (A).
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-09-19
 *
 **/

#include "../include/exageostatcore.h"


double core_dtrace (double *A, int m, int n,
        int m0, int n0, double* trace) {

    int i;
    double res = 0.0;
    for (i = 0; i < m; i++) 
    {
        res += A[i + i * m];
        trace[i] = A[i + i * m];
    }
    return res;
}

