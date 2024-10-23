/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_ddotp.c
 *
 * Calculate dot product scalar value of Z*Z. 
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/exageostatcore.h"

/***************************************************************************//**
 *
 *  core_ddopt - Calculate the dot product of the Z vector (double precision).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[in] Z
 *           The n-by-1 vector.
 *
 * @param[out] dotproduct
 *          dot product value.
 *
 * @param[in] n
 *          The number of rows in the tile A.
 *
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
double core_ddotp(double *Z, double *dotproduct, int n)
{
        return cblas_ddot(n, Z, 1, Z, 1);
}




/***************************************************************************//**
 *
 *  core_sdopt - Calculate the dot product of the Z vector (single precision).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[in] Z
 *           The m-by-1 matrix on which to calculate the determinant.
 *
 * @param[out] dotproduct
 *          dot product value.
 * @param[in] n
 *          The number of rows in the tile A.
 *
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
float core_sdotp(float *Z, float *dotproduct, int n)
{
        return cblas_sdot(n, Z, 1, Z, 1);
}
