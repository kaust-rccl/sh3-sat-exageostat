/**
 *
 * @file codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  MORSE_TASK_zplssq returns: scl * sqrt(ssq)
 *
 * with scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( 2*i )**2 * A( 2*i+1 ) )
 *                      i
 *
 * The values of A(2*i+1) are assumed to be at least unity.
 * The values of A(2*i) are assumed to be non-negative and scl is
 *
 *    scl = max( A( 2*i ) ),
 *           i
 *
 * The routine makes only one pass through the matrix A.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of couple (scale, sumsq) in the matrix A.
 *
 *  @param[in] A
 *          The 2-by-M matrix.
 *
 *  @param[out] result
 *          On exit, result contains scl * sqrt( ssq )
 *
 */
void MORSE_TASK_zplssq( const MORSE_option_t *options,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                        const MORSE_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn )
{
    struct starpu_codelet *codelet = &cl_zplssq;
    void (*callback)(void*) = options->profiling ? cl_zplssq_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(SCALESUMSQ, SCALESUMSQm, SCALESUMSQn);
    MORSE_ACCESS_RW(SCLSSQ, SCLSSQm, SCLSSQn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_R,  RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn),
        STARPU_RW, RTBLKADDR(SCLSSQ,     double, SCLSSQm,     SCLSSQn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zplssq",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq_cpu_func(void *descr[], void *cl_arg)
{
    double *SCALESUMSQ;
    double *SCLSSQ;

    SCALESUMSQ = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCLSSQ     = (double *)STARPU_MATRIX_GET_PTR(descr[1]);

    if( SCLSSQ[0] < SCALESUMSQ[0] ) {
        SCLSSQ[1] = SCALESUMSQ[1] + (SCLSSQ[1]     * (( SCLSSQ[0] / SCALESUMSQ[0] ) * ( SCLSSQ[0] / SCALESUMSQ[0] )));
        SCLSSQ[0] = SCALESUMSQ[0];
    } else {
        SCLSSQ[1] = SCLSSQ[1]     + (SCALESUMSQ[1] * (( SCALESUMSQ[0] / SCLSSQ[0] ) * ( SCALESUMSQ[0] / SCLSSQ[0] )));
    }

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq, 2, cl_zplssq_cpu_func)

void MORSE_TASK_zplssq2( const MORSE_option_t *options,
                         const MORSE_desc_t *RESULT, int RESULTm, int RESULTn )
{
    struct starpu_codelet *codelet = &cl_zplssq2;
    void (*callback)(void*) = options->profiling ? cl_zplssq2_callback : NULL;

    if ( morse_desc_islocal( RESULT, RESULTm, RESULTn ) ) {
        starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_RW, RTBLKADDR(RESULT, double, RESULTm, RESULTn),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zplssq2",
#endif
            0);
    }
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq2_cpu_func(void *descr[], void *cl_arg)
{
    double *RESULT;

    RESULT = (double *)STARPU_MATRIX_GET_PTR(descr[0]);

    RESULT[0] = RESULT[0] * sqrt( RESULT[1] );

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq2, 1, cl_zplssq2_cpu_func)
