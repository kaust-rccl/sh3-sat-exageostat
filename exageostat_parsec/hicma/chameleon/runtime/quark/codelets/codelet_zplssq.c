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
 * @brief Chameleon zplssq Quark codelet
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
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zplssq_quark(Quark *quark)
{
    double *SCALESUMSQ;
    double *SCLSSQ;

    quark_unpack_args_2( quark, SCALESUMSQ, SCLSSQ );

    if( SCLSSQ[0] < SCALESUMSQ[0] ) {
        SCLSSQ[1] = SCALESUMSQ[1] + (SCLSSQ[1]     * (( SCLSSQ[0] / SCALESUMSQ[0] ) * ( SCLSSQ[0] / SCALESUMSQ[0] )));
        SCLSSQ[0] = SCALESUMSQ[0];
    } else {
        SCLSSQ[1] = SCLSSQ[1]     + (SCALESUMSQ[1] * (( SCALESUMSQ[0] / SCLSSQ[0] ) * ( SCALESUMSQ[0] / SCLSSQ[0] )));
    }
}

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
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(opt->quark, CORE_zplssq_quark, (Quark_Task_Flags*)opt,
        sizeof(double)*2, RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INPUT,
        sizeof(double)*2, RTBLKADDR(SCLSSQ,     double, SCLSSQm,     SCLSSQn),     INOUT,
        0);
}


void CORE_zplssq2_quark(Quark *quark)
{
    double *RESULT;

    quark_unpack_args_1( quark, RESULT );
    RESULT[0] = RESULT[0] * sqrt( RESULT[1] );
}

void MORSE_TASK_zplssq2( const MORSE_option_t *options,
                         const MORSE_desc_t *RESULT, int RESULTm, int RESULTn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(opt->quark, CORE_zplssq2_quark, (Quark_Task_Flags*)opt,
        sizeof(double)*2, RTBLKADDR(RESULT, double, RESULTm, RESULTn), INOUT,
        0);
}
