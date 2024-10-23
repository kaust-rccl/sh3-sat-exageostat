/**
 *
 * @file codelet_zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zplghe_quark(Quark *quark)
{
    double bump;
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_9( quark, bump, m, n, A, lda, bigM, m0, n0, seed );
    CORE_zplghe( bump, m, n, A, lda, bigM, m0, n0, seed );
}

void MORSE_TASK_zplghe( const MORSE_option_t *options,
                        double bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLGHE;
    QUARK_Insert_Task(opt->quark, CORE_zplghe_quark, (Quark_Task_Flags*)opt,
        sizeof(double),                   &bump, VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(MORSE_Complex64_t)*lda*n, RTBLKADDR(A, MORSE_Complex64_t, Am, An),         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}
