/**
 *
 * @file codelet_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_zgetrf(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,             &m,                        sizeof(int),
        STARPU_VALUE,             &n,                        sizeof(int),
        STARPU_RW,                     RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,           &lda,                        sizeof(int),
        STARPU_VALUE,                  &IPIV,                      sizeof(int*),
        STARPU_VALUE,    &check_info,                sizeof(MORSE_bool),
        STARPU_VALUE,         &iinfo,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int *IPIV;
    MORSE_bool check_info;
    int iinfo;
    int info = 0;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &IPIV, &check_info, &iinfo);
    CORE_zgetrf( m, n, A, lda, IPIV, &info );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf, 1, cl_zgetrf_cpu_func)
