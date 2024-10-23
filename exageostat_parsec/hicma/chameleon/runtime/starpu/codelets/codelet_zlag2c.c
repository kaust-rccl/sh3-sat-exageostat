/**
 *
 * @file codelet_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void MORSE_TASK_zlag2c(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zlag2c;
    void (*callback)(void*) = options->profiling ? cl_zlag2c_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_W(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_W,         RTBLKADDR(B, MORSE_Complex32_t, Bm, Bn),
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlag2c",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlag2c_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex32_t *B;
    int ldb;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex32_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    CORE_zlag2c( m, n, A, lda, B, ldb);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

void MORSE_TASK_clag2z(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_clag2z;
    void (*callback)(void*) = options->profiling ? cl_clag2z_callback : NULL;

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( B, Bm, Bn ) )
    {
        starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,         RTBLKADDR(A, MORSE_Complex32_t, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,         RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
            STARPU_VALUE,    &ldb,               sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "clag2z",
#endif
            0);
    }
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_clag2z_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    MORSE_Complex32_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;

    A = (MORSE_Complex32_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    CORE_clag2z( m, n, A, lda, B, ldb);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlag2c, 1, cl_zlag2c_cpu_func)
/*
 * Codelet definition
 */
CODELETS_CPU(clag2z, 2, cl_clag2z_cpu_func)
