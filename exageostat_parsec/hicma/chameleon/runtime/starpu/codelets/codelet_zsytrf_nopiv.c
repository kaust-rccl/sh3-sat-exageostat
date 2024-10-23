/**
 *
 * @file codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @date 2011-10-09
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_zsytrf_nopiv(const MORSE_option_t *options,
                             MORSE_enum uplo, int n, int nb,
                             const MORSE_desc_t *A, int Am, int An, int lda,
                             int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zsytrf_nopiv;
    void (*callback)(void*) = options->profiling ? cl_zsytrf_nopiv_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        /* STARPU_SCRATCH,   options->ws_worker, */
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zsytrf_nopiv",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zsytrf_nopiv_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int iinfo;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda, &iinfo);
    CORE_zsytf2_nopiv(uplo, n, A, lda);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zsytrf_nopiv, 1, cl_zsytrf_nopiv_cpu_func)
