/**
 *
 * @file codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge StarPU codelet
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
void MORSE_TASK_zhe2ge(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)mb;
    struct starpu_codelet *codelet = &cl_zhe2ge;
    void (*callback)(void*) = options->profiling ? cl_zhe2ge_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_W(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,  &uplo,                sizeof(MORSE_enum),
        STARPU_VALUE,     &m,                        sizeof(int),
        STARPU_VALUE,     &n,                        sizeof(int),
        STARPU_R,             RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,   &lda,                        sizeof(int),
        STARPU_W,             RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
        STARPU_VALUE,   &ldb,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zhe2ge",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhe2ge_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    int M;
    int N;
    const MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    A = (const MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &LDA, &LDB);
    CORE_zhe2ge(uplo, M, N, A, LDA, B, LDB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhe2ge, 2, cl_zhe2ge_cpu_func)
