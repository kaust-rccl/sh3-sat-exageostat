/**
 *
 * @file codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"


/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zlaset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = MorseUpper: Upper part of A is set;
 *          = MorseLower: Lower part of A is set;
 *          = MorseUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
void MORSE_TASK_zlaset(const MORSE_option_t *options,
                       MORSE_enum uplo, int M, int N,
                       MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                       const MORSE_desc_t *A, int Am, int An, int LDA)
{

    struct starpu_codelet *codelet = &cl_zlaset;
    void (*callback)(void*) = options->profiling ? cl_zlaset_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,  &uplo,                sizeof(MORSE_enum),
        STARPU_VALUE,     &M,                        sizeof(int),
        STARPU_VALUE,     &N,                        sizeof(int),
        STARPU_VALUE, &alpha,         sizeof(MORSE_Complex64_t),
        STARPU_VALUE,  &beta,         sizeof(MORSE_Complex64_t),
        STARPU_W,      RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,   &LDA,                        sizeof(int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlaset",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaset_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *A;
    int LDA;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &alpha, &beta, &LDA);
    CORE_zlaset(uplo, M, N, alpha, beta, A, LDA);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlaset, 1, cl_zlaset_cpu_func)
