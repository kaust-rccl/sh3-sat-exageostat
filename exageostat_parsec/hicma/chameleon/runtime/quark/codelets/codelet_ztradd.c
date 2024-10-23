/**
 *
 * @file codelet_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2011-11-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_ztradd_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *B;
    int LDB;

    quark_unpack_args_10(quark, uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_ztradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  MORSE_TASK_ztradd adds two trapezoidal matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          = MorseUpperLower: A and B are general matrices.
 *          = MorseUpper: op(A) and B are upper trapezoidal matrices.
 *          = MorseLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = MorseNoTrans:   op(A) = A
 *          = MorseTrans:     op(A) = A'
 *          = MorseConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = MorseNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = MorseNoTrans, and k=N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
void MORSE_TASK_ztradd(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans, int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_Complex64_t beta,  const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GEADD;
    QUARK_Insert_Task(opt->quark, CORE_ztradd_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                 &uplo,  VALUE,
        sizeof(MORSE_enum),                 &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(MORSE_Complex64_t),         &alpha, VALUE,
        sizeof(MORSE_Complex64_t)*lda*n,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(MORSE_Complex64_t),         &beta,   VALUE,
        sizeof(MORSE_Complex64_t)*ldb*n,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);

    (void)nb;
}

