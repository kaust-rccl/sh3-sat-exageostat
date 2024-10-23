/**
 *
 * @file codelet_zpamm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpamm PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  ZPAMM  performs one of the matrix-matrix operations
 *
 *                    LEFT                      RIGHT
 *     OP MorseW  :  W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *     OP MorseA2 :  A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 *  where  op( V ) is one of
 *
 *     op( V ) = V   or   op( V ) = V**T   or   op( V ) = V**H,
 *
 *  A1, A2 and W are general matrices, and V is:
 *
 *        l = k: rectangle + triangle
 *        l < k: rectangle + trapezoid
 *        l = 0: rectangle
 *
 *  Size of V, both rowwise and columnwise, is:
 *
 *         ----------------------
 *          side   trans    size
 *         ----------------------
 *          left     N     M x K
 *                   T     K x M
 *          right    N     K x N
 *                   T     N x K
 *         ----------------------
 *
 *  LEFT (columnwise and rowwise):
 *
 *              |    K    |                 |         M         |
 *           _  __________   _              _______________        _
 *              |    |    |                 |             | \
 *     V:       |    |    |            V':  |_____________|___\    K
 *              |    |    | M-L             |                  |
 *           M  |    |    |                 |__________________|   _
 *              |____|    |  _
 *              \    |    |                 |    M - L    | L  |
 *                \  |    |  L
 *           _      \|____|  _
 *
 *  RIGHT (columnwise and rowwise):
 *
 *          |         K         |                   |    N    |
 *          _______________        _             _  __________   _
 *          |             | \                       |    |    |
 *     V':  |_____________|___\    N        V:      |    |    |
 *          |                  |                    |    |    | K-L
 *          |__________________|   _             K  |    |    |
 *                                                  |____|    |  _
 *          |    K - L    | L  |                    \    |    |
 *                                                    \  |    |  L
 *                                               _      \|____|  _
 *
 *  Arguments
 *  ==========
 *
 * @param[in] op
 *
 *         OP specifies which operation to perform:
 *
 *         @arg MorseW  : W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *         @arg MorseA2 : A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 * @param[in] side
 *
 *         SIDE specifies whether  op( V ) multiplies A2
 *         or W from the left or right as follows:
 *
 *         @arg MorseLeft  : multiply op( V ) from the left
 *                            OP MorseW  :  W  = A1 + op(V) * A2
 *                            OP MorseA2 :  A2 = A2 - op(V) * W
 *
 *         @arg MorseRight : multiply op( V ) from the right
 *                            OP MorseW  :  W  = A1 + A2 * op(V)
 *                            OP MorseA2 :  A2 = A2 - W * op(V)
 *
 * @param[in] storev
 *
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored in V:
 *
 *         @arg MorseColumnwise
 *         @arg MorseRowwise
 *
 * @param[in] M
 *         The number of rows of the A1, A2 and W
 *         If SIDE is MorseLeft, the number of rows of op( V )
 *
 * @param[in] N
 *         The number of columns of the A1, A2 and W
 *         If SIDE is MorseRight, the number of columns of op( V )
 *
 * @param[in] K
 *         If SIDE is MorseLeft, the number of columns of op( V )
 *         If SIDE is MorseRight, the number of rows of op( V )
 *
 * @param[in] L
 *         The size of the triangular part of V
 *
 * @param[in] A1
 *         On entry, the M-by-N tile A1.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, if OP is MorseA2 A2 is overwritten
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
 *
 * @param[in] V
 *         The matrix V as described above.
 *         If SIDE is MorseLeft : op( V ) is M-by-K
 *         If SIDE is MorseRight: op( V ) is K-by-N
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *
 * @param[in,out] W
 *         On entry, the M-by-N matrix W.
 *         On exit, W is overwritten either if OP is MorseA2 or MorseW.
 *         If OP is MorseA2, W is an input and is used as a workspace.
 *
 * @param[in] LDW
 *         The leading dimension of array WORK.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */


/**/

static inline int
CORE_zpamm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int op;
    MORSE_enum side;
    int storev;
    int M;
    int N;
    int K;
    int L;
    MORSE_Complex64_t *A1;
    int LDA1;
    MORSE_Complex64_t *A2;
    int LDA2;
    MORSE_Complex64_t *V;
    int LDV;
    MORSE_Complex64_t *W;
    int LDW;

    parsec_dtd_unpack_args(
        this_task, &op, &side, &storev, &M, &N, &K, &L, &A1, &LDA1, &A2, &LDA2, &V, &LDV, &W, &LDW );

    CORE_zpamm( op, side, storev, M, N, K, L, A1, LDA1, A2, LDA2, V, LDV, W, LDW );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void
MORSE_TASK_zpamm(const MORSE_option_t *options,
                 int op, MORSE_enum side, int storev,
                 int m, int n, int k, int l,
                 const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                 const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *W, int Wm, int Wn, int ldw)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zpamm_parsec, options->priority, "pamm",
        sizeof(int),                        &op,                VALUE,
        sizeof(MORSE_enum),                 &side,              VALUE,
        sizeof(MORSE_enum),                 &storev,            VALUE,
        sizeof(int),                        &m,                 VALUE,
        sizeof(int),                        &n,                 VALUE,
        sizeof(int),                        &k,                 VALUE,
        sizeof(int),                        &l,                 VALUE,
        PASSED_BY_REF,         RTBLKADDR( A1, MORSE_Complex64_t, A1m, A1n ), morse_parsec_get_arena_index( A1 ) | INPUT,
        sizeof(int),                        &lda1,              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A2, MORSE_Complex64_t, A2m, A2n ), morse_parsec_get_arena_index( A2 ) | INOUT | AFFINITY,
        sizeof(int),                        &lda2,              VALUE,
        PASSED_BY_REF,         RTBLKADDR( V, MORSE_Complex64_t, Vm, Vn ), morse_parsec_get_arena_index( V ) | INPUT,
        sizeof(int),                        &ldv,               VALUE,
        PASSED_BY_REF,         RTBLKADDR( W, MORSE_Complex64_t, Wm, Wn ), morse_parsec_get_arena_index( W ) | INOUT,
        sizeof(int),                        &ldw,               VALUE,
        PARSEC_DTD_ARG_END );
}
