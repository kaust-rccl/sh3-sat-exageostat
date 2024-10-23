/**
 *
 * @file codelet_zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zbuild_quark(Quark *quark)
{
    MORSE_Complex64_t *A;
    int lda;
    void *user_data;
    void (*user_build_callback)(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data) ;
    int row_min, row_max, col_min, col_max;

    quark_unpack_args_8( quark, row_min, row_max, col_min, col_max, A, lda, user_data, user_build_callback);

    user_build_callback(row_min, row_max, col_min, col_max, A, lda, user_data);
}

void MORSE_TASK_zbuild( const MORSE_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_BUILD;
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;

    QUARK_Insert_Task(opt->quark, CORE_zbuild_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                      &row_min,    VALUE,
                      sizeof(int),                      &row_max,    VALUE,
                      sizeof(int),                      &col_min,    VALUE,
                      sizeof(int),                      &col_max,    VALUE,
                      sizeof(MORSE_Complex64_t)*lda*A->nb, RTBLKADDR(A, MORSE_Complex64_t, Am, An),         OUTPUT,
                      sizeof(int),                      &lda,  VALUE,
                      sizeof(void*),                    &user_data,  VALUE,
                      sizeof(void*),                    &user_build_callback,   VALUE,
                      0);
}
