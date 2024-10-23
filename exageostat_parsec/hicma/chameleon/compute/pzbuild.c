/**
 *
 * @file pzbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  Parallel tile matrix generation
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = MorseUpperLower: All the matrix A
 *            = MorseUpper: Upper triangular part
 *            = MorseLower: Lower triangular part
 *
 * @param[in] A
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 */
void morse_pzbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback,
                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
  MORSE_context_t *morse;
  MORSE_option_t options;

  int m, n;
  int ldam;

  morse = morse_context_self();
  if (sequence->status != MORSE_SUCCESS)
    return;
  RUNTIME_options_init(&options, morse, sequence, request);

  for (m = 0; m < A->mt; m++) {
    ldam = BLKLDD(A, m);
    for (n = 0; n < A->nt; n++) {

      if ( ( uplo == MorseUpper && m <= n ) ||
           ( uplo == MorseLower && m >= n ) ||
           ( uplo == MorseUpperLower ) )
        MORSE_TASK_zbuild(
              &options,
              A(m, n), ldam,
              user_data, user_build_callback );
    }
  }

  RUNTIME_options_finalize( &options, morse);
}
