/**
 *
 * @file pzgebrd_ge2gb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgebrd_ge2gb parallel algorithm
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

void morse_pzgebrd_ge2gb(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    int k;
    int tempkm, tempkn;
    MORSE_desc_t *A1, *A2, *T1, *D1 = NULL;

    if (A->m >= A->n){
        for (k = 0; k < A->nt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

            A1 = morse_desc_submatrix(A, k*A->mb,     k*A->nb, A->m-k*A->mb, tempkn);
            A2 = morse_desc_submatrix(A, k*A->mb, (k+1)*A->nb, A->m-k*A->mb, A->n-(k+1)*A->nb);
            T1 = morse_desc_submatrix(T, k*T->mb,     k*T->nb, T->m-k*T->mb, T->nb );
            if ( D != NULL ) {
                D1 = morse_desc_submatrix(D, k*D->mb, k*D->nb, D->m-k*D->mb, tempkn);
            }

            morse_pzgeqrf( A1, T1, D1,
                           sequence, request);

            morse_pzunmqr( MorseLeft, MorseConjTrans,
                           A1, A2, T1, D1,
                           sequence, request);

            if (k+1 < A->nt){
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                A1 = morse_desc_submatrix(A,     k*A->mb, (k+1)*A->nb, tempkm,           A->n-(k+1)*A->nb);
                A2 = morse_desc_submatrix(A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb);
                T1 = morse_desc_submatrix(T,     k*T->mb, (k+1)*T->nb, T->mb,            T->n-(k+1)*T->nb);
                if ( D != NULL ) {
                    D1 = morse_desc_submatrix(D, k*D->mb, (k+1)*D->nb, tempkm,           D->n-(k+1)*D->nb);
                }

                morse_pzgelqf( A1, T1, D1,
                               sequence, request);

                morse_pzunmlq( MorseRight, MorseConjTrans,
                               A1, A2, T1, D1,
                               sequence, request);
            }
        }
    }
    else{
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

            A1 = morse_desc_submatrix(A,     k*A->mb, k*A->nb, tempkm,           A->n-k*A->nb);
            A2 = morse_desc_submatrix(A, (k+1)*A->mb, k*A->nb, A->m-(k+1)*A->mb, A->n-k*A->nb);
            T1 = morse_desc_submatrix(T,     k*T->mb, k*T->nb, T->mb,            T->n-k*T->nb);
            if ( D != NULL ) {
                D1 = morse_desc_submatrix(D, k*D->mb, k*D->nb, tempkm,           D->n-k*D->nb);
            }
            morse_pzgelqf( A1, T1, D1,
                           sequence, request);

            morse_pzunmlq( MorseRight, MorseConjTrans,
                           A1, A2, T1, D1,
                           sequence, request);

            if (k+1 < A->mt){
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                A1 = morse_desc_submatrix(A, (k+1)*A->mb,     k*A->nb, A->m-(k+1)*A->mb, tempkn);
                A2 = morse_desc_submatrix(A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb);
                T1 = morse_desc_submatrix(T, (k+1)*T->mb,     k*T->nb, T->m-(k+1)*T->mb, T->nb );
                if ( D != NULL ) {
                    D1 = morse_desc_submatrix(D, (k+1)*D->mb, k*D->nb, D->m-(k+1)*D->mb, tempkn);
                }

                morse_pzgeqrf( A1, T1, D1,
                               sequence, request);

                morse_pzunmqr( MorseLeft, MorseConjTrans,
                               A1, A2, T1, D1,
                               sequence, request);
            }
        }
    }
}
