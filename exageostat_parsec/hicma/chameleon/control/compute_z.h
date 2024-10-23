/**
 *
 * @file compute_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
/**
 *  LAPACK/Tile Descriptor accesses
 */
#define MorseDescInput  1
#define MorseDescOutput 2
#define MorseDescInout  (MorseDescInput | MorseDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
#define morse_zdesc_alloc_diag( descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );                                   \
    RUNTIME_desc_create( &(descA) );

#define morse_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    RUNTIME_desc_create( &(descA) );

/**
 *  Declarations of internal sequential functions
 */
int morse_zshift(MORSE_context_t *morse, int m, int n, MORSE_Complex64_t *A,
                  int nprob, int me, int ne, int L,
                  MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
void morse_pzbarrier_pnl2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_row2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_tl2pnl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_tl2row(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgebrd_gb2bd(MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgebrd_ge2gb(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgelqf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgemm(MORSE_enum transA, MORSE_enum transB, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetmi2(MORSE_enum idep, MORSE_enum odep, MORSE_enum storev, int m, int n, int mb, int nb, MORSE_Complex64_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, MORSE_desc_t *D, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_nopiv(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_reclap(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_rectil(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzhegst(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzhemm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzherk(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzher2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzhetrd_he2hb(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *E, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlag2c(MORSE_desc_t *A, MORSE_desc_t *SB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlange(MORSE_enum norm, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlanhe(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlascal(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaset( MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_Complex64_t beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaset2(MORSE_enum uplo, MORSE_Complex64_t alpha,                          MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaswp(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaswpc(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlauum(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzplghe(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pzplgsy(MORSE_Complex64_t bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pzplrnt(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pzpotrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzshift(int, int, int, MORSE_Complex64_t *, int *, int, int, int, MORSE_sequence_t*, MORSE_request_t*);
void morse_pzsymm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsyrk(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta,  MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsyr2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsytrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztile2band(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *descAB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztpgqrt( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pztpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pztradd(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsmpl(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungbr(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungbrrh(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqrrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglqrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungtr(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqr(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqrrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlq(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlqrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t *request );

void morse_pzgelqf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqr_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);


/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
morse_zlap2tile( MORSE_context_t *morse,
                 MORSE_desc_t *descAl, MORSE_desc_t *descAt,
                 MORSE_enum mode, MORSE_enum uplo,
                 MORSE_Complex64_t *A, int mb, int nb, int lm, int ln, int m, int n,
                 MORSE_sequence_t *seq, MORSE_request_t *req )
{
    /* Initialize the Lapack descriptor */
    *descAl = morse_desc_init_user( MorseComplexDouble, mb, nb, (mb)*(nb),
                                    lm, ln, 0, 0, m, n, 1, 1,
                                    morse_getaddr_cm, morse_getblkldd_cm, NULL  );
    descAl->mat = A;
    descAl->styp = MorseCM;

    /* Initialize the tile descriptor */
    *descAt = morse_desc_init( MorseComplexDouble, mb, nb, (mb)*(nb),
                               lm, ln, 0, 0, m, n, 1, 1 );

    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        if ( morse_desc_mat_alloc( descAt ) ) {
            morse_error( "morse_zlap2tile", "morse_desc_mat_alloc() failed");
            return MORSE_ERR_OUT_OF_RESOURCES;
        }

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & MorseDescInput ) {
            morse_pzlacpy( uplo, descAl, descAt, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_zlap2tile", "INPLACE translation not supported yet");
        descAt->mat = A;

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & MorseDescInput ) {
            /* MORSE_zgecfi_Async( lm, ln, A, MorseCM, mb, nb, */
            /*                     MorseCCRB, mb, nb, seq, req ); */
        }
        return MORSE_ERR_NOT_SUPPORTED;
    }

    return MORSE_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
morse_ztile2lap( MORSE_context_t *morse, MORSE_desc_t *descAl, MORSE_desc_t *descAt,
                 MORSE_enum mode, MORSE_enum uplo, MORSE_sequence_t *seq, MORSE_request_t *req )
{
    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        if ( mode & MorseDescOutput ) {
            morse_pzlacpy( uplo, descAt, descAl, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_ztile2lap", "INPLACE translation not supported yet");
        if ( mode & MorseDescOutput ) {
            /* MORSE_zgecfi_Async( descAl->lm, descAl->ln, descAl->mat, */
            /*                     MorseCCRB, descAl->mb, descAl->nb,   */
            /*                     MorseCM, descAl->mb, descAl->nb, seq, req ); */
        }
        return MORSE_ERR_NOT_SUPPORTED;
    }
    RUNTIME_desc_flush( descAl, seq );
    RUNTIME_desc_flush( descAt, seq );

    return MORSE_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
morse_ztile2lap_cleanup( MORSE_context_t *morse, MORSE_desc_t *descAl, MORSE_desc_t *descAt )
{
    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        morse_desc_mat_free( descAt );
    }
    RUNTIME_desc_destroy( descAl );
    RUNTIME_desc_destroy( descAt );
}
