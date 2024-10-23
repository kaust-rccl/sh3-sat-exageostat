/**
 *
 * @file runtime_zlocality.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU MORSE_Complex64_t kernel locality management
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#ifdef CHAMELEON_USE_CUDA
/* Only codelets with multiple choices are present here */
void RUNTIME_zlocality_allrestrict( uint32_t where )
{

    /* Blas 3 */
    cl_zgemm_restrict_where( where );
#if defined(PRECISION_z) || defined(PRECISION_c)
    cl_zhemm_restrict_where( where );
    cl_zher2k_restrict_where( where );
    cl_zherk_restrict_where( where );
    cl_zsytrf_nopiv_restrict_where( where );
#endif
    cl_zsymm_restrict_where( where );
    cl_zsyr2k_restrict_where( where );
    cl_zsyrk_restrict_where( where );
    cl_ztrmm_restrict_where( where );
    cl_ztrsm_restrict_where( where );

    /*
     * Lapack
     */
    /* Cholesky */
    cl_zpotrf_restrict_where( where );
    cl_zlauum_restrict_where( where );
    cl_ztrtri_restrict_where( where );

    /* LU */
    cl_zgetrf_incpiv_restrict_where( where );
    cl_zgetrf_nopiv_restrict_where( where );
    cl_zgessm_restrict_where( where );
    cl_zssssm_restrict_where( where );
    cl_ztstrf_restrict_where( where );

    /* QR */
    cl_zgeqrt_restrict_where( where );
    cl_ztsqrt_restrict_where( where );
    cl_zunmqr_restrict_where( where );
    cl_ztsmqr_restrict_where( where );

    /* QR-RH */
/*     cl_zttqrt_restrict_where( where ); */
/*     cl_zttmqr_restrict_where( where ); */

    /* LQ */
   cl_zgelqt_restrict_where( where );
   cl_ztslqt_restrict_where( where );
   cl_zunmlq_restrict_where( where );
   cl_ztsmlq_restrict_where( where );

    /* LQ-RH */
/*     cl_zttlqt_restrict_where( where ); */
/*     cl_zttmlq_restrict_where( where ); */

}

void RUNTIME_zlocality_onerestrict( MORSE_kernel_t kernel, uint32_t where )
{
    switch( kernel ) {
    /* Blas 3 */
    case MORSE_GEMM:   cl_zgemm_restrict_where( where );  break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case MORSE_HEMM:   cl_zhemm_restrict_where( where );  break;
    case MORSE_HER2K:  cl_zher2k_restrict_where( where ); break;
    case MORSE_HERK:   cl_zherk_restrict_where( where );  break;
    case MORSE_SYTRF_NOPIV: cl_zsytrf_nopiv_restrict_where( where );  break;
#endif
    case MORSE_SYMM:   cl_zhemm_restrict_where( where );  break;
    case MORSE_SYR2K:  cl_zher2k_restrict_where( where ); break;
    case MORSE_SYRK:   cl_zherk_restrict_where( where );  break;
    case MORSE_TRMM:   cl_ztrmm_restrict_where( where );  break;
    case MORSE_TRSM:   cl_ztrsm_restrict_where( where );  break;

    /*
     * Lapack
     */
    /* Cholesky */
    case MORSE_POTRF:  cl_zpotrf_restrict_where( where ); break;
    case MORSE_LAUUM:  cl_zlauum_restrict_where( where ); break;
    case MORSE_TRTRI:  cl_ztrtri_restrict_where( where ); break;

    /* LU */
    case MORSE_GETRF_INCPIV: cl_zgetrf_incpiv_restrict_where( where ); break;
    case MORSE_GETRF_NOPIV: cl_zgetrf_nopiv_restrict_where( where ); break;
    case MORSE_GESSM:  cl_zgessm_restrict_where( where ); break;
    case MORSE_SSSSM:  cl_zssssm_restrict_where( where ); break;
    case MORSE_TSTRF:  cl_ztstrf_restrict_where( where ); break;

    /* QR */
    case MORSE_GEQRT:  cl_zgeqrt_restrict_where( where ); break;
    case MORSE_UNMQR:  cl_zunmqr_restrict_where( where ); break;
    case MORSE_TSMQR:  cl_ztsmqr_restrict_where( where ); break;
    case MORSE_TSQRT:  cl_ztsqrt_restrict_where( where ); break;

    /* QR-RH */
/*     case MORSE_TTMQR:  cl_zttmqr_restrict_where( where ); break; */
/*     case MORSE_TTQRT:  cl_zttqrt_restrict_where( where ); break; */

    /* LQ */
   case MORSE_GELQT:  cl_zgelqt_restrict_where( where ); break;
   case MORSE_UNMLQ:  cl_zunmlq_restrict_where( where ); break;
   case MORSE_TSMLQ:  cl_ztsmlq_restrict_where( where ); break;
   case MORSE_TSLQT:  cl_ztslqt_restrict_where( where ); break;

    /* LQ-RH */
/*     case MORSE_TTMLQ:  cl_zttmlq_restrict_where( where ); break; */
/*     case MORSE_TTLQT:  cl_zttlqt_restrict_where( where ); break; */

    default:
      return;
    }
}

void RUNTIME_zlocality_allrestore( )
{
    /* Blas 3 */
    cl_zgemm_restore_where();
#if defined(PRECISION_z) || defined(PRECISION_c)
    cl_zhemm_restore_where();
    cl_zher2k_restore_where();
    cl_zherk_restore_where();
    cl_zsytrf_nopiv_restore_where();
#endif
    cl_zsymm_restore_where();
    cl_zsyr2k_restore_where();
    cl_zsyrk_restore_where();
    cl_ztrmm_restore_where();
    cl_ztrsm_restore_where();

    /*
     * Lapack
     */
    /* Cholesky */
    cl_zpotrf_restore_where();
    cl_zlauum_restore_where();
    cl_ztrtri_restore_where();

    /* LU incpiv */
    cl_zgetrf_incpiv_restore_where();
    cl_zgessm_restore_where();
    cl_zssssm_restore_where();
    cl_ztstrf_restore_where();

    /* QR */
    cl_zgeqrt_restore_where();
    cl_ztsqrt_restore_where();
    cl_zunmqr_restore_where();
    cl_ztsmqr_restore_where();

    /* QR-RH */
/*     cl_zttqrt_restore_where(); */
/*     cl_zttmqr_restore_where(); */

    /* LQ */
   cl_zgelqt_restore_where();
   cl_ztslqt_restore_where();
   cl_zunmlq_restore_where();
   cl_ztsmlq_restore_where();

    /* LQ-RH */
/*     cl_zttlqt_restore_where(); */
/*     cl_zttmlq_restore_where(); */

}

void RUNTIME_zlocality_onerestore( MORSE_kernel_t kernel )
{
    switch( kernel ) {
    /* Blas 3 */
    case MORSE_GEMM:   cl_zgemm_restore_where();  break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case MORSE_HEMM:   cl_zhemm_restore_where();  break;
    case MORSE_HER2K:  cl_zher2k_restore_where(); break;
    case MORSE_HERK:   cl_zherk_restore_where();  break;
    case MORSE_SYTRF_NOPIV: cl_zsytrf_nopiv_restore_where();  break;
#endif
    case MORSE_SYMM:   cl_zhemm_restore_where();  break;
    case MORSE_SYR2K:  cl_zher2k_restore_where(); break;
    case MORSE_SYRK:   cl_zherk_restore_where();  break;
    case MORSE_TRMM:   cl_ztrmm_restore_where();  break;
    case MORSE_TRSM:   cl_ztrsm_restore_where();  break;

    /*
     * Lapack
     */
    /* Cholesky */
    case MORSE_POTRF:  cl_zpotrf_restore_where(); break;
    case MORSE_LAUUM:  cl_zlauum_restore_where(); break;
    case MORSE_TRTRI:  cl_ztrtri_restore_where(); break;

    /* LU */
    case MORSE_GETRF_INCPIV: cl_zgetrf_incpiv_restore_where(); break;
    case MORSE_GETRF_NOPIV: cl_zgetrf_nopiv_restore_where(); break;
    case MORSE_GESSM:  cl_zgessm_restore_where(); break;
    case MORSE_SSSSM:  cl_zssssm_restore_where(); break;
    case MORSE_TSTRF:  cl_ztstrf_restore_where(); break;

    /* QR */
    case MORSE_GEQRT:  cl_zgeqrt_restore_where(); break;
    case MORSE_UNMQR:  cl_zunmqr_restore_where(); break;
    case MORSE_TSMQR:  cl_ztsmqr_restore_where(); break;
    case MORSE_TSQRT:  cl_ztsqrt_restore_where(); break;

    /* QR-RH */
/*     case MORSE_TTMQR:  cl_zttmqr_restore_where(); break; */
/*     case MORSE_TTQRT:  cl_zttqrt_restore_where(); break; */

    /* LQ */
   case MORSE_GELQT:  cl_zgelqt_restore_where(); break;
   case MORSE_UNMLQ:  cl_zunmlq_restore_where(); break;
   case MORSE_TSMLQ:  cl_ztsmlq_restore_where(); break;
   case MORSE_TSLQT:  cl_ztslqt_restore_where(); break;

    /* LQ-RH */
/*     case MORSE_TTMLQ:  cl_zttmlq_restore_where(); break; */
/*     case MORSE_TTLQT:  cl_zttlqt_restore_where(); break; */

    default:
      return;
    }
}
#else

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
}

void RUNTIME_zlocality_onerestrict( MORSE_kernel_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
}

void RUNTIME_zlocality_allrestore( ) {}

void RUNTIME_zlocality_onerestore( MORSE_kernel_t kernel )
{
    (void)kernel;
}

#endif
