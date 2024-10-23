/**
 *
 * @file runtime_zprofiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU MORSE_Complex64_t kernel progiling
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

void RUNTIME_zdisplay_allprofile()
{

    /* BLAS 3 */
    profiling_display_zgemm_info();
#if defined(PRECISION_z) || defined(PRECISION_c)
    profiling_display_zhemm_info();
    profiling_display_zher2k_info();
    profiling_display_zherk_info();
    profiling_display_zsytrf_nopiv_info();
#endif
    profiling_display_zsymm_info();
    profiling_display_zsyr2k_info();
    profiling_display_zsyrk_info();
    profiling_display_ztrmm_info();
    profiling_display_ztrsm_info();

    /* Lapack */
    profiling_display_zgelqt_info();
    profiling_display_zgeqrt_info();
    profiling_display_zgessm_info();
    profiling_display_zgetrf_info();
    profiling_display_zgetrf_incpiv_info();
    profiling_display_zgetrf_nopiv_info();
    profiling_display_zlauum_info();
    profiling_display_zpotrf_info();
    profiling_display_zssssm_info();
    profiling_display_ztrtri_info();
    profiling_display_ztslqt_info();
    profiling_display_ztsmqr_info();
    profiling_display_ztsqrt_info();
    profiling_display_ztstrf_info();
    profiling_display_zttlqt_info();
    profiling_display_zttmlq_info();
    profiling_display_zttmqr_info();
    profiling_display_zttqrt_info();
    profiling_display_zunmlq_info();
    profiling_display_zunmqr_info();

    profiling_display_zlange_info();
}

void RUNTIME_zdisplay_oneprofile( MORSE_kernel_t kernel )
{
    switch( kernel ) {
        /* Blas 3 */
    case MORSE_GEMM:         profiling_display_zgemm_info();         break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case MORSE_HEMM:         profiling_display_zhemm_info();         break;
    case MORSE_HER2K:        profiling_display_zher2k_info();        break;
    case MORSE_HERK:         profiling_display_zherk_info();         break;
    case MORSE_SYTRF_NOPIV:  profiling_display_zsytrf_nopiv_info();  break;
#endif
    case MORSE_SYMM:         profiling_display_zsymm_info();         break;
    case MORSE_SYR2K:        profiling_display_zsyr2k_info();        break;
    case MORSE_SYRK:         profiling_display_zsyrk_info();         break;
    case MORSE_TRMM:         profiling_display_ztrmm_info();         break; 
    case MORSE_TRSM:         profiling_display_ztrsm_info();         break;

        /* Lapack */
    case MORSE_GELQT:        profiling_display_zgelqt_info();        break;
    case MORSE_GEQRT:        profiling_display_zgeqrt_info();        break;
    case MORSE_GESSM:        profiling_display_zgessm_info();        break;
    case MORSE_GETRF:        profiling_display_zgetrf_info();        break;
    case MORSE_GETRF_INCPIV: profiling_display_zgetrf_incpiv_info(); break;
    case MORSE_GETRF_NOPIV:  profiling_display_zgetrf_nopiv_info();  break;
    case MORSE_LAUUM:        profiling_display_zlauum_info();        break;
    case MORSE_POTRF:        profiling_display_zpotrf_info();        break;
    case MORSE_SSSSM:        profiling_display_zssssm_info();        break;
    case MORSE_TRTRI:        profiling_display_ztrtri_info();        break;
    case MORSE_TSLQT:        profiling_display_ztslqt_info();        break;
    case MORSE_TSMQR:        profiling_display_ztsmqr_info();        break;
    case MORSE_TSQRT:        profiling_display_ztsqrt_info();        break;
    case MORSE_TSTRF:        profiling_display_ztstrf_info();        break;
    case MORSE_TTLQT:        profiling_display_zttlqt_info();        break;
    case MORSE_TTMLQ:        profiling_display_zttmlq_info();        break;
    case MORSE_TTMQR:        profiling_display_zttmqr_info();        break;
    case MORSE_TTQRT:        profiling_display_zttqrt_info();        break;
    case MORSE_UNMLQ:        profiling_display_zunmlq_info();        break;
    case MORSE_UNMQR:        profiling_display_zunmqr_info();        break;

    case MORSE_LANGE:        profiling_display_zlange_info();        break;

    default:
        return;
    }
}

