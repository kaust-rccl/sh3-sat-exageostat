###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
# @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 2012-07-13
#
###

# Generate the morse sources for all possible precisions
# ------------------------------------------------------

set(COREBLAS_SRCS_GENERATED "")
set(ZSRC
    core_dzasum.c
    core_zaxpy.c
    core_zgeadd.c
    core_zlascal.c
    core_zgelqt.c
    core_zgemm.c
    core_zgeqrt.c
    core_zgesplit.c
    core_zgessm.c
    core_zgessq.c
    core_zgetf2_nopiv.c
    core_zgetrf.c
    core_zgetrf_incpiv.c
    core_zgetrf_nopiv.c
    core_zhe2ge.c
    core_zherfb.c
    core_zhemm.c
    core_zher2k.c
    core_zherk.c
    core_zhessq.c
    core_zlacpy.c
    core_zlag2c.c
    core_zlange.c
    core_zlanhe.c
    core_zlansy.c
    core_zlantr.c
    core_zlaset2.c
    core_zlaset.c
    core_zlatro.c
    core_zlauum.c
    core_zpamm.c
    core_zparfb.c
    core_zpemv.c
    core_zplghe.c
    core_zplgsy.c
    core_zplrnt.c
    core_zpotrf.c
    core_zssssm.c
    core_zsymm.c
    core_zsyr2k.c
    core_zsyrk.c
    core_zsyssq.c
    core_zsytf2_nopiv.c
    core_ztplqt.c
    core_ztpqrt.c
    core_ztpmlqt.c
    core_ztpmqrt.c
    core_ztradd.c
    core_ztrasm.c
    core_ztrmm.c
    core_ztrsm.c
    core_ztrssq.c
    core_ztrtri.c
    core_ztslqt.c
    core_ztsmlq.c
    core_ztsmqr.c
    core_ztsmlq_hetra1.c
    core_ztsmqr_hetra1.c
    core_ztsqrt.c
    core_ztstrf.c
    core_zttlqt.c
    core_zttmlq.c
    core_zttmqr.c
    core_zttqrt.c
    core_zunmlq.c
    core_zunmqr.c
    )

precisions_rules_py(COREBLAS_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

set(COREBLAS_SRCS
    global.c
    ${COREBLAS_SRCS_GENERATED}
    )

# Force generation of sources
# ---------------------------
add_custom_target(coreblas_sources ALL SOURCES ${COREBLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};coreblas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(coreblas ${COREBLAS_SRCS})
list(INSERT CHAMELEON_DEP 0 -lcoreblas)
add_dependencies(coreblas coreblas_include coreblas_sources)
set_property(TARGET coreblas PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET coreblas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(coreblas
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  )

# Installation
# ------------
install(TARGETS coreblas
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###
