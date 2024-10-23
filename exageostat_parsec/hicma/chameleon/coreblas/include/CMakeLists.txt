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

# Generate header files
# ---------------------
set(COREBLAS_HDRS_GENERATED "")
set(ZHDR
  coreblas/coreblas_z.h
  coreblas/coreblas_zc.h
)
precisions_rules_py(
  COREBLAS_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  coreblas
  PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
set(COREBLAS_HDRS
    coreblas.h
    coreblas/cblas.h
    coreblas/lapacke.h
    coreblas/lapacke_config.h
    coreblas/lapacke_mangling.h
    )

# Add generated headers
# ---------------------
foreach( hdr_file ${COREBLAS_HDRS_GENERATED} )
  list(APPEND COREBLAS_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(coreblas_include ALL SOURCES ${COREBLAS_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};coreblas_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES coreblas.h
         DESTINATION include )

install( FILES ${COREBLAS_HDRS}
         DESTINATION include/coreblas )

###
### END CMakeLists.txt
###
