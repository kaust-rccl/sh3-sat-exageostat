###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Florent Pruvost
#  @date 2012-07-13
#
###

# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(CHAMELEON_HDRS_GENERATED "")
set(ZHDR
    chameleon/morse_z.h
    chameleon/morse_zc.h
    chameleon/morse_tasks_z.h
    chameleon/morse_tasks_zc.h
)

precisions_rules_py(CHAMELEON_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  "chameleon"
  PRECISIONS "z;c;d;s;zc;ds" )

configure_file("chameleon/chameleon_config.h.in"
               "chameleon/chameleon_config.h"
               @ONLY)

# Define the list of headers
# --------------------------
set(CHAMELEON_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/chameleon/chameleon_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/chameleon/morse_mangling.h
  chameleon/morse_constants.h
  chameleon/morse_kernels.h
  chameleon/morse_runtime.h
  chameleon/morse_simulate.h
  chameleon/morse_struct.h
  chameleon/morse_tasks.h
  chameleon/morse_types.h
  )

# Add generated headers
# ---------------------
foreach( hdr_file ${CHAMELEON_HDRS_GENERATED} )
  list(APPEND CHAMELEON_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(chameleon_include
  ALL SOURCES morse.h ${CHAMELEON_HDRS} )
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};chameleon_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES morse.h
         DESTINATION include )

install( FILES ${CHAMELEON_HDRS}
         DESTINATION include/chameleon )

###
### END CMakeLists.txt
###
