###
#
# @file CMakeLists.txt
#
# @copyright 2009-2015 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
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
cmake_minimum_required(VERSION 2.8)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
#set(ZHDR
#)
#
#precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
#                    PRECISIONS "s;d;c;z;ds;zc"
#                    TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_parsec.h
  )

# Add generated headers
# ---------------------
foreach( hdr_file ${RUNTIME_HDRS_GENERATED} )
  list(APPEND RUNTIME_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_parsec_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/parsec )

# Generate the Chameleon common for all possible precisions
# ---------------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC
  control/runtime_zprofiling.c
  control/runtime_zlocality.c
  )

precisions_rules_py(
  RUNTIME_COMMON_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "control")

set(RUNTIME_COMMON
  control/runtime_async.c
  control/runtime_context.c
  control/runtime_control.c
  control/runtime_descriptor.c
  control/runtime_options.c
  control/runtime_profiling.c
  ${RUNTIME_COMMON_GENERATED}
  )

# Generate the Chameleon sources for all possible precisions
# ----------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
  ${CODELETS_ZSRC}
  )

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  )

# Force generation of sources
# ---------------------------
add_custom_target(parsec_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_parsec_include;parsec_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_parsec ${RUNTIME_SRCS})
set_property(TARGET chameleon_parsec PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_parsec PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_parsec
  ${PARSEC_LIBRARIES_DEP})
if (NOT CHAMELEON_SIMULATION)
  target_link_libraries(chameleon_parsec
    coreblas)
  if(CHAMELEON_USE_CUDA)
    target_link_libraries(chameleon_parsec
      cudablas)
  endif(CHAMELEON_USE_CUDA)
endif(NOT CHAMELEON_SIMULATION)

add_dependencies(chameleon_parsec
  chameleon_include
  control_include
  runtime_parsec_include
  parsec_sources
)

if (NOT CHAMELEON_SIMULATION)
  add_dependencies(chameleon_parsec coreblas_include)
  if (CHAMELEON_USE_CUDA)
    add_dependencies(chameleon_parsec cudablas_include)
  endif()
endif()

# installation
# ------------
install(TARGETS chameleon_parsec
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###
