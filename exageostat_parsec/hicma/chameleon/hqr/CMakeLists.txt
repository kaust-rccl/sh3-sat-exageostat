###
#
#  @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 0.1.0
#  @author Mathieu Faverge
#  @date 2017-04-26
#
###
cmake_minimum_required (VERSION 3.0)
project (HQR C)

# Check if compiled independtly or within another project
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set( BUILD_SUBPROJECT OFF )
else()
  set( BUILD_SUBPROJECT ON )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
include(GenHQRPkgConfig)

# The current version number
set (HQR_VERSION_MAJOR 0)
set (HQR_VERSION_MINOR 1)
set (HQR_VERSION_MICRO 0)

set( HQR_VERSION "${HQR_VERSION_MAJOR}.${HQR_VERSION_MINOR}.${HQR_VERSION_MICRO}" )

### Misc options
if ( NOT BUILD_SUBPROJECT )
  option(BUILD_SHARED_LIBS
    "Build shared libraries" OFF)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
endif()

set(hdrs
  include/libhqr.h
  include/libhqr_internal.h
  include/libhqr_queue.h
)
set(srcs
  # Low level tree functions
  src/low_flat.c
  src/low_binary.c
  src/low_fibonacci.c
  src/low_greedy.c
  src/low_greedy1p.c
  src/low_adaptiv.c
  # High level tree functions
  src/high_flat.c
  src/high_binary.c
  src/high_fibonacci.c
  src/high_greedy.c
  # Direct access tree functions
  src/systolic.c
  src/svd.c
  src/hqr.c
  src/mtxtree.c
  # Others
  src/check.c
  src/gendot.c
  src/gensvg.c
  src/queue.c
)

include_directories(include)
add_library(hqr ${srcs})
set_property(TARGET hqr PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries( hqr m )

add_subdirectory(testings)

install(FILES
  include/libhqr.h
  DESTINATION include )

install(TARGETS hqr
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

generate_hqr_pkgconfig_file()

#-- Add a custom target to generate tags
add_custom_target (tags
  COMMAND "git ls-files | xargs etags" )
