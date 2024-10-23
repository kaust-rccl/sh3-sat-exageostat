###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  MORSE example routines
#  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
#  University of Bordeaux, Bordeaux INP
#
# @version 1.0.0
#  @author Florent Pruvost
#  @date 2016-08-23
#
###
if (NOT CHAMELEON_SCHED_STARPU)
  message(ERROR "This directory should not be included if CHAMELEON_SCHED_STARPU is not enabled")
endif()

if(CHAMELEON_SIMULATION)
  message(ERROR "This directory should not be included if CHAMELEON_SIMULATION is enabled")
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# list of sources
set(OOC_SOURCES
    out_of_core.c
   )

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_ooc)
list(APPEND libs_for_ooc
     chameleon
     ${STARPU_LIBRARIES_DEP}
)
link_directories(${STARPU_LIBRARY_DIRS})

list(APPEND libs_for_ooc
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBRARIES}
  )

link_directories(${LAPACKE_LIBRARY_DIRS})
link_directories(${LAPACK_LIBRARY_DIRS})
link_directories(${BLAS_LIBRARY_DIRS})

# message(STATUS "libs examples: ${libs_for_ooc}")
foreach(_ooc ${OOC_SOURCES})
    get_filename_component(_name_exe ${_ooc} NAME_WE)
    add_executable(${_name_exe} ${_ooc})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -undefined dynamic_lookup" )
    endif()
    target_link_libraries(${_name_exe} ${libs_for_ooc})
    install(TARGETS ${_name_exe}
            DESTINATION bin/example/out_of_core)
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
