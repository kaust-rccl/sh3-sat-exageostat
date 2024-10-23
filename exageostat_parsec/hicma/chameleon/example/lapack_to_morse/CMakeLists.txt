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
#  @date 2014-10-13
#
###
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# compilation for other sources step1 and >
set(LTM_SOURCES
    step1.c
    step2.c
    step3.c
    step4.c
    step5.c
    step6.c
    step7.c
   )

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_ltm)
list(APPEND libs_for_ltm
  chameleon
  coreblas
  ${LAPACKE_LIBRARIES}
)

# specific compilation for step0 because we potentially want to use
# multithreaded BLAS and LAPACK libraries for this step
unset(libs_for_step0)
list(APPEND libs_for_step0
  ${libs_for_ltm}
  ${CBLAS_LIBRARIES}
  )

if(LAPACK_PAR_LIBRARIES)
  list(APPEND libs_for_step0
    ${LAPACK_PAR_LIBRARIES}
    )
else()
  list(APPEND libs_for_step0
    ${LAPACK_SEQ_LIBRARIES}
    )
endif()

if(BLAS_PAR_LIBRARIES)
  list(APPEND libs_for_step0
    ${BLAS_PAR_LIBRARIES}
    )
else()
  list(APPEND libs_for_step0
    ${BLAS_SEQ_LIBRARIES}
    )
endif()

list(APPEND libs_for_step0
  ${EXTRA_LIBRARIES}
  )

list(APPEND libs_for_ltm
  ${LAPACK_SEQ_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBRARIES}
  )

link_directories(${LAPACKE_LIBRARY_DIRS})
link_directories(${LAPACK_LIBRARY_DIRS})
link_directories(${CBLAS_LIBRARY_DIRS})
link_directories(${BLAS_LIBRARY_DIRS})

# message(STATUS "libs examples: ${libs_for_ltm}")
foreach(_ltm ${LTM_SOURCES})
    get_filename_component(_name_exe ${_ltm} NAME_WE)
    add_executable(${_name_exe} ${_ltm})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_ltm})
    install(TARGETS ${_name_exe}
            DESTINATION bin/example/lapack_to_morse)
endforeach()

add_executable(step0 step0.c)
set_property(TARGET step0 PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(step0 ${libs_for_step0})
install(TARGETS step0
        DESTINATION bin/example/lapack_to_morse)

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
