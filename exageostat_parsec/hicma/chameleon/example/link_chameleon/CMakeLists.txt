### Main CMakeLists.txt for project link_chameleon

cmake_minimum_required(VERSION 2.8)
project(LINK_CHAMELEON Fortran C CXX)

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
set( MORSE_DISTRIB_DIR "" CACHE PATH "Directory of MORSE distribution")
option( MORSE_CHAMELEON_USE_QUARK "CHAMELEON to be linked with QUARK runtime" OFF)
option( MORSE_CHAMELEON_USE_STARPU "CHAMELEON to be linked with StarPU runtime (default)" ON)

# Both options cannot be ON at the same time: it is either QUARK or StarPU
if (MORSE_CHAMELEON_USE_QUARK)
    set(MORSE_CHAMELEON_USE_STARPU OFF)
elseif (MORSE_CHAMELEON_USE_STARPU)
    set(MORSE_CHAMELEON_USE_QUARK OFF)
endif()

if (MORSE_DISTRIB_DIR)
    set( MORSE_CMAKE_MODULE_DIR "${MORSE_DISTRIB_DIR}/morse_cmake/modules" CACHE PATH
         "Directory where to find MORSE CMake modules (morse_cmake/modules)")
    list(APPEND CMAKE_MODULE_PATH "${MORSE_CMAKE_MODULE_DIR}")
    list(APPEND CMAKE_MODULE_PATH "${MORSE_CMAKE_MODULE_DIR}/find")
    message(STATUS "MORSE_CMAKE_MODULE_DIR : ${MORSE_CMAKE_MODULE_DIR}/morse_cmake/modules")
    include(MorseInit)

    # Detect CHAMELEON
    if (MORSE_CHAMELEON_USE_QUARK)
        find_package(CHAMELEON COMPONENTS QUARK)
    else()
        find_package(CHAMELEON COMPONENTS STARPU MPI CUDA FXT)
    endif()
    if (CHAMELEON_FOUND)
        link_directories(${CHAMELEON_LIBRARY_DIRS_DEP})
        include_directories(${CHAMELEON_INCLUDE_DIRS_DEP})
    else()
        message(FATAL_ERROR "Chameleon libraries not found in your environment."
        "Please update your environment variables that help find it"
        "(INCLUDE, LD_LIBRARY_PATH) or give its install directory setting CHAMELEON_DIR.")
    endif()

    # link_chameleon exe
    add_executable(link_chameleon_c link_chameleon.c)
    target_link_libraries(link_chameleon_c ${CHAMELEON_LIBRARIES_DEP})
    add_executable(link_chameleon_f link_chameleon.f90)
    target_link_libraries(link_chameleon_f ${CHAMELEON_LIBRARIES_DEP})
else()
    message(STATUS "MORSE_DISTRIB_DIR is not set")
    message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
endif()

###
### END CMakeLists.txt
###
