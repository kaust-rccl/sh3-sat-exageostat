cmake_policy(SET CMP0074 NEW)
include(FindPkgConfig)
# EXAGEOSTAT depends on PaRSEC and DPLASMA 
# -------------------------------
if( EXAGEOSTAT_USE_DPLASMA )
add_definitions(-DEXAGEOSTAT_USE_DPLASMA)
# DPLASMA
pkg_search_module(DPLASMA dplasma REQUIRED)
if( DPLASMA_FOUND )
    include_directories(${DPLASMA_INCLUDE_DIRS})
    link_directories(${DPLASMA_LIBRARY_DIRS})
    list(INSERT EXAGEOSTAT_DEP 0 ${DPLASMA_LIBRARIES})
    message(STATUS "DPLASMA ${DPLASMA_VERSION} found")
else()
    message("DPLASMA is NOT found, skipping it")
    set(DPLASMA OFF)
endif( DPLASMA_FOUND )

#
#set environment variable PaRSEC_ROOT to the root of PaRSEC installation 
find_package(PaRSEC REQUIRED) 
if( NOT PaRSEC_FOUND)
    message(FATAL_ERROR "PaRSEC: NOT FOUND.") 
endif( NOT PaRSEC_FOUND )


include(ParsecCompilePTG)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma/parsec")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma/src")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma/src/include")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma/build/src/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../dplasma/src/cores")

if( EXAGEOSTAT_ON_SUMMIT )
    set( SRC
        mix_precision_internal.c
        ../../examples/examples.c
	my_dtrsm_wrapper.c
        hsdpotrf_L_wrapper.c
        hsdpotrf_L_band_wrapper.c
        spotrf_L_wrapper.c
        hsdpotrf_L_trsm_convert_wrapper.c 
        hsdpotrf_L_4taskclass_wrapper.c 
        shpotrf_L_wrapper.c
        /ccs/home/cql0536/opt/plasma/lib/liblapack.a
       )
else()
    set( SRC
        mix_precision_internal.c
        ../../examples/examples.c
	my_dtrsm_wrapper.c
        hsdpotrf_L_wrapper.c
        hsdpotrf_L_band_wrapper.c
        spotrf_L_wrapper.c
        hsdpotrf_L_trsm_convert_wrapper.c
        hsdpotrf_L_4taskclass_wrapper.c 
        shpotrf_L_wrapper.c
       )
endif( EXAGEOSTAT_ON_SUMMIT )

set( JDF_DOUBLE
  dmatrix_generation.jdf
  dZ_generation.jdf
  dmatrix_det.jdf
  get_zobs.jdf
  dmatrix_scale.jdf
  dmatrix_shift.jdf
  Lapack_to_Tile.jdf
  mse_calculation.jdf
  dmatrix_init.jdf
  dmatrix_sum.jdf
  dmatrix_dot.jdf
  dZ_sum.jdf
  dmatrix_avg.jdf
  my_dtrmv.jdf
  my_warmup.jdf
  dmatrix_set_diagonal.jdf
  my_dtrsm_LLN.jdf
  MLE_dtrace.jdf
   )

set( JDF_SINGLE
  smatrix_generation.jdf
  spotrf_L.jdf
   )

set( JDF_MIX
  ${JDF_DOUBLE}
  ${JDF_SINGLE}
  hsdpotrf_L.jdf
  hsdpotrf_L_band.jdf
  hsdpotrf_L_trsm_convert.jdf
  hsdpotrf_L_4taskclass.jdf
  shpotrf_L.jdf
  band_convert_s2d.jdf
  band_convert_d2s.jdf
  band_allocate.jdf
  band_elem_set.jdf
    )

#if( EXAGEOSTAT_USE_CUDA )
if( CUDA_FOUND )
    cuda_compile(cuda_objs1 mix_precision_internal.cu)
    cuda_compile(cuda_objs2 ../../cudacore/compute/cuda_zcmg.cu)
    cuda_compile(cuda_objs3 ../../cudacore/compute/cuda_zscale.cu)
    SET(SRC ${SRC} ${cuda_objs1} ${cuda_objs2} ${cuda_objs3})
endif( CUDA_FOUND )

set( CMAKE_C_FLAGS 
  ${CMAKE_C_FLAGS}
  -D_GNU_SOURCE
   )

#file(DOWNLOAD "https://docs.google.com/uc?export=download&id=0B6mGsgjYG40kOGdNc2cyZzBjelk" "${CMAKE_CURRENT_BINARY_DIR}/marc-sample.tar.gz" EXPECTED_HASH MD5=2978769f72250eeff401f725158587a3)
#file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1Ib7W9Y5z7nqL2t9C8bbXFxL3J0E9f34U" "${CMAKE_CURRENT_BINARY_DIR}/wind-sample.tar.gz" EXPECTED_HASH MD5=2103ba786a3b07987ef70b4a37a94e9b)

#add_custom_target( parsec_extract_tests_data ALL
#                    DEPENDS ${exe_names}
#                    COMMAND cmake -E tar zxf marc-sample.tar.gz
#                    COMMAND cmake -E tar zxf wind-sample.tar.gz
#  )


add_executable(testing_hsdpotrf ${SRC} testing_hsdpotrf.c)
if( EXAGEOSTAT_ON_SUMMIT )
    if( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
    else( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf exageostat parsec dplasma stdc++ lapacke)
    endif( CUDA_FOUND )
else()
    if( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf exageostat parsec dplasma cublas cudart cusolver stdc++)
    else( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf exageostat parsec dplasma stdc++)
    endif( CUDA_FOUND )
endif( EXAGEOSTAT_ON_SUMMIT )
target_ptg_sources(testing_hsdpotrf PRIVATE ${JDF_MIX} )


#add_executable(testing_hsdpotrf_real_1dlocs ${SRC} testing_hsdpotrf_real_1dlocs.c)
#if( EXAGEOSTAT_ON_SUMMIT )
#    if( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_1dlocs exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
#    else( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_1dlocs exageostat parsec dplasma stdc++ lapacke)
#    endif( CUDA_FOUND )
#else()
#    if( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_1dlocs exageostat parsec dplasma cublas cudart cusolver stdc++)
#    else( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_1dlocs exageostat parsec dplasma stdc++)
#    endif( CUDA_FOUND )
#endif( EXAGEOSTAT_ON_SUMMIT )
#target_ptg_sources(testing_hsdpotrf_real_1dlocs PRIVATE ${JDF_MIX} )

#add_executable(testing_hsdpotrf_real_2dlocs ${SRC} testing_hsdpotrf_real_2dlocs.c)
#if( EXAGEOSTAT_ON_SUMMIT )
#    if( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_2dlocs exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
#    else( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_2dlocs exageostat parsec dplasma stdc++ lapacke)
#    endif( CUDA_FOUND )
#else()
#    if( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_2dlocs exageostat parsec dplasma cublas cudart cusolver stdc++)
#    else( CUDA_FOUND )
#        target_link_libraries(testing_hsdpotrf_real_2dlocs exageostat parsec dplasma stdc++)
#    endif( CUDA_FOUND )
#endif( EXAGEOSTAT_ON_SUMMIT )
#target_ptg_sources(testing_hsdpotrf_real_2dlocs PRIVATE ${JDF_MIX} )

add_executable(testing_performance ${SRC} testing_performance.c)
if( EXAGEOSTAT_ON_SUMMIT )
    if( CUDA_FOUND )
        target_link_libraries(testing_performance exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
    else( CUDA_FOUND )
        target_link_libraries(testing_performance exageostat parsec dplasma stdc++ lapacke)
    endif( CUDA_FOUND )
else()
    if( CUDA_FOUND )
        target_link_libraries(testing_performance exageostat parsec dplasma cublas cudart cusolver stdc++)
    else( CUDA_FOUND )
        target_link_libraries(testing_performance exageostat parsec dplasma stdc++)
    endif( CUDA_FOUND )
endif( EXAGEOSTAT_ON_SUMMIT )
target_ptg_sources(testing_performance PRIVATE ${JDF_MIX} )

add_executable(testing_hsdpotrf_csv ${SRC} testing_hsdpotrf_csv.c)
if( EXAGEOSTAT_ON_SUMMIT )
    if( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf_csv exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
    else( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf_csv exageostat parsec dplasma stdc++ lapacke)
    endif( CUDA_FOUND )
else()
    if( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf_csv exageostat parsec dplasma cublas cudart cusolver stdc++)
    else( CUDA_FOUND )
        target_link_libraries(testing_hsdpotrf_csv exageostat parsec dplasma stdc++)
    endif( CUDA_FOUND )
endif( EXAGEOSTAT_ON_SUMMIT )
target_ptg_sources(testing_hsdpotrf_csv PRIVATE ${JDF_MIX} )

add_executable(testing_performance_band ${SRC} testing_performance_band.c)
if( EXAGEOSTAT_ON_SUMMIT )
    if( CUDA_FOUND )
        target_link_libraries(testing_performance_band exageostat parsec dplasma cublas cudart cusolver stdc++ lapacke)
    else( CUDA_FOUND )
        target_link_libraries(testing_performance_band exageostat parsec dplasma stdc++ lapacke)
    endif( CUDA_FOUND )
else()
    if( CUDA_FOUND )
        target_link_libraries(testing_performance_band exageostat parsec dplasma cublas cudart cusolver stdc++)
    else( CUDA_FOUND )
        target_link_libraries(testing_performance_band exageostat parsec dplasma stdc++)
    endif( CUDA_FOUND )
endif( EXAGEOSTAT_ON_SUMMIT )
target_ptg_sources(testing_performance_band PRIVATE ${JDF_MIX}) 


endif( EXAGEOSTAT_USE_DPLASMA )
