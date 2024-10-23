###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2017 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenPkgConfig.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.1
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 10-11-2014
#
###

###
#
# CONVERT_LIBSTYLE_TO_PKGCONFIG: convert a libraries list to follow the pkg-config style
#                                used in CLEAN_LIB_LIST
#
###
# macro(CONVERT_LIBSTYLE_TO_PKGCONFIG _liblist)
#     set(${_liblist}_CPY "${${_liblist}}")
#     set(${_liblist} "")
#     foreach(_dep ${${_liblist}_CPY})
#         if (${_dep} MATCHES "^/")
#             get_filename_component(dep_libname ${_dep} NAME)
#             get_filename_component(dep_libdir  ${_dep} DIRECTORY)
#             string(REPLACE "lib"    "" dep_libname "${dep_libname}")
#             string(REPLACE ".so"    "" dep_libname "${dep_libname}")
#             string(REPLACE ".a"     "" dep_libname "${dep_libname}")
#             string(REPLACE ".dylib" "" dep_libname "${dep_libname}")
#             string(REPLACE ".dll"   "" dep_libname "${dep_libname}")
#             list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
#         elseif(NOT ${_dep} MATCHES "^-")
#             list(APPEND ${_liblist} "-l${_dep}")
#         else()
#             list(APPEND ${_liblist} ${_dep})
#         endif()
#     endforeach()
# endmacro(CONVERT_LIBSTYLE_TO_PKGCONFIG)

###
#
# CLEAN_LIB_LIST: clean libraries lists to follow the pkg-config style
#                 used in GENERATE_PKGCONFIG_FILE
#
###
#macro(CLEAN_LIB_LIST _package)
#    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS)
#    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS_PRIVATE)
#    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
#    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
#    convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS)
#    convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS_PRIVATE)
#    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
#    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
#    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
#    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
#endmacro(CLEAN_LIB_LIST)

###
#
# GENERATE_PKGCONFIG_FILE: generate files hqr.pc
#
###
macro(GENERATE_HQR_PKGCONFIG_FILE)

  set(HQR_PKGCONFIG_LIBS "-lhqr")
  set(HQR_PKGCONFIG_LIBS_PRIVATE "-lm")
  set(HQR_PKGCONFIG_REQUIRED "")
  set(HQR_PKGCONFIG_REQUIRED_PRIVATE "")

  #clean_lib_list(HQR)

  set(_output_hqr_file "${CMAKE_BINARY_DIR}/hqr.pc")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/hqr.pc.in"
    "${_output_hqr_file}"
    @ONLY
    )
  install(
    FILES ${_output_hqr_file}
    DESTINATION lib/pkgconfig
    )

endmacro(GENERATE_HQR_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##
