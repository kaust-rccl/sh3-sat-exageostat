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

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if ( NOT BUILD_SUBPROJECT )
  set( HQR_INSTALL_DIR bin/hqr )
else()
  set( HQR_INSTALL_DIR bin )
endif()

set(TESTINGS
  draw_hqr.c
  draw_systolic.c
  testing_hqr.c
  testing_systolic.c
  )

foreach (_file ${TESTINGS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we}
    ${_file} common.c)
  target_link_libraries(${_name_we} hqr)

  install(
    TARGETS ${_name_we}
    RUNTIME DESTINATION ${HQR_INSTALL_DIR} )
  install(
    FILES   ${_file}
    DESTINATION ${HQR_INSTALL_DIR} )

endforeach()
