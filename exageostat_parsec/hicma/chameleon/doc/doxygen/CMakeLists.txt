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

cmake_minimum_required(VERSION 2.8)

# Looking for doxygen
# -------------------
find_package(Doxygen)

# Generate documentation
# ----------------------
if(DOXYGEN_EXECUTABLE)

    # Create file chameleon.dox
    # ---------------------------
    if(CHAMELEON_SCHED_STARPU)

        configure_file("conf-dev-starpu.dox.in"
                       "${CMAKE_CURRENT_BINARY_DIR}/chameleon-dev.dox"
                       @ONLY)
#        configure_file("conf-usr-starpu.dox.in"
#                       "chameleon-usr.dox"
#                       @ONLY)

    elseif(CHAMELEON_SCHED_QUARK)

        configure_file("conf-dev-quark.dox.in"
                       "${CMAKE_CURRENT_BINARY_DIR}/chameleon-dev.dox"
                        @ONLY)
#        configure_file("conf-usr-quark.dox.in"
#                       "chameleon-usr.dox"
#                       @ONLY)

    endif()

    # Message during configure
    #-------------------------
    message(STATUS "Looking for doxygen - found")

    # User's documentation
    # --------------------
    add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/out-dev
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       ARGS    ${CMAKE_CURRENT_BINARY_DIR}/chameleon-dev.dox
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/chameleon-dev.dox
                      )
#    add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/out-usr
#                       COMMAND ${DOXYGEN_EXECUTABLE}
#                       ARGS    ${CMAKE_CURRENT_BINARY_DIR}/chameleon-usr.dox
#                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/chameleon-usr.dox
#                      )
    add_custom_target(doxygen-out ALL
                      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/out-dev
                              #${CMAKE_CURRENT_BINARY_DIR}/out-usr
                              chameleon_all_sources
                              )

    # Installation for html version
    # -----------------------------
    install(DIRECTORY   ${CMAKE_CURRENT_BINARY_DIR}/out-dev/html
            DESTINATION share/chameleon/doxygen/dev)
#    install(DIRECTORY   ${CMAKE_CURRENT_BINARY_DIR}/out-usr/html
#            DESTINATION share/chameleon/doxygen/usr)

else()
    message(STATUS "Looking for doxygen - not found")

endif()

###
### END CMakeLists.txt
###
