###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



get_filename_component(PROJID ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${PROJID} LANGUAGES Fortran C)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c)

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES})

include(compiler)
set_flags()
set_defs()
set_parall()
include(dependencies)

include(lib-management)
alya_link(interpolation)

include(lib-set-install)
