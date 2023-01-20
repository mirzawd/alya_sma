###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



get_filename_component(PROJID ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${PROJID} LANGUAGES Fortran C)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c)
file(GLOB_RECURSE HEADER_FILES ${CMAKE_CURRENT_LIST_DIR}/*.inc)
include_directories(${CMAKE_CURRENT_LIST_DIR})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES} ${HEADER_FILES})

include(compiler)
set_flags()
set_defs()
set_parall()
include(dependencies)

include(lib-management)
alya_link(live-info utils ecoute log memory)

include(lib-set-install)
