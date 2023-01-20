###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



get_filename_component(ProjectId ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})

project(${ProjectId} LANGUAGES Fortran)
file(GLOB_RECURSE SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/modules/${PROJECT_NAME}/*.f90)
add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES})

include(compiler)
set_flags()
set_flags_fortran_error()
set_defs()
set_parall()
include(dependencies)

include(lib-management)
alya_link(kernel)

include(lib-set-install)
