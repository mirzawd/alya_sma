###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



get_filename_component(PROJID ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${PROJID} LANGUAGES Fortran C)

set(GITINFO_FILE "gitinfo.f90")

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${GITINFO_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/_${GITINFO_FILE}"
  COMMAND ${CMAKE_COMMAND}
    -Dlocal_dir="${CMAKE_SOURCE_DIR}/cmake"
    -Doutput_dir="${CMAKE_CURRENT_BINARY_DIR}"
    -Droot_dir="${CMAKE_SOURCE_DIR}"
    -P "${CMAKE_SOURCE_DIR}/cmake/gitinfo.cmake"
)

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} "${CMAKE_CURRENT_BINARY_DIR}/${GITINFO_FILE}")

include(compiler)
set_flags()
set_defs()
set_parall()
include(dependencies)

include(lib-set-install)
