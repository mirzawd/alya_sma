###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



get_filename_component(PROJID ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${PROJID} LANGUAGES Fortran C)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c ${CMAKE_HOME_DIRECTORY}/src/alya-core/*.f90)
add_executable(${PROJECT_NAME} ${SRC_FILES})

include(compiler)
set_flags()
set_flags_fortran_error()
set_defs()
set_parall()
include(dependencies)

include(lib-management)
alya_link(module-interface)

target_link_libraries(${PROJECT_NAME} c)

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
