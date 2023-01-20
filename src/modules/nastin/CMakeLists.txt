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

if(IS_DIRECTORY "${CMAKE_HOME_DIRECTORY}/src/private/nastin" AND NOT WITH_NASTIN_PRIVATE MATCHES "OFF")
  file(GLOB_RECURSE SRC_FILES_PRIVATE ${CMAKE_HOME_DIRECTORY}/src/private/nastin/*.f90)
  set(WITH_NASTIN_PRIVATE "ON" CACHE BOOL "Enable Nastin++" FORCE)
else()
  set(WITH_NASTIN_PRIVATE "OFF" CACHE BOOL "Enable Nastin++" FORCE)
endif()

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES} ${SRC_FILES_PRIVATE})

include(compiler)
set_flags()
set_defs()
set_parall()
include(dependencies)

if (NOT WITH_NASTIN_PRIVATE)
  add_definitions(-DNASTIN_PRIVATE_OFF)
endif()

include(lib-management)
alya_link(kernel) 

include(lib-set-install)
