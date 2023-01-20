###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



#Public libraries
add_subdirectory(alya-info)
add_subdirectory(alya-define)
add_subdirectory(alya-std)
add_subdirectory(alya-memory)
add_subdirectory(alya-time)
add_subdirectory(alya-utils)
add_subdirectory(alya-strings)
add_subdirectory(alya-arguments)
add_subdirectory(alya-io)
add_subdirectory(alya-log)
add_subdirectory(alya-ecoute)
add_subdirectory(alya-fem)
add_subdirectory(alya-live-info)
add_subdirectory(alya-parall)
if (WITH_YAML)
  add_subdirectory(alya-yaml)
endif()
add_subdirectory(alya-vecmat)
add_subdirectory(alya-maths)
add_subdirectory(alya-physics)
add_subdirectory(alya-matrix)
add_subdirectory(alya-search)
add_subdirectory(alya-htable)
add_subdirectory(alya-wrapper-metis)
add_subdirectory(alya-comm)
add_subdirectory(alya-run)
add_subdirectory(alya-mpio)
add_subdirectory(alya-options)
if (WITH_FTI)
  add_subdirectory(alya-fti)
endif()
add_subdirectory(alya-graph)
add_subdirectory(alya-soldir)
add_subdirectory(alya-solver)
add_subdirectory(alya-domain)
add_subdirectory(alya-files)
add_subdirectory(alya-interpolation)
add_subdirectory(alya-coupling)
add_subdirectory(alya-random)
add_subdirectory(alya-tab)
add_subdirectory(alya-ann)
add_subdirectory(alya-kermod)
add_subdirectory(alya-performance)
add_subdirectory(alya-proper)
add_subdirectory(alya-wrapper-petsc)

#Private libraries
if(IS_DIRECTORY "${CMAKE_HOME_DIRECTORY}/src/private/alya-proper" AND WITH_PROPER_PRIVATE)
  add_subdirectory(private/alya-proper)
  set(WITH_PROPER_PRIVATE "ON" CACHE BOOL "Disable proper-private" FORCE)
else()
  set(WITH_PROPER_PRIVATE "OFF" CACHE BOOL "Disable proper-private" FORCE)
endif()
if(IS_DIRECTORY "${CMAKE_HOME_DIRECTORY}/src/private/alya-proper-elem" AND WITH_PROPER_ELEM_PRIVATE)
  add_subdirectory(private/alya-proper-elem)
  set(WITH_PROPER_ELEM_PRIVATE "ON" CACHE BOOL "Disable proper-elem-private" FORCE)
else()
  set(WITH_PROPER_ELEM_PRIVATE "OFF" CACHE BOOL "Disable proper-elem-private" FORCE)
endif()

#Kernel (to be refactorized)
add_subdirectory(alya-kernel)

#Modules
add_subdirectory(modules)
add_subdirectory(alya-module-interface)

#Executables
if (NOT WITH_LIBRARIES_ONLY)
  add_subdirectory(alya)
  add_subdirectory(alya2pos)
endif()
