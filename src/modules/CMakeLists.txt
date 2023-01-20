###################################################################
# Copyright 2005 - 2022 Barcelona Supercomputing Center.          #
# Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE #
# for nonprofit scientific purposes only.                         #
# See companion file LICENSE.txt.                                 #
###################################################################



foreach(MOD ${MODULES})
  string(TOUPPER ${MOD} UPMOD)
  if(WITH_MODULE_${UPMOD})
    add_subdirectory(${MOD})
  endif()
endforeach()
