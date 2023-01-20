!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_extensions

#define EXT_MAX_SIZE 10

type Extensions
  character(EXT_MAX_SIZE) :: yaml = ".yaml"
end type Extensions

type(Extensions), parameter :: ext = Extensions()

end module def_extensions
