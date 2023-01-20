!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_io_config

    use def_kintyp, only: lg

    implicit none

    type IO_Configuration
        logical(lg)          :: state = .false.
    end type

    type(IO_Configuration) :: io_config

    public :: IO_Configuration, io_config

end module mod_io_config
