!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_memory_config

    use def_kintyp, only: ip, rp, lg

    implicit none

    type Memory_Configuration
        logical(lg)          :: output = .false.
        logical(lg)          :: varcount = .false.
        logical(lg)          :: high_memory = .false.
    contains
        procedure :: configure => memory_config_configure
    end type

    type(Memory_Configuration) :: memory_config

    public :: Memory_Configuration, memory_config

contains

    subroutine memory_config_configure(this)

        implicit none

        class(Memory_Configuration), intent(inout) :: this

        this%varcount = this%output .and. this%varcount

    end subroutine memory_config_configure

end module mod_memory_config
