!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_write_files

    use def_master
    implicit none

    public :: write_header, write_header_restart

contains

    subroutine write_header(nunit)

        implicit none

        integer(ip), intent(in) :: nunit
        integer(4) :: nunit4

        nunit4 = int(nunit, 4)
        write (nunit4, 1) adjustl(trim(title))

1       format(///, 5x, '   --------', /, &
                5x, '|- A L Y A: ', a, /, &
                5x, '   --------',/)
    end subroutine write_header

    subroutine write_header_restart(nunit)

        implicit none

        integer(ip), intent(in) :: nunit
        integer(4) :: nunit4

        nunit4 = int(nunit, 4)
        write (nunit4, 2)

2       format(/, 5x, '|- THIS IS A TIME RESTART RUN: CONTINUING...', //)
    end subroutine write_header_restart

end module mod_write_files
