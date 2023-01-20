!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_generic_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO generic
!> @details This module is a bridge that redirects to the sequential or the parallel
!>          versions of the MPI I/O operations
!> @}
!-----------------------------------------------------------------------

module mod_mpio_generic_io

    use def_master
    use mod_mpio_seq_generic_io
    use mod_mpio_par_generic_io

    implicit none

    private

    interface MPIO_READ
        module procedure                        MPIO_READ_INT_V,&
                                                MPIO_READ_INT_M,&
                                                MPIO_READ_REAL_V,&
                                                MPIO_READ_REAL_M
    end interface

    public :: MPIO_READ

    contains

    subroutine MPIO_READ_INT_V(buf, filename, dim)
        integer(ip), pointer, intent(inout)            :: buf(:)
        character(*),         intent(in)               :: filename
        integer(ip),          intent(inout)            :: dim
        if (ISEQUEN) then
            call mpio_seq_read(buf, filename, dim)
        else
            call mpio_par_read(buf, filename, dim)
        end if
    end subroutine

    subroutine MPIO_READ_INT_M(buf, filename, lines, columns)
        integer(ip), pointer, intent(inout)            :: buf(:,:)
        character(*),         intent(in)               :: filename
        integer(ip),          intent(inout)            :: lines, columns
        if (ISEQUEN) then
            call mpio_seq_read(buf, filename, lines, columns)
        else
            call mpio_par_read(buf, filename, lines, columns)
        end if
    end subroutine

    subroutine MPIO_READ_REAL_V(buf, filename, dim)
        real(rp), pointer, intent(inout)               :: buf(:)
        character(*),         intent(in)               :: filename
        integer(ip),          intent(inout)            :: dim
        if (ISEQUEN) then
            call mpio_seq_read(buf, filename, dim)
        else
            call mpio_par_read(buf, filename, dim)
        end if
    end subroutine

    subroutine MPIO_READ_REAL_M(buf, filename, lines, columns)
        real(rp), pointer, intent(inout)               :: buf(:,:)
        character(*),         intent(in)               :: filename
        integer(ip),          intent(inout)            :: lines, columns
        if (ISEQUEN) then
            call mpio_seq_read(buf, filename, lines, columns)
        else
            call mpio_par_read(buf, filename, lines, columns)
        end if
    end subroutine

end module
