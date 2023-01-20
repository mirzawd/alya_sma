!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_seq_log.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO log manager (sequential)
!> @details This module records information on the I/O operations (sequential)
!>          Required MPIOLOG to be defined
!>          \verbatim
!>          Output format:
!>          FIELD, IO OPERATION, PAR/SEQ, FORMAT, IO WORKING PROCESSES, ALL IO PROCESSES, FILE SIZE (B), START TIMESTAMP (s), END TIMESTAMP (s), DURATION (s)
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

 module mod_mpio_seq_log

    use def_kintyp,                     only : ip,rp,lg,r1p
    use def_master,                     only : intost, retost, namda
    use def_master,                     only : IMASTER, ISEQUEN

    implicit none

    private

    character(150)                      :: fil_timer

    integer, parameter                  :: u=19589453

    real(rp)                            :: time1, time2

    character(20)                       :: s_object, s_operation, s_format, s_par, s_workers, s_code_size
    character(25)                       :: s_file_size

    public                              :: start_timer, end_timer, write_data, openfile_timer, closefile_timer, u, fil_timer, time1, time2, s_object, s_operation, s_format, s_par, s_workers, s_code_size, s_file_size

    contains

    subroutine start_timer(object, operation, format, workers, code_size, file_size)
        character(*), intent(in)                        :: object, operation, format
        integer(ip), intent(in), optional               :: workers, code_size
        integer(8), intent(in), optional               :: file_size
#ifdef MPIOLOG
        s_object=object
        s_operation=operation
        s_format=format
        s_par='SEQ'
        s_workers='1'
        s_code_size='1'
        s_file_size='0'
        if (present(code_size)) then
            s_code_size=intost(code_size)
        end if
        if (present(file_size)) then
            write(s_file_size,*) file_size
            s_file_size = adjustl(s_file_size)
        end if
        if (IMASTER .or. ISEQUEN) then
            call cputim(time1)
        end if
#endif
    end subroutine

    subroutine end_timer()
#ifdef MPIOLOG
        if (IMASTER .or. ISEQUEN) then
            call cputim(time2)
            call openfile_timer()
            call write_data()
            call closefile_timer()
        end if
#endif
    end subroutine

    subroutine write_data()
        write(u,*) trim(s_object), ", ", trim(s_operation), ", ", trim(s_par), ", ", trim(s_format), ", ", trim(s_workers), ", ",&
             trim(s_code_size), ", ", trim(s_file_size), ", ", trim(retost(time1)), ", ", trim(retost(time2)), ", ",&
             trim(retost(time2-time1))
    end subroutine


    subroutine openfile_timer()
        integer                         :: stat
        fil_timer = adjustl(trim(namda))//".mpio.log"
        open(unit=u, iostat=stat, file=fil_timer, status='old')
        if (stat == 0) then
            close (u)
            open(u, file=fil_timer, status="old", position="append", action="write", recl=1000)
        else
            close (u)
            open(u, file=fil_timer, status="new", action="write", recl=1000)
        end if
    end subroutine

    subroutine closefile_timer()
        close(u)
    end subroutine

end module
