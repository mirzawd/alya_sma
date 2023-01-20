!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_log.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO log manager (sequential)
!> @details This module records information on the I/O operations (parallel)
!>          It extends mod_mpio_seq_log
!>          Required MPIOLOG to be defined
!>          \verbatim
!>          Output format:
!>          FIELD, IO OPERATION, PAR/SEQ, FORMAT, IO WORKING PROCESSES, ALL IO PROCESSES, FILE SIZE (B), START TIMESTAMP (s), END TIMESTAMP (s), DURATION (s)
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

 module mod_mpio_par_log

    use def_kintyp,                     only : ip,rp,lg,r1p
    use def_master,                     only : intost, retost, namda
    use def_master,                     only : IMASTER, ISEQUEN
    use mod_communications,             only : PAR_BARRIER
    use mod_mpio_seq_log

    implicit none

    private

    public                              :: par_start_timer, par_end_timer

    contains

    subroutine par_start_timer(wherein, object, operation, format, workers, code_size, file_size)
        character(*), intent(in)                        :: wherein, object, operation, format
        integer(ip), intent(in), optional               :: workers, code_size
        integer(8), intent(in), optional               :: file_size
#ifdef MPIOLOG
        s_object=object
        s_operation=operation
        s_format=format
        s_par='PAR'
        s_workers='1'
        s_code_size='1'
        s_file_size='0'
        if (present(workers)) then
            s_workers=intost(workers)
        end if
        if (present(code_size)) then
            s_code_size=intost(code_size)
        end if
        if (present(file_size)) then
            write(s_file_size,*) file_size
            s_file_size = adjustl(s_file_size)
        end if
        call PAR_BARRIER(wherein)
        if (IMASTER .or. ISEQUEN) then
            call cputim(time1)
        end if
#endif
    end subroutine

    subroutine par_end_timer(wherein)
        character(*), intent(in)                        :: wherein
#ifdef MPIOLOG
        call PAR_BARRIER(wherein)
        if (IMASTER .or. ISEQUEN) then
            call cputim(time2)
            call openfile_timer()
            call write_data()
            call closefile_timer()
        end if
#endif
    end subroutine

end module
