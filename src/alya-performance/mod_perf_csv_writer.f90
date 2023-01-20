!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

! Created by Damien Dosimont on 23/06/2020.

module mod_perf_csv_writer

    use def_kintyp, only: ip, rp
    use mod_strings, only: lower_case
    use def_master

    implicit none

    private
    character(150) :: fil_perf
    character(150) :: pname
    character(150) :: phost
    character(150) :: pdate
    integer(ip) :: pcpus

    public :: init_perf, write_perf_header, write_perf_line

contains

    subroutine init_perf()
        use mod_iofile, only: iofile_open_unit
        fil_perf = adjustl(trim(namda))//'-performance.csv' ! Unit 52
        call iofile_open_unit(lun_perf, fil_perf, 'PERFORMANCE')
        pname = title
        call date_and_time(DATE=pdate)
        call get_environment_variable("HOST", phost)
        if (phost == "") then
            phost = "unknown"
        end if
        pcpus = npart + 1

    end subroutine init_perf

    subroutine write_perf_header()
        write (lun_perf, '(a)') 'name, hostname, date, cores, module, function, type, parent_module, parent_function, parent_type, category, unit, value'
    end subroutine write_perf_header

    subroutine write_perf_line(pmodule, pfunction, ptype, punit, pvalue, pcat, ppmodule, ppfunction, pptype)
        character(len=*), intent(in) :: pmodule
        character(len=*), intent(in) :: pfunction
        character(len=*), intent(in) :: ptype
        character(len=*), intent(in) :: punit
        real(rp), intent(in) :: pvalue
        character(len=*), optional, intent(in) :: pcat
        character(len=*), optional, intent(in) :: ppmodule
        character(len=*), optional, intent(in) :: ppfunction
        character(len=*), optional, intent(in) :: pptype
        character(len=150)                     :: local_pcat
        character(len=150)                     :: local_ppmodule
        character(len=150)                     :: local_ppfunction
        character(len=150)                     :: local_pptype

        if (present(ppmodule)) then
            local_ppmodule = ppmodule
        else
            local_ppmodule = pmodule
        end if
        if (present(ppfunction)) then
            local_ppfunction = ppfunction
        else
            local_ppfunction = "all"
        end if
        if (present(pptype)) then
            local_pptype = pptype
        else
            local_pptype = ptype
        end if
        if (present(pcat)) then
            local_pcat = pcat
            if (trim(local_pcat) == "all") then
                local_pcat = "computation+communications+io"
            else if (trim(local_pcat) == "cc") then
                local_pcat = "computation+communications"
            end if
        else
            local_pcat = "computation+communications"
        end if

        write (lun_perf, 1) trim(pname), trim(phost), trim(pdate), trim(intost(pcpus)), &
            trim(lower_case(pmodule)), trim(pfunction), trim(ptype), &
            trim(lower_case(local_ppmodule)), trim(local_ppfunction), trim(local_pptype), &
            trim(local_pcat), trim(punit), pvalue

1       format(a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', a, ', ', e15.8e3)
    end subroutine write_perf_line

end module mod_perf_csv_writer
