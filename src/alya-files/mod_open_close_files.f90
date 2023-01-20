!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_open_close_files

    use def_master

    implicit none

    public :: open_data_file, close_data_file, open_output_files, open_memory_files

contains

    subroutine open_data_file
        use mod_iofile, only: iofile
        use mod_iofile, only: iofile_file_exists

        external :: par_finali

        if (ISLAVE) then
            return
        end if

        if (adjustl(trim(namda)) == 'connard') then
            call runend('CONNARD TOI-MEME')
        end if
        fil_pdata = adjustl(trim(namda))//'.dat'
        if (.not. iofile_file_exists(fil_pdata)) then
            write (6, '(a)') ''
            write (6, '(a)') 'OPENFI: WRONG PROBLEM NAME, FILE '//trim(fil_pdata)//' DOES NOT EXIST'
            write (6, '(a)') ''
            call par_finali(0_ip)
        else
            call iofile(0_ip, lun_pdata, fil_pdata, 'DATA', 'old')
        end if

    end subroutine open_data_file

    subroutine open_output_files

        use mod_iofile, only: iofile
        use mod_iofile, only: iofile_restart_run
        use mod_iofile, only: iofile_normal_run
        use mod_iofile, only: iofile_file_exists
        use mod_live_info_config, only: live_info_config
        use mod_write_files, only: write_header
        use mod_write_files, only: write_header_restart
        use mod_messages, only: messages_live
        use def_domain, only: fil_outpu_dom

        implicit none

        character(7)               :: statu
        character(11)              :: forma
        character(6)               :: posit
        character(150), save       :: fil_outpu, fil_memor, fil_livei, fil_commu
        character(150), save       :: fil_latex, fil_gnupl, fil_syste
        character(150)             :: fil_memory

        if (ISLAVE) then
            return
        end if

        fil_outpu = adjustl(trim(namda))//'.log'              ! Unit 12
        fil_memor = adjustl(trim(namda))//'.mem'              ! Unit 13
        fil_conve = adjustl(trim(namda))//'.cvg'              ! Unit 14
        fil_livei = adjustl(trim(namda))//'.liv'              ! Unit 16
        fil_postp = adjustl(trim(namda))                      ! Unit 15
        fil_rstar = adjustl(trim(namda))//'.rst'              ! Unit 17
        fil_latex = adjustl(trim(namda))//'-latex.tex'        ! Unit 18
        fil_gnupl = adjustl(trim(namda))//'-latex.plt'        ! Unit 19
        fil_commu = adjustl(trim(namda))//'.com'              ! Unit 20
        fil_binar = adjustl(trim(namda))//'.dom.bin'          ! Unit 24
        fil_syste = adjustl(trim(namda))//'-system.log'       ! Unit 28
        fil_rstib = adjustl(trim(namda))//'-IB.rst'           ! Unit 45
        fil_memory = adjustl(trim(namda))//'-memory.res'       ! Unit 34
        !
        ! Check if restart file exists
        !
        if (kfl_rstar /= 0 .and. .not. iofile_file_exists(fil_rstar)) then
            kfl_rstar = 0
            call messages_live('THIS IS NOT A RESTART', 'WARNING')
            call messages_live('THIS IS NOT A RESTART', 'WARNING')
            call messages_live('THIS IS NOT A RESTART', 'WARNING')
            call messages_live('THIS IS NOT A RESTART', 'WARNING')
            call messages_live('THIS IS NOT A RESTART', 'WARNING')
            call messages_live('KERNEL RESTART FILE DOES NOT EXIST... CONTINUING JUST IN CASE', 'WARNING')
        end if
        !
        ! Open permanent files
        !
        if (kfl_rstar == 2) then
            call iofile_restart_run(statu, forma, posit)
        else
            call iofile_normal_run(statu, forma, posit)
        end if
        call iofile(0_ip, lun_outpu, fil_outpu, 'RUN EVOLUTION', statu, forma, posit)
        call iofile(0_ip, lun_conve, fil_conve, 'CONVERGENCE HISTORY', statu, forma, posit)
        call iofile(0_ip, lun_syste, fil_syste, 'SYSTEM INFO', statu, forma, posit)
        call iofile(0_ip, lun_memory, fil_memory, 'MEMORY EVOLUTION', statu, forma, posit)
        !
        ! Open Postprocess file
        !
        if (kfl_postp_par == 1) then
            !
            ! Compose mesh and result file names
            !
            fil_outpu_dom = trim(fil_postp)//'.post.msh'
            fil_postp = trim(fil_postp)//'.post.res'

        end if
        !
        ! Open live file
        !
        if (live_info_config%lun_livei == 16) then
            if (kfl_rstar == 2) then
                call iofile(0_ip, live_info_config%lun_livei, fil_livei, 'LIVE INFORMATION', 'old', 'formatted', 'append')
            else
                call iofile(0_ip, live_info_config%lun_livei, fil_livei, 'LIVE INFORMATION', 'unknown', 'formatted')
            end if
        end if
        !
        ! Write header
        !
        if (kfl_rstar == 2) then
            call write_header_restart(lun_outpu)
            call write_header_restart(lun_syste)
        else
            call write_header(lun_outpu)
            call write_header(lun_syste)
        end if

    end subroutine open_output_files

    subroutine close_data_file

        use mod_iofile, only: iofile

        implicit none

        if (ISLAVE) then
            return
        end if

        call iofile(2_ip, lun_pdata, ' ', 'DATA')

    end subroutine close_data_file

    subroutine open_memory_files

        use mod_iofile, only: iofile
        use mod_iofile, only: iofile_append_tag
        use mod_memory_config, only: memory_config
        use mod_memory, only: lun_memor
        use mod_memory, only: lun_varcount

        implicit none

        character(150), save        :: fil_memor
        character(150)              :: fil_varcount
        !
        ! Open memory file
        !
        if (ISLAVE) then
            memory_config%output = .false.
            memory_config%varcount = .false.
        end if

        if (memory_config%output) then
            if (kfl_naked == 0) then
                call get_environment_variable('FOR013', fil_memor)
            else if (kfl_naked == 1) then
                fil_memor = adjustl(trim(namda))//'.mem'              ! Unit 13
            end if
            if (INOTSLAVE) then
                if (ISLAVE) call iofile_append_tag(fil_memor, kfl_paral)
                if (kfl_rstar == 2) then
                    call iofile(0_ip, lun_memor, fil_memor, 'MEMORY EVOLUTION', 'old', 'formatted', 'append')
                else
                    call iofile(0_ip, lun_memor, fil_memor, 'MEMORY EVOLUTION')
                end if
            end if
        end if
        if (memory_config%varcount) then
            if (kfl_naked == 0) then
                call get_environment_variable('FOR050', fil_varcount)
            else if (kfl_naked == 1) then
                fil_varcount = adjustl(trim(namda))//'.varcount'          ! Unit 50
            end if
            if (INOTSLAVE) then
                if (ISLAVE) call iofile_append_tag(fil_varcount, kfl_paral)
                if (kfl_rstar == 2) then
                    call iofile(0_ip, lun_varcount, fil_varcount, 'VARIABLE MEMORY COUNTER', 'old', 'formatted', 'append')
                else
                    call iofile(0_ip, lun_varcount, fil_varcount, 'VARIABLE MEMORY COUNTER')
                end if
            end if
        end if
    end subroutine open_memory_files

end module mod_open_close_files
