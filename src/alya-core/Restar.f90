!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    Restart_new.f90
!> @author  houzeaux
!> @date    2019-11-26
!> @brief   Restart
!> @details Restart main subroutine
!> @}
!-----------------------------------------------------------------------

subroutine Restar(itask)

    use def_master
    use mod_communications, only : PAR_BROADCAST
    use mod_messages,       only : messages_live
    use mod_moduls,         only : moduls
    use mod_restart,        only : restart_add
    use mod_restart,        only : restart_ini
    use mod_restart,        only : restart_end
    use mod_restart,        only : restart_initialize
    use mod_restart,        only : restart_finalize
    use mod_strings,        only : integer_to_string
    use mod_strings,        only : real_to_string
    use mod_performance,    only : performance_summary
    use mod_run_config,     only : run_config
    implicit none
    integer(ip), intent(in) :: itask

    read_restart  = .false.
    write_restart = .false.

    !----------------------------------------------------------------------
    !
    ! Read continue restart file
    !
    !----------------------------------------------------------------------

    if (itask == ITASK_READ_RESTART) then

        if (kfl_rstar >= 1) then
            read_restart = .true.
            kfl_reawr = 1
        end if

        if (read_restart) then

            call messages_live('READ RESTART FILES', 'START SECTION')
            !
            ! Reading restart files
            !
            do modul = 0, mmodu
                call moddef(7_ip)
            end do
            modul = 0

            call restart_initialize(itask)

            if (kfl_rstar == 2) then

                call restart_ini(itask)
                call restart_add(ittim, 'ittim')
                call restart_add(cutim, 'cutim')
                call restart_add(dtinv, 'dtinv')
                call restart_add(dtime, 'dtime')
                call restart_add(dpthe, 'dpthe')
                call restart_add(prthe, 'prthe')
                call restart_add(dtold, 'dtold')
                call restart_add(dtinv_old, 'dtinv_old')

                call restart_end(itask)

                if (cutim >= abs(timef)) kfl_gotim = 0
                if (ittim >= mitim) kfl_gotim = 0

            end if

            call Kermod(itask)
            do iblok = 1, nblok
                call moduls(itask)
            end do

            call restart_finalize(itask)
            !
            ! Close files
            !
            do modul = 0, mmodu
                call moddef(6_ip)
            end do
            modul = 0

            call messages_live('RESTARTING FROM TIME STEP '//integer_to_string(ittim)//' AND TIME '//real_to_string(cutim, '(es13.6)'))
            call messages_live('READ RESTART FILES', 'END SECTION')

        end if

    end if

    !----------------------------------------------------------------------
    !
    ! Write restart file
    !
    !----------------------------------------------------------------------

    if (itask == ITASK_WRITE_RESTART) then

        ittim_rst = ittim
        if (run_config%restart%preliminary .and. ( &
            mod(ittim, run_config%restart%preliminary_frequency) == 0 .or. &
            cutim     >= abs(timef) - 1.0e-10_rp                      .or. &
            ittim     >= mitim                                        .or. &
            kfl_timei == 0                                            .or. &
            kfl_gotim == 0                                            .or. &
            kfl_stop  /= 0 )                                          ) then
            kfl_reawr = 2
            write_restart = .true.
        end if

        if (write_restart) then

            call messages_live('WRITE RESTART FILES', 'START SECTION')
            do modul = 0, mmodu
                call moddef(8_ip)
            end do
            modul = 0

            call restart_initialize(itask)
            call restart_ini(itask)
            call restart_add(ittim, 'ittim')
            call restart_add(cutim, 'cutim')
            call restart_add(dtinv, 'dtinv')
            call restart_add(dtime, 'dtime')
            call restart_add(dpthe, 'dpthe')
            call restart_add(prthe, 'prthe')
            call restart_add(dtold, 'dtold')
            call restart_add(dtinv_old, 'dtinv_old')

            call restart_end(itask)

            call Kermod(itask)
            do iblok = 1, nblok
                call moduls(itask)
            end do

            call restart_finalize(itask)
            !
            ! Close files
            !
            do modul = 0, mmodu
                call moddef(6_ip)
            end do
            modul = 0

            call messages_live('WRITE RESTART FILES', 'END SECTION')
            !
            ! Output computing time summary
            !
            call performance_summary()

        end if

    end if
    !
    ! Close restart file
    !
    kfl_reawr = 0
    call moddef(6_ip)

end subroutine Restar
