!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

subroutine respre(itask, kfl_gores)
    !------------------------------------------------------------------------
    !****f* kernel/respre
    ! NAME
    !    respre
    ! DESCRIPTION
    !    Define if a preliminary or restart run should be carried out
    ! USES
    ! USED BY
    !    nsi_turnon
    !***
    !------------------------------------------------------------------------
    use def_parame
    use def_master
    use def_domain
    use mod_messages, only: messages_live
    use mod_memory
    use mod_run_config, only: run_config
#ifdef ALYA_FTI
    use mod_alya2fti, only: FTI_st
#endif
    implicit none
    integer(ip), intent(in)  :: itask
    integer(ip), intent(out) :: kfl_gores
    integer(ip), save        :: kfl_ptask_old, jtask

    if (itask /= 3) then
        kfl_gores = 0
        kfl_ptask_old = kfl_ptask
    end if
    jtask = abs(itask)

    select case (jtask)

    case (READ_RESTART_FILE)
        !
        ! Read from restart file
        !
        if (kfl_rstar >= 1) then
            kfl_gores = 1
            kfl_reawr = 1
            if (itask > 0) then
                kfl_ptask = 1
                if (modul == 0) then
                    call messages_live('KERNEL: READ RESTART FILE')
                else
                    call messages_live(trim(namod(modul))//': READ RESTART FILE')
                end if
#ifdef ALYA_FTI
                if (0 == FTI_st) then
#endif
                    call moddef(7_ip)
#ifdef ALYA_FTI
                end if
#endif
            end if
        end if

    case (WRITE_RESTART_FILE)
        !
        ! Write to restart file
        !

        if ( &
#ifndef ALYA_FTI
            run_config%restart%preliminary .and. ( &
#endif
            mod(ittim, run_config%restart%preliminary_frequency) == 0 .or. &
            cutim >= timef - 1.0e-10_rp .or. &
            ittim >= mitim .or. &
            kfl_timei == 0) &
#ifndef ALYA_FTI
            ) &
#endif
            then
            kfl_gores = 1
            kfl_reawr = 2
            if (itask > 0) then
                kfl_ptask = 1
                if (modul == 0) then
                    call messages_live('KERNEL: WRITE RESTART FILE')
                else
                    call messages_live(trim(namod(modul))//': WRITE RESTART FILE')
                end if
                call moddef(8_ip)
            end if
        end if

    case (3_ip)
        !
        ! Recover old values
        !
        lun_postp = lun_postp_old
        if (kfl_gores == 1) kfl_ptask = kfl_ptask_old
        kfl_reawr = 0
#ifdef ALYA_FTI
        if (FTI_st /= 1_ip) then
#endif
            call moddef(6_ip)
#ifdef ALYA_FTI
        end if
#endif
        nullify (gesca)
        nullify (gevec)

    end select

    if (kfl_gores == 1 .and. itask <= 2) then

        if (IMASTER) then
            nullify (gevec)
            nullify (gesca)
            nullify (ger3p)
        end if

        lun_postp_old = lun_postp
        if (modul > 0) lun_postp = momod(modul)%lun_rstpo

    end if

end subroutine respre
