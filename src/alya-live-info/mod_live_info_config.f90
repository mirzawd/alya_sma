!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_live_info_config

    use def_kintyp, only: ip, rp, lg
    use def_master, only: mmodu, ISLAVE

    implicit none

    integer(ip), parameter   :: MESSAGES_SECTION = -1
    integer(ip), parameter   :: MESSAGES_WARNING = -2

    type LiveInfo_Configuration
        logical(lg)          :: screen = .true.
        integer(ip)          :: lun_livei = 6_ip
        logical(lg)          :: color = .false.
        integer(ip)          :: max_section = 1000000_ip
        character(2)         :: colors(-5:mmodu)
    contains
        procedure :: configure => live_info_config_configure
    end type

    type(LiveInfo_Configuration) :: live_info_config

    public :: LiveInfo_Configuration, live_info_config, MESSAGES_SECTION, MESSAGES_WARNING

contains

    subroutine live_info_config_configure(this)

        implicit none

        class(LiveInfo_Configuration), intent(inout) :: this

        this%color = this%screen .and. this%color
        if (this%screen) then
            this%lun_livei = 6_ip
        else
            this%lun_livei = 16_ip
        end if
        if (ISLAVE) then
            this%lun_livei = 0_ip
        end if

    end subroutine live_info_config_configure

end module mod_live_info_config
