!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_fti_config

#ifdef ALYA_FTI

    use def_kintyp, only: ip, rp, lg

    implicit none

    type FTI_Configuration
        logical(lg)        :: enabled = .false.
    contains
        procedure :: broadcast => fti_config_broadcast
        procedure :: print => fti_config_print
        procedure :: set_flags => fti_set_flags
        procedure :: check_mpio => fti_check_mpio
#ifdef ALYA_YAML
        procedure :: read => fti_config_read_yaml
        procedure :: read_dat => fti_config_read_dat
#else
        procedure :: read => fti_config_read_dat
#endif

    end type

    type(FTI_Configuration) :: fti_config

    public :: FTI_Configuration, fti_config

contains

    function get_fti_yaml_file() result(yaml_file)
        use def_extensions, only: ext
        use def_master, only: namda

        character(256) :: yaml_file

        yaml_file = trim(adjustl(trim(adjustl(namda))//trim(adjustl(ext%yaml))))

    end function get_fti_yaml_file

    function exist_fti_yaml_file() result(ok)
        logical(lg) :: ok

        inquire (file=trim(adjustl(get_fti_yaml_file())), exist=ok)

    end function exist_fti_yaml_file

    subroutine fti_set_flags(this)

        use mod_alya2fti, only: FTI_usage
        use mod_alya2fti, only: FTI_DBG

        implicit none

        class(FTI_Configuration), intent(inout) :: this
        if (this%enabled) then
            FTI_usage = .true.
            FTI_DBG = 0
        end if

    end subroutine fti_set_flags

    subroutine fti_check_mpio(this)

        use mod_mpio_config, only: mpio_config

        implicit none

        class(FTI_Configuration), intent(inout) :: this

        if (mpio_config%output%restart%enabled .and. this%enabled) then
            call runend('MOD_FTI_CONFIG: MPIO RESTART AND FTI CANNOT BE ENABLED TOGETHER!')
        end if

    end subroutine fti_check_mpio

    subroutine fti_config_broadcast(this)

        use mod_exchange, only: exchange_init
        use mod_exchange, only: exchange_add
        use mod_exchange, only: exchange_end

        implicit none

        class(FTI_Configuration), intent(inout) :: this

        call exchange_init()
        call exchange_add(this%enabled)
        call exchange_end()

    end subroutine fti_config_broadcast

    subroutine fti_config_print(this)

        use def_master, only: intost
        use mod_messages, only: messages_live

        implicit none

        class(FTI_Configuration), intent(inout) :: this

        call messages_live('FTI', 'START SECTION')

        if (this%enabled) then
            call messages_live('FTI: enabled')
        end if

        call messages_live('FTI', 'END SECTION')

    end subroutine fti_config_print

#ifdef ALYA_YAML

    subroutine fti_config_read_yaml(this)
        use mod_yaml, only: AlyaYAMLFile, AlyaYAMLMap
        use def_master, only: ISLAVE, namda
        use mod_mpio_config, only: mpio_config
        use mod_messages, only: messages_live

        implicit none

        class(FTI_Configuration), intent(inout) :: this
        type(AlyaYAMLFile)                      :: yaml_file
        type(AlyaYAMLMap)                       :: map_1, map_2, map_3

        if (ISLAVE) then
            return
        end if

        if (mpio_config%output%post_process%export_only) then
            return
        end if

        if (.not. exist_fti_yaml_file()) then
            call messages_live('Reading FTI options from '//adjustl(trim(namda))//'.dat is deprecated!', 'WARNING')
            call messages_live('You should read them from '//adjustl(trim(get_fti_yaml_file()))//'!', 'WARNING')
            call this%read_dat()
            return
        end if

        call yaml_file%open(trim(adjustl(get_fti_yaml_file())))

        !.md<module>kernel
        !.md<input>case.yaml
        !.md<pos>5
        !.md<sec>
        !.md<0># FTI option
        !.md<code>yaml
        !.mdfti:

        if (yaml_file%has_map(label="fti")) then

            map_1 = yaml_file%get_map(label="fti")
            !.md<1>enabled: false # read or write restart using FTI
            call map_1%set_logical(label="enabled", option=this%enabled)
            call map_1%destroy()

        end if
        !.md</code>

        call yaml_file%close()

    end subroutine fti_config_read_yaml

#endif

    subroutine fti_config_read_dat(this)

        use mod_ecoute, only: ecoute, ecoute_set_read_unit
        use mod_messages, only: messages_live
        use def_inpout, only: words, option
        use def_master, only: ISLAVE, lun_pdata
        use mod_mpio_config, only: mpio_config

        implicit none

        class(FTI_Configuration), intent(inout) :: this

        !.md<module>kernel
        !.md<input>case.dat
        !.md<pos>5
        !.md<sec>
        !.md<0># FTI option
        !.md<code>
        !.mdFTI: ON | OFF

        if (ISLAVE) then
            return
        end if

        if (mpio_config%output%post_process%export_only) then
            return
        end if

        call ecoute_set_read_unit(lun_pdata) ! Reading file
        rewind (lun_pdata)

        do while (words(1) /= 'FTI  ')
            call ecoute('READ_FTI', STOP_END_OF_FILE=.false.)
            if (words(1) == 'EOF  ') then
                return
            end if
        end do

        if (option('FTI  ')) then
            this%enabled = .true.
        end if

        !.md</code>

    end subroutine fti_config_read_dat

#endif

end module mod_fti_config
