!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_run_config

    use def_kintyp, only: ip, rp, lg
    use def_master, only: mmodu
    use mod_memory_config, only: memory_config
    use mod_io_config, only: io_config
    use mod_live_info_config

    implicit none

    integer(ip), parameter   :: rst_none = 0_ip
    integer(ip), parameter   :: rst_initial = 1_ip
    integer(ip), parameter   :: rst_continue = 2_ip
    integer(ip), parameter   :: rst_interpolate = 3_ip

    integer(ip), parameter   :: cust_none = 0_ip
    integer(ip), parameter   :: cust_marek = 1_ip
    integer(ip), parameter   :: cust_cfwind1 = 2_ip
    integer(ip), parameter   :: cust_cfwind2 = 3_ip
    integer(ip), parameter   :: cust_cfwind0 = -2_ip

    type Run_Restart
        logical(lg)          :: init_from_cmd = .false.
        integer(ip)          :: preliminary_frequency = huge(1_ip)
        logical(lg)          :: preliminary = .false.
        integer(ip)          :: run_type = rst_none
        logical(lg)          :: append_timestep = .false.
        character(256)       :: char_run_type = "final"
    end type

    type Run_Configuration
        character(66)        :: title = ' '
        integer(ip)          :: code = 1_ip
        type(Run_Restart)    :: restart
        character(5)         :: char_customer = "none"
        integer(ip)          :: customer = cust_none
        logical(lg)          :: log = .true.
        logical(lg)          :: timeline = .false.
        logical(lg)          :: timing = .false.
        logical(lg)          :: yaml = .true.

    contains
        procedure :: initialize => run_config_initialize
        procedure :: configure => run_config_configure
        procedure :: broadcast => run_config_broadcast
        procedure :: print => run_config_print
#ifdef ALYA_YAML
        procedure :: read => run_config_read_yaml
        procedure :: read_dat => run_config_read_dat
#else
        procedure :: read => run_config_read_dat
#endif

    end type Run_Configuration

    type(Run_Configuration) :: run_config

    public :: Run_Configuration, run_config

contains

    function get_run_yaml_file() result(yaml_file)
        use def_extensions, only: ext
        use def_master, only: namda

        character(256) :: yaml_file

        yaml_file = trim(adjustl(trim(adjustl(namda))//trim(adjustl(ext%yaml))))

    end function get_run_yaml_file

    function exist_run_yaml_file() result(ok)
        logical(lg) :: ok

        inquire (file=trim(adjustl(get_run_yaml_file())), exist=ok)

    end function exist_run_yaml_file

    elemental subroutine str2int(str, int, stat)
        implicit none
        ! Arguments
        character(len=*), intent(in) :: str
        integer(ip), intent(out)     :: int
        integer(ip), intent(out)     :: stat

        read (str, *, iostat=stat) int
    end subroutine str2int

    function color_label_to_id(label) result(id)

        use def_master, only: module_name_to_id

        implicit none

        character(len=*) :: label
        integer(ip) :: id

        if (label == "section") then
            id = MESSAGES_SECTION
        else if (label == "warning") then
            id = MESSAGES_WARNING
        else
            id = module_name_to_id(label)
        end if

    end function color_label_to_id

    subroutine run_config_initialize(this)

        use def_master
        use mod_ansi_colors, only: ansi_colors_name_to_code

        implicit none

        class(Run_Configuration), intent(inout) :: this

        live_info_config%colors = ''
        live_info_config%colors(MESSAGES_SECTION) = ansi_colors_name_to_code('light blue')
        live_info_config%colors(MESSAGES_WARNING) = ansi_colors_name_to_code('light red')
        this%title = ' '

        !kfl_rstar    = 0          ! Not a restart run

    end subroutine run_config_initialize

    subroutine run_config_configure(this)

        use def_master

        implicit none

        class(Run_Configuration), intent(inout) :: this

        ! Restart

        if (this%restart%char_run_type == "final") then
            this%restart%run_type = rst_none
        else if (this%restart%char_run_type == "initial") then
            this%restart%run_type = rst_initial
        else if (this%restart%char_run_type == "continue") then
            this%restart%run_type = rst_continue
        else if (this%restart%char_run_type == "interpolate") then
            this%restart%run_type = rst_interpolate
        end if

        live_info_config%color = live_info_config%screen .and. live_info_config%color
        if (live_info_config%screen) then
            live_info_config%lun_livei = 6_ip
        else
            live_info_config%lun_livei = 16_ip
        end if

        call memory_config%configure()
        call live_info_config%configure()

        ! Set alya variables

        !title
        title = trim(adjustl(this%title))
        !code
        current_code = this%code
        !run type
        kfl_rstar = this%restart%run_type

    end subroutine run_config_configure

    subroutine run_config_broadcast(this)

        use mod_exchange, only: exchange_init
        use mod_exchange, only: exchange_add
        use mod_exchange, only: exchange_end

        implicit none

        class(Run_Configuration), intent(inout) :: this
        integer(ip)                             :: i

        call exchange_init()
        call exchange_add(this%restart%init_from_cmd)
        call exchange_add(this%restart%preliminary)
        call exchange_add(this%restart%preliminary_frequency)
        call exchange_add(this%restart%run_type)
        call exchange_add(this%restart%append_timestep)
        call exchange_add(this%restart%char_run_type)
        call exchange_add(this%char_customer)
        call exchange_add(this%customer)
        call exchange_add(this%log)
        call exchange_add(this%timeline)
        call exchange_add(this%timing)
        call exchange_add(this%yaml)
        call exchange_add(io_config%state)
        call exchange_add(memory_config%output)
        call exchange_add(memory_config%varcount)
        call exchange_add(memory_config%high_memory)
        call exchange_add(live_info_config%screen)
        call exchange_add(live_info_config%lun_livei)
        call exchange_add(live_info_config%color)
        call exchange_add(live_info_config%max_section)
        do i = -5, mmodu
            call exchange_add(live_info_config%colors(i))
        end do
        call exchange_end()

    end subroutine run_config_broadcast

    subroutine run_config_print(this)

        use def_master, only: intost, retost, namda
        use mod_messages, only: messages_live

        implicit none

        class(Run_Configuration), intent(inout) :: this

        if (.not. this%yaml) then
            call messages_live('Reading run options from '//adjustl(trim(namda))//'.dat is deprecated!', 'WARNING')
            call messages_live('You should read them from '//adjustl(trim(get_run_yaml_file()))//'!', 'WARNING')
        end if
        call messages_live('READ RUN CONFIGURATION', 'START SECTION')
        call messages_live('Title: '//trim(adjustl(this%title)))
        call messages_live('Code number: '//intost(this%code))

        if (io_config%state) call messages_live('State file: enabled')
        if (this%log) call messages_live('Log file: enabled')
        if (this%timeline) call messages_live('Timeline: enabled')
        if (this%timing) call messages_live('Timing: enabled')

        call messages_live('READ RESTART', 'START SECTION')
        call messages_live('Run type: '//trim(adjustl(this%restart%char_run_type)))
        if (this%restart%preliminary) call messages_live('Preliminary: enabled')
        if (this%restart%preliminary) call messages_live('Preliminary frequency: '//intost(this%restart%preliminary_frequency))
        if (this%restart%append_timestep) call messages_live('Append timestep to file name: enabled')
        call messages_live('READ RESTART', 'END SECTION')

        call messages_live('READ MEMORY', 'START SECTION')
        if (memory_config%output) call messages_live('Output: enabled')
        if (memory_config%varcount) call messages_live('Variable memory counter output: enabled')
        if (memory_config%high_memory) call messages_live('High memory: enabled')
        call messages_live('READ MEMORY', 'END SECTION')

        call messages_live('READ RUN CONFIGURATION', 'END SECTION')

        call messages_live('READ LIVE INFO CONFIGURATION', 'START SECTION')
        if (live_info_config%screen) then
            call messages_live('Live info: screen')
        else
            call messages_live('Live info: file')
        end if
        if (live_info_config%color) call messages_live('Color: enabled')
        call messages_live('Max sections: '//intost(live_info_config%max_section))
        call messages_live('READ LIVE INFO CONFIGURATION', 'END SECTION')

    end subroutine run_config_print

#ifdef ALYA_YAML

    subroutine run_config_read_yaml(this)
        use mod_yaml, only: AlyaYAMLFile, AlyaYAMLMap
        use def_master, only: ISLAVE
        use mod_messages, only: messages_live
        use mod_ansi_colors, only: ansi_colors_name_to_code
        use mod_strings, only: integer_to_string

        implicit none

        class(Run_Configuration), intent(inout) :: this
        type(AlyaYAMLFile)                      :: yaml_file
        type(AlyaYAMLMap)                       :: map_1, map_2
        character(256)                          :: label
        character(256)                          :: color_name
        integer(ip)                             :: color_code
        integer(ip)                             :: stat

        if (ISLAVE) then
            return
        end if

        if (.not. exist_run_yaml_file()) then
            this%yaml = .false.
            call this%read_dat()
            return
        end if

        call this%initialize()

        call yaml_file%open(trim(adjustl(get_run_yaml_file())))

        !.md<module>kernel
        !.md<input>case.yaml
        !.md<pos>1
        !.md<sec>
        !.md<0># RUN section
        !.md<code>yaml
        !.mdrun:

        if (yaml_file%has_map(label="run")) then

            map_1 = yaml_file%get_map(label="run")

            !.md<1>title: "char" # simulation title
            call map_1%set_string(label="title", option=this%title)

            !.md<1>code: 1 # code number
            call map_1%set_integer(label="code", option=this%code)

            !.md<1>restart: # restart section
            if (map_1%has_label(label="restart") .and. .not. this%restart%init_from_cmd) then

                map_2 = map_1%get_map(label="restart")

                !.md<2>type: "final" # initial | continue: run type
                call map_2%set_string(label="type", option=this%restart%char_run_type)

                !.md<2>preliminary: false # enable preliminary run
                call map_2%set_logical(label="preliminary", option=this%restart%preliminary)

                if (this%restart%preliminary) then

                    !.md<2>frequency: 1 # preliminary frequency
                    call map_2%set_integer(label="frequency", option=this%restart%preliminary_frequency)

                end if

                !.md<2>append-timestep: true # append timestep to the file name
                call map_2%set_logical(label="append-timestep", option=this%restart%append_timestep)

                call map_2%destroy()

            end if

            !.md<1>state: false # enable state file
            call map_1%set_logical(label="state", option=io_config%state)

            !.md<1>customer: "none" # marek | cfdw1 | cfdw2 | cfdw0
            call map_1%set_string(label="customer", option=this%char_customer)

            !.md<1>log: true # enable log file
            call map_1%set_logical(label="log", option=this%log)

            !.md<1>timeline: false # enable timeline
            call map_1%set_logical(label="timeline", option=this%timeline)

            !.md<1>memory: # memory section
            if (map_1%has_label(label="memory")) then

                map_2 = map_1%get_map(label="memory")

                !.md<2>output: false # output memory
                call map_2%set_logical(label="output", option=memory_config%output)

                !.md<2>variable-counter: false # output variable memory counter
                call map_2%set_logical(label="variable-counter", option=memory_config%varcount)

                !.md<2>high-memory: false # high memory usage
                call map_2%set_logical(label="high-memory", option=memory_config%high_memory)

                call map_2%destroy()

            end if

            call map_1%destroy()

        end if

        !.md</code>
        !.md<0># LIVE INFO section
        !.md<code>yaml
        !.mdlive-info:

        if (yaml_file%has_map(label="live-info")) then

            map_1 = yaml_file%get_map(label="live-info")

            !.md<1>screen: true # live information on screen (true) or in file (false)
            call map_1%set_logical(label="screen", option=live_info_config%screen)

            !.md<1>color: false # enable color
            call map_1%set_logical(label="color", option=live_info_config%color)

            !.md<1>max-section: 1000000 # enable color
            call map_1%set_integer(label="max-section", option=live_info_config%max_section)

            if (live_info_config%screen .and. live_info_config%color .and. map_1%has_label(label="colors")) then

                map_2 = map_1%get_map(label="colors")
                !.md<1>colors: # color dictionnary
                !.md<2>section: "light green"
                !.md<2>warning: "red"
                !.md<2>nastin: 41

                label = trim(adjustl(map_2%get_next_label()))
                do while (label /= "")
                    color_name = ""
                    call map_2%set_string(label=label, option=color_name)
                    call str2int(color_name, color_code, stat)
                    if (stat == 0) then
                        live_info_config%colors(color_label_to_id(label)) = integer_to_string(color_code)
                    else
                        live_info_config%colors(color_label_to_id(label)) = ansi_colors_name_to_code(color_name)
                    end if
                    label = trim(adjustl(map_2%get_next_label()))
                end do

                call map_2%destroy()

            end if

            call map_1%destroy()

        end if
        !.md</code>

        call yaml_file%close()

    end subroutine run_config_read_yaml

#endif

    subroutine run_config_read_dat(this)

        use mod_ecoute, only: ecoute, ecoute_set_read_unit
        use mod_messages, only: messages_live
        use def_inpout, only: words, option, option_not_off, exists, wname
        use def_inpout, only: getcha_long, getint
        use def_master, only: ISLAVE, lun_pdata

        implicit none

        class(Run_Configuration), intent(inout) :: this

        !.md<module>kernel
        !.md<input>case.dat
        !.md<pos>0
        !.md<sec>
        !.md<0># run section
        !.md<code>
        !.mdRUN_DATA

        if (ISLAVE) then
            return
        end if

        call ecoute_set_read_unit(lun_pdata) ! Reading file
        rewind (lun_pdata)

        call this%initialize()

        do while (words(1) /= 'RUNDA')
            call ecoute('READ_RUNDA', STOP_END_OF_FILE=.false.)
            if (words(1) == 'EOF  ') then
                call runend('RRUDAT: WRONG RUN_DATA CARD')
            end if
        end do

        do while (words(1) /= 'ENDRU')

            !.md<1>ALYA: [str]
            !.md<field>ALYA
            !.md<com>Simulation title

            if (words(1) == 'ALYA ') then
                !
                ! Read and write title
                !
                this%title = trim(adjustl(wname))

                !.md<1>CODE: [int]
                !.md<field>CODE
                !.md<com>code number
            else if (words(1) == 'CODE ') then
                !
                ! Read and write title
                !
                this%code = getint('CODE ', 1_ip, '#My code number')

                !.md<1>RUNTY: [ FINAL | INITI | CONTI | INTER ], PRELI, FREQU=[int], APPEN
                !.md<field>RUNTY
                !.md<com>code Run type. Several choices are available:
                !.md<com>
                !.md<com>    -  `FINAL`: final run; does not export restart files
                !.md<com>    -  `INITI`: restart run initial; read restart files as initial conditions
                !.md<com>    -  `CONTI`: restart run continue; read restart files to continue a run
                !.md<com>    -  `INTER`: restart run interpolate and initial; read restart files
                !.md<com>    -  `PRELI`: preliminary run option; write restart files
                !.md<com>    -  `FREQU`: preliminary frequency option
                !.md<com>    -  `APPEN`: option to append time steps to the file names
                !.md<com>
            else if (words(1) == 'RUNTY' .and. .not. this%restart%init_from_cmd) then
                !
                ! Read the type of run
                !
                if (exists('PRELI')) then

                    this%restart%preliminary = .true.
                    this%restart%preliminary_frequency = getint('FREQU', 1_ip, '#Preliminary frequency')   ! Preliminary frequency

                end if

                if (exists('RESTA')) this%restart%char_run_type = "initial"     ! Restart run: initial
                if (exists('INITI')) this%restart%char_run_type = "initial"     ! Restart run: initial
                if (exists('CONTI')) this%restart%char_run_type = "continue"    ! Restart run: continue
                if (exists('INTER')) this%restart%char_run_type = "interpolate" ! Restart run: interpolate and initial
                this%restart%append_timestep = exists('APPEN')                  ! Append time step to file names

                !.md<1>STATE: [bool]
                !.md<field>STATE
                !.md<com>State file
            else if (words(1) == 'STATE') then
                !
                ! State file
                !
                io_config%state = words(2) == 'YES  ' .or. words(2) == 'ON   '

                !.md<1>CUSTO: [ NONE | MAREK | CFDW1 | CFDW2 | CFDW0 ]
                !.md<field>CUSTO
                !.md<com>customer:
                !.md<com>
                !.md<com>    - `MAREK`: Marek
                !.md<com>    - `CFDW1`: Iberdrola: CFDWin1
                !.md<com>    - `CFDW2`: Iberdrola: CFDWin2
                !.md<com>    - `CFDW0`: Iberdrola: CFDWin1 + base field
                !.md<com>
            else if (words(1) == 'CUSTO') then
                !
                ! Read the customer
                !
                if (words(2) == 'MAREK') then
                    this%customer = cust_marek      ! Marek
                else if (words(2) == 'CFDW1') then
                    this%customer = cust_cfwind1    ! Iberdrola: CFDWind1
                else if (words(2) == 'CFDW2') then
                    this%customer = cust_cfwind2    ! Iberdrola: CFDWind2
                else if (words(2) == 'CFDW0') then
                    this%customer = cust_cfwind0    ! Iberdrola: CFDWind1 + base field
                end if

                !.md<1>LOGFI: [bool]
                !.md<field>LOGFI
                !.md<com>Enable log file
            else if (words(1) == 'LOGFI') then
                !
                ! Log file
                !
                this%log = words(2) == 'YES  ' .or. words(2) == 'ON   '

                !.md<1>TIMEL: [bool]
                !.md<field>TIMEL
                !.md<com>Enable timeline
            else if (words(1) == 'TIMEL') then
                !
                ! Output timeline
                !
                this%timeline = option('TIMEL')

                !.md<1>MEMOR: [bool], VARIA
                !.md<field>MEMOR
                !.md<com>Enable memory output. Option `VARIA`: output variable memory counter
            else if (words(1) == 'MEMOR') then
                !
                ! Memory count
                !
                if (words(2) == 'YES  ' .or. words(2) == 'ON   ') then
                    memory_config%output = .true.
                    memory_config%varcount = exists('VARIA')
                end if

                !.md<1>TIMIN: [bool]
                !.md<field>TIMIN
                !.md<com>Enable timing
            else if (words(1) == 'TIMIN') then
                this%timing = words(2) == 'YES  ' .or. words(2) == 'ON   '

                !.md<1>LOTOF: [bool]
                !.md<field>LOTOF
                !.md<com>Enable high memory consumption
            else if (words(1) == 'LOTOF') then
                memory_config%high_memory = words(2) == 'YES  ' .or. words(2) == 'ON   '

                !.md<1>LIVEI: [ SCREEN | FILE ], COLOR, SECTI=[int]
                !.md<field>LIVEI
                !.md<com>Log file is screen
            else if (words(1) == 'LIVEI') then
                !
                ! Live information
                !
                if (words(2) == 'SCREE') then
                    live_info_config%screen = .true.
                else if (words(2) == 'FILE ') then
                    live_info_config%screen = .false.
                end if
                live_info_config%color = exists('COLOR')
                if (exists('SECTI')) live_info_config%max_section = getint('SECTI', 1_ip, '#Number of sections')

                !.md<1>COLOR
                !.md<field>COLOR
                !.md<com>COLOR settings
            else if (words(1) == 'COLOR') then
                !
                ! Output colors
                !
                block

                    use def_master, only: module_name_to_id
                    use mod_ansi_colors, only: ansi_colors_name_to_code
                    use mod_strings, only: integer_to_string

                    character(200) :: name_color
                    integer(ip)    :: code_color, imodu
                    character(5)   :: mname
                    do while (words(1) /= 'ENDCO')
                        !.md<2>WARNING: red
                        if (words(1) == 'WARNI') then
                            if (trim(words(2)) /= '') then
                                name_color = getcha_long('WARNI', 'red  ', '#color', len(name_color, KIND=ip))
                                live_info_config%colors(MESSAGES_WARNING) = ansi_colors_name_to_code(name_color)
                            else
                                code_color = getint('WARNI', 1_ip, '#color')
                                live_info_config%colors(MESSAGES_WARNING) = integer_to_string(code_color)
                            end if
                            !.md<2>SECTION: light blue
                        else if (words(1) == 'SECTI') then
                            if (trim(words(2)) /= '') then
                                name_color = getcha_long('SECTI', 'red  ', '#color', len(name_color, KIND=ip))
                                live_info_config%colors(MESSAGES_SECTION) = ansi_colors_name_to_code(name_color)
                            else
                                code_color = getint('SECTI', 1_ip, '#color')
                                live_info_config%colors(MESSAGES_SECTION) = integer_to_string(code_color)
                            end if
                            !.md<2>NASTIN: 41
                        else
                            mname = words(1)
                            imodu = module_name_to_id(mname)
                            if (trim(words(2)) /= '') then
                                name_color = getcha_long(mname, 'red  ', '#color', len(name_color, KIND=ip))
                                live_info_config%colors(imodu) = ansi_colors_name_to_code(name_color)
                            else
                                code_color = getint(mname, 1_ip, '#color')
                                live_info_config%colors(imodu) = integer_to_string(code_color)
                            end if
                        end if
                        call ecoute('RRUDAT')
                    end do

                end block

            end if
            !.md<1>END_COLOR

            call ecoute('RRUDAT')

        end do

        !.mdEND_RUN_DATA
        !.md</code>

    end subroutine run_config_read_dat

end module mod_run_config
