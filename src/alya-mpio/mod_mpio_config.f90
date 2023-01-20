!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_mpio_config

    use def_kintyp, only: ip, rp, lg

    implicit none

    type MPIO_IO
        logical(lg)            :: enabled = .true.
        logical(lg)            :: force_sequential = .false.
        logical(lg)            :: parallel = .true. ! This variable is not read by the parser
    end type

    type, extends(MPIO_IO)     :: MPIO_Input
        logical(lg)            :: parallel_partition = .false.
    end type

    type, extends(MPIO_IO)     :: MPIO_PostProcess
        logical(lg)            :: light = .false.
        logical(lg)            :: export_only = .false.
    end type

    type MPIO_Output
        logical(lg)            :: enabled = .true.
        logical(lg)            :: merge = .false.
        logical(lg)            :: asynchronous = .true.
        integer(ip)            :: flush = -1_ip
        integer(ip)            :: par_min_block_size = -1_ip
        type(MPIO_PostProcess) :: post_process
        type(MPIO_IO)          :: restart
    end type

    type MPIO_Configuration
        logical(lg)            :: enabled = .false.
        type(MPIO_Input)       :: input
        type(MPIO_Output)      :: output
    contains
        procedure :: configure => mpio_config_configure
        procedure :: broadcast => mpio_config_broadcast
#ifdef MPIO_DEBUG
        procedure :: print => mpio_config_print_debug
#else
        procedure :: print => mpio_config_print
#endif
#ifdef ALYA_YAML
        procedure :: read => mpio_config_read_yaml
        procedure :: read_dat => mpio_config_read_dat
#else
        procedure :: read => mpio_config_read_dat
#endif

    end type

    type(MPIO_Configuration) :: mpio_config

    public :: MPIO_Configuration, mpio_config

contains

    function get_mpio_yaml_file() result(yaml_file)
        use def_extensions, only: ext
        use def_master, only: namda

        character(256) :: yaml_file

        yaml_file = trim(adjustl(trim(adjustl(namda))//trim(adjustl(ext%yaml))))

    end function get_mpio_yaml_file

    function exist_mpio_yaml_file() result(ok)
        logical(lg) :: ok

        inquire (file=trim(adjustl(get_mpio_yaml_file())), exist=ok)

    end function exist_mpio_yaml_file

    subroutine mpio_config_configure(this)

        use def_master, only: IPARALL, kfl_ptask
        use mod_parall, only: PAR_CODE_SIZE, PAR_PARALLEL_PARTITION
        use def_parall, only: kfl_parseq_par

        implicit none

        class(MPIO_Configuration), intent(inout) :: this

        ! Output
        this%output%enabled = this%enabled .and. this%output%enabled
        this%output%asynchronous = this%output%enabled .and. this%output%asynchronous

        ! Output/Post Process
        this%output%post_process%enabled = this%output%enabled .and. this%output%post_process%enabled
        this%output%post_process%export_only = this%output%enabled .and. this%output%post_process%export_only
        this%output%post_process%force_sequential = this%output%post_process%enabled .and. &
                    & this%output%post_process%force_sequential
        this%output%post_process%light = this%output%post_process%enabled .and. &
                    & this%output%post_process%light

        if (this%output%post_process%export_only .and. (kfl_ptask == 2_ip)) then
            call runend("MPIO MESH EXPORTING IS NOT COMPATIBLE WITH MESH RESTART")
        end if

        ! Output/Merge: enable if export only enabled
        this%output%merge = this%output%enabled .and. (this%output%merge .or. this%output%post_process%export_only)

        ! Output/Restart
        this%output%restart%enabled = this%output%enabled .and. this%output%restart%enabled
        this%output%restart%force_sequential = this%output%restart%enabled .and. &
                    & this%output%restart%force_sequential

        ! Input: disable input if export only enabled
        this%input%enabled = this%enabled .and. this%input%enabled .and. (.not. this%output%post_process%export_only) &
                    & .and. (kfl_ptask /= 2_ip)
        this%input%force_sequential = this%input%enabled .and. this%input%force_sequential
        this%input%parallel_partition = this%input%enabled .and. (kfl_parseq_par == PAR_PARALLEL_PARTITION)

        ! IO Parallel management
        this%input%parallel = this%input%enabled .and. (.not. this%input%force_sequential) .and. IPARALL .and. (PAR_CODE_SIZE > 2_ip)
        this%output%post_process%parallel = this%output%post_process%enabled .and. (.not. this%output%post_process%force_sequential) .and. &
                    & IPARALL .and. (PAR_CODE_SIZE > 2_ip)
        this%output%restart%parallel = this%output%restart%enabled .and. (.not. this%output%restart%force_sequential) .and. &
                    & IPARALL .and. (PAR_CODE_SIZE > 2_ip)

    end subroutine mpio_config_configure

    subroutine mpio_config_broadcast(this)

        use mod_exchange, only: exchange_init
        use mod_exchange, only: exchange_add
        use mod_exchange, only: exchange_end

        implicit none

        class(MPIO_Configuration), intent(inout) :: this

        call exchange_init()
        call exchange_add(this%enabled)
        call exchange_add(this%input%enabled)
        call exchange_add(this%input%force_sequential)
        call exchange_add(this%input%parallel)
        call exchange_add(this%input%parallel_partition)
        call exchange_add(this%output%enabled)
        call exchange_add(this%output%merge)
        call exchange_add(this%output%asynchronous)
        call exchange_add(this%output%flush)
        call exchange_add(this%output%par_min_block_size)
        call exchange_add(this%output%post_process%enabled)
        call exchange_add(this%output%post_process%force_sequential)
        call exchange_add(this%output%post_process%parallel)
        call exchange_add(this%output%post_process%light)
        call exchange_add(this%output%post_process%export_only)
        call exchange_add(this%output%restart%enabled)
        call exchange_add(this%output%restart%force_sequential)
        call exchange_add(this%output%restart%parallel)
        call exchange_end()

    end subroutine mpio_config_broadcast

    subroutine mpio_config_print(this)

        use def_master, only: intost
        use mod_messages, only: messages_live

        implicit none

        class(MPIO_Configuration), intent(inout) :: this

        call messages_live('READ MPIO', 'START SECTION')

        if (this%enabled) then
            call messages_live('MPIO feature: enabled')

            call messages_live('INPUT', 'START SECTION')
            if (this%input%enabled) then
                call messages_live('Reading input: enabled')
                if (this%input%force_sequential) then
                    call messages_live('Force sequential: enabled')
                end if
                if (this%input%parallel_partition) then
                    call messages_live('Parallel partition: enabled')
                end if
            else
                call messages_live('Reading input: disabled')
            end if
            call messages_live('INPUT', 'END SECTION')

            call messages_live('OUTPUT', 'START SECTION')
            if (this%output%enabled) then
                call messages_live('Writing output: enabled')
                if (this%output%merge) then
                    call messages_live('Merge subdomains: enabled')
                end if
                if (this%output%asynchronous) then
                    call messages_live('Asynchronous writing: enabled')
                    if (this%output%flush == -1_ip) then
                        call messages_live('Flush: each iteration')
                    else if (this%output%flush == 0_ip) then
                        call messages_live('Flush: each file')
                    else
                        call messages_live('Flush: each '//trim(adjustl(intost(this%output%flush)))//' files')
                    end if
                else
                    call messages_live('Asynchronous writing: disabled')
                end if
                if (this%output%par_min_block_size == -1_ip) then
                    call messages_live('Mixed sequential/parallel mode: forced to parallel')
                else if (this%output%par_min_block_size == 0_ip) then
                    call messages_live('Miwed sequential/parallel mode: forced to sequential')
                else
                    call messages_live('Mixed sequential/parallel mode: minimum block size:'&
                    &//trim(adjustl(intost(this%output%par_min_block_size)))//' MB')
                end if
                call messages_live('POST-PROCESS', 'START SECTION')
                if (this%output%post_process%enabled) then
                    call messages_live('Post-Process writing: enabled')
                    if (this%output%post_process%force_sequential) then
                        call messages_live('Force sequential: enabled')
                    end if
                    if (this%output%post_process%light) then
                        call messages_live('Light mode: enabled')
                    end if
                    if (this%output%post_process%export_only) then
                        call messages_live('Export only: enabled')
                    end if
                else
                    call messages_live('Post-Process writing: disabled')
                end if
                call messages_live('POST-PROCESS', 'END SECTION')

                call messages_live('RESTART', 'START SECTION')
                if (this%output%restart%enabled) then
                    call messages_live('Restart writing: enabled')
                    if (this%output%restart%force_sequential) then
                        call messages_live('Force sequential: enabled')
                    end if
                else
                    call messages_live('Restart writing: disabled')
                end if
                call messages_live('RESTART', 'END SECTION')
            else
                call messages_live('Writing output: disabled')
            end if
            call messages_live('OUTPUT', 'END SECTION')
        else
            call messages_live('MPIO feature: disabled')
        end if

        call messages_live('READ MPIO', 'END SECTION')

    end subroutine mpio_config_print

#ifdef MPIO_DEBUG

    subroutine mpio_config_print_debug(this)

        use def_master, only: intost
        use mod_messages, only: messages_live

        implicit none

        class(MPIO_Configuration), intent(inout) :: this

        ! ------- GLOBAL -------
        call messages_live('READ MPIO', 'START SECTION')

        ! ------- INPUT -------
        if (this%enabled) then
            call messages_live('MPIO feature: enabled')
        else
            call messages_live('MPIO feature: disabled')
        end if

        call messages_live('INPUT', 'START SECTION')

        if (this%input%enabled) then
            call messages_live('Reading input: enabled')
        else
            call messages_live('Reading input: disabled')
        end if

        if (this%input%force_sequential) then
            call messages_live('Force sequential: enabled')
        end if

        if (this%input%parallel_partition) then
            call messages_live('Parallel partition: enabled')
        end if

        call messages_live('INPUT', 'END SECTION')

        ! ------- OUTPUT -------

        call messages_live('OUTPUT', 'START SECTION')

        if (this%output%enabled) then
            call messages_live('Writing output: enabled')
        else
            call messages_live('Writing output: disabled')
        end if

        if (this%output%merge) then
            call messages_live('Merge subdomains: enabled')
        end if

        if (this%output%asynchronous) then
            call messages_live('Asynchronous writing: enabled')
        end if

        if (this%output%flush == -1_ip) then
            call messages_live('Flush: each iteration')
        else if (this%output%flush == 0_ip) then
            call messages_live('Flush: each file')
        else
            call messages_live('Flush: each '//trim(adjustl(intost(this%output%flush)))//' files')
        end if

        if (this%output%par_min_block_size == -1_ip) then
            call messages_live('Mixed sequential/parallel mode: forced to parallel')
        else if (this%output%par_min_block_size == 0_ip) then
            call messages_live('Miwed sequential/parallel mode: forced to sequential')
        else
            call messages_live('Mixed sequential/parallel mode: minimum block size:'&
            &//trim(adjustl(intost(this%output%par_min_block_size)))//' MB')
        end if

        call messages_live('POST-PROCESS', 'START SECTION')

        if (this%output%post_process%enabled) then
            call messages_live('Post-Process writing: enabled')
        else
            call messages_live('Post-Process writing: disabled')
        end if

        if (this%output%post_process%force_sequential) then
            call messages_live('Force sequential: enabled')
        end if
        if (this%output%post_process%light) then
            call messages_live('Light mode: enabled')
        end if
        if (this%output%post_process%export_only) then
            call messages_live('Export only: enabled')
        end if

        call messages_live('POST-PROCESS', 'END SECTION')

        call messages_live('RESTART', 'START SECTION')
        if (this%output%restart%enabled) then
            call messages_live('Restart writing: enabled')
        else
            call messages_live('Restart writing: disabled')
        end if
        if (this%output%restart%force_sequential) then
            call messages_live('Force sequential: enabled')
        end if
        call messages_live('RESTART', 'END SECTION')
        call messages_live('OUTPUT', 'END SECTION')

        call messages_live('READ MPIO', 'END SECTION')

    end subroutine mpio_config_print_debug

#endif

#ifdef ALYA_YAML

    subroutine mpio_config_read_yaml(this)
        use mod_yaml, only: AlyaYAMLFile, AlyaYAMLMap
        use def_master, only: ISLAVE, namda
        use mod_messages, only: messages_live

        implicit none

        class(MPIO_Configuration), intent(inout) :: this
        type(AlyaYAMLFile)                       :: yaml_file
        type(AlyaYAMLMap)                        :: map_1, map_2, map_3

        if (ISLAVE) then
            return
        end if

        if (this%output%post_process%export_only) then
            return
        end if

        if (.not. exist_mpio_yaml_file()) then
            call messages_live('Reading MPIO options from '//adjustl(trim(namda))//'.dat is deprecated!', 'WARNING')
            call messages_live('You should read them from '//adjustl(trim(get_mpio_yaml_file()))//'!', 'WARNING')
            call this%read_dat()
            return
        end if

        call yaml_file%open(trim(adjustl(get_mpio_yaml_file())))

        !.md<module>kernel
        !.md<input>case.yaml
        !.md<pos>4
        !.md<sec>
        !.md<0># MPI-IO section
        !.md<code>yaml
        !.mdmpio:

        if (yaml_file%has_map(label="mpio")) then

            map_1 = yaml_file%get_map(label="mpio")
            !------------------------------enable-------------------------------------------!
            !.md<1>enable: false # enable the MPIO format section
            call map_1%set_logical(label="enable", option=this%enabled)

            !------------------------------input--------------------------------------------!
            !.md<1>input: # input section
            if (map_1%has_label(label="input")) then

                map_2 = map_1%get_map(label="input")

                !.md<2>enable: true # read the input files using the MPIO format
                call map_2%set_logical(label="enable", option=this%input%enabled)

                !.md<2>force-sequential: false # read the input files in sequential (true) or in parallel (false)
                call map_2%set_logical(label="force-sequential", option=this%input%force_sequential)

                call map_2%destroy()

            end if

            !------------------------------output-------------------------------------------!
            !.md<1>output: # output section
            if (map_1%has_label(label="output")) then

                map_2 = map_1%get_map(label="output")

                !.md<2>enable: true # write the output files using the MPIO format
                call map_2%set_logical(label="enable", option=this%output%enabled)

                !.md<2>merge: false # merge the subdomains in the output files
                call map_2%set_logical(label="merge", option=this%output%merge)

                !.md<2>asynchronous: true # write the output files asynchronously
                call map_2%set_logical(label="asynchronous", option=this%output%asynchronous)

                !.md<2>flush: -1 # [integer] flush the asynchronous writing queue every n item (-1: flush at each iteration)
                call map_2%set_integer(label="flush", option=this%output%flush)

                !.md<2>parallel-minimum-block-size: -1 # [integer] minimum file block size per MPI process required to write in parallel (-1: write always in parallel)
                call map_2%set_integer(label="parallel-minimum-block-size", option=this%output%par_min_block_size)

                !------------------------------output: post-process-----------------------------!
                !.md<2>post-process: # post-process section
                if (map_2%has_label(label="post-process")) then

                    map_3 = map_2%get_map(label="post-process")

                    !.md<3>enable: true # write the post-process files using the MPIO format
                    call map_3%set_logical(label="enable", option=this%output%post_process%enabled)

                    !.md<3>force-sequential: false # write the post-process files in sequential (true) or in parallel (false)
                    call map_3%set_logical(label="force-sequential", option=this%output%post_process%force_sequential)

                    !.md<3>export-only: false # read the ascii input and convert it to the MPIO format without running any iteration
                    call map_3%set_logical(label="export-only", option=this%output%post_process%export_only)

                    !.md<3>light: false # do not export all the mesh files; incompatible with restart
                    call map_3%set_logical(label="light", option=this%output%post_process%light)

                    call map_3%destroy()

                end if

                !.md<2>restart: # restart section
                if (map_2%has_label("restart")) then

                    map_3 = map_2%get_map(label="restart")

                    !.md<3>enable: true # write the restart files using the MPIO format
                    call map_3%set_logical(label="enable", option=this%output%restart%enabled)

                    !.md<3>force-sequential: false # write the restart files in sequential (true) or in parallel (false)
                    call map_3%set_logical(label="force-sequential", option=this%output%restart%force_sequential)

                    call map_3%destroy()

                end if

                call map_2%destroy()

            end if

            call map_1%destroy()

        end if
        !.md</code>

        call yaml_file%close()

    end subroutine mpio_config_read_yaml

#endif

    subroutine mpio_config_read_dat(this)

        use mod_ecoute, only: ecoute, ecoute_set_read_unit
        use mod_messages, only: messages_live
        use def_inpout, only: words, option, option_not_off, exists, param
        use def_master, only: ISLAVE, lun_pdata

        implicit none

        class(MPIO_Configuration), intent(inout) :: this

        !.md<module>kernel
        !.md<input>case.dat
        !.md<pos>4
        !.md<sec>
        !.md<0># MPI-IO section
        !.md<code>
        !.mdMPI_IO: ON | OFF

        if (ISLAVE) then
            return
        end if

        if (this%output%post_process%export_only) then
            return
        end if

        call ecoute_set_read_unit(lun_pdata) ! Reading file
        rewind (lun_pdata)

        do while (words(1) /= 'MPIIO')
            call ecoute('READ_MPI_IO', STOP_END_OF_FILE=.false.)
            if (words(1) == 'EOF  ') then
                return
            end if
        end do

        if (option_not_off('MPIIO')) then

            this%enabled = .true.

            call ecoute('READ_MPI_IO')

            do while (words(1) /= 'ENDMP')

                !.md<1>GEOMETRY: ON | SEQUENTIAL | EXPORT | OFF
                !.md<field>GEOMETRY
                !.md<com>Geometry reading method.
                !.md<com>    -  `ON`:         Read the `mpio` format, in parallel if possible
                !.md<com>    -  `SEQUENTIAL`: Read the `mpio` format, in sequential
                !.md<com>    -  `EXPORT`:     Read the `ASCII` format, and export the mesh to the `mpio` format without running any step
                !.md<com>    -  `OFF`:        Read the `ASCII` format

                if (words(1) == 'GEOME') then
                    if (option_not_off('GEOME')) then
                        this%input%enabled = .true.
                        if (words(2) == 'SEQUE') then
                            this%input%force_sequential = .true.
                        else if (words(2) == 'EXPOR') then
                            this%output%enabled = .true.
                            this%output%post_process%enabled = .true.
                            this%output%post_process%export_only = .true.
                        end if
                    else
                        this%input%enabled = .false.
                    end if

                    !.md<1>MERGE: ON | OFF
                    !.md<field>MERGE
                    !.md<com>Merge the post-process subdomains

                else if (words(1) == 'MERGE') then
                    this%output%enabled = .true.

                    if (option_not_off('MERGE')) then
                        this%output%merge = .true.
                    else
                        this%output%merge = .false.
                    end if

                    !.md<1>SYNCHRONOUS: ON | OFF
                    !.md<field>SYNCHRONOUS
                    !.md<com>Force synchronous file writing, OFF by default. Use this option if the asynchronous writing fails.

                else if (words(1) == 'SYNCH') then
                    this%output%enabled = .true.

                    if (option('SYNCH')) then
                        this%output%asynchronous = .false.
                    else
                        this%output%asynchronous = .true.
                    end if

                    !.md<1>MINIMUM: (int)
                    !.md<field>MINIMUM
                    !.md<com>Minimal block size en MB to perform MPI-IO calls (sequential I/O if < this value)

                else if (words(1) == 'MINIM') then
                    this%output%enabled = .true.
                    this%output%flush = int(param(1), ip)

                    !.md<1>BUFFER | ASYNC: (int)
                    !.md<com>-  **ASYNC**: Flush asynchronous writes every [int] files.

                else if (words(1) == 'BUFFE' .or. words(1) == 'ASYNC') then
                    this%output%enabled = .true.
                    this%output%par_min_block_size = int(param(1), ip)

                    !.md<1>POSTPROCESS: ON | SEQUENTIAL | OFF, LIGHT
                    !.md<field>POSTPROCESS
                    !.md<com>Post-process writing method
                    !.md<com>    -  `ON`:         Write to the `mpio` format, in parallel if possible
                    !.md<com>    -  `SEQUENTIAL`: Write to the `mpio` format, in sequential
                    !.md<com>    -  `OFF`:        Write to the `alyabin` format
                    !.md<com>    -  `LIGHT`:      (optional, if `ON` or `SEQUENTIAL`) Only export the mesh main files necessary for post-process. Incompatible with mesh restart!

                else if (words(1) == 'POSTP') then
                    this%output%enabled = .true.

                    if (option_not_off('POSTP')) then
                        this%output%post_process%enabled = .true.
                        if (words(2) == 'SEQUE') then
                            this%output%post_process%force_sequential = .true.
                        end if
                        if (exists('LIGHT')) then
                            this%output%post_process%light = .true.
                        end if
                    else
                        this%output%post_process%enabled = .false.

                    end if

                    !.md<1>RESTART: ON | SEQUENTIAL | OFF
                    !.md<field>RESTART
                    !.md<com>Restart writing method
                    !.md<com>    -  `ON`:         Write to the `mpio` format, in parallel if possible
                    !.md<com>    -  `SEQUENTIAL`: Write to the `mpio` format, in sequential
                    !.md<com>    -  `OFF`:        Write to the `alyabin` format

                else if (words(1) == 'RESTA') then
                    this%output%enabled = .true.

                    if (option_not_off('RESTA')) then
                        this%output%restart%enabled = .true.
                        if (words(2) == 'SEQUE') then
                            this%output%restart%force_sequential = .true.
                        end if
                    else
                        this%output%restart%enabled = .false.

                    end if

                end if

                call ecoute('READ_MPI_IO')

            end do

            !.mdEND_MPI_IO
            !.md</code>

        end if

    end subroutine mpio_config_read_dat

end module mod_mpio_config
