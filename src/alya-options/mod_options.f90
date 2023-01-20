!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_options

    use def_master

    implicit none

    public :: read_options

contains

    subroutine read_options()
        use mod_arguments, only: type_options
        use mod_arguments, only: arguments_number
        use mod_arguments, only: arguments_help
        use mod_arguments, only: arguments_get
        use mod_arguments, only: arguments_find
        use mod_mpio_config, only: mpio_config
        use mod_run_config, only: run_config
        use mod_strings, only: string_to_integer

        implicit none

        integer(ip)                 :: istat
        integer(ip)                 :: iopt, nopt, ii, kk, marg
        character(100)              :: opt
        character(100), allocatable :: opt_arg(:)
        type(type_options)          :: options(11) = [ &
                                       type_options('k', 'check', 0_4, 'Check data file'), &
                                       type_options('f', 'file', 1_4, 'Specify input file arg1'), &
                                       type_options('n', 'name', 1_4, 'Dirichlet/Neumann with PLE'), &
                                       type_options('h', 'help', 0_4, 'Print help'), &
                                       type_options('r', 'read', 0_4, 'Read partition'), &
                                       type_options('w', 'write', 0_4, 'Write partition'), &
                                       type_options('e', 'export', 0_4, 'Export geometry in MPIO format'), &
                                       type_options('c', 'read-rst', 0_4, 'Read restart to continue the run'), &
                                       type_options('i', 'read-rst-init', 0_4, 'Read restart to initialize the run'), &
                                       type_options('p', 'write-rst', 0_4, 'Write restart'), &
                                       type_options('t', 'time-steps', 1_4, 'Run arg1 time steps')]

        marg = maxval(options(:)%num_arguments)
        if (marg > 0) allocate (opt_arg(marg))

        if (ISLAVE) then
            return
        end if

        nopt = arguments_number()
        iopt = 0
        if (nopt == 0) call arguments_help(options)

        do while (iopt < nopt)

            iopt = iopt + 1
            opt = arguments_get(iopt)
            ii = arguments_find(options, opt)

            if (ii == -1) then

                write (*, '(a)') ' '
                write (*, '(a)') '   WRONG COMMAND LINE OPTION: '//trim(opt)
                write (*, '(a)') ' '
                call arguments_help(options)

            else if (ii /= 0) then

                do kk = 1, options(ii)%num_arguments
                    iopt = iopt + 1
                    opt_arg(kk) = arguments_get(iopt)
                end do

                select case (options(ii)%shortopt)

                case ('f')
                    !
                    ! Problem name
                    !
                    namda = trim(opt_arg(1))

                case ('t')
                    !
                    ! Number of time steps
                    !
                    mitim = string_to_integer(opt_arg(1), istat)
                    if (istat /= 0) call runend('OPENFI: WRONG NUMBER OF TIME STEPS IN COMMAND LINE')

                case ('k')
                    !
                    ! Check data file
                    !
                    kfl_check_data_file = 1

                case ('n')

                    continue

                case ('h')
                    !
                    ! Help
                    !
                    call arguments_help(options)

                case ('p')
                    !
                    ! Write restart files: preliminary run
                    !
                    run_config%restart%preliminary = .true.
                    run_config%restart%preliminary_frequency = huge(1_ip)
                    run_config%restart%init_from_cmd = .true.

                case ('c')
                    !
                    ! Write restart files (continue run)
                    !
                    run_config%restart%char_run_type = "continue"
                    run_config%restart%init_from_cmd = .true.

                case ('i')
                    !
                    ! Write restart files (initial condition)
                    !
                    run_config%restart%char_run_type = "initial"
                    run_config%restart%init_from_cmd = .true.

                case ('w')
                    !
                    ! Write restart
                    !
                    kfl_ptask = 0

                case ('r')
                    !
                    ! Read restart
                    !
                    kfl_ptask = 2

                case ('e')
                    !
                    ! Export
                    !
                    mpio_config%enabled = .true.
                    mpio_config%output%enabled = .true.
                    mpio_config%output%post_process%enabled = .true.
                    mpio_config%output%post_process%export_only = .true.

                case default
                    !
                    ! Problem name
                    !
                    namda = trim(opt)

                end select

            else

                namda = trim(opt)

            end if

        end do
        if (marg > 0) deallocate (opt_arg)

    end subroutine read_options

end module mod_options
