!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

subroutine reapro()
    !-----------------------------------------------------------------------
    !****f* master/reapro
    ! NAME
    !    Reapro
    ! DESCRIPTION
    !    This routine starts reading data.
    ! USES
    !    inirun   to perform some initializations.
    !    openfi   to get file names and open them.
    !    rrudat   to read run data.
    !    rproda   to read general problem data.
    !    cputim
    !    Nastin
    !    Temper
    !    Codire
    !    Alefor
    ! USED BY
    !    Alya
    !***
    !-----------------------------------------------------------------------
    use def_parame
    use def_master
    use mod_options, only: read_options
    use mod_ker_timeline, only: ker_timeline_open_file
    use mod_ker_timeline, only: ker_timeline_synchronization
    use mod_messages, only: livinf
    use mod_messages, only: messages_live
    use mod_messages, only: messages_header
    use mod_output_postprocess, only: output_postprocess_allocate
    use mod_run_config, only: run_config
    use mod_mpio_config, only: mpio_config
    use mod_open_close_files, only: open_data_file
    use mod_open_close_files, only: open_output_files
    use mod_open_close_files, only: open_memory_files
#ifdef ALYA_FTI
    use mod_fti_config, only: fti_config
#endif
    implicit none

    external :: par_initialize_openacc

    !
    ! Check if process was initiated by MPI
    !
    call par_initialize_mpi()
    !
    ! Message header
    !
    call messages_header()
    !
    ! Get data file name and open it
    !
    call read_options()
    call open_data_file()
    !
    ! Read run data
    !
    !call rrudat()
    call run_config%read()
    call run_config%configure()
    call run_config%broadcast()
    !
    ! Get result file names and open them
    !
    call open_output_files()
    call ker_timeline_open_file()

    !
    ! Live information
    !
    call livinf(1_ip, ' ', zero)
    call messages_live('READ PROBLEM DATA')
    call run_config%print()
    !
    ! Checkpoint for Parall communication
    !
    call par_checkpoint()
    !
    ! Checkpoint for OpenMP
    !
    call par_initialize_omp()

    !----------------------------------------------------------------------
    !
    ! Read main data file *.dat
    !
    !----------------------------------------------------------------------
    !
    ! Read general problem data
    !
    call readat()

    !
    ! Read MPI IO data
    !
    call mpio_config%read()
    call mpio_config%configure()
    call mpio_config%broadcast()
    call mpio_config%print()

    !
    ! Read FTI data
    !
#ifdef ALYA_FTI
    call fti_config%read()
    call fti_config%broadcast()
    call fti_config%check_mpio()
    call fti_config%set_flags()
    call fti_config%print()
#endif

    !
    ! Read PETSc data
    !
#ifdef PETSC
    block
        use mod_alya2petsc, only: alya2petsc_read_configuration
        call alya2petsc_read_configuration()
    end block
#endif
    !
    ! Modules: read data
    !
    do modul = 1, mmodu
        call read_module_options()
    end do
    modul = 0
    !
    ! Parallelization
    !
    call par_turnon()
    !
    call par_initialize_openacc()
    !
    ! Exchange PETSc data
    !
#ifdef PETSC
    block
        use mod_alya2petsc, only: alya2petsc_exchange
        call alya2petsc_exchange()
    end block
#endif

    !----------------------------------------------------------------------
    !
    ! Some initializations
    !
    !----------------------------------------------------------------------
    !
    ! Define some memory output options and open files
    !
    call open_memory_files()
    !
    ! Block ordering
    !
    call modser()
    !
    ! Allocate and initialize postp type
    !
    call output_postprocess_allocate()
    !
    ! Initial variable
    !
    call inivar(1_ip)
    !
    ! Check errors
    !
    call outerr(0_ip)
    !
    ! Synchronization of timeline
    !
    call ker_timeline_synchronization()

end subroutine reapro
