!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_perf_csv

    use def_master
    use def_kermod
    use def_domain
    use def_solver
    use mod_strings, only: integer_to_string

    ! CSV writing API
    use mod_perf_csv_writer, only: init_perf
    use mod_perf_csv_writer, only: write_perf_line
    use mod_perf_csv_writer, only: write_perf_header

    ! Memory
    use def_performance, only: r_tomax
    use def_performance, only: max_mem
    use def_performance, only: ave_mem
    use def_performance, only: r_mem_max_modul
    use def_performance, only: r_mem_ave_modul

    ! Time
    use def_performance, only: cpu_total
    use def_performance, only: cpu_max_element
    use def_performance, only: cpu_ave_element
    use def_performance, only: cpu_ave_per_element
    use def_performance, only: cpu_max_boundary
    use def_performance, only: cpu_ave_boundary
    use def_performance, only: cpu_ave_per_boundary
    use def_performance, only: cpu_max_node
    use def_performance, only: cpu_ave_node
    use def_performance, only: cpu_ave_per_node
    use def_performance, only: cpu_max_particle
    use def_performance, only: cpu_ave_particle
    use def_performance, only: cpu_max_min_element
    use def_performance, only: cpu_ave_min_element
    use def_performance, only: cpu_delta_per
    use def_performance, only: cpu_max_solver
    use def_performance, only: cpu_ave_solver
    use def_performance, only: cpu_max_post
    use def_performance, only: cpu_ave_post
    use def_performance, only: cpu_ave_prope
    use def_performance, only: cpu_max_prope
    use def_performance, only: cpu_ave_rst_read
    use def_performance, only: cpu_max_rst_read
    use def_performance, only: cpu_ave_rst_write
    use def_performance, only: cpu_max_rst_write
    use def_performance, only: cpu_ave_iniunk
    use def_performance, only: cpu_max_iniunk
    use def_performance, only: cpu_ave_begrun
    use def_performance, only: cpu_max_begrun
    use def_performance, only: cpu_ave_begite
    use def_performance, only: cpu_max_begite
    use def_performance, only: cpu_ave_endite
    use def_performance, only: cpu_max_endite
    use def_performance, only: cpu_ave_othite
    use def_performance, only: cpu_max_othite
    use def_performance, only: cpu_ave_begste
    use def_performance, only: cpu_max_begste
    use def_performance, only: cpu_ave_endste
    use def_performance, only: cpu_max_endste

    ! Load balance
    use def_performance, only: lb_ave_element
    use def_performance, only: lb_ave_boundary
    use def_performance, only: lb_ave_node
    use def_performance, only: lb_ave_particle
    use def_performance, only: lb_ave_min_element
    use def_performance, only: lb_ave_solver
    use def_performance, only: lb_ave_post

    ! TALP
    use def_performance, only: lb_eff
    use def_performance, only: time_comp_ave
    use def_performance, only: time_comp_max
    use def_performance, only: time_mpi_ave
    use def_performance, only: time_mpi_max

    implicit none

    private

    real(rp) :: cpu_read
    real(rp) :: cpu_init

    public :: performance_csv

contains

    subroutine performance_csv()

        integer(ip)                   :: imodu, ivari, ii, jj
        !
        ! Write performance
        !
        if (INOTSLAVE) then
            call init_perf()
            call write_perf_header()
            !
            ! General
            !
            call write_perf_line('kernel', 'all', 'average time', 's', cpu_total, 'all', '', '', '')
            call write_perf_line('kernel', 'all', 'maximum time', 's', cpu_total, 'all', '', '', '')
#ifdef ALYA_TALP
            call write_perf_line('kernel', 'all', 'load balance', '', lb_eff(1, 0), 'cc', '', '', '')
            call write_perf_line('kernel', 'all', 'communication efficiency', '', lb_eff(2, 0), 'communications', '', '', '')
            call write_perf_line('kernel', 'all', 'parallel efficiency', '', lb_eff(3, 0), 'cc', '', '', '')
            call write_perf_line('kernel', 'all', 'average computation time', 's', time_comp_ave(0), 'computation', '', '', '')
            call write_perf_line('kernel', 'all', 'average communication time', 's', time_MPI_ave(0), 'communications', '', '', '')
            call write_perf_line('kernel', 'all', 'maximum computation time', 's', time_comp_max(0), 'computation', '', '', '')
            call write_perf_line('kernel', 'all', 'maximum communication time', 's', time_MPI_max(0), 'communications', '', '', '')
#endif
            call write_perf_line('kernel', 'all', 'master memory', 'B', r_tomax, 'memory', '', '', '')
            call write_perf_line('kernel', 'all', 'maximum memory', 'B', max_mem, 'memory', '', '', '')
            call write_perf_line('kernel', 'all', 'average memory', 'B', ave_mem, 'memory', '', '', '')
            !
            ! Reading
            !
            cpu_read = cpu_start(CPU_READ_GEO) + cpu_start(CPU_READ_SETS) + cpu_start(CPU_READ_BCS) + cpu_start(CPU_READ_FIELDS)
            ! AVG
            call write_perf_line('read', 'all', 'average time', 's', cpu_read, 'io', 'kernel')
            call write_perf_line('read', 'geometry', 'average time', 's', cpu_start(CPU_READ_GEO), 'io')
            call write_perf_line('read', 'sets', 'average time', 's', cpu_start(CPU_READ_SETS), 'io')
            call write_perf_line('read', 'boundary conditions', 'average time', 's', cpu_start(CPU_READ_BCS), 'io')
            call write_perf_line('read', 'fields', 'average time', 's', cpu_start(CPU_READ_FIELDS), 'io')
            ! MAX
            call write_perf_line('read', 'all', 'maximum time', 's', cpu_read, 'io', 'kernel')
            call write_perf_line('read', 'geometry', 'maximum time', 's', cpu_start(CPU_READ_GEO), 'io')
            call write_perf_line('read', 'sets', 'maximum time', 's', cpu_start(CPU_READ_SETS), 'io')
            call write_perf_line('read', 'boundary conditions', 'maximum time', 's', cpu_start(CPU_READ_BCS), 'io')
            call write_perf_line('read', 'fields', 'maximum time', 's', cpu_start(CPU_READ_FIELDS), 'io')
            !
            ! Initialization/Pre-process
            !
            cpu_init = cpu_start(CPU_MESH_PARTITION) + cpu_start(CPU_MESH_PARTITION) + cpu_start(CPU_CONSTRUCT_DOMAIN) + cpu_start(CPU_ADDTIONAL_ARRAYS)
            ! AVG
            call write_perf_line('init', 'all', 'average time', 's', cpu_init, 'cc', 'kernel')
            call write_perf_line('init', 'mesh partition', 'average time', 's', cpu_start(CPU_MESH_PARTITION), 'cc')
            call write_perf_line('init', 'mesh multiplication', 'average time', 's', cpu_start(CPU_MESH_MULTIPLICATION), 'cc')
            call write_perf_line('init', 'additional arrays', 'average time', 's', cpu_start(CPU_ADDTIONAL_ARRAYS), 'cc')
            ! MAX
            call write_perf_line('init', 'all', 'maximum time', 's', cpu_init, 'cc', 'kernel')
            call write_perf_line('init', 'mesh partition', 'maximum time', 's', cpu_start(CPU_MESH_PARTITION), 'cc')
            call write_perf_line('init', 'mesh multiplication', 'maximum time', 's', cpu_start(CPU_MESH_MULTIPLICATION), 'cc')
            call write_perf_line('init', 'additional arrays', 'maximum time', 's', cpu_start(CPU_ADDTIONAL_ARRAYS), 'cc')
            !
            ! Domain (construct domain)
            !
            ! AVG
            call write_perf_line('domain', 'all', 'average time', 's', cpu_start(CPU_CONSTRUCT_DOMAIN), 'cc', 'init')
            call write_perf_line('domain', 'groups', 'average time', 's', cpu_domain(CPU_GROUPS), 'cc')
            call write_perf_line('domain', 'halos', 'average time', 's', cpu_domain(CPU_HALOS), 'cc')
            call write_perf_line('domain', 'elsest', 'average time', 's', cpu_domain(CPU_ELSEST), 'cc')
            call write_perf_line('domain', 'coupling', 'average time', 's', cpu_domain(CPU_COUPLING), 'cc')
            call write_perf_line('domain', 'output', 'average time', 's', cpu_domain(CPU_OUTPUT_DOMAIN), 'cc')
            ! MAX
            call write_perf_line('domain', 'all', 'maximum time', 's', cpu_start(CPU_CONSTRUCT_DOMAIN), 'cc', 'init')
            call write_perf_line('domain', 'groups', 'maximum time', 's', cpu_domain(CPU_GROUPS), 'cc')
            call write_perf_line('domain', 'halos', 'maximum time', 's', cpu_domain(CPU_HALOS), 'cc')
            call write_perf_line('domain', 'elsest', 'maximum time', 's', cpu_domain(CPU_ELSEST), 'cc')
            call write_perf_line('domain', 'coupling', 'maximum time', 's', cpu_domain(CPU_COUPLING), 'cc')
            call write_perf_line('domain', 'output', 'maximum time', 's', cpu_domain(CPU_OUTPUT_DOMAIN), 'cc')
            !
            ! Kermod
            !
            call write_perf_line(namod(mmodu), 'properties', 'average time', 's', cpu_ave_prope, 'cc', 'kernel')
            call write_perf_line(namod(mmodu), 'properties', 'maximum time', 's', cpu_max_prope, 'cc', 'kernel')

            do imodu = 1, mmodu - 1
                if (kfl_modul(imodu) /= 0) then
                    !
                    ! General
                    !
                    call write_perf_line(namod(imodu), 'all', 'average time', 's', cpu_modul(CPU_TOTAL_MODULE, imodu), 'all', 'kernel')
                    call write_perf_line(namod(imodu), 'all', 'maximum time', 's', cpu_modul(CPU_TOTAL_MODULE, imodu), 'all', 'kernel')
#ifdef ALYA_TALP
                    call write_perf_line(namod(imodu), 'all', 'load balance', '', lb_eff(1, imodu), 'cc', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'communication efficiency', '', lb_eff(2, imodu), 'communications', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'parallel efficiency', '', lb_eff(3, imodu), 'cc', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'average computation time', 's', time_comp_ave(imodu), 'computation', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'average communication time', 's', time_MPI_ave(imodu), 'communications', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'maximum computation time', 's', time_comp_max(imodu), 'computation', '', '', '')
                    call write_perf_line(namod(imodu), 'all', 'maximum communication time', 's', time_MPI_max(imodu), 'communications', '', '', '')
#endif
                    call write_perf_line(namod(imodu), 'all', 'maximum memory', 'B', r_mem_max_modul(imodu), 'memory', 'kernel')
                    call write_perf_line(namod(imodu), 'all', 'average memory', 'B', r_mem_ave_modul(imodu), 'memory', 'kernel')
                    !
                    ! Assembly
                    !
                    if (cpu_modul(CPU_COUNT_ASSEMBLY, imodu) > 0.5_rp) then
                        call write_perf_line(namod(imodu), 'element assembly', 'average time', 's', cpu_ave_element(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'element assembly', 'maximum time', 's', cpu_max_element(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'element assembly', 'load balance', '', lb_ave_element(imodu), 'cc', '', '', '')
                        call write_perf_line(namod(imodu), 'element assembly', 'variability', '', cpu_delta_per(imodu)/100.0_rp, 'cc', '', '', '')
                        call write_perf_line(namod(imodu), 'element assembly', 'time per element', 's', cpu_ave_per_element(imodu), 'cc', '', '', '')
                    end if
                    if (cpu_modul(CPU_COUNT_BOUNDARY, imodu) > 0.5_rp) then
                        call write_perf_line(namod(imodu), 'boundary assembly', 'average time', 's', cpu_ave_boundary(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'boundary assembly', 'maximum time', 's', cpu_max_boundary(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'boundary assembly', 'load balance', '', lb_ave_boundary(imodu), 'cc', '', '', '')
                        call write_perf_line(namod(imodu), 'boundary assembly', 'time per boundary', 's', cpu_ave_per_boundary(imodu), 'cc', '', '', '')
                    end if
                    if (cpu_modul(CPU_COUNT_NODE, imodu) > 0.5_rp) then
                        call write_perf_line(namod(imodu), 'node assembly', 'average time', 's', cpu_ave_node(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'node assembly', 'maximum time', 's', cpu_max_node(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'node assembly', 'load balance', '', lb_ave_node(imodu), 'cc', '', '', '')
                        call write_perf_line(namod(imodu), 'node assembly', 'time per node', 's', cpu_ave_per_node(imodu), 'cc', '', '', '')
                    end if
                    if (cpu_modul(CPU_COUNT_PARTICLE, imodu) > 0.5_rp) then
                        call write_perf_line(namod(imodu), 'particle assembly', 'average time', 's', cpu_ave_particle(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'particle assembly', 'maximum time', 's', cpu_max_particle(imodu), 'cc')
                        call write_perf_line(namod(imodu), 'particle assembly', 'load balance', '', lb_ave_particle(imodu), 'cc', '', '', '')
                    end if
                    !
                    ! IO
                    !
                    call write_perf_line(namod(imodu), 'output', 'average time', 's', cpu_ave_post(imodu), 'io')
                    call write_perf_line(namod(imodu), 'output', 'maximum time', 's', cpu_max_post(imodu), 'io')
                    !
                    ! Begrun
                    !
                    call write_perf_line(namod(imodu), 'run beginning', 'average time', 's', cpu_ave_begrun(imodu), 'all')
                    call write_perf_line(namod(imodu), 'run beginning', 'maximum time', 's', cpu_max_begrun(imodu), 'all')
                    !
                    ! Iniunk
                    !
                    call write_perf_line(namod(imodu), 'initial solution', 'average time', 's', cpu_ave_iniunk(imodu), 'all')
                    call write_perf_line(namod(imodu), 'initial solution', 'maximum time', 's', cpu_max_iniunk(imodu), 'all')
                    !
                    ! Begste
                    !
                    call write_perf_line(namod(imodu), 'step beginning', 'average time', 's', cpu_ave_begste(imodu), 'all')
                    call write_perf_line(namod(imodu), 'step beginning', 'maximum time', 's', cpu_max_begste(imodu), 'all')
                    !
                    ! Endste
                    !
                    call write_perf_line(namod(imodu), 'step ending', 'average time', 's', cpu_ave_endste(imodu), 'all')
                    call write_perf_line(namod(imodu), 'step ending', 'maximum time', 's', cpu_max_endste(imodu), 'all')
                    !
                    ! Begite
                    !
                    call write_perf_line(namod(imodu), 'iteration beginning', 'average time', 's', cpu_ave_begite(imodu), 'all')
                    call write_perf_line(namod(imodu), 'iteration beginning', 'maximum time', 's', cpu_max_begite(imodu), 'all')
                    !
                    ! Endite
                    !
                    call write_perf_line(namod(imodu), 'iteration ending', 'average time', 's', cpu_ave_endite(imodu), 'all')
                    call write_perf_line(namod(imodu), 'iteration ending', 'maximum time', 's', cpu_max_endite(imodu), 'all')
                    !
                    ! Othite
                    !
                    call write_perf_line(namod(imodu), 'iteration other', 'average time', 's', cpu_ave_othite(imodu), 'all')
                    call write_perf_line(namod(imodu), 'iteration other', 'maximum time', 's', cpu_max_othite(imodu), 'all')
                    !
                    ! Restart
                    !
                    call write_perf_line(namod(imodu), 'restart reading', 'average time', 's', cpu_ave_rst_read(imodu), 'io')
                    call write_perf_line(namod(imodu), 'restart reading', 'maximum time', 's', cpu_max_rst_read(imodu), 'io')
                    call write_perf_line(namod(imodu), 'restart writing', 'average time', 's', cpu_ave_rst_write(imodu), 'io')
                    call write_perf_line(namod(imodu), 'restart writing', 'maximum time', 's', cpu_max_rst_write(imodu), 'io')
                    !
                    ! Solvers
                    !
                    solve_sol => momod(imodu)%solve
                    call write_perf_line(namod(imodu), 'solvers', 'average time', 's', cpu_ave_solver(imodu), 'cc')
                    call write_perf_line(namod(imodu), 'solvers', 'maximum time', 's', cpu_max_solver(imodu), 'cc')
                    !
                    if (associated(solve_sol)) then
                        do ivari = 1, size(solve_sol)
                            if (solve_sol(ivari)%kfl_algso /= -999 .and. solve_sol(ivari)%nsolv > 0) then
                                call write_perf_line(namod(imodu), 'spmv '//integer_to_string(ivari), 'average time', 's', solve_sol(ivari)%cpu_spmv(8), "cc", namod(imodu), "solver "//integer_to_string(ivari))
                                call write_perf_line(namod(imodu), 'spmv '//integer_to_string(ivari), 'maximum time', 's', solve_sol(ivari)%cpu_spmv(9), "cc", namod(imodu), "solver "//integer_to_string(ivari))
                                call write_perf_line(namod(imodu), 'dot '//integer_to_string(ivari), 'average time', 's', solve_sol(ivari)%cpu_dot(8), "cc", namod(imodu), "solver "//integer_to_string(ivari))
                                call write_perf_line(namod(imodu), 'dot '//integer_to_string(ivari), 'maximum time', 's', solve_sol(ivari)%cpu_dot(9), "cc", namod(imodu), "solver "//integer_to_string(ivari))
                                call write_perf_line(namod(imodu), 'solver '//integer_to_string(ivari), 'average time', 's', solve_sol(ivari)%cputi(1), "cc", namod(imodu), "solvers")
                                call write_perf_line(namod(imodu), 'solver '//integer_to_string(ivari), 'maximum time', 's', solve_sol(ivari)%cputi(1), "cc", namod(imodu), "solvers")
                            end if
                        end do
                    end if
                    !
                    ! Module perf counters
                    !
                    if (associated(momod(imodu)%times)) then
                        do ii = 1, size(momod(imodu)%times)
                            if (momod(imodu)%times(ii)%used) then
                                call write_perf_line(namod(imodu), trim(momod(imodu)%times(ii)%name), 'average time', 's', momod(imodu)%times(ii)%time_ave, &
                                                     momod(imodu)%times(ii)%category, namod(imodu), momod(imodu)%times(ii)%parent)
                                call write_perf_line(namod(imodu), trim(momod(imodu)%times(ii)%name), 'maximum time', 's', momod(imodu)%times(ii)%time_max, &
                                                     momod(imodu)%times(ii)%category, namod(imodu), momod(imodu)%times(ii)%parent)
                            end if
                        end do
                    end if
                end if
            end do
            imodu = 0
            if (associated(momod(imodu)%times)) then

                if (associated(times)) then
                    do ii = 1, size(times)
                        if (times(ii)%parent == 'kernel') then
                            call write_perf_line(trim(momod(imodu)%times(ii)%name), 'all', 'average time', 's', momod(imodu)%times(ii)%time_ave, &
                                                 momod(imodu)%times(ii)%category, 'kernel', 'all')
                            call write_perf_line(trim(momod(imodu)%times(ii)%name), 'all', 'maximum time', 's', momod(imodu)%times(ii)%time_max, &
                                                 momod(imodu)%times(ii)%category, 'kernel', 'all')
                            do jj = 1, size(times)
                                if (times(jj)%parent == times(ii)%name) then
                                    if (momod(imodu)%times(jj)%used) then
                                        call write_perf_line(trim(momod(imodu)%times(ii)%name), trim(momod(imodu)%times(jj)%name), 'average time', 's', momod(imodu)%times(jj)%time_ave, &
                                                             momod(imodu)%times(jj)%category, momod(imodu)%times(jj)%parent, 'all')
                                        call write_perf_line(trim(momod(imodu)%times(ii)%name), trim(momod(imodu)%times(jj)%name), 'maximum time', 's', momod(imodu)%times(jj)%time_max, &
                                                             momod(imodu)%times(jj)%category, momod(imodu)%times(jj)%parent, 'all')
                                    end if
                                end if
                            end do
                        end if
                    end do
                end if

            end if
        end if

    end subroutine performance_csv

end module mod_perf_csv
!> @}
