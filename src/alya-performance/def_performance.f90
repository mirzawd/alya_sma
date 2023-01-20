!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module def_performance

    use def_master, only: mmodu
    use def_kintyp, only: rp

    implicit none

    public

    real(rp) :: r_tomax                 ! Max memory
    real(rp) :: max_mem                 ! Should be almost=r_tomax... but taken directly from mod_memory
    real(rp) :: ave_mem                 ! Average maximum memory
    real(rp) :: r_mem_max_modul(mmodu)
    real(rp) :: r_mem_ave_modul(mmodu)

    real(rp) :: cpu_total
    real(rp) :: cpu_max_element(mmodu)
    real(rp) :: cpu_ave_element(mmodu)
    real(rp) :: lb_ave_element(mmodu)
    real(rp) :: cpu_ave_per_element(mmodu)
    real(rp) :: cpu_max_boundary(mmodu)
    real(rp) :: cpu_ave_boundary(mmodu)
    real(rp) :: lb_ave_boundary(mmodu)
    real(rp) :: cpu_ave_per_boundary(mmodu)
    real(rp) :: cpu_max_node(mmodu)
    real(rp) :: cpu_ave_node(mmodu)
    real(rp) :: lb_ave_node(mmodu)
    real(rp) :: cpu_ave_per_node(mmodu)
    real(rp) :: cpu_max_min_node(mmodu)
    real(rp) :: cpu_ave_min_node(mmodu)
    real(rp) :: lb_ave_min_node(mmodu)
    real(rp) :: cpu_delta_node_per(mmodu)

    real(rp) :: cpu_element(mmodu)
    real(rp) :: cpu_boundary(mmodu)
    real(rp) :: cpu_node(mmodu)
    real(rp) :: cpu_particle(mmodu)

    real(rp) :: cpu_max_particle(mmodu)
    real(rp) :: cpu_ave_particle(mmodu)
    real(rp) :: lb_ave_particle(mmodu)
    real(rp) :: cpu_max_min_element(mmodu)
    real(rp) :: cpu_ave_min_element(mmodu)
    real(rp) :: lb_ave_min_element(mmodu)
    real(rp) :: cpu_delta_per(mmodu)
    real(rp) :: cpu_max_solver(mmodu)
    real(rp) :: cpu_ave_solver(mmodu)
    real(rp) :: lb_ave_solver(mmodu)
    real(rp) :: cpu_max_post(mmodu)
    real(rp) :: cpu_ave_post(mmodu)
    real(rp) :: lb_ave_post(mmodu)
    real(rp) :: lb_eff(3, 0:mmodu)

    real(rp) :: time_comp_ave(0:mmodu)
    real(rp) :: time_comp_max(0:mmodu)
    real(rp) :: time_mpi_ave(0:mmodu)
    real(rp) :: time_mpi_max(0:mmodu)

    real(rp) :: cpu_ave_prope
    real(rp) :: cpu_max_prope
    real(rp) :: cpu_ave_rst_read(mmodu)
    real(rp) :: cpu_max_rst_read(mmodu)
    real(rp) :: cpu_ave_rst_write(mmodu)
    real(rp) :: cpu_max_rst_write(mmodu)
    real(rp) :: cpu_ave_begste(mmodu)
    real(rp) :: cpu_max_begste(mmodu)
    real(rp) :: cpu_ave_endste(mmodu)
    real(rp) :: cpu_max_endste(mmodu)
    real(rp) :: cpu_ave_iniunk(mmodu)
    real(rp) :: cpu_max_iniunk(mmodu)
    real(rp) :: cpu_ave_begrun(mmodu)
    real(rp) :: cpu_max_begrun(mmodu)
    real(rp) :: cpu_ave_begite(mmodu)
    real(rp) :: cpu_max_begite(mmodu)
    real(rp) :: cpu_ave_endite(mmodu)
    real(rp) :: cpu_max_endite(mmodu)
    real(rp) :: cpu_ave_doiter(mmodu)
    real(rp) :: cpu_max_doiter(mmodu)
    real(rp) :: cpu_ave_othite(mmodu)
    real(rp) :: cpu_max_othite(mmodu)

end module def_performance
!> @}
