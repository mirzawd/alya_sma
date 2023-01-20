!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Begrun
!> @{
!> @file    Begrun.f90
!> @author  Guillaume Houzeaux
!> @brief   Begin the run
!> @details Initial operations: do some computations on the domain,
!>          essentially topological operations
!>
!> @}
!-----------------------------------------------------------------------
subroutine Begrun()

  use def_kintyp,               only : ip,rp
  use def_domain,               only : meshe
  use def_kermod,               only : ndivi
  use def_master,               only : IPARALL
  use def_master,               only : cpu_start
  use def_master,               only : CPU_ADDTIONAL_ARRAYS
  use def_master,               only : ITASK_BEGRUN
  use def_master,               only : nblok
  use def_master,               only : iblok
  use mod_moduls,               only : moduls
  use mod_auto_tuning,          only : auto_tuning_SpMV_OpenMP
  use mod_messages,             only : messages_report
  use mod_communications,       only : PAR_BARRIER
  use mod_par_affinity,         only : par_affinity
  use mod_messages,             only : messages_live
  use mod_outfor,               only : outfor
  use mod_messages,             only : messages_live
  use mod_mpio_par_async_io,    only : PAR_ASYNCHRONOUS_WRITE_POP_ALL
  use mod_par_element_loop,     only : par_element_loop
  use mod_ker_detection,        only : ker_detection_open_file
  use mod_open_close_files,           only : close_data_file
  implicit none

  real(rp)    :: time1,time2

  call cputim(time1)
  !
  ! Output optimization options
  !
  call outfor(-73_ip,0_ip,' ')
  !
  ! Parallel stuffs: send some data computed from modules data
  !
  call par_sendat(7_ip)
  ! 
  ! Hybrid parallelization: this has been moved to domain
  !
  !call par_element_loop()
  !call par_boundary_loop()
!!$  BLOCK
!!$    use def_domain, only : ompss_domains
!!$    use mod_parall, only : list_elements_par
!!$    use def_domain, only : ompss_boundaries
!!$    use mod_parall, only : list_boundaries_par
!!$    use mod_parall, only : num_pack_nboun_par
!!$    use mod_parall, only : num_pack_par
!!$   call parall_openmp_adjacency_ompss_unity_test(&
!!$         ompss_domains,list_elements_par,num_pack_par,npoin,lnnod,lnods)
!!$    call parall_openmp_adjacency_ompss_unity_test(&
!!$         ompss_boundaries,list_boundaries_par,num_pack_nboun_par,npoin,lnnob,lnodb)
!!$  END BLOCK
  call outfor(-72_ip)

  !----------------------------------------------------------------------
  !
  ! Required arrays depending on modules data
  !
  !----------------------------------------------------------------------
  !
  ! List of required arrays:
  !
  call reqarr()
  !
  ! PELPO_2, LELPO_2, PELEL_2, LELEL_2 extended Graphs
  !
  call lelpo2()
  !
  ! ZDOM_*, R_DOM_*, C_DOM:*. Schur type solvers and Aii preconditioners. * = Aii, Aib, Abi, Abb
  !
  call solpre()
  !
  ! R_SYM, C_SYM: Symmetric graph if necessary
  !
  call symgra()
  !
  ! LELBF, LELFA: Global face graph and element boundary face
  !
  call lgface()
  !
  ! Element bin for neighboring elements
  !
  call elebin()

  !----------------------------------------------------------------------
  !
  ! Others
  !
  !----------------------------------------------------------------------
  !
  ! Check warnings and errors
  !
  call outerr(1_ip)
  !
  ! Read restart file
  !
  call restar(1_ip)
  !
  ! Compute groups: UNDER DEVELOPEMENT
  !
!!!call ker_groups()
  !
  ! Initialization of some variables and solvers
  !
  call inivar(2_ip)
  !
  ! Close data file/open occasional files
  !
  call close_data_file()
  call ker_detection_open_file()
  !
  ! Support geometry
  !
  call ker_submsh()
  !
  ! Auto-tuning operations
  !
  call auto_tuning_SpMV_OpenMP()
  ! 
  ! Write information
  !
  call outcpu_operations(meshe(ndivi))

  call outinf()
  !
  ! Report
  !
  call messages_report()
  !
  !call initialize coprocessing
  !
#ifdef CATA
  call coprocessorinitializewithpython("coproc.py",9)
  call messages_live('CALL INITIALIZE COPROCESSING')
#endif
  !
  ! Transient fields
  !
  call calc_kx_tran_fiel()
  !
  ! Asynchronous MPIO
  !
  if (IPARALL) then
     call PAR_ASYNCHRONOUS_WRITE_POP_ALL(0_ip)
  end if
  
  call cputim(time2)
  cpu_start(CPU_ADDTIONAL_ARRAYS) = cpu_start(CPU_ADDTIONAL_ARRAYS) + time2 - time1  
  !
  ! Call modules
  !
  call Kermod(ITASK_BEGRUN)
  do iblok = 1,nblok
     call moduls(ITASK_BEGRUN)
  end do

end subroutine Begrun
