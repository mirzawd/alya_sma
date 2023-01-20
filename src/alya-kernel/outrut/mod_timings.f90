!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup CPU_Time 
!> @{
!> @file    mod_timings.f90
!> @author  houzeaux
!> @date    2019-03-01
!> @brief   Compute assembly timings
!> @details Assembly timings
!>          
!-----------------------------------------------------------------------

module mod_timings

  use def_kintyp
  use def_master
  use def_domain
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE

  implicit none

  real(rp) :: cpu_accumulated_element(mmodu)
  real(rp) :: cpu_accumulated_boundary(mmodu)
  real(rp) :: cpu_accumulated_node(mmodu)
  real(rp) :: cpu_accumulated_particle(mmodu)
  real(rp) :: element_load_balance(mmodu)
  real(rp) :: boundary_load_balance(mmodu)
  real(rp) :: node_load_balance(mmodu)
  real(rp) :: particle_load_balance(mmodu)

  real(rp) :: cpu_element(mmodu)
  real(rp) :: cpu_boundary(mmodu)
  real(rp) :: cpu_node(mmodu)
  real(rp) :: cpu_particle(mmodu)
  real(rp) :: cpu_solvers(mmodu)

  real(rp) :: count_element(mmodu)
  real(rp) :: count_boundary(mmodu)
  real(rp) :: count_node(mmodu)
  real(rp) :: count_particle(mmodu)

  real(rp) :: time1,time2
  
  private

  public :: timings_initialization
  public :: timings_assembly
  public :: timings_node_assembly
  public :: timings_particle
  public :: timings_output
  public :: timings_doiter
  public :: timings_unity_test
  public :: timings_ini
  public :: timings_end

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-01
  !> @brief   Initialize timings
  !> @details Initializaiton of timings
  !> 
  !-----------------------------------------------------------------------

  subroutine timings_doiter(CURRENT_MODULE)

    integer(ip), intent(in), optional :: CURRENT_MODULE
    integer(ip)                       :: imod1,imod2,imodu
    real(rp)                          :: cpu_max(5)
    real(rp)                          :: cpu_ave(5)
    !
    ! Timings of solvers is continuously accumulated
    !
    if( present(CURRENT_MODULE) ) then
       imod1 = CURRENT_MODULE
       imod2 = CURRENT_MODULE
    else
       imod1 = 1
       imod2 = mmodu
    end if
    
    do imodu = imod1,imod2
       if( kfl_modul(imodu) /= 0 ) then
          cpu_solvers(imodu)  = cpu_modul(CPU_SOLVER,imodu) - cpu_solvers(imodu)
       end if
    end do
    !
    ! Output max and averages
    !
    cpu_max = 0.0_rp    
    do imodu = imod1,imod2
       if( kfl_modul(imodu) /= 0 ) then
          cpu_max(1) = cpu_element(imodu)
          cpu_max(2) = cpu_boundary(imodu)
          cpu_max(3) = cpu_node(imodu)
          cpu_max(4) = cpu_particle(imodu)
          cpu_ave    = cpu_max
          call PAR_MAX    (4_ip,cpu_max)
          call PAR_AVERAGE(4_ip,cpu_ave)
          if( INOTSLAVE ) then
             write(UNIT=momod(modul) % lun_timin,FMT=10)      &
                  ittim,itcou,cutim,                          &
                  cpu_modul(ITASK_DOITER,modul),              &
                  count_element(imodu) ,cpu_max(1),cpu_ave(1),&
                  count_boundary(imodu),cpu_max(2),cpu_ave(2),&
                  count_node(imodu)    ,cpu_max(3),cpu_ave(3),&
                  count_particle(imodu),cpu_max(4),cpu_ave(4),&
                  cpu_solvers(imodu)
          end if
       end if
    end do
    !
    ! Reset timings
    !
    cpu_element    = 0.0_rp
    cpu_boundary   = 0.0_rp
    cpu_node       = 0.0_rp
    cpu_particle   = 0.0_rp
    count_element  = 0.0_rp
    count_boundary = 0.0_rp
    count_node     = 0.0_rp
    count_particle = 0.0_rp
    do imodu = imod1,imod2
       if( kfl_modul(imodu) /= 0 ) then
          cpu_solvers(imodu) = cpu_modul(CPU_SOLVER,imodu)
       end if
    end do
    
10  format(1x,i9,1x,i9,1x,50(1x,es13.6))
    
  end subroutine timings_doiter

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-01
  !> @brief   Initialize timings
  !> @details Initializaiton of timings
  !> 
  !-----------------------------------------------------------------------

  subroutine timings_initialization()

    cpu_accumulated_element  = 0.0_rp
    cpu_accumulated_boundary = 0.0_rp
    cpu_accumulated_node     = 0.0_rp
    cpu_accumulated_particle = 0.0_rp
    element_load_balance     = 0.0_rp
    boundary_load_balance    = 0.0_rp
    node_load_balance        = 0.0_rp
    particle_load_balance    = 0.0_rp
    
    cpu_element              = 0.0_rp
    cpu_boundary             = 0.0_rp
    cpu_node                 = 0.0_rp
    cpu_particle             = 0.0_rp
    cpu_solvers              = 0.0_rp
    
    count_element            = 0.0_rp
    count_boundary           = 0.0_rp
    count_node               = 0.0_rp
    count_particle           = 0.0_rp
    
  end subroutine timings_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-01
  !> @brief   Element assembly
  !> @details Compute some element assembly timings
  !> 
  !-----------------------------------------------------------------------

  subroutine timings_assembly(time1,time2,CURRENT_MODULE,TYPE_OF_ASSEMBLY)

    real(rp),                   intent(in) :: time1
    real(rp),         optional, intent(in) :: time2
    integer(ip)     , optional, intent(in) :: CURRENT_MODULE
    character(len=*), optional, intent(in) :: TYPE_OF_ASSEMBLY
    integer(ip)                            :: imodu
    real(rp)                               :: times_max(2)
    real(rp)                               :: times_ave(2)
    logical(lg)                            :: if_element
    logical(lg)                            :: if_boundary
    logical(lg)                            :: if_node
    logical(lg)                            :: if_particle

    if( present(TYPE_OF_ASSEMBLY) ) then
       if_element  = .false.
       if_boundary = .false.
       if_node     = .false.
       if_particle = .false.
       if( index(TYPE_OF_ASSEMBLY,'ELEMENT')  /= 0 ) if_element  = .true.
       if( index(TYPE_OF_ASSEMBLY,'BOUNDARY') /= 0 ) if_boundary = .true.
       if( index(TYPE_OF_ASSEMBLY,'NODE')     /= 0 ) if_node     = .true.
       if( index(TYPE_OF_ASSEMBLY,'PARTICLE') /= 0 ) if_particle = .true.
    else
       if_element  = .true.
       if_boundary = .true.
       if_node     = .false.
       if_particle = .false.
    end if
    
    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if
    !
    ! Element assembly
    !
    if( if_element ) then
       cpu_modul(CPU_COUNT_ASSEMBLY,imodu) = cpu_modul(CPU_COUNT_ASSEMBLY,imodu) + 1.0_rp   ! k
       cpu_modul(CPU_ASSEMBLY,imodu)       = cpu_modul(CPU_ASSEMBLY,imodu)       + time1    ! sum_k ti(k)
       cpu_element(imodu)                  = cpu_element(imodu)                  + time1
       count_element(imodu)                = count_element(imodu)                + 1.0_rp
       cpu_modul(CPU_MINI_ASSEMBLY,imodu)  = min(cpu_modul(CPU_MINI_ASSEMBLY,imodu),time1)  ! min_k ti(k)
       cpu_modul(CPU_MAXI_ASSEMBLY,imodu)  = max(cpu_modul(CPU_MAXI_ASSEMBLY,imodu),time1)  ! max_k ti(k)
    end if
    !
    ! Boundary assembly
    !
    if( if_boundary ) then
       if( present(time2) ) then
          cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu) = cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu) + time2
          cpu_modul(CPU_COUNT_BOUNDARY,imodu)    = cpu_modul(CPU_COUNT_BOUNDARY,imodu)    + 1.0_rp
          cpu_boundary(imodu)                    = cpu_boundary(imodu)                    + time2
          count_boundary(imodu)                  = count_boundary(imodu)                  + 1.0_rp
       end if
    end if
    if( if_element .or. if_boundary ) then
       !
       ! Accumulate maximum times
       !
       times_max    = 0.0_rp
       times_ave    = 0.0_rp    
       times_max(1) = time1
       times_ave(1) = time1
       if( present(time2) ) then
          times_max(2) = time2
          times_ave(2) = time2
       end if    
       call PAR_MAX(2_ip,times_max)
       call PAR_AVERAGE(2_ip,times_ave)
       
       cpu_accumulated_element(imodu)  = cpu_accumulated_element(imodu)  + times_max(1)                                     ! sum_k max_i(ti)
       if( times_max(1) > zeror ) element_load_balance(imodu)  = element_load_balance(imodu)  + times_ave(1) / times_max(1) ! sum_k (max_ti/ave_ti)
       
       if( present(time2) ) then
          if( times_max(2) > zeror ) boundary_load_balance(imodu) = boundary_load_balance(imodu) + times_ave(2) / times_max(2)
          cpu_accumulated_boundary(imodu) = cpu_accumulated_boundary(imodu) + times_max(2)
       end if
    end if
    !
    ! Node assembly
    !
    if( if_node ) then
       cpu_modul(CPU_ASSEMBLY_NODE,imodu) = cpu_modul(CPU_ASSEMBLY_NODE,imodu) + time1
       cpu_modul(CPU_COUNT_NODE,imodu)    = cpu_modul(CPU_COUNT_NODE,imodu)    + 1.0_rp
       cpu_modul(CPU_MINI_NODE,imodu)     = min(cpu_modul(CPU_MINI_NODE,imodu),time1)  ! min_k ti(k)
       cpu_modul(CPU_MAXI_NODE,imodu)     = max(cpu_modul(CPU_MAXI_NODE,imodu),time1)  ! max_k ti(k)
       cpu_node(imodu)                    = cpu_node(imodu)                    + time1
       count_node(imodu)                  = count_node(imodu)                  + 1.0_rp
       times_max                          = 0.0_rp
       times_ave                          = 0.0_rp    
       times_max(1)                       = time1
       times_ave(1)                       = time1
       call PAR_MAX(1_ip,times_max)
       call PAR_AVERAGE(1_ip,times_ave)       
       cpu_accumulated_node(imodu)        = cpu_accumulated_node(imodu)  + times_max(1)                               ! sum_k max_i(ti)
       if( times_max(1) > zeror ) node_load_balance(imodu)  = node_load_balance(imodu)  + times_ave(1) / times_max(1) ! sum_k (max_ti/ave_ti)
    end if
    !
    ! Particle assembly
    !
    if( if_particle ) then
       cpu_modul(CPU_ASSEMBLY_PARTICLE,imodu) = cpu_modul(CPU_ASSEMBLY_PARTICLE,imodu) + time1
       cpu_modul(CPU_COUNT_PARTICLE,imodu)    = cpu_modul(CPU_COUNT_PARTICLE,imodu)    + 1.0_rp
       cpu_particle(imodu)                    = cpu_particle(imodu)                    + time1
       count_particle(imodu)                  = count_particle(imodu)                  + 1.0_rp
       times_max                              = 0.0_rp
       times_ave                              = 0.0_rp    
       times_max(1)                           = time1
       times_ave(1)                           = time1
       call PAR_MAX(1_ip,times_max)
       call PAR_AVERAGE(1_ip,times_ave)       
       cpu_accumulated_particle(imodu)        = cpu_accumulated_particle(imodu)  + times_max(1)                               ! sum_k max_i(ti)
       if( times_max(1) > zeror ) particle_load_balance(imodu)  = particle_load_balance(imodu)  + times_ave(1) / times_max(1) ! sum_k (max_ti/ave_ti)
    end if
        
  end subroutine timings_assembly

  subroutine timings_output(&
       CURRENT_MODULE,      &
       ! element
       cpu_element,         &  
       cpu_max_element,     &  
       cpu_ave_element,     &  
       lb_ave_element,      &
       cpu_ave_per_element, &
       cpu_max_min_element, &  
       cpu_ave_min_element, &  
       lb_ave_min_element,  &
       cpu_delta_per,       &       
       ! boundary
       cpu_boundary,        &   
       cpu_max_boundary,    &  
       cpu_ave_boundary,    &  
       lb_ave_boundary,     &    
       cpu_ave_per_boundary,&
       ! node
       cpu_node,            &  
       cpu_max_node,        &  
       cpu_ave_node,        &  
       lb_ave_node,         &    
       cpu_ave_per_node,    &       
       cpu_max_min_node,    &
       cpu_ave_min_node,    &
       lb_ave_min_node,     &
       cpu_delta_node_per,  &
       ! particle
       cpu_particle,        &  
       cpu_max_particle,    &  
       cpu_ave_particle,    &  
       lb_ave_particle,     &    
       cpu_max_solver,      &
       cpu_ave_solver,      &
       lb_ave_solver,       &
       cpu_max_post,        &
       cpu_ave_post,        &
       lb_ave_post          )

    integer(ip), intent(in)  :: CURRENT_MODULE

    real(rp),    intent(out) :: cpu_element
    real(rp),    intent(out) :: cpu_max_element
    real(rp),    intent(out) :: cpu_ave_element
    real(rp),    intent(out) :: lb_ave_element
    real(rp),    intent(out) :: cpu_ave_per_element

    real(rp),    intent(out) :: cpu_max_min_element
    real(rp),    intent(out) :: cpu_ave_min_element
    real(rp),    intent(out) :: lb_ave_min_element
    real(rp),    intent(out) :: cpu_delta_per

    real(rp),    intent(out) :: cpu_boundary
    real(rp),    intent(out) :: cpu_max_boundary
    real(rp),    intent(out) :: cpu_ave_boundary
    real(rp),    intent(out) :: lb_ave_boundary
    real(rp),    intent(out) :: cpu_ave_per_boundary

    real(rp),    intent(out) :: cpu_node
    real(rp),    intent(out) :: cpu_max_node
    real(rp),    intent(out) :: cpu_ave_node
    real(rp),    intent(out) :: lb_ave_node
    real(rp),    intent(out) :: cpu_ave_per_node

    real(rp),    intent(out) :: cpu_max_min_node
    real(rp),    intent(out) :: cpu_ave_min_node
    real(rp),    intent(out) :: lb_ave_min_node
    real(rp),    intent(out) :: cpu_delta_node_per
    
    real(rp),    intent(out) :: cpu_particle
    real(rp),    intent(out) :: cpu_max_particle
    real(rp),    intent(out) :: cpu_ave_particle
    real(rp),    intent(out) :: lb_ave_particle

    real(rp),    intent(out) :: cpu_max_solver
    real(rp),    intent(out) :: cpu_ave_solver
    real(rp),    intent(out) :: lb_ave_solver

    real(rp),    intent(out) :: cpu_max_post
    real(rp),    intent(out) :: cpu_ave_post
    real(rp),    intent(out) :: lb_ave_post

    integer(ip)              :: imodu

    cpu_max_element     = 0.0_rp
    cpu_ave_element     = 0.0_rp
    lb_ave_element      = 0.0_rp
    cpu_ave_per_element = 0.0_rp

    cpu_max_min_element = 0.0_rp
    cpu_ave_min_element = 0.0_rp
    lb_ave_min_element  = 0.0_rp
    cpu_delta_per       = 0.0_rp

    cpu_max_boundary    = 0.0_rp
    cpu_ave_boundary    = 0.0_rp
    lb_ave_boundary     = 0.0_rp
    cpu_ave_per_boundary= 0.0_rp

    cpu_max_node        = 0.0_rp
    cpu_ave_node        = 0.0_rp
    lb_ave_node         = 0.0_rp
    cpu_ave_per_node    = 0.0_rp

    cpu_max_min_node    = 0.0_rp
    cpu_ave_min_node    = 0.0_rp
    lb_ave_min_node     = 0.0_rp
    cpu_delta_node_per  = 0.0_rp
    
    cpu_max_particle    = 0.0_rp
    cpu_ave_particle    = 0.0_rp
    lb_ave_particle     = 0.0_rp
                       
    cpu_max_solver      = 0.0_rp
    cpu_ave_solver      = 0.0_rp
    lb_ave_solver       = 0.0_rp

    cpu_max_post        = 0.0_rp
    cpu_ave_post        = 0.0_rp
    lb_ave_post         = 0.0_rp
    
    imodu               = CURRENT_MODULE
    !
    ! Elements
    !
    if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
       cpu_ave_element     = cpu_modul(CPU_ASSEMBLY,imodu)
       cpu_max_element     = cpu_accumulated_element(imodu)
       cpu_element         = cpu_accumulated_element(imodu)
       lb_ave_element      = element_load_balance(imodu) / cpu_modul(CPU_COUNT_ASSEMBLY,imodu)
       cpu_ave_per_element = cpu_modul(CPU_ASSEMBLY,imodu) / ( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) * real(max(1_ip,nelem),rp) )
       call PAR_MAX(    cpu_max_element    ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_element    ,'IN MY CODE')
       call PAR_AVERAGE(lb_ave_element     ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_per_element,'IN MY CODE')
       !
       ! Theoretical elements... my duration assuming always minimum time
       !    
       cpu_ave_min_element = cpu_modul(CPU_MINI_ASSEMBLY,imodu) * cpu_modul(CPU_COUNT_ASSEMBLY,imodu)
       cpu_max_min_element = cpu_modul(CPU_MINI_ASSEMBLY,imodu) * cpu_modul(CPU_COUNT_ASSEMBLY,imodu)
       lb_ave_min_element  = lb_ave_min_element / cpu_modul(CPU_COUNT_ASSEMBLY,imodu)
       call PAR_MAX(    cpu_max_min_element,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_min_element,'IN MY CODE')
       if( cpu_max_min_element > zeror ) lb_ave_min_element  = cpu_ave_min_element / (cpu_max_min_element+zeror)
       !
       ! Variability % in element assembly
       !
       if( cpu_modul(CPU_MAXI_ASSEMBLY,imodu) < cpu_modul(CPU_MINI_ASSEMBLY,imodu) ) then
          cpu_delta_per = 0.0_rp
       else
          cpu_delta_per = &
               100.0_rp * (cpu_modul(CPU_MAXI_ASSEMBLY,imodu)-cpu_modul(CPU_MINI_ASSEMBLY,imodu)) &
               / ( (zeror+cpu_modul(CPU_ASSEMBLY,imodu))/(zeror+cpu_modul(CPU_COUNT_ASSEMBLY,imodu)))
       end if
       call PAR_MAX(cpu_delta_per,'IN MY CODE')!,'INCLUDE MASTER')
       
    end if
    !
    ! Boundaries
    !
    if( cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp ) then
       cpu_ave_boundary     = cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu)
       cpu_boundary         = cpu_accumulated_boundary(imodu)
       cpu_max_boundary     = cpu_accumulated_boundary(imodu)
       lb_ave_boundary      = boundary_load_balance(imodu) / cpu_modul(CPU_COUNT_BOUNDARY,imodu)
       cpu_ave_per_boundary = cpu_modul(CPU_ASSEMBLY_BOUNDARY,imodu) / ( cpu_modul(CPU_COUNT_BOUNDARY,imodu) * real(max(1_ip,nboun),rp) )
       call PAR_MAX(    cpu_max_boundary    ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_boundary    ,'IN MY CODE')
       call PAR_AVERAGE(lb_ave_boundary     ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_per_boundary,'IN MY CODE')
    end if
    !
    ! Nodes
    !
    if( cpu_modul(CPU_COUNT_NODE,imodu) > 0.5_rp ) then
       cpu_ave_node     = cpu_modul(CPU_ASSEMBLY_NODE,imodu)
       cpu_node         = cpu_accumulated_node(imodu)
       cpu_max_node     = cpu_accumulated_node(imodu)
       lb_ave_node      = node_load_balance(imodu) / cpu_modul(CPU_COUNT_NODE,imodu)
       cpu_ave_per_node = cpu_modul(CPU_ASSEMBLY_NODE,imodu) / ( cpu_modul(CPU_COUNT_NODE,imodu) * real(max(1_ip,npoin),rp) )
       call PAR_MAX(    cpu_max_node    ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_node    ,'IN MY CODE')
       call PAR_AVERAGE(lb_ave_node     ,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_per_node,'IN MY CODE')
       !
       ! Theoretical elements... my duration assuming always minimum time
       !    
       cpu_ave_min_node = cpu_modul(CPU_MINI_NODE,imodu) * cpu_modul(CPU_COUNT_NODE,imodu)
       cpu_max_min_node = cpu_modul(CPU_MINI_NODE,imodu) * cpu_modul(CPU_COUNT_NODE,imodu)
       lb_ave_min_node  = lb_ave_min_node / cpu_modul(CPU_COUNT_NODE,imodu)
       call PAR_MAX(    cpu_max_min_node,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_min_node,'IN MY CODE')
       if( cpu_max_min_node > zeror ) lb_ave_min_node  = cpu_ave_min_node / (cpu_max_min_node+zeror)
       !
       ! Variability % in node assembly
       !
       if( cpu_modul(CPU_MAXI_NODE,imodu) < cpu_modul(CPU_MINI_NODE,imodu) ) then
          cpu_delta_node_per = 0.0_rp
       else
          cpu_delta_node_per = &
               100.0_rp * (cpu_modul(CPU_MAXI_NODE,imodu)-cpu_modul(CPU_MINI_NODE,imodu)) &
               / ( (zeror+cpu_modul(CPU_ASSEMBLY_NODE,imodu))/(zeror+cpu_modul(CPU_COUNT_NODE,imodu)))
       end if
       call PAR_MAX(cpu_delta_node_per,'IN MY CODE')!,'INCLUDE MASTER')       
    end if
    !
    ! Particles
    !
    if( cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then
       cpu_ave_particle = cpu_modul(CPU_ASSEMBLY_PARTICLE,imodu)
       cpu_particle     = cpu_accumulated_particle(imodu)
       cpu_max_particle = cpu_accumulated_particle(imodu)
       lb_ave_particle  = particle_load_balance(imodu) / cpu_modul(CPU_COUNT_PARTICLE,imodu)
       call PAR_MAX(    cpu_max_particle,'IN MY CODE')
       call PAR_AVERAGE(cpu_ave_particle,'IN MY CODE')
       call PAR_AVERAGE(lb_ave_particle ,'IN MY CODE')
    end if
    !
    ! Solver
    !
    cpu_ave_solver = cpu_modul(CPU_SOLVER,imodu) + cpu_modul(CPU_EIGEN_SOLVER,imodu)
    cpu_max_solver = cpu_modul(CPU_SOLVER,imodu) + cpu_modul(CPU_EIGEN_SOLVER,imodu)
    call PAR_MAX(    cpu_max_solver,'IN MY CODE')
    call PAR_AVERAGE(cpu_ave_solver,'IN MY CODE')
    if( cpu_max_solver > zeror ) lb_ave_solver = cpu_ave_solver / (cpu_max_solver+zeror) 
    !
    ! Postprocess
    !
    cpu_ave_post = cpu_modul(CPU_OUTPUT,imodu)
    cpu_max_post = cpu_modul(CPU_OUTPUT,imodu)
    call PAR_MAX(    cpu_max_post,'IN MY CODE',INCLUDE_ROOT=.true.)
    call PAR_AVERAGE(cpu_ave_post,'IN MY CODE')
    if( cpu_max_post > zeror ) lb_ave_post = cpu_ave_post / (cpu_max_post+zeror) 

  end subroutine timings_output

  subroutine timings_unity_test(time_element,time_boundary)

    real(rp),  intent(out) :: time_element
    real(rp),  intent(out) :: time_boundary
    integer(ip), save      :: kk=0

    !      Matrix element assembly:   
    !           Average:                    2.42 ( 13.71 % )
    !           Maximum:                    3.50 ( 19.85 % )
    !           Load balance:               0.73
    !           ---                   
    !           Theoretical average:        1.67
    !           Theoretical maximum:        2.00
    !           Theoretical load bal.:      0.83
    !           Maximum variability:       80.00 % 

    kk = kk + 1
    
    if( kfl_paral == 1 ) then
       if(      kk == 1 ) then
          time_element  = 1.0_rp
          time_boundary = 1.0_rp
       else if( kk == 2 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       else if( kk == 3 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       else if( kk == 4 ) then
          time_element  = 1.0_rp
          time_boundary = 1.0_rp
      end if
    else if( kfl_paral == 2 ) then
       if(      kk == 1 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       else if( kk == 2 ) then
          time_element  = 1.0_rp
          time_boundary = 1.0_rp
       else if( kk == 3 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       else if( kk == 4 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       end if
    else if( kfl_paral == 3 ) then
       if(      kk == 1 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       else if( kk == 2 ) then
          time_element  = 0.25_rp
          time_boundary = 0.25_rp
       else if( kk == 3 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
      else if( kk == 4 ) then
          time_element  = 0.5_rp
          time_boundary = 0.5_rp
       end if
    end if
    
  end subroutine timings_unity_test

  subroutine timings_node_assembly(time_node,CURRENT_MODULE)

    real(rp),              intent(in) :: time_node
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(ip)                       :: imodu
    real(rp)                          :: times_max(2)
    real(rp)                          :: times_ave(2)

    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if
    !
    ! Node assembly
    !
    cpu_modul(CPU_COUNT_NODE,imodu)     = cpu_modul(CPU_COUNT_NODE,imodu)     + 1.0_rp          ! k
    cpu_modul(CPU_ASSEMBLY_NODE,imodu)  = cpu_modul(CPU_ASSEMBLY_NODE,imodu)  + time_node    ! sum_k ti
    cpu_node(imodu)                     = cpu_node(imodu)                  + time_node
    count_node(imodu)                   = count_node(imodu)                + 1.0_rp
    !cpu_modul(CPU_MINI_ASSEMBLY,imodu)  = min(cpu_modul(CPU_MINI_ASSEMBLY,imodu),time_node)  ! min_k ti
    !cpu_modul(CPU_MAXI_ASSEMBLY,imodu)  = max(cpu_modul(CPU_MAXI_ASSEMBLY,imodu),time_node)  ! max_k ti
    !
    ! Accumulate maximum times
    !
    times_max    = 0.0_rp
    times_ave    = 0.0_rp
    
    times_max(1) = time_node
    times_ave(1) = time_node
    
    call PAR_MAX(1_ip,times_max)
    call PAR_AVERAGE(1_ip,times_ave)

    cpu_accumulated_node(imodu)  = cpu_accumulated_node(imodu)  + times_max(1)                                     ! sum_k max_i(ti)
    if( times_max(1) > zeror ) node_load_balance(imodu)  = node_load_balance(imodu)  + times_ave(1) / times_max(1) ! sum_k (max_ti/ave_ti)

  end subroutine timings_node_assembly

  subroutine timings_particle(time_particle,CURRENT_MODULE)

    real(rp),              intent(in) :: time_particle
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(ip)                       :: imodu
    real(rp)                          :: times_max(2)
    real(rp)                          :: times_ave(2)

    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if
    !
    ! Particle assembly
    !
    cpu_modul(CPU_COUNT_PARTICLE,imodu)     = cpu_modul(CPU_COUNT_PARTICLE,imodu)     + 1.0_rp          ! k
    cpu_modul(CPU_ASSEMBLY_PARTICLE,imodu)  = cpu_modul(CPU_ASSEMBLY_PARTICLE,imodu)  + time_particle    ! sum_k ti
    cpu_particle(imodu)                     = cpu_particle(imodu)                  + time_particle
    count_particle(imodu)                   = count_particle(imodu)                + 1.0_rp
    !cpu_modul(CPU_MINI_ASSEMBLY,imodu)  = min(cpu_modul(CPU_MINI_ASSEMBLY,imodu),time_particle)  ! min_k ti
    !cpu_modul(CPU_MAXI_ASSEMBLY,imodu)  = max(cpu_modul(CPU_MAXI_ASSEMBLY,imodu),time_particle)  ! max_k ti
    !
    ! Accumulate maximum times
    !
    times_max    = 0.0_rp
    times_ave    = 0.0_rp
    
    times_max(1) = time_particle
    times_ave(1) = time_particle
    
    call PAR_MAX(1_ip,times_max)
    call PAR_AVERAGE(1_ip,times_ave)

    cpu_accumulated_particle(imodu)  = cpu_accumulated_particle(imodu)  + times_max(1)                                     ! sum_k max_i(ti)
    if( times_max(1) > zeror ) particle_load_balance(imodu)  = particle_load_balance(imodu)  + times_ave(1) / times_max(1) ! sum_k (max_ti/ave_ti)

  end subroutine timings_particle

    !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-14
  !> @brief   Initialize timing
  !> @details Initialize timing
  !> 
  !-----------------------------------------------------------------------
  
  subroutine timings_ini()

    call cputim(time1)
    
  end subroutine timings_ini

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-14
  !> @brief   End timing
  !> @details End timing
  !> 
  !-----------------------------------------------------------------------
  
  subroutine timings_end(itask,timec)

    integer(ip), optional, intent(in)    :: itask
    real(rp),    optional, intent(inout) :: timec
    integer(ip)                          :: jtask
    
    call cputim(time2)
    if( present(timec) ) then
       timec = timec + (time2 - time1)
    else if( present(itask) ) then
       select case( itask )
       case ( ITASK_BEGITE ) ; jtask = CPU_BEGITE
       case ( ITASK_ENDITE ) ; jtask = CPU_ENDITE
       case default ; call runend('MOD_TIMINGS: UNKNOWN TASK')
       end select

       cpu_modul(jtask,modul) = cpu_modul(jtask,modul) + (time2 - time1)

    end if
    
  end subroutine timings_end
  
end module mod_timings
!> @}
