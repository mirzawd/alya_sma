!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tem_explicit

  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_solver,         only : solver_lumped_mass_system
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MAX, PAR_MIN
  use mod_tem_entropy,    only : tem_entropy_initialization
  use mod_tem_entropy,    only : tem_entropy_memory

  implicit none

  real(rp),   pointer :: bt(:)    ! Energy RHS
  real(rp),   pointer :: tt(:)    ! Temperature
  real(rp),   pointer :: rt(:,:)  ! Energy residual
  real(rp),   pointer :: Mass(:)  ! lumped mass
  real(rp),   pointer :: Tim(:,:) ! projection dt / (rho*cp)
  real(rp)            :: dt(3)

  private

  public :: tem_explicit_solution
  public :: tem_explicit_initialization
  public :: tem_explicit_memory

contains 

  subroutine tem_explicit_initialization()
    
    nullify(rt)
    nullify(Tim)
    call tem_entropy_initialization()
    
  end subroutine tem_explicit_initialization
  
  subroutine tem_explicit_memory()

    call memory_alloca(mem_modul(1:2,modul),'RT'  ,'tem_explicit_allocate', rt,  max(1_ip,npoin),4_ip)
    call memory_alloca(mem_modul(1:2,modul),'TIM' ,'tem_explicit_allocate',Tim, max(1_ip,npoin),4_ip)

    rt  = 0.0_rp
    Tim = 0.0_rp
    dt  = 1.0_rp
    
    call tem_entropy_memory()
    
  end subroutine tem_explicit_memory


  subroutine tem_explicit_solution()
    integer(ip) :: ipoin,kfl_advec_old,kfl_timei_old, time_order
    real(rp)    :: dtinv_tmp2,dtinv_tmp,reslo
    real(rp)    :: dt01,dt012,dt12,alpha1, alpha2, alpha3, rholoc
    integer(ip), save :: iter = 0

    iter = iter + 1

    !
    ! Update inner iteration counter
    !
    itinn(modul) = itinn(modul) + 1
    ittot_tem    = ittot_tem + 1
    !-----------------------------------------------------------------
    !
    ! Allocate memory if necessary
    !
    !-----------------------------------------------------------------
    
    bt       => rhsid
    tt       => unkno
    Mass     => vmass
    

    !-----------------------------------------------------------------
    !
    ! Assemble equations
    !
    !-----------------------------------------------------------------

    ! Update boundary conditions
    !
    call tem_updbcs(ITASK_INNITE)
    !
    ! Compute temperature gradients
    !
    if( kfl_ellen_tem == -1 ) call gradie(tempe(:,1),grtem_tem)
    !
    ! If initial solution is Stokes: save original values
    !
    if( kfl_inidi_tem == 1 ) then
       kfl_advec_old = kfl_advec_tem
       kfl_timei_old = kfl_timei_tem
       dtinv_tmp2     = dtinv_tem
       kfl_advec_tem = 0
       kfl_timei_tem = 0
       dtinv_tem     = 0.0_rp
       dt            = 0.0_rp
    end if
    !
    ! Construct the system matrix and right-hand-side
    !
    if( solve(1) % kfl_algso == -2 ) then
       dtinv_tmp     = dtinv_tem
       dtinv_tem     = 0.0_rp
       kfl_timei_tem = 0
    end if

    dt(3) = dt(2)
    dt(2) = dt(1)
    dt(1) = 1.0_rp/dtinv

    do ipoin = 1,npoin
       dt_rho_cp_tem(ipoin) = 0.0_rp
    end do

    call tem_matrix()

    call solver_lumped_mass_system(1_ip,dt_rho_cp_tem)

    if( solve(1) % kfl_algso == -2 ) then
       dtinv_tem     = dtinv_tmp
       kfl_timei_tem = 1
    end if 

    !
    ! solve explicit system
    !
    do ipoin = 1,npoin
       Tim(ipoin,1) = dt_rho_cp_tem(ipoin) / dt(1)
       rt(ipoin,1) = bt(ipoin)
    end do


    call rhsmod(1_ip,rt(:,1))

    time_order = min(iter,kfl_tiacc_tem)
    if(      time_order == 1 ) then
       !
       ! 1st order
       !
       alpha1 =  1.0_rp
       alpha2 =  0.0_rp
       alpha3 =  0.0_rp

    else if( time_order == 2 ) then
       !
       ! 2nd order
       !
       dt01   = 1.0_rp
       alpha1 =  ( 1.0_rp + 0.5_rp * dt01)
       alpha2 = -0.5_rp * dt01
       alpha3 =  0.0_rp


    else if( time_order == 3 ) then
       !
       ! 3rd order
       !
       dt012  = dt(1) + dt(2) + dt(3)
       dt01   = dt(1) + dt(2)
       dt12   = dt(2) + dt(3)
       alpha1 = ( 0.25_rp*(dt012*dt01*(dt012+dt01)-dt12*dt(2)*(dt12+dt(2)))&
          &     + 1.0_rp/12.0_rp*(-dt012**3-dt01**3+dt12**3+dt(2)**3))/( dt(2)*dt12)  &
          &   / dt(1)
       alpha2 = ( 0.25_rp*(dt012*dt(1) *(dt012+dt(1))) &
          &     + 1.0_rp/12.0_rp*(-dt012**3-dt(1)**3 +dt12**3))/(-dt(2)*dt(2)) &
          &   / dt(1)
       alpha3 = ( 0.25_rp*(dt01 *dt(1) *(dt01 +dt(1)))&
          &     + 1.0_rp/12.0_rp*(-dt01**3 -dt(1)**3 +dt(2)**3))/(dt(2)*dt12) &
          &   / dt(1)

    end if

    if( INOTMASTER) then
       reslo = 0.0_rp
       do ipoin = 1,npoin
          reslo     = rt(ipoin,1) * alpha1 + rt(ipoin,2) * alpha2 + rt(ipoin,3) * alpha3
          rholoc    = Tim(ipoin,1)
          tt(ipoin) = tt(ipoin) + ( reslo*rholoc )* dt(1) / Mass(ipoin)
       end do
       !
       ! Save old values
       !
       rt(1:npoin,3)  =  rt(1:npoin,2)
       rt(1:npoin,2)  =  rt(1:npoin,1)
       Tim(1:npoin,3) = Tim(1:npoin,2)
       Tim(1:npoin,2) = Tim(1:npoin,1)

       !
       ! Dirichet bbcc 
       !
       do ipoin = 1,npoin
          if( kfl_fixno_tem(1,ipoin) > 0 ) &
             tt(ipoin) = bvess_tem(1,ipoin,1)
       end do

       !
       ! If initial solution is Stokes: recover original values
       !
       if( kfl_inidi_tem == 1 ) then
          kfl_inidi_tem = 0
          kfl_advec_tem = kfl_advec_old 
          kfl_timei_tem = kfl_timei_old 
          dtinv_tem     = dtinv_tmp2
       end if
    end if

 end subroutine tem_explicit_solution

 end module mod_tem_explicit
