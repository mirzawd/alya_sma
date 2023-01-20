!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tem_rk_explicit

  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_matrix
  use mod_tem_entropy
  use mod_solver,      only : solver_lumped_mass_system
  use mod_solver,      only : solver_explicit
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_tem_entropy, only : tem_entropy_initialization
  use mod_tem_entropy, only : tem_entropy_memory
  use mod_communications, only: PAR_MAX

  implicit none

  real(rp),   pointer :: bt(:)    ! Energy RHS
  real(rp),   pointer :: tt(:)    ! Temperature
  real(rp),   pointer :: tt0(:)   ! Temperature
  real(rp),   pointer :: rt(:,:)  ! Energy residual
  real(rp),   pointer :: Mass(:)  ! lumped mass
  real(rp),   pointer :: Tim(:,:) ! projection dt / (rho*cp)
  real(rp)            :: dt(3)

  real (rp)           :: a(5)
  real (rp)           :: b(5,5)
  real (rp)           :: sigma(5)

  type(soltyp), pointer :: solve1(:)
  type(soltyp), pointer :: solve2(:)
  real(rp),     pointer :: x(:)

  real(rp)              :: cpu_time_matrix
  
  private

  public :: tem_rk_explicit_solution
  public :: tem_rk_explicit_initialization
  public :: tem_rk_explicit_memory

contains 

  subroutine tem_rk_explicit_initialization()

    nullify(bt)    ! Energy RHS
    nullify(tt)    ! Temperature
    nullify(tt0)   ! Temperature
    nullify(rt)    ! Energy residual
    nullify(Mass)  ! lumped mass
    nullify(Tim)   ! projection dt / (rho*cp)
    call tem_entropy_initialization()
    
  end subroutine tem_rk_explicit_initialization
  
  subroutine tem_rk_explicit_memory()

    call memory_alloca(mem_modul(1:2,modul),'RT'  ,'tem_rk_explicit_allocate',rt,  max(1_ip,npoin),5_ip)
    call memory_alloca(mem_modul(1:2,modul),'TIM' ,'tem_rk_explicit_allocate',Tim, max(1_ip,npoin),4_ip)
    call memory_alloca(mem_modul(1:2,modul),'TT0' ,'tem_rk_explicit_allocate',tt0, max(1_ip,npoin))
    
    solve1 => solve(1:)
    solve2 => solve(2:)
    x      => tempe(:,1)
    
    dt = 1.0_rp

    a = 0.0_rp
    b = 0.0_rp
    sigma = 0.0_rp

    if(kfl_tiacc_tem == 1) then !RK1
       a(4) = 1.0_rp
       sigma(4) = 1.0_rp
    else if(kfl_tiacc_tem == 2) then !Heun’s method or RK2 just for validation 
       a(4) = 1.0_rp

       b(4,3) = 1.0_rp

       sigma(3) = 0.5_rp
       sigma(4) = 0.5_rp
    else if(kfl_tiacc_tem == 3) then !  3 order RK
       a(3) = 1.0_rp
       a(4) = 0.5_rp  

       b(3,2) = 1_rp
       b(4,2) = 1.0_rp/4.0_rp
       b(4,3) = 1.0_rp/4.0_rp

       sigma(2) = 1.0_rp/6.0_rp 
       sigma(3) = 1.0_rp/6.0_rp
       sigma(4) = 2.0_rp/3.0_rp
    else if(kfl_tiacc_tem == 4) then ! standard 4 order RK
       a(2) = 0.5_rp
       a(3) = 0.5_rp
       a(4) = 1.0_rp

       b(2,1) = 0.5_rp 
       b(3,2) = 0.5_rp
       b(4,3) = 1.0_rp

       sigma(1) = 1.0_rp/6.0_rp 
       sigma(2) = 1.0_rp/3.0_rp 
       sigma(3) = 1.0_rp/3.0_rp 
       sigma(4) = 1.0_rp/6.0_rp 
    end if

    call tem_entropy_memory()

  end subroutine tem_rk_explicit_memory

  subroutine tem_multi_step_fs_eval(istep)

     integer(ip) , intent(in) :: istep
     integer(ip) :: ipoin
     real(rp)    :: dtinv_tmp,time1,time2
     !real(rp)    :: time1,time2
     
     do ipoin = 1,npoin
        dt_rho_cp_tem(ipoin) = 0.0_rp
     end do

     call cputim(time1)
     call tem_matrix()
     call cputim(time2)
     cpu_time_matrix = cpu_time_matrix + time2-time1
     
     if( INOTEMPTY ) then
        if (kfl_rhs_scal_tem == 0) &
           call solver_lumped_mass_system(1_ip,dt_rho_cp_tem)

        if( solve(1) % kfl_algso == -2 ) then
           dtinv_tem     = dtinv_tmp
           kfl_timei_tem = 1
        end if 

        !
        ! Projection of 1/rhocp for non-conservative  
        !
        if ( kfl_rhs_scal_tem == 0 ) then
           do ipoin = 1,npoin
              Tim(ipoin,1)    = dt_rho_cp_tem(ipoin) / dt(1)
              rt(ipoin,istep) = bt(ipoin)
           end do
        !
        ! Projection of 1/rhocp for conservative  
        !
        else if ( kfl_rhs_scal_tem > 0 ) then
           do ipoin = 1,npoin
              Tim(ipoin,1)    = 1.0_rp
              rt(ipoin,istep) = bt(ipoin)
           end do
        endif
     end if

     call rhsmod(1_ip,rt(:,istep))

     if(istep == 4_ip) then
        call tem_multi_step_fs_solution_sij(4_ip,sigma,1.0_rp)
     else
        call tem_multi_step_fs_solution_sij(istep,b(istep+1,:),a(istep+1))
     endif

  end subroutine tem_multi_step_fs_eval

  subroutine tem_multi_step_fs_solution_sij(istep,weight1,weight2)
     integer(ip) , intent(in) :: istep
     real(rp) ,    intent(in) :: weight1(5)
     real(rp) ,    intent(in) :: weight2
     integer(ip)       :: ipoin,jstep
     real(rp) :: aux2

     do ipoin = 1,npoin
        aux2 = 0.0_rp
        do jstep=1,istep
           !aux2   = aux2 + rt(ipoin,jstep)*b(istep+1,jstep)
           aux2   = aux2 + rt(ipoin,jstep)*weight1(jstep)
        end do
        !tt(ipoin) = tt0(ipoin) + ( Tim(ipoin,1)*aux2*dt(1) ) / Mass(ipoin)
        tt(ipoin) = Tim(ipoin,1)*aux2*dt(1) 
     end do
     !
     ! Solve system 
     !
     call solver_explicit(solve1,tt,EXCHANGE=.false.,solve_consistent=solve2,x=x) 
     do ipoin = 1,npoin 
        tt(ipoin) = tt0(ipoin) + tt(ipoin)
     end do
     !
     ! Dirichet bbcc 
     !
     do ipoin = 1,npoin
        if( kfl_fixno_tem(1,ipoin) > 0 ) &
             tt(ipoin) = bvess_tem(1,ipoin,1)
     end do

  end subroutine tem_multi_step_fs_solution_sij

  subroutine tem_rk_explicit_solution()
     integer(ip) :: ipoin,kfl_advec_old,kfl_timei_old
     real(rp)    :: dtinv_tmp2,dtinv_tmp,time1,time2
     integer(ip), save :: iter = 0
 !    real(rp) :: time1,time2,cpu_time_explicit
      real(rp) :: cpu_time_explicit
     
     iter = iter + 1

     cpu_time_matrix = 0.0_rp
     call cputim(time1)
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

     !
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
     ! guardar la tt de n
     do ipoin = 1,npoin
        tt0(ipoin)   = tt(ipoin)
     end do
     kfl_entpred_tem = 0_ip

     if( kfl_entropy_tem == 1_ip ) call tem_entropy_solution()

     if(     kfl_tiacc_tem == 1) then

        !
        ! Runge Kutta stage
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(4_ip)

     ! order 2 (Heun RK2)

     else if(     kfl_tiacc_tem == 2) then

        !
        ! Runge Kutta stage
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 1_ip
        call tem_multi_step_fs_eval(3_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !
        call tem_endite(ITASK_INNITE)

        !
        ! Runge Kutta stage
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 1_ip
        call tem_multi_step_fs_eval(4_ip)

     ! order 3 (RK3 SSP)

     else if(kfl_tiacc_tem == 3) then  

        !
        ! Runge Kutta stage
        !
        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(2_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !
        call tem_endite(ITASK_INNITE)

        !
        ! Runge Kutta stage
        !
        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(3_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !
        call tem_endite(ITASK_INNITE)

        !
        ! Runge Kutta stage
        !
        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(4_ip)

     else ! order 4

        !
        ! Runge Kutta stage 1
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(1_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !

        call tem_endite(ITASK_INNITE)
        !
        ! Runge Kutta stage 2
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(2_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !
        call tem_endite(ITASK_INNITE)

        !
        ! Runge Kutta stage 3
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(3_ip)

        !
        ! Clipping + update therm(:,:,1) + update thermodynamic pressure
        !
        call tem_endite(ITASK_INNITE)

        !
        ! Runge Kutta stage 4
        !

        !kfl_entpred_tem = 1_ip
        !if(kfl_entropy_tem == 1_ip) call tem_entropy_solution()
        !kfl_entpred_tem = 0_ip
        call tem_multi_step_fs_eval(4_ip)

     end if

     if( INOTMASTER) then
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
     call cputim(time2)
     cpu_time_explicit = time2-time1 -cpu_time_matrix
     call PAR_MAX(cpu_time_explicit)
!     if( IMASTER ) print*,'time explicit= ',cpu_time_explicit
     
  end subroutine tem_rk_explicit_solution



 end module mod_tem_rk_explicit
