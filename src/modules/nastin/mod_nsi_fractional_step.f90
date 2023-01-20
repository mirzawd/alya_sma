!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_solsch.f90
!> @author  Guillaume Houzeaux
!> @date    20/09/2016
!> @brief   Fractional step
!> @details Fractional step
!> @}
!-----------------------------------------------------------------------
module mod_nsi_fractional_step

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use mod_postpr
  use mod_matrix
  use mod_memory,                      only : memory_alloca
  use mod_memory,                      only : memory_deallo
  use mod_nsi_schur_operations,        only : nsi_apuvec
  use mod_nsi_schur_operations,        only : nsi_aupvec
  use mod_nsi_schur_operations,        only : nsi_appvec
  use mod_nsi_schur_operations,        only : nsi_solini
  use mod_matrices,                    only : matrices_gradient_divergence
  use mod_matrices,                    only : matrices_laplacian
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_matrix_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_couplings,                   only : couplings_impose_dirichlet
  use mod_messages,                    only : livinf
  use mod_local_basis,                 only : local_basis_global_to_local
  use mod_local_basis,                 only : local_basis_local_to_global

  implicit none

  real(rp)            :: gamma1

  real(rp),   pointer :: vv(:)      ! Auxiliar vector
  real(rp),   pointer :: rm(:,:)    ! Residual of momentum at time
  real(rp),   pointer :: rc(:)      ! residual of continuity
  real(rp),   pointer :: divup(:)   ! div up
  real(rp),   pointer :: rt(:)      ! Total residual of momentum
  real(rp),   pointer :: rt_tmp(:)  ! Total residual of momentum
  real(rp),   pointer :: bu(:)      ! Momentum RHS
  real(rp),   pointer :: bp(:)      ! Continuity RHS
  real(rp),   pointer :: uu(:)      ! Velocity
  real(rp),   pointer :: pp(:)      ! Pressure
  real(rp),   pointer :: L(:)       ! Continuous -Laplacian = int_V grad(u).grad(v) dV
  real(rp),   pointer :: Dp(:)      ! p^{n+1} - gamma * p^n
  real(rp),   pointer :: Dpp(:)     ! p^{n+1} - p^n
  real(rp),   pointer :: Aupp(:)    ! SMVP Aup p
  real(rp),   pointer :: bp2(:)
  real(rp),   pointer :: bp1(:)

  real(rp),   pointer :: Tim(:,:)   ! Nodal projection of dt/rho
  real(rp),   pointer :: Tau(:)     ! Nodal projection of tau
  real(rp),   pointer :: Mass(:)    ! Mass matrix
  real(rp),   pointer :: Grad(:,:)  ! Gradient matrix   G_ij=-\int (div Ni) Nj
  real(rp),   pointer :: Div(:,:)   ! Divergence matrix D_ij= \int Ni (div Nj) => D=-G^t
  real(rp),   pointer :: Lapl(:)    ! Laplacian matrix

  private

  public :: nsi_fractional_step_solution
  public :: nsi_fractional_step_initialization
  public :: nsi_fractional_step_memory
  public :: nsi_fractional_step_matrices

contains

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Initialization
  !> @details Initialization
  !
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_initialization()

    nullify(vv)
    nullify(rm)
    nullify(rt)
    nullify(rc)
    nullify(divup)
    nullify(bp2)
    nullify(bp1)
    nullify(Dp)
    nullify(Dpp)
    nullify(Aupp)
    nullify(Tim)
    nullify(Tau)
    nullify(rt_tmp)

    nullify(Grad)
    nullify(Div)
    nullify(Lapl)
    nullify(L)

  end subroutine nsi_fractional_step_initialization

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Memory
  !> @details Memory
  !
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_memory()

    call memory_alloca(mem_modul(1:2,modul),'VV'    ,'mod_nsi_fractional_step',vv    ,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'RM'    ,'mod_nsi_fractional_step',rm    ,max(1_ip,ndime*npoin),4_ip) ! A discutir
    call memory_alloca(mem_modul(1:2,modul),'RT'    ,'mod_nsi_fractional_step',rt    ,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'RC'    ,'mod_nsi_fractional_step',rc    ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DIVUP' ,'mod_nsi_fractional_step',divup ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'BP2'   ,'mod_nsi_fractional_step',bp2   ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'BP1'   ,'mod_nsi_fractional_step',bp1   ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DP'    ,'mod_nsi_fractional_step',Dp    ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DPP'   ,'mod_nsi_fractional_step',Dpp   ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUPP'  ,'mod_nsi_fractional_step',Aupp  ,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'TIM'   ,'mod_nsi_fractional_step',Tim   ,max(1_ip,npoin),3_ip)
    call memory_alloca(mem_modul(1:2,modul),'TAU'   ,'mod_nsi_fractional_step',Tau   ,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'RT_TMP','mod_nsi_fractional_step',rt_tmp,max(1_ip,ndime*npoin))

  end subroutine nsi_fractional_step_memory

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Fractional step
  !> @details Compute matrices
  !
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_matrices

    !integer(ip), save :: ipass = 0

    !ipass = ipass + 1
    !if( ipass /= 1 ) return

    bu       => rhsid
    bp       => rhsid(ndbgs_nsi+1:)
    uu       => unkno
    pp       => unkno(ndbgs_nsi+1:)       
    Mass     => vmass
    gamma1   =  1.0_rp - gamma_nsi

    if( kfl_grad_div_nsi /= 0 ) call matrices_gradient_divergence(Grad,Div,DEALLOCATE_MATRICES=.true.)
    if( kfl_grad_div_nsi /= 0 ) call matrices_laplacian(Lapl,DEALLOCATE_MATRIX=.true.)
    if( kfl_grad_div_nsi /= 0 ) call nsi_dirichlet_matrix_rhs(Q=Lapl,ROTATION=.false.,DIRICHLET=.true.)
    if( kfl_grad_div_nsi /= 0 ) then
       L  => Lapl
    else
       L  => lapla_nsi
    end if

  end subroutine nsi_fractional_step_matrices

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Fractional step
  !> @details Assemble NS equations and solve them using a
  !>          fractional step method
  !
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_solution()

    use mod_matrix, only : matrix_transpose

    real(rp) :: time1

    !-----------------------------------------------------------------
    !
    ! Allocate memory if necessary
    !
    !-----------------------------------------------------------------

!!!!call nsi_fractional_step_matrices()

    !-----------------------------------------------------------------
    !
    ! Assemble equations
    !
    !-----------------------------------------------------------------

    call nsi_solini(3_ip)                                      ! Initialize solver
    if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
       call local_basis_global_to_local(kfl_fixrs_nsi,uu)      ! Global to local
    end if
    call nsi_matrix()                                          ! Assemble equation

    call livinf(56_ip,' ',modul)        
    call livinf(160_ip,' ',1_ip)
    call cputim(time1)

    !-----------------------------------------------------------------
    !
    ! Fractional step algorithm
    !
    !-----------------------------------------------------------------

    if(      kfl_press_stab_nsi == 0 ) then
       !
       ! Classical fractional step (default)
       !
       call nsi_fractional_step_solution_dt()

    else if( kfl_press_stab_nsi >= 1 ) then
       !
       ! Using tau
       !
       call nsi_fractional_step_solution_tau()

    end if

    !-----------------------------------------------------------------
    !
    ! Go back to global system
    !
    !-----------------------------------------------------------------

    call livinf(165_ip,'C',0_ip)
    if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
       call local_basis_local_to_global(kfl_fixrs_nsi,uu) ! Local to global
    end if
    call livinf(164_ip,' ',1_ip)

  end subroutine nsi_fractional_step_solution

  !-----------------------------------------------------------------------
  !>
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Fractional step
  !> @details Solution of Navier-Stokes using an explicit fractional
  !>          step method. From the following paper:
  !>
  !>          R. Codina, "Pressure Stability in Fractional Step Finite
  !>          Element Methods for Incompressible Flows",
  !>          J. Compt. Phys 170, 112-140 (2001)
  !>
  !>          Let Tim be the nodal projection of dt/rho
  !>          and Tau be the nodal projection of tau
  !>          L is minus the Laplacian: L := \int_V grad(u).grad(v) dv
  !>
  !>          1. Solve momentum explicitly:
  !>
  !>             Tim^-1 M ( u* - u^n ) = bu - Auu u^n - gamma * Aup p^n
  !>
  !>          2. Solve pressure equation:
  !>
  !>             Tim L Dp^{n+1}  = bp - Apu u* - App p^n <=>
  !>
  !>             Tim L Dp'^{n+1} = bp - Apu u* - App p^n + (gamma-1) Tim L p^n
  !>
  !>          3. Correct velocity:
  !>
  !>             u^{n+1} = u* - Tim * M^-1 Aup Dp^{n+1}
  !>
  !>          Dp^{n+1}' := p^{n+1} - p^n
  !>          Dp^{n+1}  := p^{n+1} - gamma * p^n
  !>          rm^{n+1}  := bu - Auu u^n - Aup p^{n+1}
  !>          rc^{n+1}  := bp - Apu u^{n+1} - App p^{n+1}
  !>
  !>          Subsitute Eq. (3) in (2) and (3) into (1):
  !>
  !>          Tim^-1 * M * ( u^{n+1} - u^n ) = bu - Auu u^n - Aup p^{n+1}
  !>          rc^{n+1} + ( Tim L + Apu Tim M^{-1} Aup ) Dp^{k+1} = 0
  !>
  !>          +-                 -+ +-        -+   +-                             -+
  !>          ! Tim^{-1} M   Aup  | | u^{n+1}  |   | bu - Auu u^n + Tim^{-1} M u^n |
  !>          !                   | |          | = |                               |
  !>          ! Apu           S   | | p^{n+1}  |   | bp + gamma * S * p^n          |
  !>          +-                 -+ +-        -+   +-                             -+
  !>
  !>          where S =  Tim L + Apu Tim M^{-1} Aup
  !>
  !>          To do incremental with ASGS:   put ASGS and gamma = 1
  !>          To do FS with OSS in pressure: put no stabilization and gamma = 0
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_solution_dt()

    !-----------------------------------------------------------------
    !
    ! Solve Momentum u*
    !
    !-----------------------------------------------------------------

    call nsi_fractional_step_momentum()

    !-----------------------------------------------------------------
    !
    ! Solve pressure p^{n+1} 
    !
    !-----------------------------------------------------------------

    call nsi_fractional_step_continuity()

    !-----------------------------------------------------------------
    !
    ! Correct velocity u^{n+1}
    !
    !-----------------------------------------------------------------

    call nsi_fractional_velocity_correction()

    !-----------------------------------------------------------------
    !
    ! Momentum residual ||rm||
    !
    !-----------------------------------------------------------------

    if( kfl_grad_div_nsi /= 0 )  then
       call nsi_aupvec(1_ip,Grad,pp,Aupp)                               ! rm =  rm - G p^{n+1}
    else
       call nsi_aupvec(1_ip,amatr(poaup_nsi:),pp,Aupp)                  ! rm =  rm - Aup p^{n+1}
    end if

    call rhsmod(ndime,Aupp)
    rt = rt - Aupp    

    if( INOTMASTER ) rt_tmp = rt
    call nsi_fractional_step_dirichlet(rt_tmp)
    call norm2x(ndime,rt_tmp,resin_nsi(1))                           ! || rm ||

  end subroutine nsi_fractional_step_solution_dt

  !-----------------------------------------------------------------------
  !>
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Boundary condition on velocity
  !> @details Prescribe Dirichlet boundary condition on velocity
  !>          A priori, this is useless if everything was done well
  !>          before...
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_dirichlet(xarray)

    real(rp),   intent(inout), pointer, optional :: xarray(:)
    integer(ip)                                  :: idofn,ipoin,idime

    if( present(xarray) ) then
       idofn = 0
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idofn + 1
             if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                  xarray(idofn) = 0.0_rp
          end do
       end do
    else
       idofn = 0
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idofn + 1
             if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                  uu(idofn) = bvess_nsi(idime,ipoin,1)
          end do
       end do
    end if

  end subroutine nsi_fractional_step_dirichlet
  
  !-----------------------------------------------------------------------
  !>
  !> @date    04/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Fractional step stabilized with tau
  !> @details Solve pressure equation so that final stabilization involves
  !>          tau and not dt. Only valid for gamma=0.
  !>          We want the continuity equation to be:
  !>
  !>          Apu u^{n+1} + ( Tau L + Apu Tau M^{-1} Aup ) p^{n+1} = bp
  !>
  !>          The question is: what should be S so that rewritting equation
  !>          2) as:
  !>          S p^{n+1} = bp - Apu u*
  !>          together with equation 3)
  !>          u* = u^{n+1} + Tim M^{-1} Aup p^{n+1}
  !>          gives the right stabilization. Substituting the first into the
  !>          second equation we obtain:
  !>
  !>          S p^{n+1} = bp - Apu u^{n+1} - Apu Tim M^{-1} Aup p^{n+1}
  !>
  !>          ( S + Apu T M^{-1} Aup ) p^{n+1} = bp - Apu u^{n+1}
  !>
  !>          We rthus equire that
  !>
  !>          ( S + Apu T M^{-1} Aup ) = Tau L + Apu Tau M^{-1} Aup, so that
  !>
  !>          S = Tau L + Apu (Tau-Tim) M^{-1}. Thus equation 2) is
  !>
  !>          ( L + Tau^{-1} Apu (Tau-Tim) M^{-1} Aup ) p^{n+1} = Tau^{-1} (bp - Apu u*)
  !>
  !>          Preconditioned system to solve:
  !>          -------------------------------
  !>
  !>          L^{-1} A p^{n+1} = L^{-1} b
  !>
  !>          A := L + Tau^{-1} Apu (Tau-Tim) M^{-1} Aup
  !>          b := Tau^{-1} ( bp - Apu u* )
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_solution_tau()

    real(rp)            :: resid
    real(rp)            :: pnorm

    !-----------------------------------------------------------------
    !
    ! Solve Momentum u*
    !
    !-----------------------------------------------------------------

    call nsi_fractional_step_momentum()

    !-----------------------------------------------------------------
    !
    ! Solve pressure p^{n+1}
    !
    !-----------------------------------------------------------------

    if(      kfl_press_stab_nsi == 3 ) then
       !
       ! CG
       !
       call nsi_fractional_step_continuity_tau_cg()

    else if( kfl_press_stab_nsi == 4 ) then
       !
       ! BiCGSTAB
       !
       call nsi_fractional_step_continuity_tau_bicgstab()

    else
       !
       ! Richardson and Orthomin(1)
       !
       call nsi_fractional_step_continuity_tau()

    end if
    !
    ! Residual of continuity equation: || p^{n+1}-p^{n} || / || p^{n+1} ||
    !
    call norm2x(1_ip,pp, pnorm)
    call norm2x(1_ip,Dpp,resid)
    resin_nsi(2) = resid / max(pnorm,zeror)

    !-----------------------------------------------------------------
    !
    ! Correct velocity u^{n+1}
    !
    !-----------------------------------------------------------------

    call nsi_fractional_velocity_correction()

    !-----------------------------------------------------------------
    !
    ! Momentum residual ||rm||
    !
    !-----------------------------------------------------------------

    if( kfl_grad_div_nsi /= 0 )  then
       call nsi_aupvec(1_ip,Grad,pp,Aupp)                               ! rm =  rm - G p^{n+1}
    else
       call nsi_aupvec(1_ip,amatr(poaup_nsi:),pp,Aupp)                  ! rm =  rm - Aup p^{n+1}
    end if
    call rhsmod(ndime,Aupp)
    rt = rt - Aupp
    call norm2x(ndime,rt,resin_nsi(1))                                  ! || rm ||

  end subroutine nsi_fractional_step_solution_tau

  !-----------------------------------------------------------------------
  !>
  !> @date    26/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Solve momentum equation
  !> @details Solve Momentum u*
  !>
  !>          To devise second order for du/dt = rm.
  !>          From Tyulor series:
  !>          (1) u^{n+1}  = u^{n} + rm^{n} dt^{n} + 1/2 d^2u/dt^2 (dt^n)^2 + O(dt^3)
  !>          (2) rm^{n-1} = rm^{n} - d^2u/dt^2 dt^{n-1} + O(dt^2)
  !>
  !>          From (2) we have:
  !>
  !>          d^2u/dt^2 = (rm^{n} - rm^{n-1}) / dt^{n-1} + O(dt)
  !>
  !>          Subsitute in (1):
  !>
  !>          u^{n+1} = u^{n} + rm^{n} dt^{n} + 1/2 (dt^n)^2 (rm^{n} - rm^{n-1}) / dt^{n-1} + O(dt^3)
  !>
  !>          u^{n+1} = u^{n} + rm^{n} dt^{n} ( 1 + dt^{n}/2dt^{n-1} ) - rm^{n-1} (dt^{n})^2/2dt^{n-1}
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_momentum()

    integer(ip) :: idofn,ipoin,idime
    integer(ip) :: kfl_tiacc_fs_nsi
    real(rp)    :: dt012,dt01,dt12,time1,time2
    real(rp)    :: alpha1,alpha2,alpha3

    call cputim(time1)

    if( INOTMASTER ) then
       !
       ! Projections of dt/rho and tau
       !
       Tim(:,1) = dt_rho_nsi
       Tau      = tau_nsi
       !
       ! Momentum residual rm^k
       !
       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,pp,Aupp)                           ! G p^k
       else          
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),pp,Aupp)              ! Aup p^k
       end if

       call rhsmod(ndime,Aupp)

       rm(1:ndime*npoin,1) = bu(1:ndime*npoin)                         ! rm = bu
       call rhsmod(ndime,rm(:,1))
       !
       ! Decide order
       !
       kfl_tiacc_fs_nsi = min(ittim,kfl_tiacc_nsi)

       do ipoin = 1,npoin

          if(      kfl_tiacc_fs_nsi == 1 ) then
             !
             ! 1st order
             !
             alpha1 =  1.0_rp
             alpha2 =  0.0_rp
             alpha3 =  0.0_rp

          else if( kfl_tiacc_fs_nsi == 2 ) then
             !
             ! 2nd order
             !
             dt01   =  Tim(ipoin,1) / Tim(ipoin,2)
             alpha1 =  ( 1.0_rp + 0.5_rp * dt01)
             alpha2 = -0.5_rp * dt01
             alpha3 =  0.0_rp

          else if( kfl_tiacc_fs_nsi == 3 ) then
             !
             ! 3rd order
             !
             dt012  = Tim(ipoin,1) + Tim(ipoin,2) + Tim(ipoin,3)
             dt01   = Tim(ipoin,1) + Tim(ipoin,2)
             dt12   = Tim(ipoin,2) + Tim(ipoin,3)
             alpha1 = ( 0.25_rp*(dt012*dt01*(dt012+dt01)-dt12*Tim(ipoin,2)*(dt12+Tim(ipoin,2)))&
                  &     + 1.0_rp/12.0_rp*(-dt012**3-dt01**3+dt12**3+Tim(ipoin,2)**3))/( Tim(ipoin,2)*dt12)  &
                  &   / Tim(ipoin,1)
             alpha2 = ( 0.25_rp*(dt012*Tim(ipoin,1) *(dt012+Tim(ipoin,1))) &
                  &     + 1.0_rp/12.0_rp*(-dt012**3-Tim(ipoin,1)**3 +dt12**3))/(-Tim(ipoin,2)*Tim(ipoin,2)) &
                  &   / Tim(ipoin,1)
             alpha3 = ( 0.25_rp*(dt01 *Tim(ipoin,1) *(dt01 +Tim(ipoin,1)))&
                  &     + 1.0_rp/12.0_rp*(-dt01**3 -Tim(ipoin,1)**3 +Tim(ipoin,2)**3))/( Tim(ipoin,2)*dt12) &
                  &   / Tim(ipoin,1)
          end if

          do idime = 1,ndime
             idofn         = (ipoin-1)*ndime + idime
             rt(idofn)     = rm(idofn,1) * alpha1 + rm(idofn,2) * alpha2 + rm(idofn,3) * alpha3
             rt_tmp(idofn) = ( rt(idofn) - gamma_nsi * Aupp(idofn) ) * Tim(ipoin,1) 
             !uu(idofn)   = uu(idofn) + ( rt(idofn) - gamma_nsi * Aupp(idofn) ) * Tim(ipoin,1)
          end do

       end do
       !
       ! u = u + M^-1 r
       !
       call solver_lumped_mass_system(ndime,rt_tmp,EXCHANGE=.false.) ! CAMBIO ORIOL
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn     = (ipoin-1)*ndime + idime
             uu(idofn) = uu(idofn) + rt_tmp(idofn)
          end do
       end do
       !
       ! Save old values
       !
       rm(:,3)  = rm(:,2)
       rm(:,2)  = rm(:,1)
       Tim(:,3) = Tim(:,2)
       Tim(:,2) = Tim(:,1)
       !
       ! Prescribe Dirichlet on u'^{k+1}
       !
       call nsi_fractional_step_dirichlet()                        ! Prescribe bc

    end if

    call cputim(time2)

    call livinf(165_ip,'M',0_ip)

  end subroutine nsi_fractional_step_momentum

  !-----------------------------------------------------------------------
  !>
  !> @date    26/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Solve continuity equation
  !> @details Solve pressure p^{n+1}
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_continuity()

    real(rp) :: time1,time2
    integer (ip) :: ipoin

    call cputim(time1)
    !
    ! Continuity residual rc'= Tim^-1 ( bp - Apu u* ) + L (gamma-1) p^n
    !
    if( INOTMASTER ) then
       if( kfl_grad_div_nsi /= 0 ) then
          call nsi_apuvec(1_ip,Div,uu,divup)
       else
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,divup)             ! rc' =  Apu u*
       end if
       rc = ( bp - divup ) / Tim(:,1)
       call nsi_appvec(0_ip,L,pp,rc,-gamma1)                           ! rc' <= rc' + (gamma-1) L p^n
       if( nodpr_nsi > 0 ) rc(nodpr_nsi) = 0.0_rp
    end if
    !
    ! Impose Dirichlet
    !
    if( kfl_exist_fib02_nsi == 1 .or. kfl_exist_fib20_nsi == 1 ) then
       do ipoin = 1,npoin
          if( kfl_fixpp_nsi(1,ipoin) > 0 ) then
             rc(ipoin) = 0.0_rp
          end if
       end do
    end if
    !
    ! Solve: L z = rc'
    !
    call nsi_solini(2_ip)                                              ! Initialize solver
    call solver(rc,Dpp,L,pmatr)                                        ! Solve system    
    !
    ! Actualize p^{n+1}
    ! Dp      = Dp' + (1-gamma)*p^n
    ! p^{n+1} = Dp' + p^n
    !
    if( INOTMASTER ) then
       Dp(1:npoin) = Dpp(1:npoin) + gamma1 * pp(1:npoin)
       do ipoin=1,npoin
          pp(ipoin) = Dpp(ipoin) + pp(ipoin) !- fsrot_nsi * Nu(ipoin) * divup(ipoin) / Mass(ipoin) 
       end do
       
       !call nsi_fractional_step_impose_dirichlet_pressure()            ! Prescribe BC

    end if
    !
    ! Residual norm: ||rc|| := bp - Apu u^{*} - App p^{k+1} + (gamma-1) L p^k
    !
    call norm2x(1_ip,rc,resin_nsi(2))                                  ! || rc ||

    call cputim(time2)


    call livinf(165_ip,'S',0_ip)

  end subroutine nsi_fractional_step_continuity

  !-----------------------------------------------------------------------
  !>
  !> @date    26/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Velocity correction
  !> @details Correct velocity u^{n+1}
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_velocity_correction()

    integer(ip) :: idofn,ipoin,idime
    real(rp)    :: time1,time2

    call cputim(time1)

    if( INOTMASTER ) then

       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,Dp,vv)                              ! rm =  G Dpp^{k+1}
       else
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),Dp,vv)                 ! rm =  Aup Dpp^{k+1}
       end if
       call rhsmod(ndime,vv)

       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn     = (ipoin-1)*ndime + idime
             rt_tmp(idofn) = vv(idofn) * Tim(ipoin,1)
             !uu(idofn) = uu(idofn) - vv(idofn) * Tim(ipoin,1) / Mass(ipoin) ! CAMBIO ORIOL
          end do
       end do

       call solver_lumped_mass_system(ndime,rt_tmp,EXCHANGE=.false.) ! CAMBIO ORIOL
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn     = (ipoin-1)*ndime + idime
             uu(idofn) = uu(idofn) - rt_tmp(idofn)
          end do
       end do

       call nsi_fractional_step_dirichlet()                             ! Prescribe bc

    end if

    call cputim(time2)


  end subroutine nsi_fractional_velocity_correction

  !-----------------------------------------------------------------------
  !>
  !> @date    26/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Solve pressure equation
  !> @details Solve pressure equation A x = b.
  !>          Two iteratives methods are implemented:
  !>
  !>          1. Richardson:
  !>
  !>          x^{k+1} = x^{k} + L^{-1} A r^k
  !>
  !>          2. Orthomin:
  !>
  !>          x^{k+1} = x^{k} + alpha L^{-1} r^k
  !>          r^{k+1} = r^{k} - alpha A L^{-1} r^k
  !>
  !>                      ( r^k, A L^{-1} r^k )
  !>          alpha = -----------------------------
  !>                  ( A L^{-1} r^k , A L^{-1} r^k )
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_continuity_tau()

    integer(ip)         :: kfl_recov
    integer(ip)         :: maxit,itera
    real(rp)            :: toler,resid
    real(rp)            :: denom
    real(rp),   pointer :: ww(:)            ! Working array
    real(rp),   pointer :: bb(:)            ! RHS
    real(rp),   pointer :: xx(:)            ! Solution
    real(rp),   pointer :: zz(:)            ! Preconditioned residual
    real(rp)            :: pnorm            ! || p^n ||
    real(rp)            :: bnorm            ! || b ||
    real(rp)            :: alpha            ! Orthomin relaxation factor
    real(rp)            :: invnb            ! 1 / || b ||

    call nsi_solini(2_ip)                   ! Initialize solver
    kfl_recov = solve_sol(1) % kfl_recov    ! Save solver's default
    solve_sol(1) % kfl_recov = 2            ! Solver RHS is already exchanged

    nullify(ww)
    nullify(bb)
    nullify(xx)
    nullify(zz)
    call memory_alloca(mem_modul(1:2,modul),'WW','nsi_fractional_step_solution_tau',ww,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'ZZ','nsi_fractional_step_solution_tau',zz,max(1_ip,npoin))
    !
    ! Norm of pressure || p ||
    !
    call norm2x(1_ip,pp,pnorm)
    !
    ! b = Tau^{-1} ( bp - Aup u* )
    ! r = b - A x
    !
    if( INOTMASTER ) then
       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_apuvec(1_ip,Div,uu,bb)
       else
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,bb)
       end if
       bb = ( bp - bb ) / Tau
       call rhsmod(1_ip,bb)
       xx = pp
    end if
    !
    ! || b ||
    !
    call norm2x(1_ip,bb,bnorm)
    if( bnorm <= zeror ) goto 1
    invnb = 1.0_rp / bnorm
    !
    ! r = b - A x, || r ||
    !
    if( INOTMASTER ) then
       call nsi_fractional_step_multiply(xx,rc)
       rc = bb - rc
    end if
    call norm2x(1_ip,rc,resid)
    resid = resid * invnb
    !
    ! Tolerance and max number of iterations
    !
    if( kfl_adres_nsi == 1 ) then
       toler = max(toler_nsi,adres_nsi * resid)
    else
       toler = toler_nsi
    end if
    maxit = mitri_nsi
    itera = 0
    !
    ! Solve A x = b
    !
    do while( itera < maxit .and. resid > toler )

       itera = itera + 1
       !
       ! Solve preconditioned system L z = r^k
       !
       call solver(rc,zz,L,pmatr)

       if(      kfl_press_stab_nsi == 1 ) then
          !
          ! Richardson
          !
          if( INOTMASTER ) then
             xx(1:npoin) = xx(1:npoin) + zz(1:npoin)
             call nsi_fractional_step_multiply(xx,rc)
             rc(1:npoin) = bb(1:npoin) - rc(1:npoin)
          end if

       else if( kfl_press_stab_nsi == 2 ) then
          !
          ! Orthomin(1)
          !
          call nsi_fractional_step_multiply(zz,ww)
          call prodxy(1_ip,npoin,rc,ww,alpha)
          call prodxy(1_ip,npoin,ww,ww,denom)
          !
          ! alpha = ( r^k , A L^{-1} r^k ) / ( A L^{-1} r^k , A L^{-1} r^k )
          !
          alpha = alpha / max(denom,zeror)
          !
          ! x^{k+1} = x^k + alpha * L^{-1} r^k
          ! r^{k+1} = r^k - alpha * A L^{-1} r^k
          !
          if( INOTMASTER ) then
             xx(1:npoin) = xx(1:npoin) + alpha * zz(1:npoin)
             rc(1:npoin) = rc(1:npoin) - alpha * ww(1:npoin)
          end if

       end if
       !
       ! Residual norm || r^{k+1} || / || b ||
       !
       call norm2x(1_ip,rc,resid)
       resid = resid * invnb

!!!!!write(92,*) itera,resid,bnorm  ; flush(92)

    end do

1   continue
    call livinf(165_ip,'S('//trim(intost(itera))//')',0_ip)
    !
    ! Actualize p^{n+1}
    ! Dp' = p^{n+1} - p^n
    ! Dp  = p^{n+1}
    !
    if( INOTMASTER ) then
       Dpp = xx - pp
       Dp  = xx
       pp  = xx
    end if

    call memory_deallo(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx)
    call memory_deallo(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb)
    call memory_deallo(mem_modul(1:2,modul),'WW','nsi_fractional_step_solution_tau',ww)
    call memory_deallo(mem_modul(1:2,modul),'ZZ','nsi_fractional_step_solution_tau',zz)
    !
    ! Recover solver's default
    !
    solve_sol(1) % kfl_recov = kfl_recov

  end subroutine nsi_fractional_step_continuity_tau

  !-----------------------------------------------------------------------
  !>
  !> @date    26/01/2017
  !> @author  Guillaume Houzeaux and Oriol Lehmkuhl
  !> @brief   Perform a matrix vector multiplication
  !> @details Compute y = ( L - S ) x, where
  !>           S = Tau^{-1} Apu (Tim-Tau) M^{-1} Aup
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_fractional_step_multiply(xx,yy,message)

    real(rp),    intent(in),    pointer  :: xx(:)
    real(rp),    intent(inout), pointer  :: yy(:)
    character(*),               optional :: message
    real(rp),                   pointer  :: Aupx(:,:)
    real(rp),                   pointer  :: zz(:)
    integer(ip)                          :: ipoin
    logical(lg)                          :: use_transpose_Apu
    !
    ! Check if transpose should be used: Aup <= -Apu^t
    !
    use_transpose_Apu = .false.
    if( present(message) ) then
       if( trim(message) == 'TRANSPOSE' ) then
          use_transpose_Apu = .true.
       end if
    end if

    nullify(Aupx)
    nullify(zz)
    call memory_alloca(mem_modul(1:2,modul),'AUPX','nsi_fractional_step_multiply',Aupx,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'ZZ',  'nsi_fractional_step_multiply',zz,npoin)

    if( INOTMASTER ) then
       !
       ! Aup x
       !
       if( use_transpose_Apu ) then
          call nsi_aupvec(1_ip,amatr(poapu_nsi:),xx,Aupx,-1.0_rp,'TRANSPOSE')     ! -Apu^t x
       else
          if( kfl_grad_div_nsi /= 0 )  then
             call nsi_aupvec(1_ip,Grad,xx,Aupx)                                   ! G x
          else
             call nsi_aupvec(1_ip,amatr(poaup_nsi:),xx,Aupx)                      ! Aup x
          end if
       end if
       call rhsmod(ndime,Aupx)
       !
       ! (Tau-Tim) M^{-1} Aup x
       !
       do ipoin = 1,npoin
          Aupx(:,ipoin) = (Tau(ipoin)-Tim(ipoin,1))/Mass(ipoin) * Aupx(:,ipoin)
       end do
       !
       ! y = Apu (Tau-Tim) M^{-1} Aup x ]
       ! z = L x
       !
       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_apuvec(1_ip,Div,Aupx,yy)
       else
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),Aupx,yy)
       end if

       call nsi_appvec(1_ip,L,xx,zz)
       !
       ! y = ( L + Apu (Tau-Tim) M^{-1} Aup ) x
       !
       yy = zz + yy / Tau
       call rhsmod(1_ip,yy)
    end if

    call memory_deallo(mem_modul(1:2,modul),'ZZ'  ,'nsi_fractional_step_multiply',zz)
    call memory_deallo(mem_modul(1:2,modul),'AUPX','nsi_fractional_step_multiply',Aupx)

  end subroutine nsi_fractional_step_multiply

  subroutine nsi_fractional_step_continuity_tau_bicgstab()

    integer(ip)         :: maxit,itera
    integer(ip)         :: kfl_resid
    integer(ip)         :: kfl_recov
    real(rp)            :: toler,resid
    real(rp)            :: rho,newrho,invnb
    real(rp),   pointer :: bb(:)            ! RHS
    real(rp),   pointer :: xx(:)            ! Solution
    real(rp),   pointer :: rr(:)            ! Preconditioned residual
    real(rp),   pointer :: qq(:)            !
    real(rp),   pointer :: vv(:)            !
    real(rp),   pointer :: tt(:)            !
    real(rp),   pointer :: ss(:)            !
    real(rp),   pointer :: p2(:)            !
    real(rp),   pointer :: r0(:)            ! First direction
    real(rp),   pointer :: ww(:)            ! Working array
    real(rp)            :: pnorm            ! || p^n ||
    real(rp)            :: bnorm            ! || b ||
    real(rp)            :: alpha
    real(rp)            :: beta
    real(rp)            :: omega
    real(rp)            :: raux
    !
    ! Solver initialization
    !
    call nsi_solini(2_ip)                   ! Initialize solver
    kfl_recov = solve_sol(1) % kfl_recov    ! Save solver's default
    solve_sol(1) % kfl_recov = 2            ! Solver RHS is already exchanged
    !
    ! Decide if convergence criteron is real residual (=1) or preconditioned residual (=0=
    !
    kfl_resid = 1

    nullify(bb)
    nullify(xx)
    nullify(rr)
    nullify(qq)
    nullify(vv)
    nullify(tt)
    nullify(ss)
    nullify(p2)
    nullify(r0)
    nullify(ww)
    call memory_alloca(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'RR','nsi_fractional_step_solution_tau',rr,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'QQ','nsi_fractional_step_solution_tau',qq,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'VV','nsi_fractional_step_solution_tau',vv,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'TT','nsi_fractional_step_solution_tau',tt,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'SS','nsi_fractional_step_solution_tau',ss,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'P2','nsi_fractional_step_solution_tau',p2,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'R0','nsi_fractional_step_solution_tau',r0,max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'WW','nsi_fractional_step_solution_tau',ww,max(1_ip,npoin))
    !
    ! Norm of pressure || p ||
    !
    call norm2x(1_ip,pp,pnorm)
    !
    ! b = Tau^{-1} rc' + ( L - S ) (gamma-1) p^n + L p^n
    !
    if( INOTMASTER ) then
       call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,bb)
       bb = ( bp - bb ) / Tau
       call rhsmod(1_ip,bb)
       xx = pp
    end if
    !
    ! Real residual:           INVNB = 1 / || b ||
    ! Preconditioned residual: INVNB = 1 / || M^-1 b ||
    !
    if( kfl_resid == 1 ) then
       call norm2x(1_ip,bb,bnorm)
    else
       call solver(bb,qq,L,pmatr)
       call norm2x(1_ip,qq,bnorm)
    end if
    if( bnorm < 1.0e-12_rp ) goto 1
    invnb = 1.0_rp / bnorm
    !
    ! rr = initial real residual
    !
    if( INOTMASTER ) then
       call nsi_fractional_step_multiply(xx,rr)
       rr = bb - rr
    end if
    !
    ! Initial preconditioned residual r0, with M r0 = r
    !
    call solver(rr,r0,L,pmatr)
    call prodxy(1_ip,npoin,r0,r0,newrho)
    !
    ! Initial real residual norm           RESID = || r || / || b ||
    ! Initial preconditioned residual norm RESID = || M^-1 r || / || M^-1 b ||
    !
    if( kfl_resid == 1 ) then
       call norm2x(1_ip,rr,resid)
       resid = resid * invnb
    else
       resid = sqrt(newrho) * invnb
    end if
    !
    ! Initialization
    !
    rr    = r0
    alpha = 0.0_rp
    omega = 1.0_rp
    rho   = 1.0_rp
    p2    = 0.0_rp
    qq    = 0.0_rp
    !
    ! Tolerance and max number of iterations
    !
    if( kfl_adres_nsi == 1 ) then
       toler = max(toler_nsi,adres_nsi * resid)
       maxit = mitri_nsi
    else
       toler = toler_nsi
       maxit = mitri_nsi
    end if
    itera = 0
    !
    ! Solve A Dp' = b
    !
    do while( itera < maxit .and. resid > toler )

       itera = itera + 1
       !
       ! beta = (rho^k/rho^{k-1})*(alpha/omega)
       !
       if( abs(rho) == 0.0_rp .or. abs(omega) == 0.0_rp ) then
          !print*,'RHO,OMEGA=',rho,omega
          goto 1
       else
          beta = (newrho/rho) * (alpha/omega)
       end if
       !
       ! p^{k+1} = r^{k} + beta*(p^k - omega*q^k)
       !
       p2(1:npoin) = rr(1:npoin) + beta * (p2(1:npoin) - omega * qq(1:npoin))
       !
       ! L q^{k+1} = A ( R^{-1} p^{k+1} )
       !
       call nsi_fractional_step_multiply(p2,vv)
       call solver(vv,qq,L,pmatr)
       !
       ! alpha = rho^k / <r0,q^{k+1}>
       !
       call prodxy(1_ip,npoin,r0,qq,raux)

       if( raux == 0.0_rp ) goto 1
       alpha = newrho / raux
       if( abs(alpha) <= zeror ) goto 1
       !
       ! s = r^k - alpha*q^{k+1}
       !
       ss(1:npoin) = rr(1:npoin) - alpha * qq(1:npoin)
       !
       ! L t = A ( R^{-1} s )
       !
       call nsi_fractional_step_multiply(ss,vv)
       call solver(vv,tt,L,pmatr)
       !
       ! omega = <t,s> / <t,t>
       !
       call prodts(1_ip,npoin,tt,ss,raux,omega)
       if( abs(raux) == 0.0_rp ) then
          goto 1
       end if
       omega = omega / raux
       !
       ! x^{k+1} = x^k + alpha*p^{k+1} + omega*s
       ! r^{k+1} = s - omega*t
       !
       xx(1:npoin) = xx(1:npoin) + alpha * p2(1:npoin) + omega * ss(1:npoin)
       rr(1:npoin) = ss(1:npoin) - omega * tt(1:npoin)
       !
       ! NEWRHO = rho^k = <r0,r^k> and RESID = || r^{k+1} ||
       !
       rho = newrho

       if( kfl_resid == 1 ) then
          call prodxy(1_ip,npoin,rr,r0,newrho)
          call nsi_fractional_step_multiply(xx,ww)
          ww(1:npoin) = bb(1:npoin) - ww(1:npoin)
          call norm2x(1_ip,ww,resid)
          resid = resid * invnb
       else
          call prodts(1_ip,npoin,rr,r0,resid,newrho)
          resid = sqrt(resid) * invnb
       end if

       write(92,*) itera,resid,bnorm  ; flush(92)

    end do

1   continue
    call livinf(165_ip,'S('//trim(intost(itera))//')',0_ip)
    !
    ! Actualize p^{n+1}
    ! Dp      = Dp' + (1-gamma)*p^n
    ! p^{n+1} = Dp' + p^n
    !
99  continue
    if( INOTMASTER ) then
       Dpp(1:npoin) = xx(1:npoin) - pp(1:npoin)
       Dp(1:npoin)  = xx(1:npoin)
       pp(1:npoin)  = xx(1:npoin)
    end if
    !
    ! Recover solver's default
    !
    solve_sol(1) % kfl_recov = kfl_recov

    call memory_deallo(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx)
    call memory_deallo(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb)
    call memory_deallo(mem_modul(1:2,modul),'RR','nsi_fractional_step_solution_tau',rr)
    call memory_deallo(mem_modul(1:2,modul),'QQ','nsi_fractional_step_solution_tau',qq)
    call memory_deallo(mem_modul(1:2,modul),'VV','nsi_fractional_step_solution_tau',vv)
    call memory_deallo(mem_modul(1:2,modul),'TT','nsi_fractional_step_solution_tau',tt)
    call memory_deallo(mem_modul(1:2,modul),'SS','nsi_fractional_step_solution_tau',ss)
    call memory_deallo(mem_modul(1:2,modul),'P2','nsi_fractional_step_solution_tau',p2)
    call memory_deallo(mem_modul(1:2,modul),'R0','nsi_fractional_step_solution_tau',r0)
    call memory_deallo(mem_modul(1:2,modul),'WW','nsi_fractional_step_solution_tau',ww)

  end subroutine nsi_fractional_step_continuity_tau_bicgstab

  subroutine nsi_fractional_step_continuity_tau_cg()

    integer(ip)         :: maxit,itera
    real(rp)            :: toler,resid
    real(rp)            :: invnb
    real(rp)            :: rho,newrho
    real(rp),   pointer :: bb(:)            ! RHS
    real(rp),   pointer :: xx(:)            ! Solution
    real(rp),   pointer :: zz(:)            !
    real(rp),   pointer :: rr(:)            !
    real(rp),   pointer :: qq(:)            !
    real(rp),   pointer :: vv(:)            !
    real(rp)            :: bnorm            ! || b ||
    real(rp)            :: alpha            ! Orthomin relaxation factor
    real(rp)            :: beta             ! Orthomin relaxation factor
    integer(ip)         :: imeth
    !
    ! Decide if we sovle for Dp' (=0) or p (=1)
    !
    imeth = 1

    call nsi_solini(2_ip)

    nullify(bb)
    nullify(xx)
    nullify(zz)
    nullify(rr)
    nullify(qq)
    nullify(vv)
    call memory_alloca(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb,npoin)
    call memory_alloca(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx,npoin)
    call memory_alloca(mem_modul(1:2,modul),'ZZ','nsi_fractional_step_solution_tau',zz,npoin)
    call memory_alloca(mem_modul(1:2,modul),'RR','nsi_fractional_step_solution_tau',rr,npoin)
    call memory_alloca(mem_modul(1:2,modul),'QQ','nsi_fractional_step_solution_tau',qq,npoin)
    call memory_alloca(mem_modul(1:2,modul),'VV','nsi_fractional_step_solution_tau',vv,npoin)
    !
    ! b = Tau^{-1} rc' + ( L - S ) (gamma-1) p^n
    !
    if( INOTMASTER ) then
       if( imeth == 0 ) then
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,rc)
          rc = ( bp - rc ) / Tau
          call nsi_fractional_step_multiply(pp,bb)
          bb = rc - bb
          xx = Dpp
       else if( imeth == 1 ) then
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,bb)
          bb = ( bp - bb ) / Tau
          xx = pp
       end if
    end if
    !
    ! INVNB = 1 / || M^-1 b ||
    !
    call solver(bb,qq,L,pmatr)
    call prodxy(1_ip,npoin,bb,qq,bnorm)
    bnorm = sqrt(bnorm)
    if( bnorm < 1.0e-12_rp ) goto 1
    invnb = 1.0_rp / bnorm
    !
    ! r0 = initial residual
    !
    if( INOTMASTER ) then
       call nsi_fractional_step_multiply(xx,rr,'TRANSPOSE')
       rr = bb - rr
    end if
    !
    ! M z0 = r0
    ! p0 = z0
    !
    call solver(rr,zz,L,pmatr)
    call prodxy(1_ip,npoin,rr,zz,rho)
    resid = sqrt(rho) * invnb
    qq = zz
    !
    ! Tolerance and max number of iterations
    !
    if( kfl_adres_nsi == 1 ) then
       toler = max(toler_nsi,adres_nsi * resid)
    else
       toler = toler_nsi
    end if
    maxit = mitri_nsi
    itera = 0
    !
    ! Solve A Dp' = b
    !
    do while( itera < maxit .and. resid > toler )

       itera = itera + 1
       !
       ! v = A p^k
       !
       call nsi_fractional_step_multiply(qq,vv,'TRANSPOSE')
       !
       ! alpha = (r^k,z^k) / (p^k,v)
       !
       call prodxy(1_ip,npoin,qq,vv,alpha)
       alpha = rho / alpha
       !
       ! x^{k+1} = x^k + alpha * p^k
       ! r^{k+1} = r^k - alpha * A p^k
       !
       if( INOTMASTER ) then
          xx(1:npoin) = xx(1:npoin) + alpha * qq(1:npoin)
          rr(1:npoin) = rr(1:npoin) - alpha * vv(1:npoin)
       end if
       !
       ! M z^{k+1} = r^{k+1}
       !
       call solver(rr,zz,L,pmatr)
       !
       ! beta = (z^{k+1},r^{k+1}) / (z^{k},r^{k})
       !
       call prodxy(1_ip,npoin,rr,zz,newrho)
       beta = newrho / rho
       rho  = newrho
       !
       ! q^{k+1} = z^{k+1} + beta * p^k
       !
       if( INOTMASTER ) then
          qq(1:npoin) = zz(1:npoin) + beta * qq(1:npoin)
       end if
       !
       ! Residual = || r || / || b ||
       !
       call norm2x(1_ip,rr,resid)
       resid = resid * invnb

!!!!write(92,*) itera,resid,bnorm,toler  ; flush(92)

    end do
1   continue
    call livinf(165_ip,'S('//trim(intost(itera))//')',0_ip)
    !
    ! Actualize p^{n+1}
    ! Dp      = Dp' + (1-gamma)*p^n
    ! p^{n+1} = Dp' + p^n
    !
    if( INOTMASTER ) then
       if( imeth == 1 ) then
          Dpp(1:npoin) = xx(1:npoin) - pp(1:npoin)
          Dp(1:npoin)  = xx(1:npoin) - gamma_nsi * pp(1:npoin)
          pp(1:npoin)  = xx(1:npoin)
       else
          Dpp(1:npoin) = xx(1:npoin)
          Dp(1:npoin)  = xx(1:npoin) + gamma1 * pp(1:npoin)
          pp(1:npoin)  = Dpp(1:npoin) + pp(1:npoin)
       end if
    end if

    call memory_deallo(mem_modul(1:2,modul),'XX','nsi_fractional_step_solution_tau',xx)
    call memory_deallo(mem_modul(1:2,modul),'BB','nsi_fractional_step_solution_tau',bb)
    call memory_deallo(mem_modul(1:2,modul),'ZZ','nsi_fractional_step_solution_tau',zz)
    call memory_deallo(mem_modul(1:2,modul),'RR','nsi_fractional_step_solution_tau',rr)
    call memory_deallo(mem_modul(1:2,modul),'QQ','nsi_fractional_step_solution_tau',qq)
    call memory_deallo(mem_modul(1:2,modul),'VV','nsi_fractional_step_solution_tau',vv)

  end subroutine nsi_fractional_step_continuity_tau_cg

  subroutine nsi_fractional_step_transpose(nbnodes,ndof1,ndof2,ia,ja,aa_in,aa_out)

    integer(ip), intent(in)              :: nbnodes
    integer(ip), intent(in)              :: ndof1
    integer(ip), intent(in)              :: ndof2
    integer(ip), intent(in)              :: ia(*)
    integer(ip), intent(in)              :: ja(*)
    real(rp)                             :: aa_in(ndof1,ndof2,*)
    real(rp),    intent(inout), optional :: aa_out(ndof1,ndof2,*)
    integer(ip)                          :: ii,jj,iz,jz

    do ii = 1,nbnodes
       do iz = ia(ii),ia(ii+1)-1
          jj = ja(iz)
          loop_jz: do jz = ia(jj),ia(jj+1)-1
             if( ja(jz) == ii ) then
                if( maxval(abs(aa_out(:,:,jz)+aa_in(:,:,iz))) > 1.0e-12_rp ) then
                   print*,'marde=',maxval(abs(aa_out(:,:,jz)+aa_in(:,:,iz)))
                end if
                exit loop_jz
             end if
          end do loop_jz
       end do
    end do

  end subroutine nsi_fractional_step_transpose

end module mod_nsi_fractional_step
