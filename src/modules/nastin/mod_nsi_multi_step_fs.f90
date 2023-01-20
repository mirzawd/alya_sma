!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_multi_step_fs.f90
!> @author  Guillaume Houzeaux
!> @date    20/09/2016
!> @brief   Fractional step
!> @details Fractional step
!> @} 
!-----------------------------------------------------------------------
module mod_nsi_multi_step_fs

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use def_kermod,                      only : kfl_noslw_ker
  use def_kermod,                      only : gasco
  use mod_ker_proper,                  only : ker_proper
  use mod_memory,                      only : memory_alloca
  use mod_memory,                      only : memory_deallo
  use mod_memory,                      only : memory_size
  use mod_nsi_schur_operations
  use mod_matrices,                    only : matrices_all
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_matrix_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_messages,                    only : livinf
  use mod_solver,                      only : solver_solve
  use mod_solver,                      only : solver_exchange
  use mod_solver,                      only : solver_explicit
  use mod_local_basis,                 only : local_basis_global_to_local
  use mod_local_basis,                 only : local_basis_local_to_global
  use mod_nsi_corio_stab,              only : nsi_corio_stab_matrix
  implicit none
  private

  real(rp)               :: gamma1

  real(rp),   pointer    :: vv(:)       ! auxiliar vector
  real(rp),   pointer    :: xx(:)       ! auxiliar vector
  real(rp),   pointer    :: aux(:)      ! auxiliar vector
  real(rp),   pointer    :: uu0(:)      ! auxiliar vector
  real(rp),   pointer    :: rm(:,:)     ! residual of momentum at time 
  real(rp),   pointer    :: rc(:)       ! residual of continuity 
  real(rp),   pointer    :: divup(:)    ! div of up 
  real(rp),   pointer    :: rt(:)       ! Total residual of momentum 
  real(rp),   pointer    :: rt_tmp(:)   ! Total residual of momentum
  real(rp),   pointer    :: rf_tmp(:)   ! idem rt_tmp but just gravity and Boussinesq components 
!  real(rp),   pointer    :: rhsid_gravb(:)   ! rhsid only gravity and Boussinesq components -- lo mande a def_nastin
  real(rp),   pointer    :: proj_gravb(:)    ! projection of gravity and Boussinesq components
  real(rp),   pointer    :: divprgravb(:)    ! div of projection of gravity and  Boussinesq

  
  real(rp),   pointer    :: bu(:)         
  real(rp),   pointer    :: bp(:)
  real(rp),   pointer    :: uu(:)
  real(rp),   pointer    :: pp(:)
  
  real(rp),   contiguous, pointer    :: Q(:)
  real(rp),   pointer    :: Dp(:)
  real(rp),   pointer    :: Aupp(:)
  
  real(rp),   pointer    :: Tim(:)      ! Nodal projection of dt/rho
  real(rp),   pointer    :: Mass(:)     ! Mass matrix
  real(rp),   pointer    :: Grad(:,:)   ! Gradient matrix   G_ij=-\int (div Ni) Nj
  real(rp),   pointer    :: Div(:,:)    ! Divergence matrix D_ij= \int Ni (div Nj) => D=-G^t
  real(rp),   contiguous, pointer    :: Lapl(:)     ! Laplacian matrix

  real(rp)               :: c3 
  real(rp)               :: a(5)
  real(rp)               :: b(5,5)
  real(rp)               :: sigma(5)
  integer(ip), parameter :: conservative = 0
  integer(ip), parameter :: n_moin_steps = 4   ! Number of steps as Moin would do it. That is, without the enhanced extrapolation proposed by Capuano for the case kfl_fscon_nsi==0 

  real(rp),   pointer    :: Nua(:)      ! Nodal projection of mu/rho

  real(rp)               :: rho(1)
  real(rp)               :: dt_inv0
  real(rp)               :: dt_inv00
  real(rp)               :: dt_inv000
  integer(ip)            :: time_iter
  integer(ip)            :: dummi
 
  type(soltyp), pointer  :: solve_mome(:)
  type(soltyp), pointer  :: solve_mass(:)
  
  public :: nsi_multi_step_fs_solution
  public :: nsi_multi_step_fs_initialization
  public :: nsi_multi_step_fs_memory
  public :: nsi_multi_step_fs_matrices

  public :: Grad
  
contains 

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Initializaation 
  !> @details Initializaation
  !
  !-----------------------------------------------------------------------

  subroutine nsi_multi_step_fs_initialization()
    
    nullify(vv)
    nullify(xx)
    nullify(aux)
    nullify(uu0)
    nullify(rm)
    nullify(rc)
    nullify(divup)
    nullify(rt)
    nullify(rt_tmp)
    if(kfl_fsgrb_nsi==1) then
       nullify(rf_tmp)
       nullify(rhsid_gravb)
       nullify(proj_gravb)
       nullify(divprgravb)
    end if

    nullify(bu)
    nullify(bp)
    nullify(uu)
    nullify(pp)

    nullify(Q)
    nullify(Dp)
    nullify(Aupp)
        
    nullify(Tim)
    nullify(Mass)
    nullify(Grad)
    nullify(Div)
    nullify(Lapl)
    nullify(Scorio_nsi)

    nullify(Nua)    

    time_iter = 1

  end subroutine nsi_multi_step_fs_initialization
  
  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Compute some matrices
  !> @details Compute some matrices
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_matrices()
    
    if( kfl_grad_div_nsi /= 0 ) then
       solve_sol => solve(2:)
       call matrices_all                (Lapl,Grad,Div,DEALLOCATE_MATRIX=.true.)
       call nsi_dirichlet_matrix_rhs    (Q=Lapl,ROTATION=.false.,DIRICHLET=.true.)
       call nsi_dirichlet_matrix_rhs    (Aup=Grad,Apu=Div,ROTATION=.true.,DIRICHLET=.false.)
       Q => Lapl
    else
       Q => lapla_nsi
    end if
    
    bu     => rhsid
    bp     => rhsid(ndbgs_nsi+1:)
    uu     => unkno
    pp     => unkno(ndbgs_nsi+1:)   
    Tim    => dt_rho_nsi
    Mass   => vmass
    !
    ! Coriolis stabilisation
    !
    if( kfl_stab_corio_nsi > 0_ip ) call nsi_corio_stab_matrix(Scorio_nsi,DEALLOCATE_MATRICES=.true.)
    
  end subroutine nsi_multi_step_fs_matrices

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate memory
  !> @details Allocate memory and define module constants
  !
  !-----------------------------------------------------------------------

  subroutine nsi_multi_step_fs_memory()

    integer(ip) :: ipoin
    
    gamma1 =  1.0_rp - gamma_nsi

    call memory_alloca(mem_modul(1:2,modul),'VV'    ,'mod_nsi_multi_step_fs',vv,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'XX'    ,'mod_nsi_multi_step_fs',xx,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'UU0'   ,'mod_nsi_multi_step_fs',uu0,   max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'RM'    ,'mod_nsi_multi_step_fs',rm,    max(1_ip,ndime*npoin),5_ip) 
    call memory_alloca(mem_modul(1:2,modul),'RT'    ,'mod_nsi_multi_step_fs',rt,    max(1_ip,ndime*npoin)) 
    call memory_alloca(mem_modul(1:2,modul),'RTTMP' ,'mod_nsi_multi_step_fs',rt_tmp,max(1_ip,ndime*npoin))
    if(kfl_fsgrb_nsi==1) then
       call memory_alloca(mem_modul(1:2,modul),'RHSID_GRAVB' ,'mod_nsi_multi_step_fs',rhsid_gravb,max(1_ip,ndime*npoin))
       call memory_alloca(mem_modul(1:2,modul),'PROJ_GRAVB'  ,'mod_nsi_multi_step_fs',proj_gravb, max(1_ip,ndime*npoin))
       call memory_alloca(mem_modul(1:2,modul),'DIVPRGRAVB'  ,'mod_nsi_multi_step_fs',divprgravb, max(1_ip,npoin))
       call memory_alloca(mem_modul(1:2,modul),'RF_TMP'      ,'mod_nsi_multi_step_fs',rf_tmp,     max(1_ip,ndime*npoin))
    end if

    call memory_alloca(mem_modul(1:2,modul),'RC'    ,'mod_nsi_multi_step_fs',rc,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUX'   ,'mod_nsi_multi_step_fs',aux,   max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DIVUP' ,'mod_nsi_multi_step_fs',divup, max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DP'    ,'mod_nsi_multi_step_fs',Dp,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUPP'  ,'mod_nsi_multi_step_fs',Aupp,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'NUA'   ,'mod_nsi_multi_step_fs',Nua,   max(1_ip,npoin))

    solve_mome => solve(NSI_SOLVER_MOMENTUM:)
    solve_mass => solve(NSI_SOLVER_CONSISTENT_MASS:)
    
    do ipoin = 1,npoin
       Nua(ipoin) = 1.0_rp
    end do
    
    a     = 0.0_rp
    b     = 0.0_rp
    sigma = 0.0_rp

    if(kfl_tiacc_nsi == 2) then !Heun’s method or RK2 just for validation

       a(4)     = 1.0_rp
       
       b(4,3)   = 1.0_rp

       sigma(3) = 0.5_rp
       sigma(4) = 0.5_rp

    else if(kfl_tiacc_nsi == 4) then ! standard 4 order RK

       a(2)     = 0.5_rp
       a(3)     = 0.5_rp
       a(4)     = 1.0_rp

       b(2,1)   = 0.5_rp 
       b(3,2)   = 0.5_rp
       b(4,3)   = 1.0_rp

       sigma(1) = 1.0_rp/6.0_rp 
       sigma(2) = 1.0_rp/3.0_rp 
       sigma(3) = 1.0_rp/3.0_rp 
       sigma(4) = 1.0_rp/6.0_rp 

    else ! energy preserving 3 order RK 

       if(conservative == 1) then 

          c3       = 1.0_rp/4.0_rp

          a(2)     = (c3 - 1.0_rp)/(4.0_rp*c3 - 3.0_rp)
          a(3)     = c3
          a(4)     = 1.0_rp

          b(2,1)   = (c3-1.0_rp)/(4.0_rp*c3-3.0_rp) 
          b(3,1)   = c3 - ((2_rp*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))/(2.0_rp*(c3-1.0_rp)) 
          b(4,1)   = -1.0_rp * ((2.0_rp*c3-1.0_rp)**2/(2.0_rp*(c3-1.0_rp)*(4.0_rp*c3-3.0_rp)))  

          b(3,2)   = ((2_rp*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))/(2.0_rp*(c3-1.0_rp)) 
          b(4,2)   = (6.0_rp*c3**2-8.0_rp*c3+3.0_rp)/(2_rp*(c3-1.0_rp)*(2_rp*c3-1.0_rp))

          b(4,3)   = (c3-1.0_rp)/((2*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))

          sigma(1) = -(1.0_rp)/(12.0_rp*(c3-1.0_rp)) 
          sigma(2) = ((4.0_rp*c3-3.0_rp)**2)/(12.0_rp*(c3-1.0_rp)*(2.0_rp*c3-1.0_rp)) 
          sigma(3) = (-1.0_rp)/(12.0_rp*(c3-1.0_rp)*(2.0_rp*c3-1.0_rp))
          sigma(4) = (4.0_rp*c3-3.0_rp)/(12.0_rp*(c3-1.0_rp))

       else

          !  a(3) = 0.5_rp
          !  a(4) = 1.0_rp  

          !  b(3,2) = 0.5_rp
          !  b(4,2) = -1.0_rp
          !  b(4,3) = 2.0_rp

          !  sigma(2) = 1.0_rp/6.0_rp 
          !  sigma(3) = 2.0_rp/3.0_rp
          !  sigma(4) = 1.0_rp/6.0_rp

          a(3)     = 1.0_rp
          a(4)     = 0.5_rp  

          b(3,2)   = 1_rp
          b(4,2)   = 1.0_rp/4.0_rp
          b(4,3)   = 1.0_rp/4.0_rp

          sigma(2) = 1.0_rp/6.0_rp 
          sigma(3) = 1.0_rp/6.0_rp
          sigma(4) = 2.0_rp/3.0_rp

          ! a(3) = 8.0_rp/15.0_rp
          ! a(4) = 2.0_rp/3.0_rp  

          ! b(3,2) = 8.0_rp/15.0_rp
          ! b(4,2) = 1.0_rp/4.0_rp
          ! b(4,3) = 5.0_rp/12.0_rp

          ! sigma(2) = 1.0_rp/4.0_rp 
          ! sigma(3) = 0.0_rp
          ! sigma(4) = 3.0_rp/4.0_rp

       end if

    end if

  end subroutine nsi_multi_step_fs_memory

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl and houzeaux
  !> @date    2020-03-21
  !> @brief   Multi-step FS
  !> @details Multistep fractional step method from:
  !>          Explicit Runge–Kutta schemes for incompressible
  !>          flow with improved energy-conservation properties
  !>          F.Capuano, G.Coppola, L.Rández, L.deLuca, Journal of
  !>          ComputationalPhysics 328: 86–94 (2017)
  !>
  !>          Schemme 3p5q(4) from the paper 3o time, 5o energy and 4 steps
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution()

    integer(ip) :: idofn,ipoin,idime
       
     !
     ! Save alues at old step
     ! 
     do ipoin = 1,npoin
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime + idime
           uu0(idofn)   = uu(idofn)
        end do
     end do
     
     if(kfl_tiacc_nsi == 2) then
        !
        ! Order 2
        !
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
        
     else if(kfl_tiacc_nsi == 3) then  
        !
        ! Order 3
        !
        if (conservative == 1) then 
           call nsi_multi_step_fs_eval(1_ip)
           call nsi_multi_step_fs_eval(2_ip)
           call nsi_multi_step_fs_eval(3_ip)
           call nsi_multi_step_fs_eval(4_ip)
        else 
           call nsi_multi_step_fs_eval(2_ip)
           call nsi_multi_step_fs_eval(3_ip)
           call nsi_multi_step_fs_eval(4_ip)           
        end if
        
     else 
        !
        ! Order 4
        !
        call nsi_multi_step_fs_eval(1_ip)
        call nsi_multi_step_fs_eval(2_ip)
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
        
     end if

     dt_inv000 = dt_inv00 
     dt_inv00  = dt_inv0 
     dt_inv0   = 1.0_rp/dtinv_nsi
     time_iter = time_iter + 1

  end subroutine nsi_multi_step_fs_solution

  !-----------------------------------------------------------------------
  !> 
  !> @author  lemhkuhl and houzeaux
  !> @date    2020-03-20
  !> @brief   Assemble equations
  !> @details Assemble equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_eval(istep)

    integer(ip) , intent(in) :: istep
    integer(ip)              :: ipoin
    real(rp)                 :: dummr
    !
    ! Obtain projection for Coriolis stabilization
    !
    call times(5) % ini()
    if( kfl_stab_corio_nsi > 0_ip ) then
       call nsi_apuvec(1_ip,Scorio_nsi,uu,prdivcor_nsi)                       ! prdivcor =  Scorio * u
       call rhsmod(1_ip,prdivcor_nsi)
       !
       !       call solver_explicit(solve_mome,prdivcor_nsi,EXCHANGE=.false.)    estoy aca es increible como se ha complicado hacer una pryeccion!!! !,solve_consistent=solve_mass,x=proj_gravb)   ! Obtain proj_grav
       !       no si si esto sirve. como distingue si es un vector de ndime*npoin o npoin
       !
       !  mejor lo hago a la antigua -- borrowed from nsi_solsgs  -- I understand no BC is needed     --- ojo para proj_gravb si use lo de solver 
       do ipoin = 1,npoin
          dummr                    = 1.0_rp / vmass(ipoin)
          prdivcor_nsi(ipoin)      = prdivcor_nsi(ipoin) * dummr 
       end do    
    end if

    !-----------------------------------------------------------------
    !
    ! Go back to local system
    !
    !-----------------------------------------------------------------
    
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local
    call times(5) % add()
    
    !-----------------------------------------------------------------
    !
    ! Assemble equations
    !
    !-----------------------------------------------------------------

    call nsi_matrix()
    if(kfl_assem_nsi >39_ip .and. kfl_paral/=0_ip) then   ! for master de debug version marked conflict I guess in master
       dt_rho_nsi(1:npoin) = 1.0_rp/(dtinv_nsi *  densi_aux)  ! por ahoar así de trucho 
       mass_rho_nsi(:,1)   =  densi_aux * vmass
    end if
    
    !-----------------------------------------------------------------
    !
    ! Compute density
    !
    !-----------------------------------------------------------------

    if( kfl_regim_nsi == 3 .and. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 ) then
       !
       ! Low Mach
       !
       do ipoin = 1,npoin
          if (kfl_lookg_nsi > 0) then
             call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
             Nua(ipoin)  = rho(1)
          else
             Nua(ipoin)  = (prthe(1)/gasco) * (wmean(ipoin,1)/tempe(ipoin,1))
          endif
       end do
       !Nua = 0.5*(3.0_rp*Nu-Nu0)

    else if( kfl_regim_nsi == 3 .and. (time_iter /= 1) .and. kfl_coupl(ID_NASTIN,ID_CHEMIC) == 0 ) then
       
       do ipoin = 1,npoin
          Nua(ipoin)  = (prthe(1)/gasco) * (1.0_rp/tempe(ipoin,1))
       end do

    else if( kfl_surte_nsi == 2_ip ) then

       Nua = 1.0_rp

    else
       !
       ! Incompressible
       !
       do ipoin = 1,npoin
          Nua(ipoin) = mass_rho_nsi(ipoin,1)
       end do
       call solver_lumped_mass_system(1_ip,Nua,EXCHANGE=.false.)
    
    end if

    call livinf(56_ip,' ',modul)    
    call livinf(160_ip,' ',1_ip)    

    if(istep == 4_ip) then
       call nsi_multi_step_fs_solution_final()
    else
       call nsi_multi_step_fs_solution_sj(istep)
    endif

    !-----------------------------------------------------------------
    !
    ! Go back to global system and actualize velocity and pressure,
    ! required for next assembly
    !
    !-----------------------------------------------------------------

    call livinf(165_ip,'C',0_ip)

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)
 
    call livinf(164_ip,' ',1_ip) 
 
  end subroutine nsi_multi_step_fs_eval  

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Momentum
  !> @details Intermediate solution of the momentum equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution_sj(istep)

    integer(ip) , intent(in) :: istep
    integer(ip)              :: idofn,ipoin,idime,jstep
    real(rp)                 :: aux2

    call times(5) % ini()
    !
    ! Momentum residual rm^k
    !
    if( INOTEMPTY ) then
       if( kfl_fscon_nsi == 0 ) then 
          if( ittim > n_moin_steps ) then 
             do ipoin = 1,npoin
                do idime = 1,ndime
                   idofn = (ipoin-1)*ndime + idime
                   rm(idofn,istep) = bu(idofn)
                end do
                aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*a(istep+1)*(press(ipoin,3)-press(ipoin,4)) !Capuano et al
             end do
          else
             do ipoin = 1,npoin
                do idime = 1,ndime
                   idofn = (ipoin-1)*ndime + idime
                   rm(idofn,istep) = bu(idofn)
                end do
                aux(ipoin) = pp(ipoin)  ! pressure solved from Laplacian
             end do
          endif
       else
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                rm(idofn,istep) = bu(idofn)
                if (kfl_fsgrb_nsi == 1_ip ) rf_tmp(idofn) = rhsid_gravb(idofn)
             end do
             aux(ipoin) = 0.0_rp
          end do
       end if
       !
       ! vv = GRAD*press
       !
       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,aux,vv)           
       else
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),aux,vv) 
       end if

       call rhsmod(ndime,rm(:,istep))
       if(kfl_fsgrb_nsi==1) call rhsmod(ndime,rf_tmp)
       call rhsmod(ndime,vv)
       
       if( kfl_regim_nsi==3 .and.  kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
          !
          ! Temporal term for skew-symmetric
          !
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn           = (ipoin-1)*ndime + idime
                aux2            = 0.5_rp * drhodt_nsi(ipoin) * uu(idofn)
                rm(idofn,istep) = rm(idofn,istep) - aux2
             end do
          end do
       endif

       if (kfl_surte_nsi == 2_ip) then
          do ipoin = 1,npoin
             call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                aux2 = 0.0_rp
                do jstep = 1,istep
                   aux2 = aux2 + rm(idofn,jstep)*b(istep+1,jstep)
                end do
                rt_tmp(idofn) = (aux2 - a(istep+1)*vv(idofn)/rho(1))
             end do
          end do
       else 
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                aux2 = 0.0_rp
                do jstep = 1,istep
                   aux2 = aux2 + rm(idofn,jstep)*b(istep+1,jstep)
                end do
                rt_tmp(idofn) = (aux2 - a(istep+1)*vv(idofn))
             end do
          end do
       end if
    end if
    call times(5) % add()
    !
    ! Solve system
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)    
    !
    ! Update solution
    !
    call times(5) % ini()
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + (rt_tmp(idofn)/dtinv_nsi)/Nua(ipoin)       ! u = u0 + r * dt/rho
       end do
    end do
    call nsi_multi_step_fs_dirichlet()                                        ! Prescribe Drichlet bc on  u'^{k+1}
    call times(5) % add()

    call livinf(165_ip,'M',0_ip)

    if( kfl_fscon_nsi == 1 ) then
       call nsi_multi_step_fs_pressure_and_correction(a(istep+1),0_ip)
    end if

  end subroutine nsi_multi_step_fs_solution_sj

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Momentum
  !> @details Fimal solution of the momentum equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution_final()

    integer(ip) :: idofn,ipoin,idime
    real(rp)    :: aux

    call times(5) % ini()

    if( INOTEMPTY ) then
       !
       ! Momentum residual rm^k
       !
       do idofn = 1,ndime*npoin                                        ! rm = bu 
          rm(idofn,4) = bu(idofn)
          if(kfl_fsgrb_nsi==1) rf_tmp(idofn) = rhsid_gravb(idofn)
       end do
       call rhsmod(ndime,rm(:,4))
       if(kfl_fsgrb_nsi==1) call rhsmod(ndime,rf_tmp)     ! not sure why this is needed
       !
       ! Treat last step:
       !
       if( kfl_regim_nsi==3 .and.  kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
          !
          ! Temporal term for skew-symmetric
          !
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn           = (ipoin-1)*ndime + idime
                aux             = 0.5_rp * drhodt_nsi(ipoin) * uu(idofn)
                rm(idofn,4)     = rm(idofn,4) - aux
             end do
          end do
       endif

       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn         = (ipoin-1)*ndime + idime
             aux           = rm(idofn,1)*sigma(1) + rm(idofn,2)*sigma(2) + rm(idofn,3)*sigma(3) + rm(idofn,4)*sigma(4)
             rt(idofn)     = aux
             rt_tmp(idofn) = aux
          end do
       end do
    end if
    call times(5) % add()

    !
    ! Solve
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    
    call times(5) % ini()
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + (rt_tmp(idofn)/dtinv_nsi)/Nua(ipoin) ! ORIOL
       end do
    end do

    call nsi_multi_step_fs_dirichlet()                            ! Prescribe Dirichlet bc on u'^{k+1}

    call livinf(165_ip,'M',0_ip)
    call times(5) % add()

    call nsi_multi_step_fs_pressure_and_correction(1.0_rp,1_ip)

  end subroutine nsi_multi_step_fs_solution_final

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Impose Dirichlet
  !> @details Impose Dirichlet conditions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_dirichlet(xarray)

    real(rp),   intent(inout), pointer, optional :: xarray(:)
    integer(ip)                                  :: idofn,ipoin,idime

    if( present(xarray) ) then
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idime + ndime * (ipoin -1)
             if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
                vafor_nsi(idime,ipoin) = xarray(idofn) ! actualy in the multistep case we could think of averaging all the substeps
                xarray(idofn) = 0.0_rp
             end if
          end do
       end do
    else
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idime + ndime * (ipoin -1)
             if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                  uu(idofn) = bvess_nsi(idime,ipoin,1)
          end do
       end do
    end if

  end subroutine nsi_multi_step_fs_dirichlet

  subroutine nsi_multi_step_fs_dirichlet2(xarray)   ! it was getting too complicate to generalise nsi_multi_step_fs_dirichlet

    real(rp),   intent(inout), pointer           :: xarray(:)
    integer(ip)                                  :: idofn,ipoin,idime

    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn = idime + ndime * (ipoin -1)
          if( kfl_fixno_nsi(idime,ipoin) > 0 ) xarray(idofn) = 0.0_rp
       end do
    end do

  end subroutine nsi_multi_step_fs_dirichlet2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Pressure and correction
  !> @details Solve pressure equation and correct velocity
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_pressure_and_correction(fact_substep,eval_rc)

    real(rp),    intent(in) :: fact_substep
    integer(ip), intent(in) :: eval_rc
    integer(ip)             :: idofn,ipoin,idime

    !-----------------------------------------------------------------
    !
    ! Solve pressure
    !
    !-----------------------------------------------------------------

    if( INOTEMPTY ) then
       !
       ! Scale solution
       !
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = (ipoin-1)*ndime + idime
             uu(idofn) = uu(idofn)*Nua(ipoin)
          end do
       end do
       !
       ! Continuity residual rc'^{k+1}
       !        
       call times(5) % ini()
       if( kfl_grad_div_nsi /= 0 ) then
          call nsi_apuvec(1_ip,Div,uu,divup)
       else
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,divup)             ! rc' =  Apu u*
       end if
       call times(5) % add()
       
       if(kfl_fsgrb_nsi==1) then
          !
          ! Div of projection of gravity and Boussinesq 
          ! x is the initial guess - I set it idem to the rhs
          !
          !          proj_gravb = rhsid_gravb
          proj_gravb = rf_tmp
          call solver_explicit(solve_mome,rf_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=proj_gravb)   ! Obtain proj_gravb
          call times(5) % ini()
          proj_gravb = rf_tmp   ! Put the output into proj_gravb
          call nsi_multi_step_fs_dirichlet2(proj_gravb)   ! This supposes that we are working on local axes here or that no local axes are used -- CHECK!!!! -- well the previous lines for divup also suppose that
                                                          ! Probably this step is not needed
          if( kfl_grad_div_nsi /= 0 ) then
             call nsi_apuvec(1_ip,Div,proj_gravb,divprgravb)
          else
             call nsi_apuvec(1_ip,amatr(poapu_nsi:),proj_gravb,divprgravb)     
          end if
          !
          ! Residual
          ! Beware, contrary to our emac paper and Ramon's FS in JCP, Q is just gradQ*gradP - without the minus sign!!! - Thus we need to change sign to all of the terms 
          ! In bp I have the part for FS with gravity with logical signs  
          !
          do ipoin = 1,npoin
             rc(ipoin) = (-bp(ipoin) - divup(ipoin) + (1.0_rp/dtinv_nsi) * divprgravb(ipoin))*(dtinv_nsi/fact_substep)
          end do
          call times(5) % add()
       else     ! normal behaviour 
          call times(5) % ini()
          do ipoin = 1,npoin
             rc(ipoin) = (bp(ipoin) - divup(ipoin) )*(dtinv_nsi/fact_substep)
          end do
          call times(5) % add()
       end if

       call times(5) % ini()
       call nsi_appvec(0_ip,Q,pp,rc,-1.0_rp)                           ! rc <= rc - Q p^k
       call rhsmod(1_ip,rc)
       !
       ! Add mass source from spray
       ! mass_sink is already exchanged, thus it's after rhsmod
       !
       if (associated(mass_sink)) then
          do ipoin = 1,npoin
             rc(ipoin) = rc(ipoin) +  mass_sink(ipoin)*(dtinv_nsi/fact_substep)
          end do
       endif
       !
       ! Low mach: drho/dt from old values.
       ! DRHODT is assumed to be already assembled
       !
       if( kfl_regim_nsi == 3 ) then
          do ipoin = 1,npoin
             rc(ipoin) = rc(ipoin) - drhodt_nsi(ipoin)*(dtinv_nsi/fact_substep)
          end do
       end if

       do ipoin = 1,npoin
          Dp(ipoin) = 0.0_rp
       end do
       call times(5) % add()

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
    ! Solve: Q z = rc'^k
    !
    call nsi_solini(2_ip)                                                         ! Initialize solver
    
    call solver_solve(solve_sol,Q,rc,Dp,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)   ! Solve system Q Dp = rc^k

    call times(5) % ini()
    !
    ! Dp      = Dp' + p^k
    ! p^{k+1} = Dp  +  p^k
    !
    do ipoin = 1,npoin
       Dp(ipoin) = Dp(ipoin) + pp(ipoin)  
       pp(ipoin) = Dp(ipoin) !- fsrot_nsi*Nu(ipoin)*divup(ipoin)
    end do
    
    if(eval_rc == 1) call norm2x(1_ip,rc,resin_nsi(2))                                  ! || rc ||

    call times(5) % add()
    call livinf(165_ip,'S',0_ip)

    !-----------------------------------------------------------------
    !
    ! Correct velocity
    !
    !-----------------------------------------------------------------
    
    call times(5) % ini()
    if( INOTMASTER ) then

       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,Dp,vv)                              ! rm =  G Dpp^{k+1}
       else
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),Dp,vv)                 ! rm =  Aup Dpp^{k+1}
       end if
       call rhsmod(ndime,vv)
       !
       ! Obtain vafor without the contribution from the pressure
       !
#ifdef VAFOR_NO_PR
       rt_tmp = rt
       call nsi_multi_step_fs_dirichlet(rt_tmp)    ! This obtains vafor - done here the master will not enter but it does not matter
#endif
       if (kfl_surte_nsi == 2_ip) then
          do ipoin = 1,npoin
             call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                rt_tmp(idofn) = fact_substep*vv(idofn)/rho(1)
                rt(idofn) = rt(idofn) - vv(idofn) 
             end do
          end do
       else 
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                rt_tmp(idofn) = fact_substep*vv(idofn)
                rt(idofn) = rt(idofn) - vv(idofn) 
             end do
          end do
       end if
    end if
    call times(5) % add()
    !
    ! Solve
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    call times(5) % ini()
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = (uu(idofn) - rt_tmp(idofn)/dtinv_nsi ) /Nua(ipoin)
       end do
    end do
    call nsi_multi_step_fs_dirichlet()                            ! Prescribe bc 

!   if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
       do idofn = 1,npoin*ndime
          rt_tmp(idofn) = rt(idofn)
       end do
#ifndef VAFOR_NO_PR
       call nsi_multi_step_fs_dirichlet(rt_tmp)
#endif
       call norm2x(ndime,rt_tmp,resin_nsi(1))                        ! || rm ||
!   else
!      call norm2x(ndime,rt,resin_nsi(1))                            ! || rm ||       
!   end if
    call times(5) % add()

  end subroutine nsi_multi_step_fs_pressure_and_correction


end module mod_nsi_multi_step_fs
