!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    def_nastin.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Nastin definitions
!> @details All definitions of Nastin
!-----------------------------------------------------------------------

module def_nastin_aux

  use def_kintyp
  use def_nastin_parameters, only : mmsgs_nsi 
  use def_coupli,            only : typ_color_coupling
  !
  ! Nastin types
  !
  type nsimat
     integer(ip)           :: kfl_exist
     real(rp), pointer     :: Auu(:,:,:)
     real(rp), pointer     :: Aup(:,:)
     real(rp), pointer     :: bu(:)
  end type nsimat
  !
  ! Boundary conditions
  !
  integer(ip), pointer     :: &
       kfl_fixno_nsi(:,:),    &      ! Nodal fixity 
       kfl_fixpr_nsi(:,:),    &      ! Nodal fixity for the pressure Schur complement
       kfl_fixpp_nsi(:,:),    &      ! Nodal fixity for the pressure 
       kfl_fixbo_nsi(:),      &      ! Element boundary fixity
       kfl_fixrs_nsi(:),      &      ! Reference system for the BV
       kfl_funno_nsi(:),      &      ! Functions for node bc
       kfl_funbo_nsi(:),      &      ! Functions for node bc
       kfl_funtn_nsi(:),      &      ! Function type on nodes
       kfl_funtb_nsi(:),      &      ! Function type on nodes
       kfl_wlawf_nsi(:),      &      ! Flag to identify if a point is part of a wall law boundary
       lexlo_nsi(:,:),        &      ! List for the exchange location for wall law
       ielem_wel(:),          &
       lbpse(:)                      ! List of boundary sets passed to nodes
  real(rp),    pointer     :: &
       bvess_nsi(:,:,:),      &      ! Essential velocity bc values
       bpess_nsi(:,:,:),      &      ! Essential bc values for pressure
       bvnat_nsi(:,:,:),      &      ! Natural bc values
       skcos_nsi(:,:,:),      &      ! Cosine matrices of NS skew systems
       velel_nsi(:,:),        &      ! Velocity at the exchange location for wall law  
       massb_nsi(:),          &      ! Boundary surface
       notra_nsi(:,:),        &      ! Traction on the boundary nodes for postprocess
       avntr_nsi(:,:),        &      ! Time average traction on the boundary nodes - variational calculation 
       avgtr_nsi(:,:),        &      ! Time average traction on the boundary nodes - using velocity gradients
       shape_wel(:,:),        &      ! shape functions associated to excange location (for implicit)
       btrac_nsi(:,:),        &      ! Traction on the boundary nodes from auxiliary RANS simulation
       tracr_nsi(:,:),        &      ! Traction calculated in auxiliary RANS simulation
       tluav_nsi(:,:)                ! Average velocity for two-layer coupling
  !
  ! Dimensions
  !
  integer(ip)              :: &
       ndofn_nsi(2),          &      ! # of d.o.f. of the NSI problem
       ndof2_nsi(2),          &      ! ndofn_nsi*ndofn_nsi
       ncomp_nsi,             &      ! Number of components of the velocity (NSI)
       nprev_nsi,             &      ! Previous time step or iteration
       nunkn_nsi(2),          &      ! # of unknonws ndofn*npoin  
       nevat_nsi,             &      ! Element matrix dim.=(ndime+1)*mnode
       nzsol_nsi,             &      ! Matrix size (per d.o.f.)
       nzmat_nsi(2),          &      ! Matrix size
       nzrhs_nsi(2),          &      ! RHS size
       nzpre_nsi(2),          &      ! Preconditioner size
       lperp_nsi(8),          &      ! List of periodic prescribed pressure 
       kfl_perip_nsi,         &      ! If pressure is prescribed on periodic nodes
       kfl_dodem_nsi,         &
       mflow_nsi                     ! Max # flow rates
  !
  ! Internal variables
  !
  logical(lg)              ::   &
       NSI_MONOLITHIC,          &    ! Monolithic algorithm
       NSI_SCHUR_COMPLEMENT,    &    ! Schur complement algorithm
       NSI_FRACTIONAL_STEP,     &    ! Fractional step
       NSI_SEMI_IMPLICIT,       &    ! Semi implicit
       NSI_ASSEMBLY_CONVECTIVE, &    ! Convective term is assembled
       NSI_ASSEMBLY_VISCOUS          ! Viscous term is assembled
  integer(ip)              :: &
       ittot_nsi,             &      ! Total number of iteration
       kfl_resid_nsi,         &      ! If velocity residual is required for post.
       kfl_grvis_nsi,         &      ! If velocity gradients exist
       kfl_goite_nsi,         &      ! Keep iterating
       kfl_rmom2_nsi,         &      ! Off-diagonal part of momentum operator exists
       kfl_p1ve2_nsi,         &      ! Off-diagonal part of momentum test function exists
       ndbgs_nsi,             &      ! Number dof for BGS
       kfl_stead_nsi,         &      ! Steady-state has been reached 
       kfl_tiaor_nsi,         &      ! Original time accuracy
       kfl_sgste_nsi,         &      ! Temperature subgrid scale considered
       kfl_autom_nsi,         &      ! Automatic boundaries
       ivari_nsi,             &      ! Equation being solved (momentum and/or continuity)
       iteqn_nsi(2),          &      ! Internal iterations for momentum+continuity
       itbgs_nsi,             &      ! BGS iteration number
       nbdfp_nsi,             &      ! Number of terms in the temporal derivative
       kfl_exist_fixi7_nsi,   &      ! exists nodal fixity of type 7 
       kfl_exist_fib20_nsi,   &      ! exists boundary fixity of type 20 
       kfl_exist_fib02_nsi,   &      ! exists boundary fixity of type 2
       itsta_nsi(mmsgs_nsi)          ! Statistics sgs
  real(rp)                 :: &
       dtinv_nsi,             &      ! 1/dt , from vers 772 theta is now longer included in dtinv_nsi
       dtsgs_nsi,             &      ! 1/(theta'*dt)
       err01_nsi(2),          &      ! L1 error u
       err02_nsi(2),          &      ! L2 error u
       err0i_nsi(2),          &      ! Linf error u
       err11_nsi(2),          &      ! L1 error grad(u)
       err12_nsi(2),          &      ! L2 error grad(u)
       err1i_nsi(2),          &      ! Linf error grad(u)
       corio_nsi,             &      ! Coriolis force
       pabdf_nsi(10),         &      ! BDF parameters, actually now (vers 772) we will extend it also for CN 
       rgsve_nsi,             &      ! residual BGS velocity
       rgspr_nsi,             &      ! residual BGS pressure
       resin_nsi(2),          &      ! Algebraic inner residual
       resou_nsi(2),          &      ! Algebraic outer residual
       resss_nsi(2),          &      ! Algebraic steady state residual
       reinf_nsi(2),          &      ! Algebraic Linf residual 
       tamin_nsi,             &      ! Min tau
       tamax_nsi,             &      ! Max tau
       vemin_nsi,             &      ! Min velocity
       vemax_nsi,             &      ! Max velocity
       prmin_nsi,             &      ! Min pressure
       prmax_nsi,             &      ! Max pressure
       pcoef_nsi,             &      ! Pressure coefficient 1-R/Cp
       relpa_nsi(2),          &      ! Relaxation parameter
       cpu_ass_sol_nsi(4),    &      ! CPU time assembly and solver at each iteration
       cputi_assembly_nsi(10),&      ! COU time for element assembly
       gamth_nsi,             &      ! gamma=Cp/(Cp-R)
       xmass_nsi,             &      ! Low-Mach: Mass computed from state equation
       actav_nsi,             &      ! Accumulated time for averaging
       difve_nsi,             &      ! Velocity residual w/r reference solution
       difpr_nsi,             &      ! Pressure residual w/r reference solution
       vinvt_nsi(4),          &      ! for lowmac, = integ(1/T)
       hydro_nsi,             &      ! Hydrostatic z-plane
       dtmax_nsi,             &      ! for local time step stores the maximum time step
       rmsgs_nsi,             &      ! Maximum subgrid scale residual
       resgs_nsi(2),          &      ! Subgrid scale residual (numerator/denominator)
       resis_nsi(2,mmsgs_nsi),&      ! Subgrid scale inner residual
       press_lev_nsi                 ! Pressure level
  real(rp),    pointer     :: &    
       veold_nsi(:,:),        &      ! Velocity for residual
       gradv_nsi(:,:),        &      ! velocity gradient (postprocess)
       unk2n_nsi(:,:),        &      ! Nastin second variables (pressure or density)
       dunkn_nsi(:,:),        &      ! Delta velocity
       dunkp_nsi(:),          &      ! Delta pressure
       avvel_nsi(:,:),        &      ! Average velocity
       avve2_nsi(:,:),        &      ! Average velocity**2
       avvxy_nsi(:,:),        &      ! Average vx*vy
       avpre_nsi(:),          &      ! Average pressure
       avpr2_nsi(:),          &      ! Average pressure**2
       avtan_nsi(:,:),        &      ! Average tangential force
       avmut_nsi(:),          &      ! Average turbulent viscosity
       avstx_nsi(:,:),        &      ! Average stress mu_t * grad(u)
       avsty_nsi(:,:),        &      ! Average stress mu_t * grad(v)
       avstz_nsi(:,:),        &      ! Average stress mu_t * grad(w)
       avmos_nsi(:,:),        &      ! Average momentum source from spray
       av_mass_flux_nsi(:,:), &      ! Average mass flux
       av_mom_flux_diag_nsi(:,:), &  ! Average momentum flux, diagonal terms rho*Vx*Vx, rho*Vy*Vy, rho*Vz*Vz
       av_mom_flux_off_nsi(:,:), &   ! Average momentum flux, off-diagonal terms rho*Vx*Vy, rho*Vy*Vz, rho*Vz*Vx
       vmaxp_nsi(:,:),        &      ! maximum nodewise velocity module during a certain widnow  
       vminp_nsi(:,:),        &      ! minimum nodewise velocity module during a certain widnow  
       vavep_nsi(:,:),        &      ! average nodewise velocity module during a certain widnow  
       pinde_nsi(:,:),        &      ! nodewise pulsatiltiy index during a certain window  
       envel_nsi(:,:),        &      ! Ensemble velocity
       enve2_nsi(:,:),        &      ! Ensemble velocity**2
       envxy_nsi(:,:),        &      ! Ensemble vx*vy
       enpre_nsi(:),          &      ! Ensemble pressure
       enpr2_nsi(:),          &      ! Ensemble pressure**2
       entan_nsi(:,:),        &      ! Ensemble tangential force
       enmut_nsi(:),          &      ! Ensemble turbulent viscosity
       enstx_nsi(:,:),        &      ! Ensemble stress mu_t * grad(u)
       ensty_nsi(:,:),        &      ! Ensemble stress mu_t * grad(v)
       enstz_nsi(:,:),        &      ! Ensemble stress mu_t * grad(w)
       resch_nsi(:),          &      ! Schur complement residual
       remom_nsi(:,:),        &      ! Momentum residual
       prope_nsi(:,:),        &      ! Smoothed fluid Properties
       norle_nsi(:,:),        &      ! Normal to the zero Level Set 
       curle_nsi(:),          &      ! Curvature
       outflow_mass(:),       &      ! Outflow mass
       dt_rho_nsi(:),         &      ! Projection of dt / rho
       mass_rho_nsi(:,:),     &      ! M rho (not exchanged) saved for last timestep
       tau_nsi(:),            &      ! Projection of tau
       bubble_nsi(:),         &      ! Pressure bubble
       bubble_aqq_nsi(:),     &      ! Bubble matrix Aqq
       bubble_aqu_nsi(:,:),   &      ! Bubble matrix Aqu
       bubble_aqp_nsi(:,:),   &      ! Bubble matrix Aqp 
       bubble_bq_nsi(:),      &      ! Bubble RHS
       lagra_nsi(:,:,:),      &      ! Lagrange multiplier velocity
       tauib_nsi(:,:,:),      &      ! Lagrange multiplier tau
       vafor_nsi(:,:),        &      ! variational force
       avvaf_nsi(:,:),        &      ! time averaged variational force
       bupor_nsi(:,:),        &      ! forces due to porous media
       rhsid_gravb(:),        &      ! rhsid only gravity and Boussinesq components
       drhodt_nsi(:),         &      ! dt / rho for explicit
       prdivcor_nsi(:),       &      ! Projection of the divergence of the Coriolis term
       Scorio_nsi(:,:),       &      ! Matrix for projection of the divergence of the Coriolis term
       fsifo_nsi(:,:),        &      ! FSI body force received from solid in IB coupling for deformable bodies
       vefix(:,:)                    ! Velocity interp/projected from solid in EFECT

  !$acc declare create(mass_rho_nsi)  
  !$acc declare create(dt_rho_nsi)  
  
  type(r3p),   pointer     :: &
       turmu_nsi(:)                  ! LES turbulent viscosity
  real(rp)                 :: &
       porfo_nsi(3)                  ! Don not know what
  type(nsimat),   pointer  :: &
       intfo_nsi(:)                  ! Internal force
  type(typ_color_coupling) :: &
       wallcoupling                  ! Wall exchange
  !
  ! Solver
  ! 
  integer(ip), pointer     :: &
       kfl_fixno_div_nsi(:,:)        ! Nodal fixity for divergence free correction 
  integer(ip)              :: &
       nmauu_nsi,             &      ! Size of Auu
       nmaup_nsi,             &      ! Size of Aup
       nmapu_nsi,             &      ! Size of Apu
       nmapp_nsi,             &      ! Size of App
       poauu_nsi,             &      ! Pointer to Auu
       poaup_nsi,             &      ! Pointer to Aup
       poapu_nsi,             &      ! Pointer to Apu
       poapp_nsi,             &      ! Pointer to App
       nschu_nsi,             &      ! # Schur complement solves
       nmome_nsi                     ! # Momentum solves
  real(rp),    contiguous,    &
       &       pointer     :: &
       lapla_nsi(:)                  ! Laplacian matrix
       
  real(rp),    pointer     :: &
       Auu_nsi(:),            &      ! Auu
       Aup_nsi(:),            &      ! Aup
       Apu_nsi(:),            &      ! Apu
       App_nsi(:),            &      ! App
       amatr_nsi(:),          &      ! Linear matrix
       visco_nsi(:),          &      ! Viscous matrix
       cmama_nsi(:),          &      ! Consistent_mass_ matrix
       deltp_nsi(:),          &      ! Schur complement: Dp (used for mass correction)
       vepro_nsi(:,:),        &      ! Velocity projection
       grpro_nsi(:,:),        &      ! Pressure gradient projection
       prpro_nsi(:),          &      ! Pressure projection
       vepr2_nsi(:,:),        &      ! Velocity projection
       grpr2_nsi(:,:),        &      ! Pressure gradient projection
       prpr2_nsi(:)                  ! Pressure projection
  type(r1p),   pointer     :: &
       hydro_density_nsi(:)          ! Hydrostatic density
  !
  ! File names
  !
  character(150)           :: &
       fil_conve_nsi,         &      ! Convergence file name
       fil_dynin_nsi,         &      ! Input dynamic model
       fil_dynou_nsi,         &      ! Output dynamic model
       fil_windt_nsi                 ! Wind turbine output file
  
end module def_nastin_aux
!> @}
