!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_elmres.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute the residual of the momentum and continuity equations
!> @details Read physical data
!>    Compute the residual of the momentum and continuity equations:
!>    rho/(dt*theta)*u + sig*u - mu*Lap(u) 
!>    rho*(uc.grad)u is computed in nsi_elmsgs as uc can depend on the
!>    subgrid scale
!>
!>    Three forms of the viscous term are considered.
!>    1. Laplacian form:  div[mu*grad(u)]
!>    2. Divergence form: div[2*mu*eps(u)]
!>    3. Complete form:   div[2*mu*epsi(u)]
!>    eps(u)  = 1/2[grad(u)+grad(u)^t]
!>    eps'(u) = 1/2[grad(u)+grad(u)^t]-1/3*div(u)I
!>    
!>    <------ Laplacian
!>    <------------------------ Divergence
!>    <------------------------------------------------- Complete
!>    -mu*(d^2ui/dxj^2)-mu*(d^2uj/dxj*dxi)+2/3*mu*d(div(u))/dxi
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_elmres(                                              &
     pnode,pgaus,plapl,gpsha,gpcar,gphes,gpgvi,gpden,gpvis,         &
     gppor,gptem,gpsgs,elvel,elpre,elvep,elprp,elgrp,               &
     eltem,elmsh,elcod,elnor,elcur,elbub,elwmean,hleng,chale,gpvel, &
     gpgpr,rmomu,rmom2,rcont,gprhs,gprhc,gplap,gpadv,gpvep,gpprp,   &
     gpgrp,gpdhy,gpmsh,gpgve,gpnor,gpcur,gpfle,ielem,               &
     gprh2,gppre,gprhs_sgs,dtinv_loc,gpgde,gpsha_bub,densi)

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_coupl,zeror,&
       &                      ID_NASTIN,ID_ALEFOR,ID_TEMPER,         &
       &                      ID_CHEMIC,tesgs,kfl_lumped
  use def_kermod, only     :  thicl,gasco,kfl_adj_prob
  use def_domain, only     :  ndime,ntens,mnode
  use def_nastin, only     :  grnor_nsi,gravi_nsi,     &
       &                      kfl_sgsti_nsi,bougr_nsi,boube_nsi,     &
       &                      boutr_nsi,kfl_cotem_nsi,kfl_advec_nsi, &
       &                      dtsgs_nsi,kfl_stabi_nsi,kfl_sgsco_nsi, &
       &                      kfl_hydro_gravity_nsi,                 &
       &                      nbdfp_nsi,pabdf_nsi,bemol_nsi,         &
       &                      corio_nsi,facca_nsi,fvela_nsi,         &
       &                      kfl_regim_nsi,                         &
       &                      kfl_prthe_nsi,                         &
       &                      faccl_nsi,frotc_nsi,fvins_nsi,         &
       &                      kfl_surte_nsi,surte_nsi,kfl_rmom2_nsi, &
       &                      kfl_linea_nsi,centr_nsi,gravb_nsi,     &
       &                      hydro_density_nsi,kfl_confi_nsi,       &
       &                      kfl_sgsli_nsi,kfl_penal_nsi

!  use def_nastin, only     :  lowtr_nsi,penal_nsi
  use def_master, only     :  prthe
  use def_parame, only     :  pi

  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,plapl,ielem
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gphes(ntens,mnode,pgaus) 
  real(rp),    intent(in)  :: gpgvi(ndime,pgaus)
  real(rp),    intent(in)  :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)  :: gppor(pgaus)
  real(rp),    intent(out) :: gptem(pgaus, nbdfp_nsi)
  real(rp),    intent(in)  :: gpsgs(ndime,pgaus,*)
  real(rp),    intent(in)  :: gpgde(ndime,pgaus)
  real(rp),    intent(in)  :: elvel(ndime,pnode,nbdfp_nsi)
  real(rp),    intent(in)  :: elpre(pnode)
  real(rp),    intent(in)  :: elvep(ndime,pnode)
  real(rp),    intent(in)  :: elprp(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(in)  :: eltem(pnode,nbdfp_nsi) 
  real(rp),    intent(in)  :: elwmean(pnode,nbdfp_nsi) 
  real(rp),    intent(in)  :: elmsh(ndime,pnode) 
  real(rp),    intent(in)  :: elcod(ndime,pnode) 
  real(rp),    intent(in)  :: elnor(ndime,pnode) 
  real(rp),    intent(in)  :: elcur(pnode) 
  real(rp),    intent(in)  :: elbub                      ! bubble
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: chale(2)    
  real(rp),    intent(in)  :: dtinv_loc
  real(rp),    intent(out) :: gpvel(ndime,pgaus,nbdfp_nsi)
  real(rp),    intent(out) :: gpgpr(ndime,pgaus)
  real(rp),    intent(out) :: rmomu(pnode,pgaus)
  real(rp),    intent(out) :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(out) :: rcont(ndime,pnode,pgaus)
  real(rp),    intent(out) :: gprhs(ndime,pgaus)
  real(rp),    intent(out) :: gprhs_sgs(ndime,pgaus)
  real(rp),    intent(out) :: gprhc(pgaus)
  real(rp),    intent(out) :: gprh2(pgaus)
  real(rp),    intent(out) :: gplap(pnode,pgaus)
  real(rp),    intent(out) :: gpadv(ndime,pgaus)
  real(rp),    intent(out) :: gpvep(ndime,pgaus)
  real(rp),    intent(out) :: gpprp(pgaus)
  real(rp),    intent(out) :: gpgrp(ndime,pgaus)  
  real(rp),    intent(out) :: gpdhy(pgaus)
  real(rp),    intent(out) :: gpmsh(ndime,pgaus)
  real(rp),    intent(out) :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(out) :: gpnor(ndime,pgaus)         ! LS normal
  real(rp),    intent(out) :: gpcur(pgaus)               ! LS curvature
  real(rp),    intent(in)  :: gpfle(pgaus)               ! LS - already calculated in nsi_proper
  real(rp),    intent(out) :: gppre(pgaus)               ! pressure
  real(rp),    intent(out) :: gpsha_bub(pgaus)           ! bubble shape function
  real(rp),    intent(out) :: densi(pgaus,nbdfp_nsi)
  integer(ip)              :: idime,inode,igaus,jdime
  integer(ip)              :: itime
  real(rp)                 :: fact0,fact1,fact2,fact3
  real(rp)                 :: w1
  real(rp)                 :: gpcod(3),alpha(3),centf(3) ! Axes motion
  real(rp)                 :: gpgrt(pgaus),gpdet(pgaus)  ! Low mach intermmediate vars
  real(rp)                 :: xvis2
  real(rp)                 :: sgste(pgaus,nbdfp_nsi)
  real(rp)                 :: gpwmean(pgaus,nbdfp_nsi),w_time
  real(rp)                 :: dtinv_res
#ifdef matiaslma
  real(rp)                 :: gpgte(ndime,pgaus)
#endif

  !----------------------------------------------------------------------
  !
  ! Time step
  !
  !----------------------------------------------------------------------

  if( kfl_lumped == 2 ) then
     dtinv_res = 0.0_rp
  else
     dtinv_res = dtinv_loc
  end if

  gprhs     = 0.0_rp
  gprhs_sgs = 0.0_rp
  gprhc     = 0.0_rp
  gprh2     = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! Gauss point values
  !
  !----------------------------------------------------------------------
  !
  ! GPVEL, GPGPR: u, grad(p)
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpvel(idime,igaus,1) = 0.0_rp
        gpvel(idime,igaus,2) = 0.0_rp
        do itime=3,nbdfp_nsi
           gpvel(idime,igaus,itime) = 0.0_rp
        end do
        gpgpr(idime,igaus)   = 0.0_rp

     end do
     do itime =1, nbdfp_nsi
        sgste(igaus, itime) = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpvel(idime,igaus,1) = gpvel(idime,igaus,1) + elvel(idime,inode,1) * gpsha(inode,igaus)
           gpvel(idime,igaus,2) = gpvel(idime,igaus,2) + elvel(idime,inode,2) * gpsha(inode,igaus)
           do itime=3,nbdfp_nsi
              gpvel(idime,igaus,itime) = gpvel(idime,igaus,itime) + elvel(idime,inode,itime) * gpsha(inode,igaus)
           end do
           gpgpr(idime,igaus)   = gpgpr(idime,igaus)   + elpre(inode) * gpcar(idime,inode,igaus)
        end do
     end do
  end do
  !
  ! Pressure
  !
  if( kfl_penal_nsi /= 0 .or. ( kfl_confi_nsi == 1 .and. kfl_regim_nsi == 3 ) ) then
     gppre = 0.0_rp
     do igaus = 1,pgaus
        do inode = 1,pnode
           gppre(igaus) = gppre(igaus) + elpre(inode) * gpsha(inode,igaus)
        end do
     end do
  end if
  !
  ! GPADV: Advection velocity
  !
  if( kfl_advec_nsi /= 0 ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpadv(idime,igaus) = gpvel(idime,igaus,1)
        end do
     end do
  else
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpadv(idime,igaus) = 0.0_rp
        end do
     end do
  end if
  !
  ! ALE
  !
  if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then ! TESTEO 
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpmsh(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           do idime = 1,ndime
              gpmsh(idime,igaus) = gpmsh(idime,igaus) + elmsh(idime,inode) * gpsha(inode,igaus)
           end do
        end do
        do idime = 1,ndime
           gpadv(idime,igaus) = gpadv(idime,igaus) - gpmsh(idime,igaus)
        end do
     end do
  end if
  !
  ! Convection tracking: uc = u + u'
  !
  if( kfl_sgsco_nsi >= 1 ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpadv(idime,igaus) = gpadv(idime,igaus) + gpsgs(idime,igaus,1)
        end do
     end do
  end if
  !
  ! Temperature and wmean
  !
  if( kfl_cotem_nsi /= 0 .or. kfl_regim_nsi == 3) then

     do igaus = 1,pgaus
        do itime =1,nbdfp_nsi
           gptem(igaus,itime)   = 0.0_rp
           gpwmean(igaus,itime) = 0.0_rp
           do inode = 1,pnode
              gptem(igaus,itime)   = gptem(igaus,itime)   + eltem(inode,itime)   * gpsha(inode,igaus)
              gpwmean(igaus,itime) = gpwmean(igaus,itime) + elwmean(inode,itime) * gpsha(inode,igaus)
           end do
        end do
     end do

     if (kfl_regim_nsi == 3) then   ! loads density in last time steps
        if(associated(tesgs)) then
           do igaus = 1,pgaus
              do itime =1,nbdfp_nsi  
                 sgste(igaus, itime) =  tesgs(ielem)%a(1,igaus,itime)
                 densi(igaus, itime) =  prthe(itime) * gpwmean(igaus,itime) /gasco /(gptem(igaus,itime)+sgste(igaus, itime))
              end do
           end do
        else
           do igaus = 1,pgaus
              do itime =1,nbdfp_nsi
                 densi(igaus, itime) =  prthe(itime) * gpwmean(igaus,itime) /gasco /gptem(igaus,itime)
              end do
           end do
        end if
     end if
  end if
  !
  ! Surface tension
  !
  if( kfl_surte_nsi /= 0 ) then
     do igaus = 1,pgaus
        gpcur(igaus) = 0.0_rp
        do idime = 1,ndime
           gpnor(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           gpcur(igaus) = gpcur(igaus) + elcur(inode) * gpsha(inode,igaus)
           do idime = 1,ndime
              gpnor(idime,igaus) = gpnor(idime,igaus) + elnor(idime,inode) * gpsha(inode,igaus)
           end do
        end do
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Projections at Gauss points
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus
     do idime = 1,ndime
        gpvep(idime,igaus) = 0.0_rp
     end do
     gpprp(igaus) = 0.0_rp
  end do

  if( kfl_stabi_nsi /= 0 .and. kfl_stabi_nsi /= -1 ) then

     do igaus = 1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              gpvep(idime,igaus) = gpvep(idime,igaus) + gpsha(inode,igaus) * elvep(idime,inode)
           end do
           gpprp(igaus) = gpprp(igaus) + gpsha(inode,igaus) * elprp(inode)
        end do
     end do

     if( kfl_stabi_nsi == 2 ) then
        do igaus = 1,pgaus
           do idime = 1,ndime
              gpgrp(idime,igaus) = 0.0_rp
           end do
        end do
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 gpgrp(idime,igaus) = gpgrp(idime,igaus) + gpsha(inode,igaus) * elgrp(idime,inode)
              end do
           end do
        end do
     end if

  end if


  if( kfl_stabi_nsi == 2 .or. kfl_stabi_nsi == -1 ) then

     !----------------------------------------------------------------------
     !
     ! Split OSS: only initialize RMOMU
     !
     !----------------------------------------------------------------------

     do igaus = 1,pgaus
        do inode = 1,pnode
           rmomu(inode,igaus) = 0.0_rp
        end do
     end do

  else

     !----------------------------------------------------------------------
     !
     ! ASGS and full OSS: residual Gauss point values
     !
     !----------------------------------------------------------------------

     if( plapl == 1 .and. ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              gplap(inode,igaus) =        &
                   + gphes(1,inode,igaus) & 
                   + gphes(2,inode,igaus)
           end do
        end do
     else if( plapl == 1 .and. ndime == 3 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              gplap(inode,igaus) =        &
                   + gphes(1,inode,igaus) &
                   + gphes(2,inode,igaus) &
                   + gphes(3,inode,igaus)
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              gplap(inode,igaus) = 0.0_rp
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! RMOMU, RCONT: Momentum and continuity residuals
     !
     !----------------------------------------------------------------------
     !
     ! rho/(dt*theta)*u + sig*u - mu*d^2u/dxk^2
     !
     do igaus = 1,pgaus
        fact1 = dtinv_res * pabdf_nsi(1) * gpden(igaus) + gppor(igaus)
        do inode = 1,pnode
           rmomu(inode,igaus) = fact1 * gpsha(inode,igaus) - gpvis(igaus) * gplap(inode,igaus)
        end do
     end do
     !
     ! Continuity and Momentum turbulent viscosity: -grad(mu).grad(u)
     !
     if(ndime==2) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              rcont(1,inode,igaus) = gpcar(1,inode,igaus)
              rcont(2,inode,igaus) = gpcar(2,inode,igaus)
              rmomu(inode,igaus)   = rmomu(inode,igaus) - gpgvi(1,igaus) * gpcar(1,inode,igaus)
              rmomu(inode,igaus)   = rmomu(inode,igaus) - gpgvi(2,igaus) * gpcar(2,inode,igaus)
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              rcont(1,inode,igaus) = gpcar(1,inode,igaus)
              rcont(2,inode,igaus) = gpcar(2,inode,igaus)
              rcont(3,inode,igaus) = gpcar(3,inode,igaus)
              rmomu(inode,igaus)   = rmomu(inode,igaus) - gpgvi(1,igaus) * gpcar(1,inode,igaus)
              rmomu(inode,igaus)   = rmomu(inode,igaus) - gpgvi(2,igaus) * gpcar(2,inode,igaus)
              rmomu(inode,igaus)   = rmomu(inode,igaus) - gpgvi(3,igaus) * gpcar(3,inode,igaus)
           end do
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! GPRHS
  !
  ! 1. ASGS + Full OSS: 
  !    RHS = rho*dt*u^n + rho*g - rho*beta*g*(T-Tr)
  !    RHS_SGS = rho/dt*u'n + proje
  ! 2. Split OSS: only force term without temporal one
  !    RHS = rho*g - rho*beta*g*(T-Tr)
  !
  !----------------------------------------------------------------------
  !
  ! GPDHY: Hydrostatic density
  !
  if( kfl_hydro_gravity_nsi /= 0 ) then
     !
     ! rho of the hydrostatic state when couplign with level set
     !
     gpdhy(1:pgaus) = hydro_density_nsi(ielem) % a(1:pgaus)

  else if (kfl_regim_nsi == 3 .and. kfl_prthe_nsi==0 ) then 
     !
     ! Low-Mach regime with open flow, constant thermodynamic pressure
     ! 
     if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 .or. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 ) then
        w1 = 0.0289_rp ! Molecular weight of air for reference
     else
        w1 = 1.0_rp
     endif
     ! gpdhy(1:pgaus) =  prthe(1)*w1/(lowtr_nsi * gasco)
     gpdhy = 0.0_rp
  else
     gpdhy = 0.0_rp
  end if
  !
  ! GPGVE: grad(u)
  !
  if( ( fvins_nsi > 0.9_rp ) .or. bemol_nsi > 0.001_rp .or. kfl_linea_nsi == 2  .or. kfl_adj_prob == 1.or.kfl_sgsli_nsi==0) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           do jdime = 1,ndime
              gpgve(jdime,idime,igaus) = 0.0_rp
           end do
        end do
        do inode = 1,pnode
           do idime = 1,ndime
              do jdime = 1,ndime
                 gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                      + gpcar(jdime,inode,igaus) * elvel(idime,inode,1)
              end do
           end do
        end do
     end do
  end if

  if( ndime == 2 ) then
     !
     ! Time derivative: rho/dt*u^n 
     !
     if( kfl_stabi_nsi /= 2 .and. kfl_stabi_nsi /= -1 ) then
        do itime = 2, nbdfp_nsi
           do igaus = 1,pgaus
              fact1 = gpden(igaus) * dtinv_res * pabdf_nsi(itime)
              gprhs(1,igaus) = gprhs(1,igaus) - fact1 * gpvel(1,igaus,itime)  
              gprhs(2,igaus) = gprhs(2,igaus) - fact1 * gpvel(2,igaus,itime)
           end do
        end do
     end if
     !
     ! Gravity: rho*g
     !
     do igaus = 1,pgaus
        fact2 = ( gpden(igaus) - gpdhy(igaus) ) * grnor_nsi
        gprhs(1,igaus) = gprhs(1,igaus) + fact2 * gravi_nsi(1)
        gprhs(2,igaus) = gprhs(2,igaus) + fact2 * gravi_nsi(2)
     end do
     !
     ! Boussinesq coupling: -rho*beta*g*(T-Tr)
     !
     if( kfl_cotem_nsi == 1) then 
        fact1 = bougr_nsi*boube_nsi
        do igaus = 1,pgaus
           fact2 = gpden(igaus) * fact1 * ( gptem(igaus,1) - boutr_nsi )
           gprhs(1,igaus) = gprhs(1,igaus) - fact2*gravb_nsi(1)
           gprhs(2,igaus) = gprhs(2,igaus) - fact2*gravb_nsi(2)
        end do
     end if
     !
     ! Surface tension: sigma*curvature*regularized_delta*normal
     !
     if( kfl_surte_nsi == 1) then
        do igaus = 1,pgaus
           if (abs(gpfle(igaus)) < thicl) then
              fact1 = (0.5_rp/thicl) * (1.0_rp + cos ( pi * gpfle(igaus) / thicl ) )
              fact2 = surte_nsi * gpcur(igaus) * fact1
              gprhs(1,igaus) = gprhs(1,igaus) - fact2 * gpnor(1,igaus)
              gprhs(2,igaus) = gprhs(2,igaus) - fact2 * gpnor(2,igaus)
           end if
        end do
     end if
     !
     ! Subgrid scale time derivative: rho/dt*u'^n
     !
     if( kfl_sgsti_nsi == 1 .and. kfl_stabi_nsi /= 2 .and. kfl_stabi_nsi /= -1 ) then
        if (kfl_regim_nsi==3) then
           do itime =2, 2
              do igaus = 1,pgaus              
                 gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + densi(igaus, itime)* gpsgs(1,igaus,itime)*dtsgs_nsi
                 gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + densi(igaus, itime)* gpsgs(2,igaus,itime)*dtsgs_nsi
              end do
           end do
        else
           do igaus = 1,pgaus
              fact1 = gpden(igaus) * dtsgs_nsi    
              gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + fact1 * gpsgs(1,igaus,2) ! + gpvep(1,igaus) / gpst1(igaus)
              gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + fact1 * gpsgs(2,igaus,2) ! + gpvep(2,igaus) / gpst1(igaus)
           end do
        end if
     end if

     !
     !!FER Note: Variable viscosity is now assembled on the left side
     !

  else
     !
     ! Time derivative: rho/dt*u^n 
     !
     if( kfl_stabi_nsi /= 2 .and. kfl_stabi_nsi /= -1 ) then
        do itime = 2, nbdfp_nsi
           do igaus = 1,pgaus
              fact1 = gpden(igaus) * dtinv_res * pabdf_nsi(itime)
              gprhs(1,igaus) = gprhs(1,igaus) - fact1 * gpvel(1,igaus,itime)  
              gprhs(2,igaus) = gprhs(2,igaus) - fact1 * gpvel(2,igaus,itime)
              gprhs(3,igaus) = gprhs(3,igaus) - fact1 * gpvel(3,igaus,itime)
           end do
        end do
     end if
     !
     ! Gravity: rho*g
     !
     do igaus = 1,pgaus 
        fact2 = ( gpden(igaus) - gpdhy(igaus) ) * grnor_nsi
        gprhs(1,igaus) = gprhs(1,igaus) + fact2 * gravi_nsi(1) 
        gprhs(2,igaus) = gprhs(2,igaus) + fact2 * gravi_nsi(2) 
        gprhs(3,igaus) = gprhs(3,igaus) + fact2 * gravi_nsi(3) 
     end do
     !
     ! Boussinesq coupling: -rho*beta*g*(T-Tr)
     !
     if( kfl_cotem_nsi == 1) then
        fact1 = bougr_nsi * boube_nsi
        do igaus = 1,pgaus
           fact2 = gpden(igaus) * fact1 * ( gptem(igaus,1) - boutr_nsi )
           gprhs(1,igaus) = gprhs(1,igaus) - fact2 * gravb_nsi(1)
           gprhs(2,igaus) = gprhs(2,igaus) - fact2 * gravb_nsi(2)
           gprhs(3,igaus) = gprhs(3,igaus) - fact2 * gravb_nsi(3)
        end do
     end if
     !
     ! Surface tension: sigma*curvature*regularized_delta*normal
     !
     if( kfl_surte_nsi == 1) then
        do igaus = 1,pgaus
           if ( abs(gpfle(igaus)) < thicl ) then
              fact1 = (0.5_rp/thicl) * (1.0_rp + cos ( pi * gpfle(igaus) / thicl ) )
              fact2 = surte_nsi * gpcur(igaus) * fact1
              gprhs(1,igaus) = gprhs(1,igaus) - fact2 * gpnor(1,igaus)
              gprhs(2,igaus) = gprhs(2,igaus) - fact2 * gpnor(2,igaus)
              gprhs(3,igaus) = gprhs(3,igaus) - fact2 * gpnor(3,igaus)
           end if
        end do
     end if
     !
     ! Subgrid scale time derivative: rho/dt*u'n
     !
     if( kfl_sgsti_nsi == 1 .and. kfl_stabi_nsi /= 2 ) then
        if (kfl_regim_nsi==3) then      
           do itime =2, 2
              do igaus = 1,pgaus
                 gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + densi(igaus, itime)*gpsgs(1,igaus,itime)*dtsgs_nsi
                 gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + densi(igaus, itime)*gpsgs(2,igaus,itime)*dtsgs_nsi
                 gprhs_sgs(3,igaus) = gprhs_sgs(3,igaus) + densi(igaus, itime)*gpsgs(3,igaus,itime)*dtsgs_nsi
              end do
           end do
        else
           do igaus = 1,pgaus
              fact1 = gpden(igaus) * dtsgs_nsi
              gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + fact1 * gpsgs(1,igaus,2) !+ gpvep(1,igaus)/ gpst1(igaus)
              gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + fact1 * gpsgs(2,igaus,2) !+ gpvep(2,igaus)/ gpst1(igaus)
              gprhs_sgs(3,igaus) = gprhs_sgs(3,igaus) + fact1 * gpsgs(3,igaus,2) !+ gpvep(3,igaus)/ gpst1(igaus)      
           end do
        end if
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! RMOM2
  !
  ! Off-diagonal momentum operator
  !
  !----------------------------------------------------------------------

  if( kfl_rmom2_nsi /= 0 ) rmom2 = 0.0_rp
  !
  ! Coriolis force
  !
  ! 2*rho*(w x u)
  ! x-equation: w x u = wy*uz - wz*uy
  ! y-equation: w x u = wz*ux - wx*uz
  ! z-equation: w x u = wx*uy - wy*ux
  !
  !
  if( corio_nsi > zeror ) then
     if( ndime == 2 ) then
        do igaus = 1,pgaus           
           fact3 = 2.0_rp * gpden(igaus) * fvela_nsi(3)
           do inode = 1,pnode 
              rmom2(1,2,inode,igaus) = -fact3 * gpsha(inode,igaus)     ! -wz*uy 
              rmom2(2,1,inode,igaus) =  fact3 * gpsha(inode,igaus)     !  wz*ux
           end do
        end do
     else
        do igaus = 1,pgaus   
           fact0 = 2.0_rp * gpden(igaus) 
           fact1 = fact0  * fvela_nsi(1)
           fact2 = fact0  * fvela_nsi(2)
           fact3 = fact0  * fvela_nsi(3)
           do inode = 1,pnode
              rmom2(1,2,inode,igaus) = - fact3 * gpsha(inode,igaus)  ! -wz*uy  
              rmom2(1,3,inode,igaus) =   fact2 * gpsha(inode,igaus)  !  wy*uz  
              rmom2(2,1,inode,igaus) =   fact3 * gpsha(inode,igaus)  !  wz*ux  
              rmom2(2,3,inode,igaus) = - fact1 * gpsha(inode,igaus)  ! -wx*uz  
              rmom2(3,1,inode,igaus) = - fact2 * gpsha(inode,igaus)  ! -wy*ux  
              rmom2(3,2,inode,igaus) =   fact1 * gpsha(inode,igaus)  !  wx*uy                  
           end do
        end do
     end if
  end if
  !
  !!FER CAREFUL THIS DOES NOT TAKE INTO ACCOUNT CHEMIC
  ! Viscous term: only terms from divergence and complete forms
  !
  ! - ( div[-2 mu eps(u)] , v ) = - d/dxj ( -2*mu* ( dui/dxj + duj/dxi - 2/3*mu div(u) delta_ij )
  !
  ! Laplacian  form: A=0, B=0, eps(u) = 1/2 grad(u)
  ! Divergence form: A=1, B=0, eps(u) = 1/2 ( grad(u) + grad(u)^t )
  ! Complete   form: A=1, B=1, eps(u) = 1/2 ( grad(u) + grad(u)^t ) - 1/3 div(u) I
  !
  ! + ( mu d^2ui/dxj^2 , vi )               (6)        already computed
  ! + ( dmu/dxj dui/dxj , vi )              (7)        already computed
  ! + A * ( mu  d^2uj/dxidxj , vi )         (8)        divergence                           
  ! + A * ( dmu/dxj duj/dxi , vi )          (9)        divergence
  ! - 2/3 * B * ( dmu/dxi (div u) , vi )   (10)        complete
  ! - 2/3 * B * ( mu d(div u)/dxi , vi )   (11)        complete  
  !       
  if( fvins_nsi > 0.9_rp ) then

     if( fvins_nsi > 1.9_rp ) then
        xvis2 = 2.0_rp / 3.0_rp
     else
        xvis2 = 0.0_rp
     end if

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact0 = xvis2 * gpvis(igaus) - gpvis(igaus)
           fact1 = xvis2 * gpgvi(1,igaus) 
           fact2 = xvis2 * gpgvi(2,igaus) 
           do inode = 1,pnode
              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact0 * gphes(1,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dx^2 + 2/3 * mu * d^2ux/dx^2
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact0 * gphes(3,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dxdy + 2/3 * mu * d^2uy/dxdy
              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact0 * gphes(3,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdy + 2/3 * mu * d^2ux/dxdy
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact0 * gphes(2,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dy^2 + 2/3 * mu * d^2uy/dy^2

              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) - gpgvi(1,igaus) * gpcar(1,inode,igaus)   ! (9)         - dmu/dx * dux/dx
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) - gpgvi(2,igaus) * gpcar(1,inode,igaus)   ! (9)         - dmu/dy * duy/dx
              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) - gpgvi(1,igaus) * gpcar(2,inode,igaus)   ! (9)         - dmu/dx * dux/dy
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) - gpgvi(2,igaus) * gpcar(2,inode,igaus)   ! (9)         - dmu/dy * dvy/dy

              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact1 * gpcar(1,inode,igaus)            ! (10)        + 2/3 * dmu/dx * dux/dx 
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact1 * gpcar(2,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duy/dy
              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact2 * gpcar(1,inode,igaus)            ! (10)        + 2/3 * dmu/dy * dux/dx
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact2 * gpcar(2,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duy/dy         
           end do
        end do

     else if( ndime == 3 ) then
        do igaus = 1,pgaus
           fact0 = xvis2 * gpvis(igaus) - gpvis(igaus)
           fact1 = xvis2 * gpgvi(1,igaus) 
           fact2 = xvis2 * gpgvi(2,igaus) 
           fact3 = xvis2 * gpgvi(3,igaus) 
           do inode = 1,pnode
              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact0 * gphes(1,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dx^2 + 2/3 * mu * d^2ux/dx^2
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact0 * gphes(4,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dxdy + 2/3 * mu * d^2uy/dxdy
              rmom2(1,3,inode,igaus) = rmom2(1,3,inode,igaus) + fact0 * gphes(5,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dxdz + 2/3 * mu * d^2uz/dxdz

              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact0 * gphes(4,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdy + 2/3 * mu * d^2ux/dxdy
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact0 * gphes(2,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dy^2 + 2/3 * mu * d^2uy/dy^2
              rmom2(2,3,inode,igaus) = rmom2(2,3,inode,igaus) + fact0 * gphes(6,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dydz + 2/3 * mu * d^2uz/dydz  

              rmom2(3,1,inode,igaus) = rmom2(3,1,inode,igaus) + fact0 * gphes(5,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdz + 2/3 * mu * d^2ux/dxdz  
              rmom2(3,2,inode,igaus) = rmom2(3,2,inode,igaus) + fact0 * gphes(6,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dydz + 2/3 * mu * d^2uy/dydz            
              rmom2(3,3,inode,igaus) = rmom2(3,3,inode,igaus) + fact0 * gphes(3,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dz^2 + 2/3 * mu * d^2uz/dz^2  

              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) - gpgvi(1,igaus) * gpcar(1,inode,igaus)   ! (9)         - dmu/dx * dux/dx
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) - gpgvi(2,igaus) * gpcar(1,inode,igaus)   ! (9)         - dmu/dy * duy/dx
              rmom2(1,3,inode,igaus) = rmom2(1,3,inode,igaus) - gpgvi(3,igaus) * gpcar(1,inode,igaus)   ! (9)         - dmu/dz * duz/dx

              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) - gpgvi(1,igaus) * gpcar(2,inode,igaus)   ! (9)         - dmu/dx * dux/dy
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) - gpgvi(2,igaus) * gpcar(2,inode,igaus)   ! (9)         - dmu/dy * duy/dy
              rmom2(2,3,inode,igaus) = rmom2(2,3,inode,igaus) - gpgvi(3,igaus) * gpcar(2,inode,igaus)   ! (9)         - dmu/dz * duz/dy  

              rmom2(3,1,inode,igaus) = rmom2(3,1,inode,igaus) - gpgvi(1,igaus) * gpcar(3,inode,igaus)   ! (9)         - dmu/dx * dux/dz
              rmom2(3,2,inode,igaus) = rmom2(3,2,inode,igaus) - gpgvi(2,igaus) * gpcar(3,inode,igaus)   ! (9)         - dmu/dy * duy/dz           
              rmom2(3,3,inode,igaus) = rmom2(3,3,inode,igaus) - gpgvi(3,igaus) * gpcar(3,inode,igaus)   ! (9)         - dmu/dz * duz/dz

              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact1 * gpcar(1,inode,igaus)            ! (10)        + 2/3 * dmu/dx * dux/dx 
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact1 * gpcar(2,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duy/dy
              rmom2(1,3,inode,igaus) = rmom2(1,3,inode,igaus) + fact1 * gpcar(3,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duz/dz

              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact2 * gpcar(1,inode,igaus)            ! (10)        + 2/3 * dmu/dy * dux/dx
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact2 * gpcar(2,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duy/dy         
              rmom2(2,3,inode,igaus) = rmom2(2,3,inode,igaus) + fact2 * gpcar(3,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duz/dz         

              rmom2(3,1,inode,igaus) = rmom2(3,1,inode,igaus) + fact3 * gpcar(1,inode,igaus)            ! (10)        + 2/3 * dmu/dz * dux/dx
              rmom2(3,2,inode,igaus) = rmom2(3,2,inode,igaus) + fact3 * gpcar(2,inode,igaus)            ! (10)        + 2/3 * dmu/dz * duy/dy         
              rmom2(3,3,inode,igaus) = rmom2(3,3,inode,igaus) + fact3 * gpcar(3,inode,igaus)            ! (10)        + 2/3 * dmu/dz * duz/dz         

           end do
        end do

     end if

  end if

  !
  ! Adjoint Newton term: rho*( uj d(adv)/d(xj) , v ) 
  !
  if( kfl_adj_prob == 1 ) then
     if( ndime == 2 ) then
        do igaus = 1,pgaus           
           do inode = 1,pnode
              fact0 = gpden(igaus) * gpsha(inode,igaus)
              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact0 * gpgve(1,1,igaus)   ! rho * ux * d(adv_x)/dx
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact0 * gpgve(2,1,igaus)   ! rho * uy * d(adv_x)/dy
              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact0 * gpgve(1,2,igaus)   ! rho * ux * d(adv_y)/dx
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact0 * gpgve(2,2,igaus)   ! rho * uy * d(adv_y)/dy
           end do
        end do
     else
        do igaus = 1,pgaus           
           do inode = 1,pnode
              fact0 = gpden(igaus) * gpsha(inode,igaus)
              rmom2(1,1,inode,igaus) = rmom2(1,1,inode,igaus) + fact0 * gpgve(1,1,igaus)   ! rho * ux * d(adv_x)/dx
              rmom2(1,2,inode,igaus) = rmom2(1,2,inode,igaus) + fact0 * gpgve(2,1,igaus)   ! rho * uy * d(adv_x)/dy
              rmom2(1,3,inode,igaus) = rmom2(1,3,inode,igaus) + fact0 * gpgve(3,1,igaus)   ! rho * uz * d(adv_x)/dz
              rmom2(2,1,inode,igaus) = rmom2(2,1,inode,igaus) + fact0 * gpgve(1,2,igaus)   ! rho * ux * d(adv_y)/dx
              rmom2(2,2,inode,igaus) = rmom2(2,2,inode,igaus) + fact0 * gpgve(2,2,igaus)   ! rho * uy * d(adv_y)/dy
              rmom2(2,3,inode,igaus) = rmom2(2,3,inode,igaus) + fact0 * gpgve(3,2,igaus)   ! rho * uz * d(adv_y)/dz
              rmom2(3,1,inode,igaus) = rmom2(3,1,inode,igaus) + fact0 * gpgve(1,3,igaus)   ! rho * ux * d(adv_z)/dx
              rmom2(3,2,inode,igaus) = rmom2(3,2,inode,igaus) + fact0 * gpgve(2,3,igaus)   ! rho * uy * d(adv_z)/dy
              rmom2(3,3,inode,igaus) = rmom2(3,3,inode,igaus) + fact0 * gpgve(3,3,igaus)   ! rho * uz * d(adv_z)/dz
           end do
        end do
     end if
  end if
  
  !----------------------------------------------------------------------
  !
  ! GPRHS
  !
  ! Rotation term:            f = f - rho * ( w x w x r + dw/dt x r )
  ! Linear acceleration term: f = f - rho * a
  !
  !----------------------------------------------------------------------

  if( corio_nsi > zeror ) then

     do igaus = 1,pgaus
        !
        ! Gauss point coordinates wrt center of rotation
        !
        gpcod(1) = 0.0_rp
        gpcod(2) = 0.0_rp
        gpcod(3) = 0.0_rp
        do idime = 1,ndime
           do inode = 1,pnode
              gpcod(idime) = gpcod(idime) + gpsha(inode,igaus) * elcod(idime,inode)
           end do
           gpcod(idime) = gpcod(idime) - frotc_nsi(idime)
        end do
        !
        ! Angular acceleration dw/dt x r
        ! 
        call vecpro(facca_nsi,gpcod,alpha,3_ip)
        !
        ! Centrifugal force w x (w x r)
        !           
        if( frotc_nsi(1) < 1.0e10_rp ) then
           call vecpro(fvela_nsi,gpcod,centf,3_ip)
           call vecpro(fvela_nsi,centf,centf,3_ip)
           centf = centf * centr_nsi
        else
           centf = 0.0_rp
        end if
        !
        ! Total force: rho * [ - w x (w x r) - dw/dt x r - a ]
        !
        do idime = 1,ndime
           gprhs(idime,igaus) = gprhs(idime,igaus) - gpden(igaus) &
                &               * (   centf(idime)                &     ! w x (w x r)
                &                   + alpha(idime)                &     ! dw/dt x r
                &                   + faccl_nsi(idime)            )     ! a
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Iterative penalization
  !
  !----------------------------------------------------------------------

  !if( kfl_penal_nsi == 2 ) then
  !   gprhc(1:pgaus) = gprhc(1:pgaus) + penal_nsi * gppre(1:pgaus)
  !end if

  !----------------------------------------------------------------------
  !
  ! Low-Mach regime: We calculate RHS of continuity equation
  !
  !----------------------------------------------------------------------

  !
  ! Integrate by parts mass conservation equation
  !
!!DMM-LM #ifdef matiaslma 
!!DMM-LM #endif

  if( kfl_regim_nsi == 3 ) then 
     !
     ! rho = p^th / RT
     ! 
#ifdef matiaslma 
     ! integrate by pats mass conservation equation
     !  gprhc2 = (vel grad rho + dtrho )/rho
     ! as the divergence term is assembled integrating by parts, we need to 
     ! only take account for time derivative of density
     !
     ! q div(u) = - [ 1/rho * d(rho)/dt + u.grad(rho)/rho ]
     !
     do igaus = 1,pgaus
        ! add time derivative of density
        fact0 = gpden(igaus) *pabdf_nsi(1)
        do itime = 2, nbdfp_nsi           
           fact0 = fact0 + pabdf_nsi(itime)* densi(igaus,itime)
        end do
        fact0 = fact0 * dtinv_res

        gprhc(igaus) = - fact0/gpden(igaus)
        gprh2(igaus)=  0.0_rp
        gppre(igaus) = 0.0_rp
        do inode =1, pnode
           gppre(igaus)   = gppre(igaus)   + elpre(inode) * gpsha(inode,igaus)
        end do

        do idime = 1,ndime
           gpgte(idime, igaus) =0.0_rp
           do inode = 1, pnode
              gpgte(idime,igaus)  = gpgte(idime,igaus)   + eltem(inode,1) * gpcar(idime,inode,igaus)
           end do
           gprh2(igaus) = gprh2(igaus) +  gpadv(idime,igaus)*gpgte(idime, igaus)
        end do
        gprh2(igaus)= gprh2(igaus)/gptem(igaus,1) - fact0/gpden(igaus)               

     end do

#else

     do igaus = 1,pgaus
        gpdet(igaus) = 0.0_rp ! time derivative of density
        gpgrt(igaus) = 0.0_rp ! vel grad rho 
        fact0 = gpden(igaus) *pabdf_nsi(1)
        do itime = 2, nbdfp_nsi
           if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
             w_time = gpwmean(igaus,itime)
           else
             w_time = 1.0_rp 
           endif
           fact0 = fact0 + pabdf_nsi(itime)* densi(igaus,itime)
        end do
        !
        ! gpdet = Delta rho / Delta t = rho^n+1 - rho^n / dt
        !
        gpdet(igaus) = fact0*dtinv_res
        do idime = 1,ndime
          gpgrt(igaus)  = gpgrt(igaus) +  gpadv(idime,igaus) * gpgde(idime,igaus)
        end do
        gprhc(igaus) = ( - gpdet(igaus) - gpgrt(igaus) ) / gpden(igaus)
        gprh2(igaus) = gprhc(igaus)

     enddo

#endif
  endif

end subroutine nsi_elmres


