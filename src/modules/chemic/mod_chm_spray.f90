!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_spray

  implicit none

  private

  public :: chm_elmpre_spray
  public :: chm_elmprc_spray
  public :: chm_rhodt_spray

contains

   subroutine chm_elmpre_spray(ielem,pnode,pgaus,elvel,gpcar,gpcon,gpden,gptur,hleng,gpsigma)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_elmpre_spray
     ! NAME
     !    chm_elmpre_spray
     ! DESCRIPTION
     !    Compute gauss point terms for spray model
     ! USES
     ! USED BY
     !    chm_elmcfi
     !***
     !-----------------------------------------------------------------------
     use def_kintyp, only      :  ip,rp
     use def_domain, only      :  ndime,mnode
     use def_chemic, only      :  nclas_chm,surf_tension_chm
     use def_chemic, only      :  sigma_gp_chm,sigma0_gp_chm,d32_gp_chm
     use def_kermod, only      :  turmu_ker

     implicit none

     integer(ip),  intent(in)  :: ielem
     integer(ip),  intent(in)  :: pnode
     integer(ip),  intent(in)  :: pgaus
     real(rp),     intent(in)  :: elvel(ndime,pnode)
     real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
     real(rp),     intent(in)  :: gpcon(pgaus,nclas_chm)
     real(rp),     intent(in)  :: gpden(pgaus)
     real(rp),     intent(in)  :: gptur(pgaus)
     real(rp),     intent(in)  :: hleng(3)
     real(rp),     intent(out) :: gpsigma(pgaus)

     real(rp)                  :: gpgve(ndime,ndime,pgaus)
     real(rp)                  :: delta,delta2,rdelta,rdelta2
     real(rp)                  :: seci4,rtau_turb,k_sgs,sigma_crit
     real(rp)                  :: phi_clip
     integer(ip)               :: igaus,idime,jdime,inode
     integer(ip)               :: iclas_phi
     integer(ip)               :: iclas_sigma

     !
     ! Define liquid volume fraction and surface density variables
     !
     iclas_phi   = nclas_chm - 1
     iclas_sigma = nclas_chm

     !
     ! Initialization
     !
     sigma_gp_chm(ielem) % a  = 0.0_rp
     sigma0_gp_chm(ielem) % a = 0.0_rp
     d32_gp_chm(ielem) % a    = 0.0_rp

     do igaus = 1,pgaus

        phi_clip = min( 1.0_rp,max(0.0_rp,gpcon(igaus,iclas_phi)) )

        !
        ! Compute length scale from filter size
        !
        if (ndime == 2 ) then
           delta   = ( hleng(1) * hleng(2) )**(0.5_rp)
        else
           delta   = ( hleng(1) * hleng(2) * hleng(3) )**(0.3333333_rp)
        endif

        delta2  = delta*delta
        rdelta  = 1.0_rp / delta
        rdelta2 = 1.0_rp / delta2

        !
        ! Compute Sigma_min = alpha / Delta * sqrt (phi * (1 - phi))
        !    Chesnel et al., (2009) Int. J. Multifase Flow
        !
        sigma0_gp_chm(ielem) % a(igaus,1,1) = 2.4_rp * rdelta * sqrt( min(0.25_rp,max(0.0_rp, (phi_clip * (1.0_rp - phi_clip)) )) )

        !
        ! Compute Sigma = Sigma_min + Sigma'
        !
        sigma_gp_chm(ielem) % a(igaus,1,1) = max(0.0_rp,sigma0_gp_chm(ielem) % a(igaus,1,1) + gpcon(igaus,iclas_sigma))

        !
        ! Compute Sauter mean diameter d32 ( scaled by phi_L*(1-phi_L)/0.25 )
        !
        if ( sigma_gp_chm(ielem) % a(igaus,1,1) > 0.0_rp ) then
           d32_gp_chm(ielem) % a(igaus,1,1) = 4.0_rp * phi_clip * (1.0_rp - phi_clip)&
               * max(0.0_rp,6.0_rp * phi_clip / sigma_gp_chm(ielem) % a(igaus,1,1))
        else
           d32_gp_chm(ielem) % a(igaus,1,1) = 0.0_rp
        end if

        !
        ! LES
        !
        if(turmu_ker % kfl_exist /= 0_ip) then
           !
           ! Compute strain rate square |S_ij|*|S_ij|
           !
           gpgve = 0.0_rp
           do inode = 1,pnode
              do idime = 1,ndime
                 do jdime = 1,ndime
                    gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                 end do
              end do
           end do

           seci4 = 0.0_rp
           do idime = 1,ndime                     ! |S_ij|*|S_ij|
              do jdime = 1,ndime
                 seci4 = seci4 + 0.5_rp * gpgve(idime,jdime,igaus) * (gpgve(idime,jdime,igaus) + gpgve(jdime,idime,igaus))
              end do
           end do

           !
           ! Compute inverse of turbulent time scale (tau_b ~= tau_t): 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
           !
           rtau_turb = ( 3.24_rp * gptur(igaus) * seci4 * rdelta2 ) ** (0.33333_rp)

           !
           ! Compute turbulent kinetic energy, k = eps_sgs * tau_t
           !
           if ( rtau_turb == 0.0_rp ) then
              gpsigma(igaus) = 0.0_rp

           else
              k_sgs          = gptur(igaus) * abs( seci4 ) / rtau_turb

              !
              ! Compute production/destruction of surface area by mean and shear flow, turbulence and interactions
              !
              sigma_crit = max(0.0_rp,sigma0_gp_chm(ielem) % a(igaus,1,1) + phi_clip * (1.0_rp - phi_clip) * gpden(igaus) * k_sgs&
                  / surf_tension_chm)

              if ( sigma_crit == 0.0_rp ) then
                 gpsigma(igaus) = 0.0_rp
              else
                 gpsigma(igaus) = sigma_gp_chm(ielem) % a(igaus,1,1) * rtau_turb * (1.0_rp - sigma_gp_chm(ielem) % a(igaus,1,1)&
                     / sigma_crit)
              endif

           end if
        !
        ! Laminar
        !
        else

           gpsigma(igaus) = sigma0_gp_chm(ielem) % a(igaus,1,1)

        end if

     end do

   end subroutine chm_elmpre_spray

   subroutine chm_elmprc_spray(&
         iclas,pgaus,gptur,gpsigma,gpden,gpdif,gprhs)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_elmprc_spray
     ! NAME
     !    chm_elmprc_spray
     ! DESCRIPTION
     !    Compute spray terms for assembly in ADR equation
     ! USES
     ! USED BY
     !    chm_element_operations
     !***
     !-----------------------------------------------------------------------
     use def_kintyp, only      :  ip,rp
     use def_chemic, only      :  kfl_spray_chm
     use def_chemic, only      :  diffu_chm
     use def_chemic, only      :  nclas_chm
     use def_kermod, only      :  turmu_ker

     implicit none

     integer(ip),  intent(in)  :: iclas
     integer(ip),  intent(in)  :: pgaus
     real(rp),     intent(in)  :: gptur(pgaus)
     real(rp),     intent(in)  :: gpsigma(pgaus)
     real(rp),     intent(out) :: gpden(pgaus)
     real(rp),     intent(out) :: gpdif(pgaus,nclas_chm)
     real(rp),     intent(out) :: gprhs(pgaus)

     !
     ! Solve for non-density weighted transport
     !
     if (iclas == nclas_chm ) then
        gpden(1:pgaus) = 1.0_rp
     end if

     !
     ! LES calculation
     !
     if (turmu_ker % kfl_exist /= 0_ip .and. kfl_spray_chm == 1) then
        gpdif(1:pgaus,iclas) = gptur(1:pgaus) / diffu_chm(1,1)

     !
     ! Level set calculation or laminar
     !
     else
        gpdif(1:pgaus,iclas) = 0.0_rp
     end if

     !
     !  Add production/destruction of surface area by interactions Sigma_int to the
     !  liquid-gas interface density
     !
     if (iclas == nclas_chm ) &
         gprhs(1:pgaus) = gprhs(1:pgaus) + gpsigma(1:pgaus)


   end subroutine chm_elmprc_spray

   subroutine chm_rhodt_spray(pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,dt_rho_chm_loc,dt_chm_loc)
     !------------------------------------------------------------------------
     ! NAME
     !    chm_rhodt_spray
     ! DESCRIPTION
     !    Projection of rho/dt for explicit time step for the gas phase
     !    and 1/dt for the liquid phase
     ! USES
     ! USED BY
     !    chm_element_operations
     !***
     !------------------------------------------------------------------------

     use def_master,      only : dtinv
     use def_kintyp,      only : ip, rp
     implicit none

     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: pgaus
     integer(ip), intent(in)  :: porde
     integer(ip), intent(in)  :: lnods_loc(pnode)
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: gpvol(pgaus)
     real(rp),    intent(in)  :: gpden(pgaus)
     real(rp),    intent(out) :: dt_rho_chm_loc(*)
     real(rp),    intent(out) :: dt_chm_loc(*)

     integer(ip)             :: inode,jnode,ipoin,igaus
     real(rp)                :: eldtrho(pnode),eldt(pnode)
     real(rp)                :: fact1,fact2
     real(rp)                :: elmat1(pnode,pnode),elmat2(pnode,pnode)
     real(rp)                :: trace1,elmass1
     real(rp)                :: trace2,elmass2


     if( porde == 1 ) then
        !
        ! Element assembly
        !
        eldtrho = 0.0_rp
        eldt    = 0.0_rp
        do igaus = 1,pgaus
           fact1 = gpvol(igaus) / ( gpden(igaus)  * dtinv )
           fact2 = gpvol(igaus) / dtinv
           do inode = 1,pnode
              eldtrho(inode) = eldtrho(inode) + gpsha(inode,igaus) * fact1
              eldt(inode)    = eldt(inode)    + gpsha(inode,igaus) * fact2
           end do
        end do
        !
        ! Nodal projection
        !
        do inode = 1,pnode
           ipoin = lnods_loc(inode)
           dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + eldtrho(inode)
           dt_chm_loc(ipoin)     = dt_chm_loc(ipoin)     + eldt(inode)
        end do
     else
        !
        ! Element assembly
        !
        eldtrho = 0.0_rp
        eldt    = 0.0_rp

        elmat1 = 0.0_rp
        elmat2 = 0.0_rp

        do igaus=1,pgaus
           do inode=1,pnode
              fact1 = gpvol(igaus)*gpsha(inode,igaus)/(gpden(igaus)*dtinv)
              fact2 = gpvol(igaus)*gpsha(inode,igaus)/dtinv
              do jnode=1,pnode
                 elmat1(inode,jnode)=elmat1(inode,jnode) + fact1*gpsha(jnode,igaus)
                 elmat2(inode,jnode)=elmat2(inode,jnode) + fact2*gpsha(jnode,igaus)
              end do
           end do
        end do

        trace1  = 0.0_rp
        trace2  = 0.0_rp
        elmass1 = 0.0_rp
        elmass2 = 0.0_rp

        do inode = 1,pnode
           trace1 = trace1 + elmat1(inode,inode)
           trace2 = trace2 + elmat2(inode,inode)
           do jnode = 1,pnode
              elmass1 = elmass1 + elmat1(inode,jnode)
              elmass2 = elmass2 + elmat2(inode,jnode)
           end do
        end do

        !
        ! Nodal projection
        !
        do inode = 1,pnode
           ipoin = lnods_loc(inode)
           dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + elmat1(inode,inode)*(elmass1/trace1)
           dt_chm_loc(ipoin)     = dt_chm_loc(ipoin)     + elmat2(inode,inode)*(elmass2/trace2)
        end do
     end if


   end subroutine chm_rhodt_spray

 end module mod_chm_spray

