!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_hydrostatic.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Hydrostatic stuffs
!> @details Hydrostatic stuffs. 
!>
!>          +---------------------------+ <= HYDRO_NSI: p=0
!>          |                           |
!>          |            air            |
!>          |                           |
!>          +---------------------------+ <= HEIHY_NSI
!>          |           water           |
!>          +---------------------------+
!>
!>
!>          The hydrostatic state (p_hyd, rho_hyd) appears at three
!>          different places:            
!>
!>          1. Equation:         (...) + grad(p) = rho*g - rho_hyd*g
!>          2. Traction:         sig.n  = -p_hyd n
!>          3. Initial pressure: p_init = p_hyd
!>
!>          - The hydrostatic state (p_hyd, rho_hyd) is computed at 
!>            INIUNK. It is actualized in BEGITE if 
!>            kfl_update_hydro_nsi = ITASK_BEGITE 
!>          - The hydrostatic state can be computed analytically using
!>            the interface height HEIHY_NSI or using a PDE
!>          - The interface height HEIHY_NSI can be prescribed or computed
!>            automatically using the boundary code 12
!>          - In both cases, the zero pressure level is imposed to
!>            HYDRO_NSI
!>          - In 1), rho_hyd is saved if kfl_hydro_gravity_nsi /= 0
!>          - In 2), the hydrostatic state is imposed as a Neumann 
!>            condition for BPESS_NSI(:) when boundary code 12 is imposed
!>          - In 3), the hydrostatic state is imposed as initial 
!>            condition for the pressure PRESS(:,NPREV_NSI)
!>
!>
!-----------------------------------------------------------------------
module mod_nsi_hydrostatic
  use def_master
  use def_parame,         only : pi
  use def_elmtyp,         only : ELEXT
  use def_domain
  use def_nastin
  use def_solver
  use def_kermod,         only : thicl,kfl_prope,densi_ker
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_messages,       only : livinf
  implicit none 
  private

  public :: nsi_hydrostatic_pressure

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/10/2015
  !> @brief   Compute hydrostatic pressure
  !> @details Compute hydrostatic pressure
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_hydrostatic_pressure(itask)

    integer(ip), intent(in) :: itask
    integer(ip)             :: kfl_auxhh,dummi
    integer(ip)             :: pelty,pnode,pgaus
    integer(ip)             :: ielem,inode,ipoin
    real(rp)                :: gpsha(mnode,mgaus)         ! N
    real(rp)                :: gpder(ndime,mnode,mgaus)   ! dN/dsi          
    real(rp)                :: gpcar(ndime,mnode,mgaus)   ! dN/dxi
    real(rp)                :: gphes(ntens,mnode,mgaus)   ! d2N/dxidxj
    real(rp)                :: gpvol(mgaus)               ! |J|w
    real(rp)                :: gpden(mgaus)               ! rho
    real(rp)                :: elcod(ndime,mnode)

    if( kfl_hydro_nsi /= 0 .and. ( itask == ITASK_INIUNK .or. ( itask == ITASK_BEGITE .and. kfl_update_hydro_nsi == ITASK_BEGITE ) ) ) then

       !----------------------------------------------------------------------
       !
       ! NSI_HYDRO_DENSITY: Compute hydrostatic density
       !
       !----------------------------------------------------------------------

       if( kfl_hydro_gravity_nsi /= 0 .and. INOTMASTER ) then
          do ielem = 1,nelem
             pelty = abs(ltype(ielem))
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             do inode = 1,pnode
                elcod(1:ndime,inode) = coord(1:ndime,lnods(inode,ielem))
             end do
             call elmca2(&
                  pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                  elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpsha,&
                  gpder,gpcar,gphes,ielem)
             call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
             hydro_density_nsi(ielem) % a(1:pgaus) = gpden(1:pgaus)
          end do
       end if

       !----------------------------------------------------------------------
       !
       ! HYDRO_NSI: Look for zero pressure height 
       !
       !----------------------------------------------------------------------

       if( itask == ITASK_INIUNK ) then
          hydro_nsi = -huge(1.0_rp)
          if( INOTMASTER ) then
             if( kfl_confi_nsi == 1 ) then
                if( nodpr_nsi > 0 ) then
                   hydro_nsi = coord(ndime,nodpr_nsi)
                end if
             else
                do ipoin = 1,npoin
                   hydro_nsi = max(hydro_nsi,coord(ndime,ipoin))
                end do
             end if
          end if
          call PAR_MAX(hydro_nsi,'IN MY ZONE')
       end if

       !----------------------------------------------------------------------
       !
       ! HEIHY_NSI: Compute hydrostatic interface height
       !
       !----------------------------------------------------------------------

       if( kfl_hydro_interface_nsi == 1 ) then
          call nsi_hydrostatic_interface_height() 
       end if

       !----------------------------------------------------------------------
       !
       ! UNKNO: Compute hydrostatic pressure
       !
       !----------------------------------------------------------------------

       if( kfl_hydro_nsi == NSI_ANALYTICAL_HYDROSTATIC_PRESSURE ) then
          !
          ! Analytical hydrostatic pressure
          !
          if( INOTMASTER ) call nsi_analytical_pressure(unkno)

       else if( kfl_hydro_nsi == NSI_PDE_HYDROSTATIC_PRESSURE ) then
          !
          ! Initialize solver
          !
          solve_sol => solve(5:)
          call inisol()


          !allocate(kfl_fixno_hydro(1,npoin))
          !allocate(bvess_hydro(1,npoin))

          !
          ! PDE: Assemble matrix and RHS
          !
          kfl_auxhh     = kfl_ellen_nsi
          kfl_ellen_nsi = 1                                     ! Maximum length
          if( INOTMASTER ) call nsi_hydrostatic_elmoperations()
          kfl_ellen_nsi = kfl_auxhh                             ! Recover length
          !
          ! Solve system
          !
          call nsi_updunk(2000_ip)                              ! UNKNO <= PRESS(:,:,1)

          !nullify(solve_sol(1) % vecto_gs) 
          !allocate(solve_sol(1) % vecto_gs(ndime,npoin))
          !do ipoin = 1,npoin
          !   solve_sol(1) % vecto_gs(1,ipoin) =  0.0_rp
          !   solve_sol(1) % vecto_gs(2,ipoin) = -1.0_rp
          !end do
          !solve_sol(1) % kfl_fixno => kfl_fixno_hydro
          !solve_sol(1) % bvess     => bvess_hydro

          call solver(rhsid,unkno,amatr,pmatr)                  ! Solve system

          !deallocate(solve_sol(1) % vecto_gs)
          !call runend('O.K.!')

       else

          call runend('NSI_HYDROSTATIC_PRESSURE: NOT AVAILABLE')

       end if

       !----------------------------------------------------------------------
       !
       ! BPESS_NSI: Update hydrostatic pressure
       !
       !----------------------------------------------------------------------

       do ipoin = 1,npoin
          bpess_nsi(1,ipoin,1) = unkno(ipoin) 
       end do

       !----------------------------------------------------------------------
       !
       ! PRESS: Compute initial pressure 
       !
       !----------------------------------------------------------------------

       if( itask == ITASK_INIUNK .and. kfl_inipr_nsi == 2 ) then
          do ipoin = 1,npoin
             press(ipoin,1) = unkno(ipoin) 
          end do
       end if

    end if

  end subroutine nsi_hydrostatic_pressure

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/10/2015
  !> @brief   Compute hydrostatic interface height
  !> @details Compute hydrostatic interface height
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_hydrostatic_interface_height()

    integer(ip) :: ipoin,jpoin,iboun
    real(rp)    :: fi,fj,xx(3),t
    real(rp)    :: heihy_nsi_old

    if( associated(fleve) ) then

       heihy_nsi_old =  heihy_nsi 
       heihy_nsi     = -huge(1.0_rp) 

       if( INOTMASTER ) then
          if( ndime == 3 ) then
             call runend('NOT CODED')
          else
             do iboun = 1,nboun
                if( kfl_fixbo_nsi(iboun) == 12 ) then 
                   ipoin = lnodb(1,iboun)
                   jpoin = lnodb(2,iboun)
                   fi    = fleve(ipoin,1)
                   fj    = fleve(jpoin,1)
                   if( fi*fj < 0.0_rp ) then
                      t = -fi / ( fj - fi + zeror)
                      if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                         xx(1:ndime) = coord(1:ndime,ipoin) + t * ( coord(1:ndime,jpoin) - coord(1:ndime,ipoin) )
                      end if
                      heihy_nsi = max(heihy_nsi,xx(ndime))
                   end if
                end if
             end do
          end if
       end if
       call PAR_MAX(heihy_nsi,'IN MY CODE')
       routp = heihy_nsi
       call livinf(-16_ip,'WATER HEIGHT',0_ip)

    end if

  end subroutine nsi_hydrostatic_interface_height

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/10/2015
  !> @brief   Compute analytical hydrostatic pressure
  !> @details Compute analytical hydrostatic pressure
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_analytical_pressure(hydro_pressure)

    real(rp),    intent(out)   :: hydro_pressure(*)
    integer(ip)                :: ipoin,kount,idime,kaxhy,imate
    real(rp)                   :: z,p,rhoa,rhof,zint,pint
    !
    ! Water and air densities
    !
    imate = 1
    if( nmate > 1 ) then
       call runend('NSI_ANALYTICAL_PRESSURE: DO SOMETHING WHEN HAVING MORE THAN 1 MATERIAL')
    end if
    if( kfl_colev_nsi == 0 ) then
       rhof = densi_ker % rlaws(1,imate)
       rhoa = densi_ker % rlaws(1,imate)
    else
       rhof = densi_ker % rlaws(1,imate)
       rhoa = densi_ker % rlaws(2,imate)
    end if
    !
    ! Height
    !
    zint  = heihy_nsi
    kount = 0
    do idime = 1,ndime
       if( abs(gravi_nsi(idime)) > (10.0_rp*zeror) ) then
          kount = kount + 1
          kaxhy = idime
       end if
    end do
    if( kount /= 1 ) call runend('NSI_ANALYTICAL_PRESSURE: is only ready for gravity oriented in one of the coordinate axis') 
    !
    !
    ! Hydrostatic pressure profile
    !
    !                         z_int
    !      +--------------------.-------+--> HYDRO_NSI, P = 0 (or VALPR_NSI)
    !      |                    .      - 
    !      |      air           .    -  
    !      |                    .  -   
    !      +      interface     + ----> P_int
    !      |                   -
    !      |                  - 
    !      |      water      -  
    !      |                -
    !      |               -
    !
    if( thicl < 1.0e-10_rp ) then

       pint = gravi_nsi(kaxhy) * rhoa * grnor_nsi * (zint-hydro_nsi)

       do ipoin = 1,npoin
          z = coord(kaxhy,ipoin)
          if( z > zint ) then
             p = gravi_nsi(kaxhy) * rhoa * grnor_nsi * (z-hydro_nsi)        ! Air
          else 
             p = pint + gravi_nsi(kaxhy) * rhof * grnor_nsi * (z-zint)      ! Water
          end if
          hydro_pressure(ipoin) = p + valpr_nsi
       end do

    else
       do ipoin = 1,npoin
          z = coord(kaxhy,ipoin)
          if( z >= zint+thicl ) then
             p = gravi_nsi(kaxhy) * rhoa*grnor_nsi*(z-hydro_nsi)
          else if( (z<zint+thicl) .and. (z>zint-thicl) ) then
             p = gravi_nsi(kaxhy) * rhoa*grnor_nsi*(zint+thicl-hydro_nsi)
             p = p + gravi_nsi(kaxhy) * grnor_nsi*( (rhoa+(0.5_rp*(rhof-rhoa)*(1.0_rp+(zint/thicl))))*z  &
                  - (0.25_rp/thicl)*(rhof-rhoa)*z*z    &
                  +(0.5_rp*thicl/(pi*pi))*(rhof-rhoa)*cos((zint-z)*pi/thicl)     &
                  -(rhoa+(0.5_rp*(rhof-rhoa)*(1.0_rp+(zint/thicl))))*(zint+thicl)  &
                  + (0.25_rp/thicl)*(rhof-rhoa)*(zint+thicl)*(zint+thicl)       &
                  +(0.5_rp*thicl/(pi*pi))*(rhof-rhoa)  )
          else
             p = gravi_nsi(kaxhy) * rhoa*grnor_nsi*(zint+thicl-hydro_nsi)
             p = p + gravi_nsi(kaxhy)*grnor_nsi*( (rhoa+(0.5_rp*(rhof-rhoa)*(1.0_rp+(zint/thicl))))*(zint-thicl)  &
                  - (0.25_rp/thicl)*(rhof-rhoa)*(zint-thicl)*(zint-thicl)    &
                  +(0.5_rp*thicl/(pi*pi))*(rhof-rhoa)*cos((zint-(zint-thicl))*pi/thicl)     &
                  -(rhoa+(0.5_rp*(rhof-rhoa)*(1.0_rp+(zint/thicl))))*(zint+thicl)  &
                  + (0.25_rp/thicl)*(rhof-rhoa)*(zint+thicl)*(zint+thicl)       &
                  +(0.5_rp*thicl/(pi*pi))*(rhof-rhoa)  )
             p = p + gravi_nsi(kaxhy)*rhof*grnor_nsi*(z-(zint-thicl))
          end if
          hydro_pressure(ipoin) = p + valpr_nsi
       end do
    end if

  end subroutine nsi_analytical_pressure

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/10/2015
  !> @brief   Assemble system to solve hydrostatic pressure
  !> @details Aseemble equation 
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_hydrostatic_elmoperations()
    real(rp)    :: elmat(mnode,mnode)
    real(rp)    :: elrhs(mnode)
    real(rp)    :: elvel(ndime,mnode)
    real(rp)    :: elcod(ndime,mnode)
    real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
    real(rp)    :: gphes(ntens,mnode,mgaus)              ! d2N/dxidxj
    real(rp)    :: gpvol(mgaus)                          ! w*|J|, |J|
    real(rp)    :: gpden(mgaus)                          ! rho
    real(rp)    :: tragl(9),chave(6),chale(2),hleng(3)   ! Stabilization
    integer(ip) :: ielem,pnode,pgaus,idime,inode,ipoin
    integer(ip) :: pelty,plapl,porde,ptopo,dummi

    elements: do ielem = 1,nelem
       pelty = ltype(ielem)
       if( pelty > 0 ) then
          !
          ! Element properties and dimensions
          !
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)
          plapl = 0
          porde = lorde(pelty)
          ptopo = ltopo(pelty)
          !
          ! Gather operations: ELVEL, ELCOD
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elvel(idime,inode) = veloc(idime,ipoin,1)
                elcod(idime,inode) = coord(idime,ipoin)
             end do
          end do
          !
          ! Cartesian derivatives, Hessian matrix and volume: GPCAR, PGVOL
          !
          call elmcar(&
               pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
               elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
               gphes,ielem)
          !
          ! HLENG and TRAGL at center of gravity and CHALE
          !
          call elmlen(&
               ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
               hnatu(pelty),hleng)
          call elmchl(&
               tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
               porde,hnatu(pelty),kfl_advec_nsi,kfl_ellen_nsi)
          !
          ! Density
          !
          call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape)
          !
          ! Assemble elemental matrix
          !
          call nsi_hydrostatic_element_matrix(&
               pgaus,pnode,lnods(:,ielem),lelch(ielem),gpvol,&
               gpcar,elmar(pelty)%shape,hleng,gpden,&
               elmat,elrhs)
          !
          ! Assembly: AMATR and RHSID
          !
          call assrhs(&
               1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
          call assmat(&
               solve(5)%ndofn,pnode,pnode,solve(5)%nunkn,&
               solve(5)%kfl_algso,ielem,lnods(1,ielem),elmat,amatr)

       end if

    end do elements

  end subroutine nsi_hydrostatic_elmoperations

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/10/2015
  !> @brief   Assemble system to solve hydrostatic pressure
  !> @details 
  !>    Assemble the elemental matrix from Gauss point contributions
  !>    Solve the following advection problem with source term:
  !>
  !>         grad(p) = rho g    =>
  !>    g  . grad(p) = rho g^2  =>
  !>    g' . grad(p) = rho g
  !>
  !>    where g' is the unit gravity vector.
  !>
  !>    ( rho g - g'.grad(p) , v + tau*g'.grad(v) ) = 0
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_hydrostatic_element_matrix(&
       pgaus,pnode,lnods,lelch,gpvol,gpcar,&
       gpsha,hleng,gpden,elmat,elrhs)

    integer(ip), intent(in)  :: pgaus                    !< # Gauss points
    integer(ip), intent(in)  :: pnode                    !< # nodes
    integer(ip), intent(in)  :: lnods(pnode)             !< # element connectivity
    integer(ip), intent(in)  :: lelch                    !< # element characteristic
    real(rp),    intent(in)  :: gpvol(pgaus)             !< Element volume
    real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus) !< Shape function Cartesian derivatives
    real(rp),    intent(in)  :: gpsha(pnode,pgaus)       !< Shape function
    real(rp),    intent(in)  :: hleng(ndime)             !< Element length
    real(rp),    intent(in)  :: gpden(pgaus)             !< Density
    real(rp),    intent(out) :: elmat(pnode,pnode)       !< Matrix
    real(rp),    intent(out) :: elrhs(pnode)             !< RHS
    integer(ip)              :: igaus,inode,jnode
    integer(ip)              :: ipoin,ibopo,idime
    real(rp)                 :: fact1,fact2,tau
    real(rp)                 :: ggrap(mnode),vtest,dummr(2)

    !----------------------------------------------------------------------
    !
    ! Compute LHS and RHS
    !
    !----------------------------------------------------------------------

    elmat = 0.0_rp
    elrhs = 0.0_rp
    ggrap = 0.0_rp

    !tau = chale(1)
    tau = hleng(1)

    do igaus = 1,pgaus
       !
       ! g' . grad(Ni)
       !
       do inode = 1,pnode
          do idime = 1,ndime
             ggrap(inode) = ggrap(inode) + gravi_nsi(idime) * gpcar(idime,inode,igaus)
          end do
       end do

       fact2 = gpden(igaus) * grnor_nsi 

       do inode = 1,pnode 
          vtest        = ( gpsha(inode,igaus) + 2.0_rp*tau*ggrap(inode) )*gpvol(igaus)
          elrhs(inode) = elrhs(inode) + fact2 * vtest
          do jnode = 1,pnode
             elmat(inode,jnode) = elmat(inode,jnode) + vtest * ggrap(jnode)
          end do
       end do

    end do

    !----------------------------------------------------------------------
    !
    ! Extension elements
    !
    !----------------------------------------------------------------------

    if( lelch == ELEXT ) then
       call elmext(&
            4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
            dummr,elrhs,dummr)
    end if

    !----------------------------------------------------------------------
    !
    ! Prescribe boundary conditions at z-plane z=hydro_nsi
    !
    !----------------------------------------------------------------------

    do inode = 1,pnode
       ipoin = lnods(inode)
       ibopo = lpoty(ipoin)
       if(ibopo /= 0) then
          if( abs(coord(ndime,ipoin)-hydro_nsi) < 1.0e-10_rp ) then

             fact1 = elmat(inode,inode)
             if( abs(fact1) < 1.0e-10_rp ) fact1 = 1.0_rp
             elrhs(inode) = 0.0_rp
             do jnode = 1,pnode
                elmat(jnode,inode) = 0.0_rp
                elmat(inode,jnode) = 0.0_rp              
             end do
             elmat(inode,inode) = fact1
          end if
       end if
    end do

  end subroutine nsi_hydrostatic_element_matrix

end module mod_nsi_hydrostatic

!> @} 
