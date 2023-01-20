!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_efect
!> @author  David Oks 
!> @date    2022-01-12
!> @brief   Embedded Finite Element Coupling Technique (EFECT) for FSI
!> @details Embedded Finite Element Coupling Technique (EFECT) for FSI
!> @}
!-----------------------------------------------------------------------

module mod_nsi_efect

  use def_kintyp,         only :  ip, rp
  use def_master,         only :  INOTMASTER, veloc, press
  use def_master,         only :  kfl_paral
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use def_coupli,         only :  coupling_type
  use def_nastin,         only :  kfl_fixno_nsi
  use def_domain

  implicit none
  private

  public                                            :: &
       nsi_compute_penalty_coefficient,                &
       nsi_compute_penalty_force,                      &
       nsi_element_penalty_force,                      &
       nsi_cauchy_stress_tensor,                       &
       nsi_set_fringe_node_bcs,                        &
       nsi_extrapolate_velocity_from_nearest_neighbor, &
       nsi_algebraic_reaction_force,                   &
       nsi_normalize_indicator_field


contains

  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Penalty coefficients for PIFEM
  !> @details Compute penalty coefficients from LHS diagonal
  !>          Kp_{a,i} = LHS_{a,i,a,i} for a=1,...,npoin; i=1,...,ndime
  !
  !-----------------------------------------------------------------------
  subroutine nsi_compute_penalty_coefficient(Auu, kpena)
  
  
    implicit none
  
    real(rp), intent(in)    :: Auu(ndime,ndime,*)
    real(rp), intent(inout) :: kpena(ndime,*)
    integer(ip)             :: ipoin, jpoin, idime, izdom
  
    ! Extract diagonal of Auu LHS matrix
    if (INOTMASTER) then
       kpena(1:ndime,1:npoin) = 0.0_rp
       do ipoin = 1,npoin
          do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
             jpoin  = c_dom(izdom)
             if (ipoin == jpoin) then
                do idime = 1,ndime
                   if (kfl_fixno_nsi(idime,ipoin) /= 1) then
                      kpena(idime,ipoin) = Auu(idime,idime,izdom)
                   end if
                end do
             end if
          end do
       end do
    end if
  
    ! Sum values at interface nodes
    call PAR_INTERFACE_NODE_EXCHANGE(ndime,kpena,'SUM','IN MY CODE')
  
  end subroutine nsi_compute_penalty_coefficient
  
  
  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Penalty force for PIFEM coupling
  !> @details Compute velocity penalty force for PIFEM coupling
  !>          f_p = k * (u^f - u^s);  k = alpha / h
  !>          where:  u^f = fluid velocity
  !>                  u^s = solid velocity
  !
  !-----------------------------------------------------------------------
  subroutine nsi_compute_penalty_force(force, icoup, kpena, indic, dvelo)
  
    implicit none
  
    real(rp),    intent(out) :: force(ndime,npoin)
    integer(ip), intent(in)  :: icoup
    real(rp),    intent(in)  :: kpena(ndime,*)
    real(rp),    intent(in)  :: indic(npoin)
    real(rp),    intent(in)  :: dvelo(ndime,*)
    real(rp),    parameter   :: alpha = 0.3_rp
    integer(ip)              :: pelty, pnode, pdime, pgaus, plapl, porde
    integer(ip)              :: ielem

    ! Integrate force as source term
    !
    ! b_u += M_ij * alpha / h_j * (u^f_j - u^s_j)
    !
    if ( INOTMASTER ) then
  
       ! Initialize force to zero
       force(1:ndime,1:npoin) = 0.0_rp

       ! Element loop
       do ielem=1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pdime = ldime(pelty)
              pgaus = ngaus(pelty)
              plapl = llapl(pelty)
              porde = lorde(pelty)
          end if

          ! Compute element contribution to penalty force
          call nsi_element_penalty_force(&
               force, kpena, indic, dvelo, pelty, pnode, pdime, pgaus, plapl, porde, ielem)
  
       end do
    end if
  
    call PAR_INTERFACE_NODE_EXCHANGE(ndime,force,'SUM','IN MY CODE')
  
  end subroutine nsi_compute_penalty_force
  
  
  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Element contribution to penalty force 
  !> @details Compute element contribution to penaty force for PIFEM
  !>          coupling. 
  !
  !-----------------------------------------------------------------------
  subroutine nsi_element_penalty_force(force, kpena, indic, dvelo, pelty, pnode, pdime, pgaus, plapl, porde, ielem)
  
    use mod_elmgeo,              only : element_type
    use mod_elmgeo,              only : elmgeo_element_length
    use mod_elmgeo,              only : elmgeo_element_characteristic_length
    use mod_element_integration, only : element_shape_function_derivatives_jacobian
  
    implicit none
  
    real(rp),    intent(inout) :: force(ndime,npoin)
    real(rp),    intent(in)    :: kpena(ndime,*)
    real(rp),    intent(in)    :: indic(npoin)
    real(rp),    intent(in)    :: dvelo(ndime,*)
    integer(ip), intent(in)    :: pelty
    integer(ip), intent(in)    :: pnode 
    integer(ip), intent(in)    :: pdime
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: plapl
    integer(ip), intent(in)    :: porde 
    integer(ip), intent(in)    :: ielem
    real(rp)                   :: elcod(ndime,pnode)
    real(rp)                   :: gpvol(pgaus)
    real(rp)                   :: gpsha(pnode,pgaus)
    real(rp)                   :: gpder(ndime,pnode,pgaus)
    real(rp)                   :: gpcar(ndime,mnode,pgaus)
    real(rp)                   :: gphes(ntens,mnode,pgaus)
    real(rp)                   :: xjacm(ndime,ndime)
    real(rp)                   :: elforce(ndime,pnode)
    real(rp)                   :: eldvelo(ndime,pnode)
    real(rp)                   :: elkpena(ndime,pnode)
    real(rp)                   :: elindic(pnode)
    real(rp)                   :: gpkpena(ndime,pgaus)
    real(rp)                   :: gpindic(pgaus)
    real(rp)                   :: tragl(ndime,ndime) 
    real(rp)                   :: hleng(ndime) 
    real(rp)                   :: chave(ndime,2) 
    real(rp)                   :: chale(2) 
    real(rp)                   :: gpdet
    real(rp)                   :: gpmas
    real(rp),    parameter     :: alpha = 150.0_rp
    integer(ip)                :: igaus, ipoin, inode, jnode, idime

    ! Only add contribution if it's not a "hole" element
    if(ltype(ielem) > 0_ip) then

       ! Gather
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          eldvelo(1:ndime,inode) = dvelo(1:ndime,ipoin)
          elindic(inode)         = indic(ipoin)
          elkpena(1:ndime,inode) = kpena(1:ndime,ipoin)
          elcod(1:ndime,inode)   = coord(1:ndime,ipoin)
       end do

       ! Get shape functions
       call element_shape_function_derivatives_jacobian(&
            pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
            elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
            gpder,gpcar,gphes,ielem)

       ! Element characteristic length
       hleng = 0.0_rp
       chave = 0.0_rp
       chale = 0.0_rp
       call elmgeo_element_characteristic_length(&
            ndime,pnode,elmar(pelty) % dercg(:,:),elcod,hleng,element_type(pelty) % natural_length,tragl)
       call elmgeo_element_length(&
            ndime,pnode,porde,tragl,hleng,elcod,eldvelo,chave,chale,&
            element_type(pelty) % natural_length,0_ip,2_ip)

       ! Interpolate to gauss points
       gpindic = 0.0_rp
       gpkpena = 0.0_rp
       do igaus = 1,pgaus
          do inode = 1,pnode
             gpindic(igaus) = gpindic(igaus) &
                            + gpsha(inode,igaus) * elindic(inode)
             do idime = 1,ndime
                gpkpena(idime,igaus) = gpkpena(idime,igaus) &
                                     + gpsha(inode,igaus) * elkpena(idime,inode)
             end do
          end do
       end do

       ! Assemble element contribution to penalty force
       elforce = 0.0_rp
       do igaus = 1,pgaus
          ! Get jacobian determinant of element
          call jacdet(&
               pdime,pnode,elcod,elmar(pelty) % deriv(1,1,igaus),&
               xjacm,gpdet)

          do inode = 1,pnode
             do jnode = 1,pnode
                gpmas = gpsha(inode,igaus) * gpsha(jnode,igaus) * gpvol(igaus)
                do idime = 1,ndime
                   elforce(idime,inode) = elforce(idime,inode) &
                                        + alpha &
                                        * gpmas &
                                        * eldvelo(idime,jnode) &
                                        * gpindic(igaus)
     !                                   * gpkpena(idime,igaus) &
     !                                   / chale(1)
                end do
             end do
          end do
       end do

       ! Scatter
       do inode=1,pnode
          ipoin = lnods(inode,ielem)
          force(1:ndime,ipoin) = force(1:ndime,ipoin) + elforce(1:ndime,inode)
       end do

    end if
  
  end subroutine nsi_element_penalty_force
  
  
  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Compute fluid Cauchy stress tensor
  !> @details Compute fluid Cauchy stress tensor
  !>          sigma_ij = -p*delta_ij + mu*(du_i/dx_j + du_jj/dx_i)
  !
  !-----------------------------------------------------------------------
  subroutine nsi_cauchy_stress_tensor( stres )
    
    use mod_memory,     only : memory_alloca, memory_deallo
    use def_master,     only : modul, mem_modul, ITER_K
    use mod_ker_proper, only : ker_proper
    use mod_gradie
  
    implicit none
  
    real(rp), intent(inout) :: stres(ntens,*)
    real(rp), pointer       :: gradv(:,:)
    real(rp), pointer       :: viscd(:)
    integer(ip)             :: ltens(ndime,ndime) 
    integer(ip)             :: ipoin, itens, idime, dummi
  
    nullify( gradv )
    nullify( viscd )
  
    if ( INOTMASTER ) then
       !
       ! Get tensor indices
       !
       if ( ndime == 2 ) then
          ltens(1,1) = 1
          ltens(2,2) = 2
          ltens(1,2) = 3
          ltens(2,1) = 3
       else if ( ndime == 3 ) then
          ltens(1,1) = 1
          ltens(2,2) = 2
          ltens(1,2) = 3
          ltens(2,1) = 3
          ltens(3,3) = 4
          ltens(1,3) = 5
          ltens(3,1) = 5
          ltens(2,3) = 6
          ltens(3,2) = 6
       end if
       !
       ! Allocate memory
       !
       call memory_alloca(mem_modul(1:2,modul),'GRADV','nsi_plugin',gradv,ntens,npoin)
       call memory_alloca(mem_modul(1:2,modul),'VISCD','nsi_plugin',viscd,npoin)
       !
       ! Compute velocity gradients at nodes
       !
       call gradie(veloc(1:ndime,1:npoin,1),gradv)
       !
       ! Get dynamic viscosity at nodes
       !
       call ker_proper('VISCO','NPOIN',dummi,dummi,viscd)
       !
       ! Compute Cauchy stress tensor at nodes
       !
       stres(1:ntens,1:npoin) = 0.0_rp
       do ipoin = 1,npoin
          !
          ! Pressure stresses (only isotropic componnets of tensor)
          !
          do idime = 1,ndime
             itens = ltens(idime,idime)
             stres(itens,ipoin) = stres(itens,ipoin) - press(ipoin,ITER_K)
          end do
          !
          ! Viscous stresses
          !
          stres(1:ntens,ipoin) = stres(1:ntens,ipoin) + viscd(ipoin)*gradv(1:ntens,ipoin)
       end do
  
    end if
  
    if(associated(gradv)) call memory_deallo(mem_modul(1:2,modul),'GRADV','nsi_plugin',gradv)
    if(associated(viscd)) call memory_deallo(mem_modul(1:2,modul),'VISCD','nsi_plugin',viscd)
  
  end subroutine nsi_cauchy_stress_tensor
  
  
  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Reset node fixities to initial state
  !> @details Reset node fixities to initial state
  !
  !-----------------------------------------------------------------------
  subroutine nsi_reset_node_fixities()

    implicit none

    integer(ip)         ::  ipoin, idime 

    do ipoin = 1,npoin
       do idime = 1,ndime
          if ( kfl_fixno_nsi(idime,ipoin) == 99 ) then
             kfl_fixno_nsi(idime,ipoin) = 0
          end if
       end do 
    end do

  end subroutine nsi_reset_node_fixities


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Extrapolate velocity from nearest neighbor
  !> @details Set previous time-step velocity of freshly uncovered node as
  !>          velocity of nearest "dry" neighbor node. "Dry" means that
  !>          it's not a wetnode in the immersed-type couplings.
  !
  !-----------------------------------------------------------------------
  subroutine nsi_extrapolate_velocity_from_nearest_neighbor(icoup)

    use def_master,    only :  iblok, TIME_N
    use def_coupli,    only :  coupling_driver_iteration
    use mod_immersed,  only :  cou_find_nearest_neighbor

    implicit none

    integer(ip), intent(in) :: icoup
    integer(ip)             :: ipoin, jpoin, kpoin

    if (coupling_driver_iteration(iblok) == 1_ip) then
       do kpoin = 1,coupling_type(icoup) % wet % npoin_fresh
          ipoin = coupling_type(icoup) % wet % lpoin_fresh(kpoin)

          ! Find nearest dry node
          call cou_find_nearest_neighbor(jpoin, ipoin, icoup, .false.)

          ! Copy its last timestep velocity
          veloc(1:ndime,ipoin,TIME_N) = veloc(1:ndime,jpoin,TIME_N)
       end do
    end if

  end subroutine nsi_extrapolate_velocity_from_nearest_neighbor


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Set boundary conditions of fringe nodes
  !> @details Set boundary conditions of fringe nodes
  !
  !-----------------------------------------------------------------------
  subroutine nsi_set_fringe_node_bcs(icoup, velin)

    use mod_memory, only : memory_alloca, memory_deallo
    use def_master, only : modul, mem_modul

    implicit none

    integer(ip),          intent(in) :: icoup
    real(rp),   optional, intent(in) :: velin(ndime,npoin)
    real(rp),   pointer              :: veloc_fix(:,:)

    nullify(veloc_fix)

    call memory_alloca(mem_modul(1:2,modul),'VELOC_FIX','mod_nsi_efect',veloc_fix,ndime,max(1_ip,npoin))

    if( present(velin) ) then
       veloc_fix(1:ndime,1:npoin) = velin(1:ndime,1:npoin)
    else
       veloc_fix = 0.0_rp
    end if

    if ( INOTMASTER ) then
       !
       ! Reset node fixities
       !
       call nsi_reset_node_fixities()
       !
       ! Map velocities on uncovered wet nodes 
       !
       !! call nsi_map_velocity_in_uncovered_wetnodes(icoup, veloc_fix) ! ### debugging ###
       !
       ! Map velocities in interior wet nodes 
       !
       call nsi_map_velocity_in_interior_wetnodes(icoup, veloc_fix)
       !
       ! Set fixity of fringe nodes to 1 and Dirichlet value to zero.
       !
       call nsi_set_fringe_node_fixities(icoup, veloc_fix)
       !
    end if

    if( associated(veloc_fix) )  call memory_deallo(mem_modul(1:2,modul),'VELOC_FIX','mod_nsi_efect',veloc_fix)

  end subroutine nsi_set_fringe_node_bcs


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Set fixity of fringe nodes
  !> @details Set fixity of fringe nodes
  !
  !-----------------------------------------------------------------------
  subroutine nsi_set_fringe_node_fixities(icoup, velin)

    use def_master, only : ITER_K
    use def_nastin, only : bvess_nsi

    implicit none

    integer(ip), intent(in) :: icoup
    real(rp),    intent(in) :: velin(ndime,npoin)
    integer(ip)             :: ipoin, idime

    do ipoin = 1,npoin
       if (coupling_type(icoup) % wet % kfl_fringe_wetnodes(ipoin) == 1_ip) then
          do idime = 1,ndime
             if (kfl_fixno_nsi(idime,ipoin) == 0_ip) then
                kfl_fixno_nsi(idime,ipoin)    = 99_ip
                bvess_nsi(idime,ipoin,ITER_K) = velin(idime,ipoin)
                veloc(idime,ipoin,ITER_K)     = velin(idime,ipoin)
             end if
          end do
       end if
    end do

  end subroutine nsi_set_fringe_node_fixities


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Algebraic momentum reaction forces
  !> @details Compute algebraic momentum reaction forces
  !
  !-----------------------------------------------------------------------
  subroutine nsi_algebraic_reaction_force(icoup, ru)

    use def_nastin,               only :  intfo_nsi
    use mod_nsi_algebraic_forces, only :  nsi_algebraic_forces_allocate, nsi_algebraic_forces

    implicit none

    integer(ip), intent(in)    :: icoup
    real(rp),    intent(inout) :: ru(ndime,*)
    integer(ip)                :: ipoin, jpoin, izdom, jzdom, idime, jdime, jcoup

    ! Get index of mirror coupling
    jcoup = coupling_type(icoup) % mirror_coupling
    !
    ! Compute algebraic reaction forces on Dirichlet nodes
    !
    if (INOTMASTER) then
       ru(1:ndime,1:npoin) = 0.0_rp
       do ipoin = 1,npoin
          do idime = 1,ndime
             if (kfl_fixno_nsi(idime,ipoin) > 0) then
                !
                ! Equal residual to RHS: ru = bu
                !
                ru(idime,ipoin) = intfo_nsi(ipoin) % bu(idime)
                jzdom = 0
                do izdom = r_dom(ipoin), r_dom(ipoin+1) - 1
                   jpoin = c_dom(izdom)
                   !
                   ! Only use data from neighbors which are free
                   !
                   if (coupling_type(jcoup) % wet % kfl_fringe_wetnodes(ipoin) > 0_ip) then
                      jzdom = jzdom + 1
                      !
                      ! Compute residual: ru = ru - Aup*p - Auu*u
                      !
                      do jdime = 1,ndime
                         ru(idime,ipoin) = ru(idime,ipoin) - intfo_nsi(ipoin)%Auu(jdime,idime,jzdom)*veloc(jdime,jpoin,1)
                      end do
                      ru(idime,ipoin) = ru(idime,ipoin) - intfo_nsi(ipoin)%Aup(idime,jzdom)*press(jpoin,1)
                   end if
                end do
             end if
          end do
       end do
       !
       ! Assemble contributions to global vector: ru
       !
       call rhsmod(ndime, ru)
    end if

  end subroutine nsi_algebraic_reaction_force


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Normalize indicator field
  !> @details Normalize indicator field used as smoothly varying levelset
  !>          to identify where overlapping solid is located. The
  !>          normalization is computed by dividing by the max value.
  !
  !-----------------------------------------------------------------------
  subroutine nsi_normalize_indicator_field(indic)

    use mod_communications, only :  PAR_MAX

    implicit none

    real(rp), intent(inout) :: indic(*)
    real(rp)                :: max_indicator

    max_indicator = 0.0_rp
    if ( INOTMASTER ) max_indicator = maxval( indic(1:npoin) )
    call PAR_MAX(max_indicator, 'IN MY CODE')
    if ( INOTMASTER ) indic(1:npoin) = indic(1:npoin) / max_indicator

  end subroutine nsi_normalize_indicator_field


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Map solid velocity to uncovered fluid wet nodes 
  !> @details Map solid velocity to last timestep velocity of uncovered
  !>          fluid wet nodes so that they have a historic value in order
  !>          to compute acceleration when uncovered
  !
  !-----------------------------------------------------------------------
  subroutine nsi_map_velocity_in_uncovered_wetnodes(icoup, velin)

    use def_master,    only :  iblok, TIME_N, ITER_K
    use def_coupli,    only :  coupling_driver_iteration

    implicit none

    integer(ip), intent(in) :: icoup
    real(rp),    intent(in) :: velin(ndime,npoin )
    integer(ip)             :: ipoin, kpoin

    if (coupling_driver_iteration(iblok) == 1_ip) then
       do kpoin = 1,coupling_type(icoup) % wet % npoin_fresh
          ipoin = coupling_type(icoup) % wet % lpoin_fresh(kpoin)
          !
          ! Set last time-step velocity to that interpolated from solid
          !
          veloc(1:ndime,ipoin,ITER_K) = velin(1:ndime,ipoin)
       end do
    end if

  end subroutine nsi_map_velocity_in_uncovered_wetnodes


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Map solid velocity to interior wet nodes
  !> @details Map solid velocity to interior wet nodes
  !
  !-----------------------------------------------------------------------
  subroutine nsi_map_velocity_in_interior_wetnodes(icoup, velin)

    use def_master,    only :  iblok, TIME_N, ITER_K
    use def_coupli,    only :  coupling_driver_iteration

    implicit none

    integer(ip), intent(in) :: icoup
    real(rp),    intent(in) :: velin(ndime,npoin)
    integer(ip)             :: ipoin

    if (coupling_driver_iteration(iblok) == 1_ip) then
       do ipoin = 1,npoin
          if (abs(coupling_type(icoup) % wet % kfl_fringe_wetnodes(ipoin)) == 1_ip) then
             !
             ! Set last-timestep velocity to that interpolated from solid
             !
             veloc(1:ndime,ipoin,TIME_N) = velin(1:ndime,ipoin)
          end if
       end do
    end if

  end subroutine nsi_map_velocity_in_interior_wetnodes


!  !-----------------------------------------------------------------------
!  !
!  !> @author  David Oks
!  !> @brief   Cancel stress tensor on overlapping domain
!  !> @details Set Cauchy stress tensor to zero on fluid nodes overlapped
!  !>          by an immersed domain. This is used in the PIFEM coupling.
!  !
!  !-----------------------------------------------------------------------
!  subroutine nsi_cancel_in_overlapping_region(icoup, ntens, varia)
!  
!    implicit none
!  
!    integer(ip), intent(in)    :: icoup
!    integer(ip), intent(in)    :: ntens
!    real(rp),    intent(inout) :: varia(ntens,*)
!    integer(ip)                :: ipoin, jcoup
!  
!    ! Get mirror coupling index
!    jcoup = coupling_type(icoup) % mirror_coupling
!  
!    if ( INOTMASTER ) then
!       !
!       ! Set Cauchy to zero in interior of overlapping region
!       !
!       do ipoin = 1,npoin
!          if (coupling_type(jcoup) % wet % kfl_fringe_wetnodes(ipoin) < 0_ip) then
!             varia(1:ntens,ipoin) = 0.0_rp
!          end if
!       end do
!  
!    end if
!  
!  end subroutine nsi_cancel_in_overlapping_region
  
  
!  !-----------------------------------------------------------------------
!  !
!  !> @author  David Oks
!  !> @brief   Compute tractions from Cauchy stress tensor
!  !> @details Compute boundary tractions by integrating Cauchy stress
!  !>          tensor along boundary, this is to compare with the algebraic
!  !>          form used in the ALEFOR FSI coupling.
!  !
!  !-----------------------------------------------------------------------
!  subroutine nsi_compute_tractions_from_cauchy(stres, tract)
!  
!    use mod_bouder
!  
!    implicit none
!  
!    real(rp),    intent(in)  :: stres(ntens,*)
!    real(rp),    intent(out) :: tract(ndime,*)
!    real(rp)                 :: baloc(ndime,ndime)
!    real(rp)                 :: bocod(ndime,mnodb)
!    real(rp)                 :: elcod(ndime,mnode)
!    real(rp)                 :: streb(ndime,ndime)
!    real(rp)                 :: tractb(ndime)
!    real(rp)                 :: eucta
!    integer(ip)              :: ltens(ndime,ndime) 
!    integer(ip)              :: pblty, pnodb, pgaub, ielem, pelty, pnode 
!    integer(ip)              :: iboun, igaub, inodb, itens, inode, icoup 
!    integer(ip)              :: ipoin, kpoin, idime, jdime, imate, ibopo
!  
!    tract(1:ndime,1:npoin) = 0.0_rp
!  
!    if( INOTMASTER ) then
!       !
!       ! Translate tensor to matrix indices
!       !
!       if ( ndime == 2 ) then
!          ltens(1,1) = 1_ip 
!          ltens(2,2) = 2_ip 
!          ltens(1,2) = 3_ip
!          ltens(2,1) = 3_ip
!       else if ( ndime == 3 ) then
!          ltens(1,1) = 1_ip
!          ltens(2,2) = 2_ip
!          ltens(1,2) = 3_ip
!          ltens(2,1) = 3_ip
!          ltens(3,3) = 4_ip
!          ltens(1,3) = 5_ip
!          ltens(3,1) = 5_ip
!          ltens(2,3) = 6_ip
!          ltens(3,2) = 6_ip
!       end if
!       !
!       ! Compute fluid tractions at solid boundary 
!       !
!       boundaries: do iboun = 1,nboun
!  
!          pblty = ltypb(iboun)
!          pnodb = nnode(pblty)
!          pgaub = ngaus(pblty)
!          ielem = lelbo(iboun)
!          pelty = ltype(ielem)
!          pnode = nnode(pelty)
!          streb = 0.0_rp
!  
!          gauss_points: do igaub = 1,pgaub
!             !
!             ! Interpolate stresses to solid boundary
!             !
!             do inodb = 1,pnodb
!                ipoin = lnodb(inodb,iboun)
!                do idime = 1,ndime
!                   bocod(idime,inodb) = coord(idime,ipoin)
!                   do jdime = 1,ndime
!                      itens = ltens(idime,jdime)
!                      streb(idime,jdime) = streb(idime,jdime) &
!                                         + stres(itens,ipoin) &
!                                         * elmar(pblty) % shape(inodb,igaub)
!                   end do
!                end do 
!             end do
!             !
!             ! Compute normal to surface at gauss point
!             !
!             do inode = 1,pnode
!                ipoin = lnods(inode,ielem)
!                do idime = 1,ndime
!                   elcod(idime,inode) = coord(idime,ipoin)
!                end do
!             end do
!             call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)
!             call chenor(pnode,baloc,bocod,elcod)
!             !
!             ! Calculate tractions ( t = sigma . n ) at gauss points
!             !
!             tractb = 0.0_rp
!             do idime = 1,ndime
!                do jdime = 1,ndime
!                   tractb(idime) = tractb(idime) - streb(idime,jdime) * baloc(jdime,ndime)
!                end do
!             end do
!             !
!             ! Spread tractions to nodes 
!             !
!             do inodb = 1,pnodb
!                ipoin = lnodb(inodb,iboun)
!                do idime = 1,ndime
!                   tract(idime,ipoin) = tract(idime,ipoin) &
!                                      + tractb(idime) * elmar(pblty) % shape(inodb,igaub)
!                end do
!             end do
!  
!          end do gauss_points
!  
!       end do boundaries
!  
!    end if
!  
!  end subroutine nsi_compute_tractions_from_cauchy
  
  
end module mod_nsi_efect
!> @}
