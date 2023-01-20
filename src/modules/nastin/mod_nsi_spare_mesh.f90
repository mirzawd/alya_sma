!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_spare_mesh_source.f90
!> @author  houzeaux
!> @date    2022-01-26
!> @brief   Source term
!> @details Integrate a source term coming from a spare mesh
!-----------------------------------------------------------------------

module mod_nsi_spare_mesh

  use def_kintyp,                        only : ip,rp
  use def_domain,                        only : mnode,mgaus,ndime
  use def_domain,                        only : meshe
  use def_kermod,                        only : num_spare_meshes,ndivi
  use def_kermod,                        only : spare_meshes
  use def_kintyp_mesh_basic,             only : mesh_type_basic
  use mod_elmgeo,                        only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo,                        only : element_type
  use def_master,                        only : INOTMASTER
  use def_master,                        only : amatr,rhsid,solve
  use def_master,                        only : mem_modul,modul
  use def_domain,                        only : c_dom,r_dom,nzdom
  use mod_memory_basic,                  only : memory_alloca
  use mod_memory_basic,                  only : memory_deallo
  use mod_solver,                        only : solver_diagonal
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_kintyp_spare_mesh,             only : typ_spare_mesh
  use def_nastin

  implicit none
  private
  
  real(rp),         allocatable :: shapf(:)
  real(rp),         allocatable :: elcod(:,:)
  real(rp),         allocatable :: gpcar(:,:,:)
  real(rp),         allocatable :: xjaci(:,:,:)
  real(rp),         allocatable :: gpvol(:)
  real(rp),         allocatable :: gpdet(:)
    
  public :: nsi_spare_mesh_dirichlet

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Loop over elements
  !> @details Impose source along spare mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_spare_mesh_dirichlet()

    integer(ip) :: imesh,pdime
    integer(ip) :: ii,mn,md,mg

    if( INOTMASTER .and. num_spare_meshes > 0 .and. associated(tscod_nsi) ) then
       
       if( associated(tscod_nsi(1) % l) ) then

          do ii = 1,size(tscod_nsi(1) % l)

             imesh =  tscod_nsi(1) % l(ii) % lcode(1)

             if( imesh <= num_spare_meshes ) then
                pdime =  maxval(element_type(spare_meshes(imesh) % mesh % ltype(:)) % dimensions)
                mg    =  maxval(spare_meshes(imesh) % mesh % quad(:) % ngaus)
                mn    =  maxval(element_type(spare_meshes(imesh) % mesh % ltype(:)) % number_nodes)
                md    =  pdime
                allocate(shapf(mn))
                allocate(elcod(md,mn))
                allocate(gpcar(md,mn,mg))
                allocate(xjaci(md,md,mg))
                allocate(gpvol(mg))
                allocate(gpdet(mg))

                select case ( pdime )
                case ( 1_ip ) ; call nsi_spare_mesh_dirichlet_1D(spare_meshes(imesh),spare_meshes(imesh) % mesh,tscod_nsi(1) % l(1) % bvess)
                end select

                if( allocated(shapf) ) deallocate(shapf)
                if( allocated(elcod) ) deallocate(elcod)
                if( allocated(gpcar) ) deallocate(gpcar)
                if( allocated(xjaci) ) deallocate(xjaci)
                if( allocated(gpvol) ) deallocate(gpvol)
                if( allocated(gpdet) ) deallocate(gpdet)

             end if

          end do

       end if
    end if

  end subroutine nsi_spare_mesh_dirichlet

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Source is imposed along a 1D mesh
  !> @details Source is imposed along a 1D mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_spare_mesh_dirichlet_1D(spare,mesh,bvess)

    type(typ_spare_mesh),  intent(in)         :: spare
    type(mesh_type_basic), intent(in)         :: mesh
    real(rp),              intent(in)         :: bvess(:)
    integer(ip)                               :: pnode,jnode,igaus,jelem
    integer(ip)                               :: pgaus,pdime,jpoin,pelty
    integer(ip)                               :: ipoin1,ipoin2,ielem,jdime,nz
    integer(ip)                               :: knode,kpoin,kcolu,izdom,np,nd
    real(rp)                                  :: xfact
    integer(ip),                  pointer     :: fringe(:)
    real(rp),         contiguous, pointer     :: Auu(:,:,:)
    real(rp),         contiguous, pointer     :: Aup(:,:)
    real(rp),         contiguous, pointer     :: App(:)
    real(rp),         contiguous, pointer     :: Apu(:,:)
    real(rp),         contiguous, pointer     :: bu(:,:) 
    real(rp),                     pointer     :: diag(:,:)

    nullify(diag,fringe)
    !
    ! Diagonal
    !
    call memory_alloca(mem_modul(1:2,modul),'DIAG','nsi_spare_mesh_dirichlet_1D',diag,meshe(ndivi) % ndime,meshe(ndivi) % npoin)
    call solver_diagonal(solve,amatr(poauu_nsi:),diag)
    !
    ! Initialize
    !
    pdime =  1
    np    =  meshe(ndivi) % npoin
    nd    =  meshe(ndivi) % ndime
    nz    =  nzdom

    Auu(1:nd,1:nd,1:nz) => amatr(poauu_nsi:)
    Aup(1:nd,1:nz)      => amatr(poaup_nsi:)
    Apu(1:nd,1:nz)      => amatr(poapu_nsi:)
    App(1:nz)           => amatr(poapp_nsi:)
    bu(1:nd,1:np)       => rhsid(1:)
    !
    ! Penalty or weak method
    !
    call memory_alloca(mem_modul(1:2,modul),'FRINGE','nsi_spare_mesh_dirichlet_1D',fringe,meshe(ndivi) % npoin)
    
    if( kfl_immersed_nsi == 0 ) then
       xfact = fact_immersed_nsi
    else
       xfact = 1.0_rp
    end if

    do ielem = 1,mesh % nelem
       if( associated(spare % eleme(ielem) % l) ) then            
          pelty = mesh % ltype(ielem)
          pgaus = mesh % quad(pelty) % ngaus
          do igaus = 1,pgaus           
             jelem = spare % eleme(ielem) % l(igaus)
             if( jelem /= 0 ) then                   
                do jnode = 1,meshe(ndivi) % lnnod(jelem)
                   jpoin = meshe(ndivi) % lnods(jnode,jelem)
                   fringe(jpoin) = 1
                end do
             end if
          end do
       end if
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(fringe,'MAX')
    do jpoin = 1,meshe(ndivi) % npoin
       if( spare % dista(jpoin) < 0.0_rp .and. fringe(jpoin) /= 1 ) then
          fringe(jpoin) = -1
       end if
    end do
    do jpoin = 1,meshe(ndivi) % npoin
       if( fringe(jpoin) == 1 ) then
          if( kfl_immersed_nsi == 1 ) then
             bu(:,jpoin) = 0.0_rp
             do izdom = r_dom(jpoin),r_dom(jpoin+1)-1
                Auu(:,:,izdom) = 0.0_rp
                Aup(:,izdom)   = 0.0_rp
             end do
          end if
          do izdom = r_dom(jpoin),r_dom(jpoin+1)-1
             if( fringe(c_dom(izdom)) == -1 ) then
          !      Auu(:,:,izdom) = 0.0_rp
          !      Aup(:,izdom)   = 0.0_rp
          !      Apu(:,izdom)   = 0.0_rp
          !      App(izdom)     = 0.0_rp
             end if
          end do
       end if
    end do
    call memory_deallo(mem_modul(1:2,modul),'FRINGE','nsi_spare_mesh_dirichlet_1D',fringe)
    !
    ! Assemble source term in HEAT_SOURCE
    !
    do ielem = 1,mesh % nelem

       if( associated(spare % eleme(ielem) % l) ) then

          pelty      = mesh % ltype(ielem)
          pgaus      = mesh % quad(pelty) % ngaus
          pnode      = element_type(pelty) % number_nodes
          ipoin1     = mesh % lnods(1,ielem)
          ipoin2     = mesh % lnods(2,ielem)
          elcod(1,1) = 0.0_rp
          elcod(1,2) = sqrt(dot_product(&
               mesh % coord(:,ipoin1)-mesh % coord(:,ipoin2),&
               mesh % coord(:,ipoin1)-mesh % coord(:,ipoin2)))

          call elmgeo_cartesian_derivatives_jacobian(pdime,pnode,pnode,pgaus,&
               elcod,mesh % iso(pelty) % deriv,xjaci,gpcar,gpdet)
          gpvol(1:pgaus) = mesh % quad(pelty) % weigp(1:pgaus) * gpdet(1:pgaus)
          !
          ! \int_S1D q vi ds (vi: shape function in background mesh, S1D: BAR3D element)
          !
          gauss_points: do igaus = 1,pgaus

             jelem = spare % eleme(ielem) % l(igaus)
             if( jelem /= 0 ) then
                do jnode = 1,meshe(ndivi) % lnnod(jelem)
                   jpoin = meshe(ndivi) % lnods(jnode,jelem)
                   do jdime = 1,nd
                      bu(jdime,jpoin) = bu(jdime,jpoin) &
                           + xfact * diag(jdime,jpoin) &
                           * spare % shapf(ielem) % a (jnode,igaus)   &
                           * gpvol(igaus) * bvess(jdime)
                   end do
                   do knode = 1,meshe(ndivi) % lnnod(jelem)
                      kpoin = meshe(ndivi) % lnods(knode,jelem)
                      izdom = r_dom(jpoin)
                      kcolu = c_dom(izdom)
                      do while( kcolu /= kpoin .and. izdom < r_dom(jpoin+1)-1 )
                         izdom = izdom + 1
                         kcolu = c_dom(izdom)
                      end do
                      if( kcolu == kpoin ) then
                         do jdime = 1,nd                            
                            Auu(jdime,jdime,izdom) = Auu(jdime,jdime,izdom) &
                                 + xfact * diag(jdime,jpoin) &
                                 * spare % shapf(ielem) % a (jnode,igaus)   &
                                 * spare % shapf(ielem) % a (knode,igaus)   &
                                 * gpvol(igaus)
                         end do
                      end if
                   end do
                end do
             end if

          end do gauss_points

       end if

    end do

    call memory_deallo(mem_modul(1:2,modul),'DIAG','nsi_spare_mesh_dirichlet_1D',diag)

  end subroutine nsi_spare_mesh_dirichlet_1D

end module mod_nsi_spare_mesh
!> @}
