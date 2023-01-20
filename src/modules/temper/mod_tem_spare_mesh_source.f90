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

module mod_tem_spare_mesh_source

  use def_kintyp,                only : ip,rp
  use def_kermod,                only : spare_meshes
  use def_domain,                only : nnode,mnode,mgaus,ndime
  use def_domain,                only : meshe,ngaus
  use def_kermod,                only : ndivi
  use def_kintyp_mesh_basic,     only : mesh_type_basic
  use mod_elmgeo,                only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo,                only : element_type
  use def_master,                only : INOTMASTER
  use def_kintyp_spare_mesh,     only : typ_spare_mesh
  use def_temper

  implicit none
  private
   
  real(rp),         allocatable :: shapf(:)
  real(rp),         allocatable :: elcod(:,:)
  real(rp),         allocatable :: gpcar(:,:,:)
  real(rp),         allocatable :: xjaci(:,:,:)
  real(rp),         allocatable :: gpvol(:)
  real(rp),         allocatable :: gpdet(:)
    
  public :: tem_spare_mesh_source

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Loop over elements
  !> @details Impose source along spare mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine tem_spare_mesh_source()

    integer(ip)                    :: imesh,pdime,mn,md,mg
    type(mesh_type_basic), pointer :: mesh

    if( kfl_sourc_tem == SOURCE_TERM_SPARE_MESH .and. INOTMASTER ) then

       imesh =  kfl_sonum_tem
       mesh  => spare_meshes(imesh) % mesh
       pdime =  maxval(element_type(spare_meshes(imesh) % mesh % ltype(1:spare_meshes(imesh) % mesh % nelem)) % dimensions)
       mg    =  maxval(spare_meshes(imesh) % mesh % quad(:) % ngaus)
       mn    =  maxval(element_type(spare_meshes(imesh) % mesh % ltype(:)) % number_nodes)
       md    =  pdime

       if( md > 0 ) then
          allocate(shapf(mn))
          allocate(elcod(md,mn))
          allocate(gpcar(md,mn,mg))
          allocate(xjaci(md,md,mg))
          allocate(gpvol(mg))
          allocate(gpdet(mg))
       end if
       
       select case ( pdime )
       case ( 0_ip ) ; call tem_spare_mesh_source_0D(spare_meshes(imesh),spare_meshes(imesh) % mesh)
       case ( 1_ip ) ; call tem_spare_mesh_source_1D(spare_meshes(imesh),spare_meshes(imesh) % mesh)
       end select
       
       if( allocated(shapf) ) deallocate(shapf)
       if( allocated(elcod) ) deallocate(elcod)
       if( allocated(gpcar) ) deallocate(gpcar)
       if( allocated(xjaci) ) deallocate(xjaci)
       if( allocated(gpvol) ) deallocate(gpvol)
       if( allocated(gpdet) ) deallocate(gpdet)
    end if

  end subroutine tem_spare_mesh_source

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Source is imposed along a 1D mesh
  !> @details Source is imposed along a 1D mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tem_spare_mesh_source_1D(spare,mesh)

    type(typ_spare_mesh),  intent(in) :: spare
    type(mesh_type_basic), intent(in) :: mesh
    integer(ip)                       :: pnode,jnode,jpoin,kgaus
    integer(ip)                       :: igaus,jelem,pelty
    integer(ip)                       :: ipoin,pgaus,pdime
    integer(ip)                       :: ipoin1,ipoin2,ielem
    !
    ! Initialize
    !
    pdime = 1
    do ipoin = 1,meshe(ndivi) % npoin
       heat_source(1,ipoin) = 0.0_rp
    end do
    !
    ! Assemble source term in HEAT_SOURCE
    !
    kgaus = 0

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
             kgaus = kgaus + 1

             jelem = spare % eleme(ielem) % l(igaus)
             if( jelem /= 0 ) then

                do jnode = 1,meshe(ndivi) % lnnod(jelem)
                   jpoin = meshe(ndivi) % lnods(jnode,jelem)
                   heat_source(1,jpoin) = heat_source(1,jpoin)  &
                        + mesh % tags(1) % values(ielem)        &
                        * spare % shapf(ielem) % a (jnode,igaus)   &
                        * gpvol(igaus)
                end do
             end if

          end do gauss_points

       else

          pelty = mesh % ltype(ielem)
          pgaus = mesh % quad(pelty) % ngaus
          kgaus = kgaus + pgaus

       end if
    end do    
    
  end subroutine tem_spare_mesh_source_1D

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Source is imposed along a 1D mesh
  !> @details Source is imposed along a 1D mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tem_spare_mesh_source_0D(spare,mesh)

    type(typ_spare_mesh),  intent(in) :: spare
    type(mesh_type_basic), intent(in) :: mesh
    integer(ip)                       :: jnode,jpoin,kgaus
    integer(ip)                       :: jelem,ipoin
    integer(ip)                       :: ielem
    !
    ! Initialize
    !
    do ipoin = 1,meshe(ndivi) % npoin
       heat_source(1,ipoin) = 0.0_rp
    end do
    !
    ! Assemble source term in HEAT_SOURCE
    !
    kgaus = 0

    do ielem = 1,mesh % nelem

       if( associated(spare % eleme(ielem) % l) ) then
          !
          ! \int_S0D q vi delta(x-x0) ds = q vi (vi: shape function in background mesh, S1D: POI3D element)
          !
          kgaus = kgaus + 1
          jelem = spare % eleme(ielem) % l(1)
          if( jelem /= 0 ) then             
             do jnode = 1,meshe(ndivi) % lnnod(jelem)
                jpoin = meshe(ndivi) % lnods(jnode,jelem)
                heat_source(1,jpoin) = heat_source(1,jpoin)  &
                     + mesh % tags(1) % values(ielem)        &
                     * spare % shapf(ielem) % a (jnode,1)
             end do
          end if

       else

          kgaus = kgaus + 1

       end if
       
    end do    

  end subroutine tem_spare_mesh_source_0D

end module mod_tem_spare_mesh_source
!> @}
