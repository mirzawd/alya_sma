!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_witness.f90
!> @author  houzeaux
!> @date    2020-01-28
!> @brief   Witness 
!> @details Witness points and geometries
!-----------------------------------------------------------------------

module mod_witness

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_elmtyp
  use mod_elmgeo,               only : element_type
  use mod_elmgeo,               only : elmgeo_natural_coordinates
  use mod_communications,       only : PAR_ALL_TO_ALL_ARRAY_OPERATION
  use mod_communications,       only : PAR_SUM
  use mod_communications,       only : PAR_MAX
  use mod_elsest,               only : elsest_host_element
  use mod_messages,             only : messages_live
  use mod_outfor,               only : outfor
  use mod_memory,               only : memory_alloca
  use mod_memory,               only : memory_deallo
  use mod_memory,               only : memory_resize
  use mod_memory,               only : memory_size
  use mod_memory,               only : memory_copy
  use mod_domain,               only : domain_memory_allocate
  use mod_domain,               only : domain_memory_deallocate
  use mod_maths,                only : maths_in_box 
  use mod_maths,                only : maths_vector_from_new_basis
  use mod_maths,                only : maths_vector_to_new_basis
  use mod_maths,                only : maths_local_orthonormal_basis
  use mod_maths_basic,          only : maths_inverse_permutation
  use mod_maths_geometry,       only : maths_point_plane_distance
  use mod_maths_geometry,       only : maths_point_line_distance
  use mod_mesh_type_basic,      only : mesh_type_basic_parallel
  use mod_mesh_type_basic,      only : mesh_type_basic_copy_boundary_mesh
  use mod_mesh_type_basic,      only : mesh_type_basic_output
  use mod_mesh_type,            only : mesh_number_boundary_nodes
  use mod_mesh_type,            only : mesh_coherent_boundary_mesh
  use mod_postpr_on_basic_mesh, only : postpr_mesh_file_list
  use mod_strings,              only : integer_to_string
  use mod_surface_mesh,         only : surface_mesh_coarsening
  use mod_kdtree,               only : kdtree_construct
  use mod_kdtree,               only : typ_kdtree
  use mod_kdtree,               only : kdtree_nearest_boundary
  use mod_kdtree,               only : kdtree_deallocate
  use mod_kdtree,               only : kdtree_initialize
  use mod_exchange,             only : exchange_init
  use mod_exchange,             only : exchange_add
  use mod_exchange,             only : exchange_end
  implicit none

  private

  character(11), parameter :: vacal='mod_witness'
  
  integer(ip), parameter :: &
       WITNESS_INSTANTANEAOUS =  0, &
       WITNESS_AVERAGE        =  1, &
       WITNESS_ACCUMULATE     =  2, &
       WITNESS_SPHERE         =  0, &
       WITNESS_BOX            =  1, &
       WITNESS_RING           =  2, &
       WITNESS_ELEMENT_SET    =  3, & 
       WITNESS_BOUNDARY       =  4, & 
       WITNESS_MATERIAL       =  5, &
       WITNESS_PLANE          =  6, &
       WITNESS_BOUNDARY_SET   =  7, &
       WITNESS_BAR3D          =  8, &
       WITNESS_3D_CIRCLE      =  9, &
       WITNESS_3D_RING        = 10, &
       WITNESS_IN_BOX         = 11
  
  public :: witness_point
  public :: witness_geometry_initialization
  public :: witness_parall
  public :: witness_mesh_initialization
  public :: witness_in_geometry
  public :: witness_ring_transform
  public :: witness_volume_geometry
  public :: witness_geometry_averaging_ini
  public :: witness_geometry_averaging_end
  public :: witness_point_averaging_ini
  public :: witness_point_averaging_end
  public :: witness_value_ini
  public :: witness_value_increment_denom
  public :: witness_value_end
  public :: witness_mesh_create

  public :: witness_mesh_plane ! should be removed
  
  public :: WITNESS_INSTANTANEAOUS
  public :: WITNESS_AVERAGE
  public :: WITNESS_ACCUMULATE
  public :: WITNESS_SPHERE
  public :: WITNESS_BOX
  public :: WITNESS_RING
  public :: WITNESS_ELEMENT_SET
  public :: WITNESS_BOUNDARY
  public :: WITNESS_MATERIAL
  public :: WITNESS_PLANE
  public :: WITNESS_BOUNDARY_SET
  public :: WITNESS_BAR3D
  public :: WITNESS_3D_CIRCLE
  public :: WITNESS_3D_RING
  public :: WITNESS_IN_BOX
    
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Averaging
  !> @details Averaging and desavearging when doing output
  !> 
  !-----------------------------------------------------------------------
  
  subroutine witness_geometry_averaging_ini()

    integer(ip) :: ivari,iwitg

    do ivari = 1,nvarg
       if( postp(1) % npp_witng(ivari) > 0 ) then
          do iwitg = 1,nwitg
             if( postp(1) % witng_kfl_time(ivari) == WITNESS_AVERAGE .and. postp(1) % witng_denom(ivari,iwitg) > 0.0_rp ) then
                postp(1) % witng(ivari,iwitg) = postp(1) % witng(ivari,iwitg) / postp(1) % witng_denom(ivari,iwitg)
              end if
          end do
       end if
    end do
    
  end subroutine witness_geometry_averaging_ini
  
  subroutine witness_point_averaging_ini()

    integer(ip) :: ivari,iwitg

    do ivari = 1,nvarw
       if( postp(1) % npp_witne(ivari) > 0 ) then
          do iwitg = 1,nwitn
             if( postp(1) % witne_kfl_time(ivari) == WITNESS_AVERAGE .and. postp(1) % witne_average(ivari) > 0.0_rp ) then
                postp(1) % witne(ivari,iwitg) = postp(1) % witne(ivari,iwitg) / postp(1) % witne_average(ivari)
              end if
          end do
       end if
    end do
    
  end subroutine witness_point_averaging_ini
  
  subroutine witness_geometry_averaging_end()
    
    integer(ip) :: ivari,iwitg

    do ivari = 1,nvarg
       if( postp(1) % npp_witng(ivari) > 0 ) then
          do iwitg = 1,nwitg
             if( postp(1) % witng_kfl_time(ivari) == WITNESS_AVERAGE .and. postp(1) % witng_denom(ivari,iwitg) > 0.0_rp ) then
                postp(1) % witng(ivari,iwitg) = postp(1) % witng(ivari,iwitg) * postp(1) % witng_denom(ivari,iwitg)
             end if
          end do
       end if
    end do
    
  end subroutine witness_geometry_averaging_end
  
  subroutine witness_point_averaging_end()
    
    integer(ip) :: ivari,iwitg

    do ivari = 1,nvarw
       if( postp(1) % npp_witne(ivari) > 0 ) then
          do iwitg = 1,nwitn
             if( postp(1) % witne_kfl_time(ivari) == WITNESS_AVERAGE .and. postp(1) % witne_average(ivari) > 0.0_rp ) then
                postp(1) % witne(ivari,iwitg) = postp(1) % witne(ivari,iwitg) * postp(1) % witne_average(ivari)
             end if
          end do
       end if
    end do
    
  end subroutine witness_point_averaging_end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Initialize
  !> @details Initialize a value according to the strategy
  !>
  !-----------------------------------------------------------------------

  subroutine witness_value_ini(postp_in)

    type(typos), pointer, intent(inout) :: postp_in(:)
    integer(ip)                         :: iwitg,ivari    
    real(rp)                            :: overt

    do ivari = 1,nvarg
       if( postp_in(1) % npp_witng(ivari) > 0 ) then
          do iwitg = 1,nwitg
             if(  postp_in(1) % witng_kfl_time(ivari) == WITNESS_INSTANTANEAOUS ) then
                !
                ! Instantaneous values
                !
                postp_in(1) % witng_deldenom(ivari,iwitg)   = 1.0_rp
                postp_in(1) % witng_delperio(ivari,iwitg)   = 1.0_rp
                postp_in(1) % witng(ivari,iwitg)            = 0.0_rp
              
             else if( postp_in(1) % witng_kfl_time(ivari) == WITNESS_AVERAGE .or. & 
                &     postp_in(1) % witng_kfl_time(ivari) == WITNESS_ACCUMULATE ) then
                !
                ! Average values: step in denominator, can be overwritten by modules
                !
                overt                                       = postp_in(1) % witng_perio_count(ivari,iwitg) - postp_in(1) % witng_perio_max(ivari)
                postp_in(1) % witng_deldenom(ivari,iwitg)   = dtime - max(overt,0.0_rp)
                postp_in(1) % witng_delperio(ivari,iwitg)   = dtime - max(overt,0.0_rp)
             end if
          end do
       end if
    end do
     do ivari = 1,nvarw
       if( postp_in(1) % npp_witne(ivari) > 0 ) then
          do iwitg = 1,nwitn
             if(  postp_in(1) % witne_kfl_time(ivari) == WITNESS_INSTANTANEAOUS ) then
                !
                ! Instantaneous values
                !
                postp_in(1) % witne_dt(ivari)      = 1.0_rp
                postp_in(1) % witne(ivari,iwitg)   = 0.0_rp
                
             else if( postp_in(1) % witne_kfl_time(ivari) == WITNESS_AVERAGE .or. & 
                &     postp_in(1) % witne_kfl_time(ivari) == WITNESS_ACCUMULATE ) then
                !
                ! Average values
                !
                overt                              = postp_in(1) % witne_average(ivari) - postp_in(1) % witne_period(ivari)
                postp_in(1) % witne_dt(ivari)      = dtime - max(overt,0.0_rp)
             end if
          end do
       end if
    end do
   
  end subroutine witness_value_ini

  !-----------------------------------------------------------------------
  !> 
  !> @author  aboth
  !> @date    2020-03-06
  !> @brief   Increment
  !> @details Increment averaging values
  !>
  !-----------------------------------------------------------------------

  subroutine witness_value_increment_denom(postp_in)
    type(typos), pointer, intent(inout) :: postp_in(:)
    integer(ip)                         :: iwitg,ivari    

    do ivari = 1,nvarg
       if( postp_in(1) % npp_witng(ivari) > 0 ) then
          do iwitg = 1,nwitg
             if( postp_in(1) % witng_kfl_time(ivari) == WITNESS_AVERAGE .or. &
                &postp_in(1) % witng_kfl_time(ivari) == WITNESS_ACCUMULATE) then
                !
                ! Average values
                ! 
                postp_in(1) % witng_perio_count(ivari,iwitg)   = postp_in(1) % witng_perio_count(ivari,iwitg)  + postp_in(1) % witng_delperio(ivari,iwitg)
                postp_in(1) % witng_denom(ivari,iwitg)         = postp_in(1) % witng_denom(ivari,iwitg)        + postp_in(1) % witng_deldenom(ivari,iwitg)
             end if
          end do
       end if
    end do
     do ivari = 1,nvarw
       if( postp_in(1) % npp_witne(ivari) > 0 ) then
          do iwitg = 1,nwitn
             if( postp_in(1) % witne_kfl_time(ivari) == WITNESS_AVERAGE .or. &
                &postp_in(1) % witne_kfl_time(ivari) == WITNESS_ACCUMULATE) then
                !
                ! Average values
                !
                postp_in(1) % witne_average(ivari) = postp_in(1) % witne_average(ivari) + postp_in(1) % witne_dt(ivari)
             end if
          end do
       end if
    end do

  end subroutine witness_value_increment_denom
  
  subroutine witness_value_end(postp_in)

    type(typos), pointer, intent(inout) :: postp_in(:)
    integer(ip)                         :: iwitg,ivari    

    do ivari = 1,nvarg
       if( postp_in(1) % npp_witng(ivari) > 0 ) then
          do iwitg = 1,nwitg
             if(  postp_in(1) % witng_kfl_time(ivari) == WITNESS_INSTANTANEAOUS ) then
                !
                ! Instantaneous values
                !
                postp_in(1) % witng_perio_count(ivari,iwitg)   = 0.0_rp
                postp_in(1) % witng_denom(ivari,iwitg)         = 0.0_rp
                postp_in(1) % witng(ivari,iwitg)               = 0.0_rp
             else if( postp_in(1) % witng_kfl_time(ivari) == WITNESS_AVERAGE .or. &
                &     postp_in(1) % witng_kfl_time(ivari) == WITNESS_ACCUMULATE) then
                !
                ! Average or accumulated values
                !
                if( postp_in(1) % witng_perio_count(ivari,iwitg) > postp_in(1) % witng_perio_max(ivari) ) then
                   postp_in(1) % witng_perio_count(ivari,iwitg)   = 0.0_rp
                   postp_in(1) % witng_denom(ivari,iwitg)         = 0.0_rp
                   postp_in(1) % witng(ivari,iwitg)               = 0.0_rp
                end if
             end if
          end do
       end if
    end do
    do ivari = 1,nvarw
       if( postp_in(1) % npp_witne(ivari) > 0 ) then
          do iwitg = 1,nwitn
             if(  postp_in(1) % witne_kfl_time(ivari) == WITNESS_INSTANTANEAOUS ) then
                !
                ! Instantaneous values
                !
                postp_in(1) % witne_average(ivari) = 0.0_rp
                postp_in(1) % witne(ivari,iwitg)   = 0.0_rp
             else if( postp_in(1) % witne_kfl_time(ivari) == WITNESS_AVERAGE .or. &
                &     postp_in(1) % witne_kfl_time(ivari) == WITNESS_ACCUMULATE) then
                !
                ! Average values
                !
                if( postp_in(1) % witne_average(ivari) > postp_in(1) % witne_period(ivari) ) then
                   postp_in(1) % witne_average(ivari) = 0.0_rp
                   postp_in(1) % witne(ivari,iwitg)   = 0.0_rp
                end if
             end if
          end do
       end if
    end do
    
  end subroutine witness_value_end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Check geometry
  !> @details Chekc if a point is inside a witness geometry
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function witness_in_geometry(xcoor,kfl_geometry,param)

    real(rp),    intent(in) :: xcoor(ndime)
    integer(ip), intent(in) :: kfl_geometry
    real(rp),    intent(in) :: param(*)
    real(rp)                :: rad,center(3),dist
    real(rp)                :: dr,dn,xpol(3)

    witness_in_geometry = .false.
    
    if( kfl_geometry == WITNESS_SPHERE ) then
       
       rad         = param(1)
       center(1:3) = param(2:4)
       dist        = dot_product(xcoor(1:ndime)-center(1:ndime),xcoor(1:ndime)-center(1:ndime))      
       if( dist <= rad**2 ) witness_in_geometry = .true.
       
    else if( kfl_geometry == WITNESS_BOX .or. kfl_geometry == WITNESS_IN_BOX ) then
       !
       ! 6 parameters: lower corner, higher corner
       !
       witness_in_geometry = maths_in_box(ndime,xcoor,param(1:ndime),param(ndime+1:2_ip*ndime))
      
    else if( kfl_geometry == WITNESS_RING ) then
       !
       ! 9 parameters: 
       !
       center(1:3)       = param(1:3)
       rad               = param(7)
       dr                = param(8)
       dn                = param(9)
       call witness_ring_transform(xcoor,param,xpol=xpol)
       !
       ! Test axial and radial distance:
       !
       if ((abs(xpol(1)) <= dn/2.0_rp) .and.( abs(xpol(2)-rad) <= dr/2.0_rp ) ) then
          witness_in_geometry = .true.
       endif
       
    end if
    
  end function witness_in_geometry
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  aboth
  !> @date    2020-04-03
  !> @brief   Geometric witness ring transformation
  !> @details Transform coordinates to align with ring directions
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_ring_transform(xcoor,param,basis,xpol,vec)
     real(rp),           intent(in)    :: xcoor(ndime)
     real(rp),           intent(in)    :: param(*)
     real(rp), optional, intent(out)   :: basis(ndime, ndime)
     real(rp), optional, intent(out)   :: xpol(ndime)
     real(rp), optional, intent(inout) :: vec(ndime)
     real(rp)                          :: rad,center(3)
     real(rp)                          :: basis_loc(ndime, ndime),xcart(3),xpol_loc(3)
!     real(rp)                          :: rcoord(3)
     real(rp)                          :: rlen
     !
     ! Parameters: 
     !
     center(1:3)           = param(1:3)
     basis_loc(1:ndime,1)  = param(4:3+ndime)
     rad                   = param(7)
     !
     ! First basis computation: only knowing the normal
     !
     call maths_local_orthonormal_basis(ndime, basis_loc)
     !
     ! Transform to local coordinates to know normal:
     !
     xcart = 0.0_rp
     xcart(1:ndime) = xcoor(1:ndime) - center(1:ndime)
     xpol_loc(1:ndime)  = xcart(1:ndime)
     call maths_vector_to_new_basis(ndime,basis_loc,xpol_loc)
     !
     ! Now it is possible to update the basis 
     !
     rlen = sqrt(dot_product(xpol_loc(2:ndime),xpol_loc(2:ndime)))
     if (rlen > 1.0e-14_rp) then
        !
        ! Second component of basis is radial
        !
        !basis_loc(1:ndime,2)  = rcoord(1:ndime) / rlen
        basis_loc(1:ndime,2)  = (xcart(1:ndime) -  basis_loc(1:ndime,1) * xpol_loc(1))
        basis_loc(1:ndime,2)  = basis_loc(1:ndime,2) /  sqrt(dot_product(basis_loc(1:ndime,2),basis_loc(1:ndime,2)))
        if (ndime >= 3) then
           !
           ! Third component of basis is tangential
           !
           basis_loc(1,3) = basis_loc(2,1)*basis_loc(3,2) - basis_loc(3,1)*basis_loc(2,2)
           basis_loc(2,3) = basis_loc(3,1)*basis_loc(1,2) - basis_loc(1,1)*basis_loc(3,2)
           basis_loc(3,3) = basis_loc(1,1)*basis_loc(2,2) - basis_loc(2,1)*basis_loc(1,2)
        endif
     endif
       
     !
     ! Optional output
     ! 
     if (present(basis)) basis = basis_loc
     if (present(xpol)) then
        xpol = xcart
        call maths_vector_to_new_basis(ndime,basis_loc,xpol)
     endif
     if (present(vec)) then
        call maths_vector_to_new_basis(ndime,basis_loc,vec)
     endif

  end subroutine witness_ring_transform
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  aboth
  !> @date    2020-04-03
  !> @brief   Witness geometry volume
  !> @details Calculate volume
  !> 
  !-----------------------------------------------------------------------

  real(rp) function witness_volume_geometry(kfl_geometry,param)

    integer(ip), intent(in) :: kfl_geometry
    real(rp),    intent(in) :: param(*)
    real(rp)                :: rad,dbox(3)
    real(rp)                :: dr,dn

    witness_volume_geometry = 0.0_rp
    
    if( kfl_geometry == WITNESS_SPHERE ) then
       rad                      = param(1)
       witness_volume_geometry  = 4.0_rp* oneo3 * rad**3 * pi 
    else if( kfl_geometry == WITNESS_BOX ) then
       !
       ! 6 parameters: lower corner, higher corner
       !
       dbox = 0.0_rp
       dbox(1:ndime) = abs( param(ndime+1:2_ip*ndime) - param(1:ndime) )
       witness_volume_geometry  = dbox(1)
       if (ndime >= 2) witness_volume_geometry  =  witness_volume_geometry * dbox(2)
       if (ndime >= 3) witness_volume_geometry  =  witness_volume_geometry * dbox(3)
    else if( kfl_geometry == WITNESS_RING ) then
       !
       ! 9 parameters: 
       !
       rad               = param(7)
       dr                = param(8)
       dn                = param(9)
       witness_volume_geometry  = dn * dr * 2.0_rp * rad * pi 
    end if
    
  end function witness_volume_geometry

   !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness point parallelization
  !> @details Broadcast data of witness points
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_parall()

    call witness_point_parall()
    call witness_geometry_parall()
    call witness_mesh_parall()    

  end subroutine witness_parall
  
 !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness point parallelization
  !> @details Broadcast data of witness points
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_point_parall()
    
    if( nwitn /= 0 ) then
       if( ISLAVE  ) call domain_memory_allocate  ('COWIT_ORIGI') 
       call exchange_init()
       call exchange_add(cowit_origi)
       call exchange_end()
       if( IMASTER ) call domain_memory_deallocate('COWIT_ORIGI') 
    end if
 
  end subroutine witness_point_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness geometry parallelization
  !> @details Broadcast data of witness geometries
  !> 
  !-----------------------------------------------------------------------
  
  subroutine witness_geometry_parall()

    integer(ip) :: ii,kk

    if( nwitg /= 0 ) then
       if( ISLAVE ) call domain_memory_allocate('GEWIT')
       call exchange_init()
       do ii = 1,nwitg
          call exchange_add(gewit(ii) % kfl_geometry)
          do kk = 1,size(gewit(ii) % param,KIND=ip)
             call exchange_add(gewit(ii) % param(kk))
          end do
       end do
       call exchange_end()
    end if
     
  end subroutine witness_geometry_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness geometry parallelization
  !> @details Broadcast data of witness geometries
  !> 
  !-----------------------------------------------------------------------
  
  subroutine witness_mesh_parall()

    integer(ip) :: ii,kk

    if( nwith /= 0 ) then
       if( ISLAVE ) call domain_memory_allocate('WITNESS_MESH')
       call exchange_init()
       do ii = 1,nwith           
          call exchange_add(witness_mesh(ii) % name)
          call exchange_add(witness_mesh(ii) % geom % kfl_geometry)
          do kk = 1,size(witness_mesh(ii) % geom % param,KIND=ip)
             call exchange_add(witness_mesh(ii) % geom % param(kk))
          end do
       end do
       call exchange_end()
    end if
    
  end subroutine witness_mesh_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness geometry parallelization
  !> @details Broadcast data of witness geometries
  !> 
  !-----------------------------------------------------------------------
  
  subroutine witness_mesh_initialization()
 
    integer(ip) :: iwitg

    do iwitg = 1,nwith
       call witness_mesh(iwitg) % init()
    end do
    
  end subroutine witness_mesh_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Initialize witness geometry
  !> @details Initialize witness geometry
  !> 
  !-----------------------------------------------------------------------
  
  subroutine witness_geometry_initialization(witness_in)

    type(witness_geo), optional, intent(out) :: witness_in
    integer(ip)                             :: ii

    if( present(witness_in) ) then
       witness_in % kfl_geometry = 0                      ! No type 
       witness_in % param        = 0.0_rp                 ! Geometry parameters
    else
       if( associated(gewit) .and. nwitg > 0 ) then      
          do ii = 1,size(gewit)
             gewit(ii) % kfl_geometry = WITNESS_SPHERE         ! Geometry type 
             gewit(ii) % param        = 0.0_rp                 ! Geometry parameters
          end do
       end if
    end if
    
  end subroutine witness_geometry_initialization
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-28
  !> @brief   Witness points
  !> @details Get witness point information
  !>          LEWIT ... Host element
  !>          SHWIT ... Shape function in host element
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_point()

    integer(ip)          :: iwitn,kwitn,idime,inode,ierro
    integer(ip)          :: ielem,ipoin,pnode
    real(rp)             :: coloc(3),dNds(ndime,mnode),elcod(ndime,mnode)
    real(rp)             :: xjaci(9),xjacm(9),gpdet,dista
    integer(ip), pointer :: lnwit_tmp(:)

    if( nwitn_all > 0 ) then

       nwitn = nwitn_all
       nullify( lnwit_tmp )
       call messages_live('COMPUTE WITNESS POINT INFORMATION')
       !
       ! Find host element: LEWIT
       !
       if( INOTMASTER ) then

          call domain_memory_allocate('WITNESS POINTS')

          if( INOTEMPTY ) then
             do iwitn = 1,nwitn
                call elsest_host_element(&
                     ielse,relse,1_ip,meshe(ndivi),cowit_origi(:,iwitn),lewit(iwitn),&
                     shwit(:,iwitn),dNds,coloc,dista)
                if( lewit(iwitn) < 1 ) then
                   lewit(iwitn) = 0 
                else
                   ielem = lewit(iwitn)
                   pnode = lnnod(ielem)
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                   end do
                   call jacobi(&
                        ndime,pnode,elcod,dNds,&
                        xjacm,xjaci,dewit(1,1,iwitn),gpdet) 
                end if

             end do
          end if

       else

          call memory_alloca(memor_dom,'LEWIT','witnes',lewit,nwitn)

       end if
       !
       ! LNWIT: Allocate memory
       !
       call domain_memory_allocate('LNWIT')
       !
       ! NWITN: Count my own witness points
       !
       call PAR_ALL_TO_ALL_ARRAY_OPERATION(nwitn,lewit,'CHOOSE ONLY ONE')
       if( INOTMASTER ) then  
          kwitn = 0
          do iwitn = 1,nwitn
             if( lewit(iwitn) /= 0 ) kwitn = kwitn + 1
          end do
       end if
       !
       ! Renumber my witness nodes
       !
       if( ISLAVE ) then
          call domain_memory_allocate('COWIT')
          kwitn = 0
          do iwitn = 1,nwitn
             if( lewit(iwitn) /= 0 ) then
                kwitn        = kwitn + 1
                lewit(kwitn) = lewit(iwitn)
                lnwit(kwitn) = iwitn
                do inode = 1,mnode
                   shwit(inode,kwitn) = shwit(inode,iwitn)
                   do idime = 1,ndime
                      dewit(idime,inode,kwitn) = dewit(idime,inode,iwitn)
                   end do
                end do
                do idime = 1,ndime
                   cowit(idime,kwitn) = cowit_origi(idime,iwitn)
                end do
             end if
          end do
          nwitn = kwitn
       end if
       !
       ! Prepare all gather to replace call Parall(45_ip)
       !
       !if( IPARALL ) then
       !   allocate(nwitn_par(0:npart))
       !   if( IMASTER ) then
       !      nwitn_par = 0
       !      do iwitn = 1,nwitn
       !         ipart = lewit(iwitn)
       !         if( ipart /= 0 ) then
       !            nwitn_par(ipart) = nwitn_par(ipart) + 1
       !         end if
       !      end do
       !   else
       !      nwitn_par = 0
       !      allocate( lnwit_tmp(nwitn) )
       !      lnwit_tmp(1:nwitn) = lnwit(1:nwitn)
       !   end if
       !   call PAR_GATHERV(lnwit_tmp,lnwit,nwitn_par,'IN MY CODE')
       !   if( ISLAVE ) then
       !      deallocate( lnwit_tmp )
       !      deallocate( nwitn_par )
       !   end if
       !end if
       !
       ! Master gets numbering LNWIT
       !
       call par_comset(4_ip)
       !
       ! Detect a possible problem
       !
       ierro = 0
       if( IMASTER ) then
          if( mwitn /= nwitn ) ierro = mwitn - nwitn
          if (ierro < 0) ierro= -ierro
       else if( ISEQUEN ) then
          if( kwitn /= nwitn ) ierro = kwitn - nwitn
          if (ierro < 0) ierro= -ierro
       end if
       call outfor(50_ip,lun_outpu,'')
       ioutp(2) = 0
       if( ierro > 0 ) then
          call messages_live(&
               'WARNING: '//trim(intost(ierro))//' WITNESS POINTS ARE LOST: ZERO RESULTS WILL APPEAR','WARNING')
          if( IMASTER ) then
             call memgen(1_ip,nwitn,0_ip)
             do iwitn = 1,nwitn
                if( lnwit(iwitn) /= 0 ) gisca(lnwit(iwitn)) = 1
             end do
             do iwitn = 1,nwitn
                if( gisca(iwitn) == 0 )then
                   ioutp(1) = iwitn
                   call outfor(51_ip,lun_outpu,'')
                   ioutp(2) = ioutp(2) + 1
                end if
             end do
             call memgen(3_ip,nwitn,0_ip)
          else if( ISEQUEN ) then        
             do iwitn = 1,nwitn
                if( lewit(iwitn) == 0 )then
                   ioutp(1) = iwitn
                   call outfor(51_ip,lun_outpu,'')
                   ioutp(2) = ioutp(2) + 1
                end if
             end do
          end if
       else
          call outfor(52_ip,lun_outpu,'')
       end if

    end if

  end subroutine witness_point

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-03
  !> @brief   Volume integral
  !> @details Compute integral
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_geometry_integral()

    integer(ip) :: pnode,iwitg,pelty,pgaus
    integer(ip) :: ix,iy,iz,nn(3),nz
    real(rp)    :: bbmin(3),bbmax(3)
    real(rp)    :: xx(3),delta(3),elcod(ndime,8)
!    real(rp)    :: gpvol(8)
!    real(rp)    :: gpder(ndime,mnode,8)
!    real(rp)    :: gpcar(ndime,mnode,8)
    
    nn=10
    if( ndime == 2 ) then
       pelty = QUA04
       pnode = 4
       pgaus = 4
       nz    = 1
    else
       pelty = HEX08
       pnode = 8
       pgaus = 8
       nz    = nn(3)
    end if
    
    do iwitg = 1,nwitg
       bbmin(1:3) = gewit(iwitg) % param(1:3)
       bbmax(1:3) = gewit(iwitg) % param(1:3)
       delta(1:3) = (bbmax(1:3) - bbmin(1:3))/real(nn(1:3),rp)
       do iz = 1,nz
          xx(3) = bbmin(3) + real(iz-1,rp)*delta(3)
          do iy = 1,nn(2)
             xx(2) = bbmin(2) + real(iy-1,rp)*delta(2)
             do ix = 1,nn(1)
                xx(1) = bbmin(1) + real(ix-1,rp)*delta(1)
                elcod(1:ndime,1) = xx(1:ndime)
                elcod(1,2)       = xx(1) + delta(1)
                elcod(2,2)       = xx(2)
                elcod(1,3)       = xx(1) + delta(1)
                elcod(2,3)       = xx(2) + delta(2)
                elcod(1,4)       = xx(1) 
                elcod(2,4)       = xx(2) + delta(2)
                if( ndime == 3 ) then
                   elcod(3,1:4)   = xx(3)
                   elcod(3,5:8)   = xx(3) + delta(3)
                   elcod(1:2,5:8) = elcod(1:2,1:4)
                end if
                !call element_shape_function_derivatives_jacobian(&
                !     pnode,pgaus,0_ip,weigp,elmar(pelty) % shapf,elmar(pelty) % deriv,&
                !     elmar(pelty) % heslo,elcod,gpvol,gpsha,gpder,gpcar)
             end do
          end do
       end do
    end do
    
  end subroutine witness_geometry_integral

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-03
  !> @brief   Create witness mesh
  !> @details Create witness mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_create()

    use def_master, only : kfl_paral
    use mod_parall, only : PAR_COMM_MY_CODE

    logical(lg), pointer  :: lmask(:)
    integer(ip)           :: ii,ielem,inode,pnode,ipoin,ichek,iboun
    integer(ip)           :: nelem_rad,nelem_cir
    real(rp)              :: xx(3)
    real(rp),    pointer  :: r1(:)
    
    if( nwith > 0 ) then
       
       call messages_live('CREATE AND OUTPUT WITNESS MESHES','START SECTION')
       nullify(lmask)

       do ii = 1,nwith

          call messages_live('WITNESS MESH '//trim(witness_mesh(ii) % name))
          call witness_mesh(ii) % mesh % init()

          if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_BOUNDARY ) then
             !
             ! Boundary mesh and coarse boundary mesh
             !
             ipoin = int(witness_mesh(ii) % geom % param(1),ip)
             if( ipoin == 0_ip ) then
                call mesh_type_basic_copy_boundary_mesh(witness_mesh(ii) % mesh,meshe(ndivi))
             else
                call witness_coarse_boundary(witness_mesh(ii) % mesh,meshe(ndivi),witness_mesh(ii) % geom % param,witness_mesh(ii) % inte)
             end if

          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_BOUNDARY_SET ) then
             !
             ! Boundary set
             !
             call memory_alloca(memor_dom,'LMASK','mod_witness',lmask,nboun)
             do iboun = 1,nboun
                if( lbset(iboun) == int(witness_mesh(ii) % geom % param(1),ip) ) lmask(iboun) = .true.
             end do
             call mesh_type_basic_copy_boundary_mesh(witness_mesh(ii) % mesh,meshe(ndivi),LMASK=lmask)
             call memory_deallo(memor_dom,'LMASK','mod_witness',lmask)

          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_PLANE ) then
             !
             ! Boundary mesh
             !
             call witness_mesh_plane(witness_mesh(ii) % mesh,meshe(ndivi),witness_mesh(ii) % geom % param,witness_mesh(ii) % inte) 
             
          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_BOX ) then
             !
             ! Cartesian mesh in a box
             !
             call witness_mesh_box(witness_mesh(ii) % mesh,meshe(ndivi),witness_mesh(ii) % geom % param,witness_mesh(ii) % inte) 

          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_BAR3D ) then
             !
             ! 3D mesh from BAR3D mesh
             !
             nelem_rad =  int(witness_mesh(ii) % geom % param(1),ip)
             nelem_cir =  int(witness_mesh(ii) % geom % param(2),ip)
             r1        => areas(:,1)
             call witness_mesh(ii) % mesh % mesh_from_BAR3D(meshe(ndivi),nelem_rad,nelem_cir,r1,witness_mesh(ii) % perm)
             
          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_3D_CIRCLE ) then
             !
             ! 3D circle mesh
             !
             call witness_mesh_3Dcircle(witness_mesh(ii) % mesh,meshe(ndivi),witness_mesh(ii) % geom % param,witness_mesh(ii) % inte) 
              
          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_3D_RING ) then
             !
             ! 3D ring mesh
             !
             call witness_mesh_3Dring(witness_mesh(ii) % mesh,meshe(ndivi),witness_mesh(ii) % geom % param,witness_mesh(ii) % inte) 
              
          else if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_IN_BOX ) then
             !
             ! Meshes with mask
             !
             call memory_alloca(memor_dom,'LMASK','mod_witness',lmask,nelem)             
             do ielem = 1,nelem
                xx    = 0.0_rp
                pnode = element_type(abs(ltype(ielem))) % number_nodes
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   xx(1:ndime) = xx(1:ndime) + coord(1:ndime,ipoin)
                end do
                xx = xx / real(pnode,rp)
                lmask(ielem) = witness_in_geometry(xx,witness_mesh(ii) % geom % kfl_geometry,witness_mesh(ii) % geom % param)
             end do
             call witness_mesh(ii) % mesh % extract(meshe(ndivi),lmask)
             call memory_deallo(memor_dom,'LMASK','mod_witness',lmask)

          else
             !
             ! Meshes with mask
             !
             call memory_alloca(memor_dom,'LMASK','mod_witness',lmask,nelem)
             if( witness_mesh(ii) % geom % kfl_geometry == WITNESS_ELEMENT_SET ) then
                do ielem = 1,nelem
                   if( leset(ielem) == int(witness_mesh(ii) % geom % param(1),ip) ) lmask(ielem) = .true.
                end do
             else
                do ielem = 1,nelem
                   xx    = 0.0_rp
                   pnode = element_type(abs(ltype(ielem))) % number_nodes
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      xx(1:ndime) = xx(1:ndime) + coord(1:ndime,ipoin)
                   end do
                   xx = xx / real(pnode,rp)
                   lmask(ielem) = witness_in_geometry(xx,witness_mesh(ii) % geom % kfl_geometry,witness_mesh(ii) % geom % param)
                end do
             end if
             call witness_mesh(ii) % mesh % extract(meshe(ndivi),lmask)
             call memory_deallo(memor_dom,'LMASK','mod_witness',lmask)

          end if
          !
          ! Name of the mesh
          !
          witness_mesh(ii) % mesh % name = witness_mesh(ii) % name
          !
          ! Check validity of the mesh
          !
          ichek = witness_mesh(ii) % mesh % nelem
          call PAR_SUM(ichek)
          if( ichek == 0 ) then
             !the warning was pointless
             !the problem is that if you start saving the witness mesh late in the simulation, 
             !you will spend several hours of compute time until alya dies with an error  
             call runend('WITNESS MESH '//trim(witness_mesh(ii) % mesh % name)//' IS EMPTY!')
          else
             !
             ! Parallelizaiton and output
             !
             call mesh_type_basic_parallel(witness_mesh(ii) % mesh,PAR_COMM_MY_CODE,kfl_paral)
             !call mesh_type_basic_parallel_numbering(witness_mesh(ii) % mesh)
             call mesh_type_basic_output(witness_mesh(ii) % mesh,MESH_OUTPUT=kfl_wimes/=0)
             call postpr_mesh_file_list(witness_mesh(ii) % mesh) 
          end if
       end do
       call messages_live('CREATE AND OUTPUT WITNESS MESHES','END SECTION')
       
    end if

  end subroutine witness_mesh_create

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a plane cartesian mesh
  !> @details Create a plane cartesian mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_plane(mesh_out,mesh_in,param,inte)

    type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    real(rp),                        intent(in)    :: param(:)
    type(typ_interp),                intent(inout) :: inte
    real(rp)                                       :: a,b,c,d
    integer(ip)                                    :: pdime
    integer(ip)                                    :: ipoin
    real(rp),              pointer                 :: dista(:)
  
    nullify(dista)
    
    pdime = mesh_in % ndime
    a     = param(1)
    b     = param(2)
    c     = param(3)
    d     = param(pdime+1)
    !
    ! Distance function
    !   
    call memory_alloca(memor_dom,'DISTA',vacal,dista,mesh_in % npoin)
    if( pdime == 2 ) then
       do ipoin = 1,mesh_in % npoin
          dista(ipoin) = maths_point_line_distance(mesh_in % coord(1:pdime,ipoin),a,b,c) 
       end do
    else
       do ipoin = 1,mesh_in % npoin
          dista(ipoin) = maths_point_plane_distance(mesh_in % coord(1:pdime,ipoin),a,b,c,d) 
       end do
    end if
    !
    ! Create mesh
    !
    call mesh_out % cut_level(mesh_in,dista,memor_dom)
    call memory_deallo(memor_dom,'DISTA',vacal,dista)   
    !
    ! Compute interpolation array
    !
    call witness_mesh_interpolation(mesh_out,mesh_in,inte)

  end subroutine witness_mesh_plane

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a plane cartesian mesh
  !> @details Create a plane cartesian mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_box(mesh_out,mesh_in,param,inte)

   type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    real(rp),                        intent(in)    :: param(:)
    type(typ_interp),                intent(inout) :: inte
    integer(ip)                                    :: knode,pnode,ielem
    integer(ip)                                    :: inode,ipoin,nn,kpoin
    logical(lg),           pointer                 :: mask(:)
    integer(ip),           pointer                 :: lelem_cpy(:)
    integer(ip),           pointer                 :: lemin(:)
    real(rp),              pointer                 :: shapf_cpy(:,:)
    type(mesh_type_basic)                          :: mesh_tmp
    integer(ip)                                    :: boxes(3)
    real(rp)                                       :: wimin(3)
    real(rp)                                       :: wimax(3)

    nullify(mask,lelem_cpy,shapf_cpy,lemin)
    
    boxes = (/int(param(7:9),ip)/)
    if( ndime == 2 ) then
       wimin(1:2) =(/param(1:2)/)
       wimax(1:2) =(/param(3:4)/)
    else
       wimin(1:3) =(/param(1:3)/)
       wimax(1:3) =(/param(4:6)/)
    end if    
    if( INOTMASTER ) then
       call mesh_tmp % init('WITNESS MESH')
       call mesh_tmp % mesh_3Dbox(ndime,boxes,wimin,wimax)
    end if
    !
    ! Compute interpolation array
    !
    call witness_mesh_interpolation(mesh_tmp,mesh_in,inte,OUTSIDE_VALUE=-1_ip)
    !
    ! Check repeated points
    !
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'LEMIN',vacal,lemin,mesh_tmp % npoin)
       do ipoin = 1,mesh_tmp % npoin
          if( inte % lelem(ipoin) > 0 ) lemin(ipoin) = kfl_paral
       end do
       call PAR_MAX(lemin,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
       do ipoin = 1,mesh_tmp % npoin
          if( lemin(ipoin) /= kfl_paral ) inte % lelem(ipoin) = -1
       end do
       call memory_deallo(memor_dom,'LEMIN',vacal,lemin)
    end if
    !
    ! Limit mesh to element owners
    !
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'MASK',vacal,mask,mesh_tmp % nelem)
       do ielem = 1,mesh_tmp % nelem
          mask(ielem) = .true.
          pnode       = element_type(mesh_tmp % ltype(ielem)) % number_nodes
          knode       = 0
          do inode = 1,pnode
             ipoin = mesh_tmp % lnods(inode,ielem)
             if( inte % lelem(ipoin) == -1 ) knode = knode + 1
          end do
          if( knode > 0 ) mask(ielem) = .false.
       end do
       
       call mesh_out % extract(mesh_tmp,mask)
       call mesh_tmp % deallo ()
       call memory_deallo     (memor_dom,'MASK',vacal,mask)

       nullify(mesh_out % parent) ! To tell parent is original mesh
       
       nn = memory_size(inte % shapf,1_ip)
       call memory_copy  (memor_dom,'INTE % LELEM',vacal,inte % lelem,lelem_cpy)
       call memory_copy  (memor_dom,'INTE % SHAPF',vacal,inte % shapf,shapf_cpy)
       call memory_deallo(memor_dom,'INTE % LELEM',vacal,inte % lelem)
       call memory_deallo(memor_dom,'INTE % SHAPF',vacal,inte % shapf)
       call memory_alloca(memor_dom,'INTE % LELEM',vacal,inte % lelem,   mesh_out % npoin)
       call memory_alloca(memor_dom,'INTE % SHAPF',vacal,inte % shapf,nn,mesh_out % npoin)
             
       do ipoin = 1,mesh_out % npoin
          kpoin                   = mesh_out % permn(ipoin)
          inte % lelem(ipoin)     = lelem_cpy(kpoin)
          inte % shapf(:,ipoin)   = shapf_cpy(:,kpoin)
          mesh_out % permn(ipoin) = ipoin
       end do
       do ielem = 1,mesh_out % nelem
          mesh_out % perme(ielem) = ielem
       end do
       
       call memory_deallo(memor_dom,'INTE % LELEM',vacal,lelem_cpy)
       call memory_deallo(memor_dom,'INTE % SHAPF',vacal,shapf_cpy)
       
    end if
    
  end subroutine witness_mesh_box

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a 3D circular mesh
  !> @details Create a p3D circular mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_interpolation(mesh_out,mesh_in,inte,OUTSIDE_VALUE)
    
    type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    type(typ_interp),                intent(inout) :: inte
    integer(ip),           optional, intent(in)    :: OUTSIDE_VALUE
    integer(ip)                                    :: ipoin,ifoun
    real(rp)                                       :: toler
    real(rp)                                       :: coloc(3)
    real(rp)                                       :: dNds(mesh_in % ndime,mesh_in % mnode)
    
    ifoun = 0
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'INTE % LELEM',vacal,inte % lelem,mesh_out % npoin)
       call memory_alloca(memor_dom,'INTE % SHAPF',vacal,inte % shapf,mesh_in  % mnode,mesh_out % npoin)
       do ipoin = 1,mesh_out % npoin
          call elsest_host_element(&
               ielse,relse,1_ip,mesh_in,mesh_out % coord(:,ipoin),inte % lelem(ipoin),&
               inte % shapf(:,ipoin),dNds,coloc,toler)
          if( inte % lelem(ipoin) == 0 ) then
             ifoun = ifoun + 1
             if( present(OUTSIDE_VALUE) ) then
                inte % lelem(ipoin) = OUTSIDE_VALUE
             else
                inte % lelem(ipoin) = 1
             end if
             inte % shapf(:,ipoin) = 0.0_rp
          end if
       end do
    end if
    call PAR_SUM(ifoun)
    if( ifoun /= 0 ) call messages_live('SOME WITNESS MESH POINT DO NOT HAVE ELEMENTS= '//integer_to_string(ifoun),'WARNING')

  end subroutine witness_mesh_interpolation

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a 3D circular mesh
  !> @details Create a p3D circular mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_3Dring(mesh_out,mesh_in,param,inte)

    type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    real(rp),                        intent(in)    :: param(:)
    type(typ_interp),                intent(inout) :: inte
    integer(ip)                                    :: nelem_rad,nelem_cir
    integer(ip)                                    :: knode,pnode,ielem
    integer(ip)                                    :: inode,ipoin,nn,kpoin
    real(rp)                                       :: center(3)
    real(rp)                                       :: normal(3)
    real(rp)                                       :: radius
    logical(lg),           pointer                 :: mask(:)
    integer(ip),           pointer                 :: lelem_cpy(:)
    integer(ip),           pointer                 :: lemin(:)
    real(rp),              pointer                 :: shapf_cpy(:,:)
    type(mesh_type_basic)                          :: mesh_tmp
    
    nullify(mask,lelem_cpy,shapf_cpy,lemin)

    nelem_rad   = int(param(8),ip)
    nelem_cir   = int(param(9),ip)
    radius      = param(1)
    center(1:3) = param(2:4)
    normal(1:3) = param(5:7)
    !
    ! Create mesh
    !
    if( INOTMASTER ) then
       call mesh_tmp % init('WITNESS MESH')
       call mesh_tmp % mesh_3Dring(center,normal,radius,nelem_rad,nelem_cir)
    end if
    !
    ! Compute interpolation array
    !
    call witness_mesh_interpolation(mesh_tmp,mesh_in,inte,OUTSIDE_VALUE=-1_ip)
    !
    ! Check repeated points
    !
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'LEMIN',vacal,lemin,mesh_tmp % npoin)
       do ipoin = 1,mesh_tmp % npoin
          if( inte % lelem(ipoin) > 0 ) lemin(ipoin) = kfl_paral
       end do
       call PAR_MAX(lemin,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
       do ipoin = 1,mesh_tmp % npoin
          if( lemin(ipoin) /= kfl_paral ) inte % lelem(ipoin) = -1
       end do
       call memory_deallo(memor_dom,'LEMIN',vacal,lemin)
    end if
    !
    ! Limit mesh to element owners
    !
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'MASK',vacal,mask,mesh_tmp % nelem)
       do ielem = 1,mesh_tmp % nelem
          mask(ielem) = .true.
          pnode = element_type(mesh_tmp % ltype(ielem)) % number_nodes
          knode = 0
          do inode = 1,pnode
             ipoin = mesh_tmp % lnods(inode,ielem)
             if( inte % lelem(ipoin) == -1 ) knode = knode + 1
          end do
          if( knode > 0 ) mask(ielem) = .false.
       end do
       call mesh_out % extract(mesh_tmp,mask)

       call mesh_tmp % deallo ()
       call memory_deallo     (memor_dom,'MASK',vacal,mask)
       
       nullify(mesh_out % parent ) ! To tell parent is original mesh

       nn = memory_size(inte % shapf,1_ip)
       call memory_copy  (memor_dom,'INTE % LELEM',vacal,inte % lelem,lelem_cpy)
       call memory_copy  (memor_dom,'INTE % SHAPF',vacal,inte % shapf,shapf_cpy)
       call memory_deallo(memor_dom,'INTE % LELEM',vacal,inte % lelem)
       call memory_deallo(memor_dom,'INTE % SHAPF',vacal,inte % shapf)
       call memory_alloca(memor_dom,'INTE % LELEM',vacal,inte % lelem,   mesh_out % npoin)
       call memory_alloca(memor_dom,'INTE % SHAPF',vacal,inte % shapf,nn,mesh_out % npoin)
       
       do ipoin = 1,mesh_out % npoin
          kpoin                   = mesh_out % permn(ipoin)
          inte % lelem(ipoin)     = lelem_cpy(kpoin)
          inte % shapf(:,ipoin)   = shapf_cpy(:,kpoin)
          mesh_out % permn(ipoin) = ipoin
       end do
       do ielem = 1,mesh_out % nelem
          mesh_out % perme(ielem) = ielem
       end do
       
       call memory_deallo(memor_dom,'INTE % LELEM',vacal,lelem_cpy)
       call memory_deallo(memor_dom,'INTE % SHAPF',vacal,shapf_cpy)

    end if
    
  end subroutine witness_mesh_3Dring
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a 3D circular mesh
  !> @details Create a p3D circular mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_mesh_3Dcircle(mesh_out,mesh_in,param,inte)

    type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    real(rp),                        intent(in)    :: param(:)
    type(typ_interp),                intent(inout) :: inte
    integer(ip)                                    :: nelem_rad,nelem_cir
    integer(ip)                                    :: knode,pnode,ielem
    integer(ip)                                    :: inode,ipoin,nn,kpoin
    real(rp)                                       :: center(3)
    real(rp)                                       :: normal(3)
    real(rp)                                       :: radius
    logical(lg),           pointer                 :: mask(:)
    integer(ip),           pointer                 :: lelem_cpy(:)
    real(rp),              pointer                 :: shapf_cpy(:,:)
    type(mesh_type_basic)                          :: mesh_tmp

    nullify(mask,lelem_cpy,shapf_cpy)
    
    nelem_rad   = int(param(8),ip)
    nelem_cir   = int(param(9),ip)
    radius      = param(1)
    center(1:3) = param(2:4)
    normal(1:3) = param(5:7)
    !
    ! Create mesh
    !
    if( INOTMASTER ) then
       call mesh_tmp % init('WITNESS MESH')
       call mesh_tmp % mesh_3Dcircle(center,normal,radius,nelem_rad,nelem_cir)
    end if
    !
    ! Compute interpolation array
    !
    call witness_mesh_interpolation(mesh_tmp,mesh_in,inte,OUTSIDE_VALUE=-1_ip)
    !
    ! Limit mesh to element owners
    !
    if( INOTMASTER ) then
       call memory_alloca(memor_dom,'MASK',vacal,mask,mesh_tmp % nelem)
       do ielem = 1,mesh_tmp % nelem
          mask(ielem) = .true.
          pnode = element_type(mesh_tmp % ltype(ielem)) % number_nodes
          knode = 0
          do inode = 1,pnode
             ipoin = mesh_tmp % lnods(inode,ielem)
             if( inte % lelem(ipoin) == -1 ) knode = knode + 1
          end do
          if( knode > 0 ) mask(ielem) = .false.
       end do
       call mesh_out % extract(mesh_tmp,mask)
       call mesh_tmp % deallo ()
       call memory_deallo     (memor_dom,'MASK',vacal,mask)
       nullify(mesh_out % parent ) ! To tell parent is original mesh

       nn = memory_size(inte % shapf,1_ip)
       call memory_copy  (memor_dom,'LELEM_CPY'   ,vacal,inte % lelem,lelem_cpy)
       call memory_copy  (memor_dom,'SHAPF_CPY'   ,vacal,inte % shapf,shapf_cpy)
       call memory_deallo(memor_dom,'INTE % LELEM',vacal,inte % lelem)
       call memory_deallo(memor_dom,'INTE % SHAPF',vacal,inte % shapf)
       call memory_alloca(memor_dom,'INTE % LELEM',vacal,inte % lelem,   mesh_out % npoin)
       call memory_alloca(memor_dom,'INTE % SHAPF',vacal,inte % shapf,nn,mesh_out % npoin)
       
       do ipoin = 1,mesh_out % npoin
          kpoin                   = mesh_out % permn(ipoin)
          inte % lelem(ipoin)     = lelem_cpy(kpoin)
          inte % shapf(:,ipoin)   = shapf_cpy(:,kpoin)
          mesh_out % permn(ipoin) = ipoin
       end do
       do ielem = 1,mesh_out % nelem
          mesh_out % perme(ielem) = ielem
       end do
       
       call memory_deallo(memor_dom,'LELEM_CPY'   ,vacal,lelem_cpy)
       call memory_deallo(memor_dom,'SHAPF_CPY'   ,vacal,shapf_cpy)
    end if
        
  end subroutine witness_mesh_3Dcircle

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Create a plane cartesian mesh
  !> @details Create a plane cartesian mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine witness_coarse_boundary(mesh_out,mesh_in,param,inte)

    type(mesh_type_basic),           intent(inout) :: mesh_out
    type(mesh_type),                 intent(in)    :: mesh_in
    real(rp),                        intent(in)    :: param(:)
    type(typ_interp),                intent(inout) :: inte
    integer(ip)                                    :: pelty,pnode
    integer(ip)                                    :: ifoun,ipoin,inode
    integer(ip)                                    :: ielem,iboun,kfoun
    integer(ip)                                    :: kpoin
    real(rp)                                       :: coloc(3),dNds(ndime,mnode)
    integer(ip)                                    :: npoin_boun(2)
    real(rp)                                       :: rpoin,toler
    type(mesh_type_basic)                          :: mesh_bou

    type(typ_kdtree)                               :: kdtree
    real(rp)                                       :: elcod(mesh_in % ndime,mnode)
    real(rp)                                       :: dista,proje(3)
    real(rp),              pointer                 :: coord_host(:,:)

    nullify(coord_host)
    !
    ! Create basic mesh MESH_BOU from boundary mesh of MESH_IN
    !
    call mesh_coherent_boundary_mesh(mesh_in,mesh_bou,INCLUDE_HALOS=.false.)
    npoin_boun = mesh_number_boundary_nodes(meshe(ndivi))
    rpoin      = param(1) 
    ipoin      = int(real(npoin_boun(1),rp)*rpoin/real(max(npoin_boun(2),1_ip),rp),ip)
    call surface_mesh_coarsening(mesh_bou,mesh_out,'npoin_coarse',ipoin)
    mesh_out % name = mesh_in % name
    call mesh_bou % deallo()

    if( mesh_out % npoin > 0 .and. INOTMASTER ) then
       !
       ! Allocate interpolation arrays
       !
       call memory_alloca(memor_dom,'INTE % LELEM',vacal,inte % lelem,mesh_out % npoin)
       call memory_alloca(memor_dom,'INTE % SHAPF',vacal,inte % shapf,mesh_in  % mnode,mesh_out % npoin)
       !
       ! Find minimum distance using kdtree
       !
       call memory_alloca(memor_dom,'COORD_HOST',vacal,coord_host,mesh_out % ndime,mesh_out % npoin)
       call kdtree_initialize(kdtree)
       call kdtree_construct(mesh_in % nboun,mesh_in % npoin,mesh_in % lnodb,mesh_in % ltypb,mesh_in % coord,kdtree)
       do ipoin = 1,mesh_out % npoin
          call kdtree_nearest_boundary(mesh_out % coord(:,ipoin),kdtree,iboun,dista,proje)
          coord_host(1:mesh_out % ndime,ipoin) = proje(1:mesh_out % ndime)
          inte % lelem(ipoin)                  = mesh_in % lelbo(iboun)
       end do
       call kdtree_deallocate(kdtree)
       !
       ! Compute interpolation array
       !
       ifoun = 0
       if( 1 == 2 ) then
          do ipoin = 1,mesh_out % npoin
             call elsest_host_element(&
                  ielse,relse,1_ip,meshe(ndivi),mesh_out % coord(:,ipoin),inte % lelem(ipoin),&
                  inte % shapf(:,ipoin),dNds,coloc,toler)
             if( inte % lelem(ipoin) == 0 ) then
                ifoun = ifoun + 1
                inte % lelem(ipoin)   = 1
                inte % shapf(:,ipoin) = 0.0_rp
             end if
          end do
       else
          do ipoin = 1,mesh_out % npoin
             ielem = inte % lelem(ipoin) 
             pelty = abs(mesh_in % ltype(ielem))
             pnode = element_type(pelty) % number_nodes
             do inode = 1,pnode
                kpoin = mesh_in % lnods(inode,ielem)
                elcod(1:mesh_in % ndime,inode) = mesh_in % coord(1:mesh_in % ndime,kpoin)
             end do
             call elmgeo_natural_coordinates(                                     &
                  mesh_in % ndime,pelty,pnode,elcod,inte % shapf(:,ipoin),dNds,   &
                  coord_host(1:mesh_out % ndime,ipoin),coloc,kfoun,toler_opt=0.01_rp)
             if( kfoun == 0 ) then
                ifoun                 = ifoun + 1
                inte % lelem(ipoin)   = 1
                inte % shapf(:,ipoin) = 0.0_rp
             end if
          end do
       end if
    end if

    call PAR_SUM(ifoun)

    if( ifoun /= 0 ) call messages_live('SOME WITNESS MESH POINT DO NOT HAVE ELEMENTS= '//integer_to_string(ifoun),'WARNING')
    !
    ! Deallocate
    !
    call memory_deallo(memor_dom,'COORD_HOST',vacal,coord_host)

  end subroutine witness_coarse_boundary

end module mod_witness
!> @}
