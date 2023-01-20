!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    setlbe.f90
!> @author  houzeaux
!> @date    2020-10-26
!> @brief   This routine constucts boundary arrays
!> @details This routine constucts boundary arrays:
!>          LBOEL(MNODB,NBOUN) ... Boundary/Element connectivity
!>          LELBO(NBOUN) ......... Boundary/ element 
!>          LTYPB(NBOUN) ......... Type of boundary
!> @} 
!-----------------------------------------------------------------------

subroutine setlbe()

  use def_kintyp,    only : ip
  use def_domain,    only : nboun
  use def_domain,    only : nelem
  use def_domain,    only : lboel
  use def_domain,    only : lnods
  use def_domain,    only : lnodb
  use def_domain,    only : ltype
  use def_domain,    only : lelbo
  use mod_elmgeo,    only : elmgeo_order_boundary_nodes
  use mod_domain,    only : domain_memory_allocate
  use mod_mesh_type, only : mesh_type_update_last_mesh
  use mod_strings,   only : integer_to_string
  implicit none 
  integer(ip) :: ielem,pelty,iboun,istat
  !
  ! Calculate only the boundary/element nodal connectivity
  !
  call domain_memory_allocate('LBOEL')
  do iboun = 1,nboun
     ielem = lelbo(iboun)
     if( ielem < 1 .or. ielem > nelem ) then
        call runend('SETLBE: BOUNDARY/ELEMENT CONNECTIVITY IS WRONG FOR BOUNDARY '//integer_to_string(iboun)//'= '//integer_to_string(lelbo(iboun)))
     end if
     pelty = abs(ltype(ielem))
     call elmgeo_order_boundary_nodes(pelty,lnods(:,ielem),lnodb(:,iboun),lboel(:,iboun),istat)
  end do
  !
  ! Repoint to account for lboel
  !
  call mesh_type_update_last_mesh() 

end subroutine setlbe
