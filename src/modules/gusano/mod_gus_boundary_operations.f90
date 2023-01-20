!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 !-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gus_boundary_operations.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Gusano assembly of Neumann conditions
!> @details Assembly
!>
!-----------------------------------------------------------------------

module mod_gus_boundary_operations

  use def_kintyp_basic
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_output
  use def_gusano
  implicit none
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-23
  !> @brief   Impose Neumann condition
  !> @details Impose Neumann condition
  !> 
  !-----------------------------------------------------------------------

  subroutine gus_boundary_operations()

    integer(ip) :: iboun,ielem,ipoin

    do iboun = 1,nboun
       ielem = lelbo(iboun)
       if( kfl_fixbo_gus(iboun) == 2 ) then
          !
          ! Prescribe pressure
          !
          ipoin = lnodb(1,iboun)
          rhsid((ipoin-1)*2+1) = rhsid((ipoin-1)*2+1) - exn1d_gus(ipoin) * bvnat_gus(1,iboun,1)
       end if
    end do

  end subroutine gus_boundary_operations
  
end module mod_gus_boundary_operations
!> @}

