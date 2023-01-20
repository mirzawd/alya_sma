!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin 
!> @{
!> @file    ker_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!>          1. Allocate a minimum memory so that all ranks can enter 
!>             COU_INTERPOLATE_NODAL_VALUES without blowing up
!>             (see ker_membcs.f90 as an example)
!> @}
!------------------------------------------------------------------------
subroutine ker_plugin(icoup)

  use def_kintyp,    only :  ip,rp
  use def_coupli,    only :  coupling_type
  use def_coupli,    only :  UNKNOWN
  use def_coupli,    only :  RESIDUAL
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_matrix,    only :  matrix_initialize
  use mod_parall,    only :  I_AM_IN_COLOR
  !
  ! Possible variables => 
  ! 
  use def_domain,    only :  walld
  use def_master,    only :  INOTMASTER
  use def_master,    only :  gesca
  use def_domain,    only :  npoin
  !use mod_projec
  implicit none
  integer(ip) :: ipoin
  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup    !< Coupling number
  character(5)            :: variable
  integer(ip)             :: color_target

  variable     = coupling_type(icoup) % variable 
  color_target = coupling_type(icoup) % color_target
  !
  ! Wall distance 
  ! 
  if( variable == 'WALLD' ) then 

     if( INOTMASTER .and. I_AM_IN_COLOR(color_target) ) then
        call memgen(0_ip,npoin,0_ip)
     else
        allocate(gesca(1))
     end if

     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,gesca,walld)     

     if( INOTMASTER .and. I_AM_IN_COLOR(color_target) ) then
        if( .not. associated(walld) ) call runend('KER_PLUGIN: WALLD NOT ALLOCATED')
        do ipoin = 1,npoin
           walld(ipoin) = gesca(ipoin)
        end do
        call memgen(2_ip,npoin,0_ip) 
     else
        deallocate(gesca)
     end if
  end if

end subroutine ker_plugin

