!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin 
!> @{
!> @file    gus_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!>          1. Allocate a minimum memory so that all ranks can enter 
!>             COU_INTERPOLATE_NODAL_VALUES without blowing up
!>             (see gus_membcs.f90 as an example)
!> @}
!------------------------------------------------------------------------
subroutine gus_plugin(icoup)

  use def_kintyp,                only : ip,rp
  use def_master,                only : press
!  use def_master,                only : lninv_loc
  use def_domain,                only : kfl_codno
  use def_domain,                only : npoin
  use def_coupli,                only : coupling_type
  use def_coupli,                only : SCALAR
  use mod_communications_global, only : PAR_BROADCAST
  use mod_communications_global, only : PAR_MAX
  use mod_couplings,             only : cou_residual_scalar
  !
  ! Possible variables => 
  ! 
  use def_gusano,                only :  bvess_gus
  !
  ! <= end coupling variables
  !
  implicit none
  integer(ip), intent(in) :: icoup    !< Coupling number
  character(5)            :: variable
  integer(ip)             :: ipoin
  real(rp)                :: flowr_nsi
  real(rp)                :: press_gus

  variable = coupling_type(icoup) % variable 
   !if(kfl_paral==0) print*,'GUS ENTRA: ',icoup
 !
  ! Flow rate 
  ! 
  if( variable == 'FLOWR' .and. coupling_type(icoup) % what == SCALAR ) then
     flowr_nsi = -huge(1.0_rp)
     call PAR_MAX(flowr_nsi,'IN THE WORLD')
!!!     call PAR_BROADCAST(flowr_nsi,'IN THE WORLD')
   !  if(kfl_paral==0) print*,'GUS RECV FLOWR: ',icoup,flowr_nsi

     do ipoin = 1,npoin
        if( kfl_codno(1,ipoin) == coupling_type(icoup) % where_number ) then           
           bvess_gus(1,ipoin,1) = flowr_nsi
           ! print*,'impongo el Q en nodo',icoup,lninv_loc(ipoin),flowr_nsi  
          exit
        end if
     end do
  end if
  !
  ! Pressure
  ! 
  if( variable == 'PRESS' .and. coupling_type(icoup) % what == SCALAR ) then
     press_gus = -huge(1.0_rp)
     do ipoin = 1,npoin
        if( kfl_codno(1,ipoin) == coupling_type(icoup) % where_number ) then           
           press_gus = press(ipoin,1)
          !  print*,'calculo una presion=',icoup,lninv_loc(ipoin),press_gus
          exit
        end if
     end do
     call PAR_MAX(press_gus,'IN THE WORLD')
     !if(kfl_paral==0) print*,'GUS SEND PRESS: ',icoup,press_gus

!!!!     call PAR_BROADCAST(press_gus,'IN THE WORLD')
  end if

end subroutine gus_plugin

