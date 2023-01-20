!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_begset.f90
!> @author  Guillaume Houzeaux
!> @date    13/12/2019
!> @brief   Begin a time step
!> @details Begin a time step of the coupling
!> @} 
!-----------------------------------------------------------------------

subroutine cou_begste()

  use def_kintyp,    only : ip,rp
  use def_master,    only : itti2,nblok
  use def_domain,    only : voave,ndime
  use def_coupli,    only : kfl_absolute_cou,toler_absolute_cou
  use def_coupli,    only : mcoup
  use def_coupli,    only : coupling_driver_iteration
  use mod_couplings, only : COU_TEMPORAL_PREDICTOR
  implicit none
  integer(ip) :: icoup
  !
  ! Temporal predction for zonal coupling (only call it after two time steps and before the modules)
  !
  if( itti2 > 2 ) then 
     do icoup = 1_ip, mcoup
        call COU_TEMPORAL_PREDICTOR(icoup)
     end do
  end if
  !
  ! Compute the absolute tolerance aas a function of the mesh size
  !
  if( kfl_absolute_cou < 0 .and. voave > 0.0_rp ) then
     toler_absolute_cou = abs(real(kfl_absolute_cou,rp)) * voave**(1.0_rp/real(ndime,rp))
  end if
  !
  ! Coupling: Put counters to zero
  ! 
  if( mcoup > 0 ) then
     coupling_driver_iteration(1:nblok) = 0
  end if

end subroutine cou_begste
