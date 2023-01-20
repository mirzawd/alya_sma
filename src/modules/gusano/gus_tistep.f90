!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_tistep.f90
!> @author  houzeaux
!> @date    2020-10-26
!> @brief   Sets time step
!> @details This routine sets the time step DTINV_GUS
!> @} 
!-----------------------------------------------------------------------

subroutine gus_tistep

  use def_parame
  use def_master
  use def_gusano
  use mod_outfor,                only : outfor
  use mod_communications_global, only : par_max
  implicit none

  if( kfl_timco /= 2 ) then
     !
     ! Global time step
     !
     dtinv_gus = dtinv
     if( kfl_timei_gus == 0 ) dtinv_gus = 0.0_rp 
     if( kfl_stead_gus == 1 ) dtinv_gus = 0.0_rp

  else
     !
     ! local time step
     !
     call PAR_MAX(dtmax_gus)
     
  end if

  routp(1) = dtcri_gus
  ioutp(1) = kfl_timei_gus
  ioutp(2) = kfl_stead_gus
  !
  ! minimum time step when using local time step
  !
  routp(2) = safet_gus*dtcri_gus 
  !
  ! Look for maximum  (time step) over subdomains
  !
  routp(3) = dtmax_gus
  call outfor(8_ip,lun_outpu,' ')

end subroutine gus_tistep
