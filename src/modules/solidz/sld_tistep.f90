!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_tistep.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine sets the time step
!> @details This routine sets the time step
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_tistep

  use def_kintyp, only : ip, rp
  use def_master, only : kfl_timco, dtinv
  use def_master, only : routp, ioutp, lun_outpu
  use mod_outfor, only : outfor
  use def_solidz, only : kfl_timei_sld, kfl_stead_sld
  use def_solidz, only : dtcri_sld, dtinv_sld, safet_sld
  use def_solidz, only : SLD_STATIC_PROBLEM

  implicit none

  if ( kfl_timco /= 2_ip ) then
     !
     ! Prescribed or from critical
     !
     dtinv_sld = dtinv
     if ( kfl_timei_sld == SLD_STATIC_PROBLEM ) dtinv_sld = 0.0_rp
     if ( kfl_stead_sld == 1_ip ) dtinv_sld = 0.0_rp
     !
     ! If dt is prescribed, compute corresponding safety factor
     !
     if ( kfl_timco == 0_ip ) then
        if ( dtcri_sld /= 0.0_rp .and. dtinv_sld /= 0.0_rp ) then
           safet_sld = 1.0_rp/(dtcri_sld*dtinv_sld)
        else
           safet_sld = 1.0_rp
        end if
     end if
     
  end if

  !
  ! Actualize time integration parameters
  !
  routp(1) = dtcri_sld
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp

  ioutp(1) = kfl_timei_sld
  ioutp(2) = kfl_stead_sld


   call outfor(8_ip,lun_outpu,' ')

end subroutine sld_tistep
