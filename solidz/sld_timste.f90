!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_timste.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine prepares for a new time step
!> @details Calculation of the time step and updates on the safety factor
!>          when is necessary.
!> @}
!-----------------------------------------------------------------------

subroutine sld_timste

  use def_kintyp,  only : ip, rp
  use def_master,  only : ittim, dtinv
  use def_master,  only : kfl_timco
  use def_solidz,  only : kfl_stead_sld, kfl_rigid_sld, kfl_savdt_sld
  use def_solidz,  only : safet_sld, safex_sld, safma_sld
  use def_solidz,  only : kfl_safet_table_sld
  use def_solidz,  only : safet_table_sld, resid_sld, nisaf_sld
  use def_solidz,  only : dtcri_sld, dtinv_sld
  
  implicit none

  integer(ip)          :: iauxi

  if( kfl_rigid_sld == 0_ip ) then

     !-------------------------------------------------------------------
     !
     ! Deformable body
     !
     !-------------------------------------------------------------------
     !
     ! Safety factor updates (do not modify at initial or user-defined time step)
     !
     if( ittim > nisaf_sld ) safet_sld = min(safet_sld*safex_sld, safma_sld)
     !
     ! Safety factor using a discrete table
     !
     if( kfl_safet_table_sld > 0 ) then
        safet_sld = safet_table_sld(1,1)
        if (ittim > 1) then
           do iauxi= 1,kfl_safet_table_sld
              if (resid_sld < safet_table_sld(2,iauxi)) then
                 safet_sld = safet_table_sld(1,iauxi)
                 cycle
              end if
           end do
        end if
     end if
     !
     ! Time step size
     !
     if( kfl_stead_sld /= 1_ip ) then
        if( kfl_savdt_sld == 1_ip ) then
           if( dtcri_sld /= 0.0_rp ) dtinv_sld = 1.0_rp/(dtcri_sld*safet_sld)
           if( kfl_timco == 1 )      dtinv     = max(dtinv,dtinv_sld)
        else
           call sld_updtss()
        end if
     end if
     
  else if( kfl_rigid_sld == 1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Rigid body
     !
     !-------------------------------------------------------------------
     !
     ! Time step size
     !
     if( kfl_stead_sld /= 1_ip ) then
        dtcri_sld = 1.0_rp
        dtinv_sld = 1.0_rp
        if( kfl_timco == 1 ) dtinv = max(dtinv,dtinv_sld)
     end if
     
  end if

end subroutine sld_timste

