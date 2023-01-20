!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_updtss.f90
!> @author  houzeaux
!> @date    2020-10-21
!> @brief   Compute time step
!> @details Compute time step
!> @} 
!-----------------------------------------------------------------------

subroutine gus_updtss()

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  implicit none
  real(rp) :: dtmin
  
  if( kfl_timei_gus /= 0 ) then
     !
     ! Critical time step
     !
     dtmin = 1.0_rp
     !
     ! Update time step
     !
     dtcri_gus = dtmin
     if( dtcri_gus /= 0.0_rp ) dtinv_gus = 1.0_rp/(dtcri_gus*safet_gus)
     if( kfl_timco == 1      ) dtinv = max(dtinv,dtinv_gus)
  end if

end subroutine gus_updtss
