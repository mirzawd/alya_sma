!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_endste.f90
!> @author  houzeaux
!> @date    2020-10-21
!> @brief   End of a time step
!> @details End of a time step
!> @} 
!-----------------------------------------------------------------------

subroutine gus_endste()

  use def_parame
  use def_master
  use def_gusano
  
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  !
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if( kfl_stead_gus == 0 .and. kfl_timei_gus == 1 ) then
     call gus_cvgunk(ITASK_ENDSTE)
     call gus_updunk(ITASK_ENDSTE)
  end if
  !
  ! If not steady, go on
  !
  if(kfl_stead_gus == 0 .and. kfl_timei_gus == 1 .and. kfl_conve(modul) == 1 ) kfl_gotim = 1

end subroutine gus_endste
