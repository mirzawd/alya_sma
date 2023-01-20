!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_endste.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   End a time step
!> @details End a time step
!> @}
!------------------------------------------------------------------------

subroutine neu_endste()

  use def_parame
  use def_master
  use def_neutro
  implicit none
  !
  ! Compute convergence residual of the time evolution 
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  !if( kfl_stead_neu == 0 .and. ADR_neu % kfl_time_integration == 1 ) then
  !   call neu_cvgunk(3_ip)
  !   call neu_updunk(5_ip)
  !end if
  !
  ! Write restart file
  !
  !call neu_restar(2_ip)
  !
  ! If not steady, go on
  !
  !if( kfl_stead_neu == 0 .and. ADR_neu % kfl_time_integration == 1 .and. kfl_conve(modul) == 1 ) kfl_gotim = 1

end subroutine neu_endste
