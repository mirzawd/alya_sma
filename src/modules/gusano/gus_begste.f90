!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_begste.f90
!> @author  houzeaux
!> @date    2020-10-21
!> @brief   Begin of a time step
!> @details Begin of a time step
!> @} 
!-----------------------------------------------------------------------

subroutine gus_begste()

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  implicit none
  
  !if(kfl_stead_tem/=1) then     
     !
     ! Initial guess fo the temperature: T(n,0,*) <-- T(n-1,*,*).
     !
     call gus_updunk(ITASK_BEGSTE)
     !
     ! Update boundary conditions
     !
     call gus_updbcs(ITASK_BEGSTE)

  !end if
    
end subroutine gus_begste

