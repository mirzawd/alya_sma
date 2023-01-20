!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_doiter.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Solve
!> @details Solve Gusano
!> @} 
!-----------------------------------------------------------------------

subroutine gus_doiter()

  use def_master
  use def_gusano
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  
  implicit none

  if( kfl_stead_gus == 0 ) then

     call timings_ini()
     call gus_begite()
     call timings_end(ITASK_BEGITE)
     
     do while( kfl_goite_gus == 1 )
        call gus_solite()
        call gus_endite(ITASK_ENDINN)
     end do

     call timings_ini()
     call gus_endite(ITASK_ENDITE)
     call timings_end(ITASK_ENDITE)

  end if
  
end subroutine gus_doiter
