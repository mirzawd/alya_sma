!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_endite.f90
!> @author  houzeaux
!> @date    2020-10-22
!> @brief   End of iterations
!> @details End of inner and outer iterations
!> @} 
!-----------------------------------------------------------------------

subroutine gus_endite(itask)

  use def_kintyp_basic, only : ip
  use def_master,       only : ITASK_ENDINN
  use def_master,       only : ITASK_ENDITE
  use def_master,       only : itinn
  use def_master,       only : modul
  use mod_messages
  
  implicit none
  
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( ITASK_ENDINN )

     call gus_cvgunk(ITASK_ENDINN)
     call gus_updunk(ITASK_ENDINN)
     call gus_cvgunk(0_ip)   

  case ( ITASK_ENDITE )

     call gus_cvgunk(ITASK_ENDITE)
     call gus_updunk(ITASK_ENDITE)
     call livinf(16_ip,' ',itinn(modul))
     
  end select
  
end subroutine gus_endite
