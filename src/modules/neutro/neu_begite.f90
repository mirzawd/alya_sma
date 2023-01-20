!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_begite.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Begin inner iterations
!> @details Begin inner iterations
!> @}
!------------------------------------------------------------------------

 subroutine neu_begite()

  use def_parame
  use def_master
  use def_domain
  use def_neutro
  use mod_messages, only : livinf
  implicit none
  external :: neu_updunk
  !
  ! Initializations
  !
  call livinf(15_ip,' ',modul)
  kfl_goite_neu = 1
  itinn(modul)  = 0
  !if( momod(modul) % miinn == 0 ) kfl_goite_neu = 0
  
  if( miinn_neu == 0 ) kfl_goite_neu = 0
  
  !if( itcou == 1 ) call neu_tistep()
  !
  ! Obtain the initial guess for inner iterations
  !
  call neu_updunk(2_ip) 

end subroutine neu_begite



