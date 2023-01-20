!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_begite.f90
!> @author  houzeaux
!> @date    2020-10-22
!> @brief   End of iterations
!> @details End of inner and outer iterations
!> @} 
!-----------------------------------------------------------------------

subroutine gus_begite()
  
  use def_kintyp_basic, only : ip
  use def_gusano,       only : GUS_BEND_OFF
  use def_gusano,       only : kfl_goite_gus
  use def_gusano,       only : kfl_bendi_gus
  use def_gusano,       only : densi_gus
  use def_gusano,       only : visco_gus
  use mod_ker_proper,   only : ker_proper
  use def_master
  use mod_messages
  
  implicit none
  integer(ip) :: dummi
  
  call livinf(15_ip,' ',modul)
  !
  ! Initializations
  !
  kfl_goite_gus = 1
  itinn(modul)  = 0
  if( momod(modul) % miinn == 0 ) kfl_goite_gus = 0
  if( itcou == 1 ) call gus_tistep()
  !
  ! Boundary conditions
  !
  call gus_updbcs(ITASK_BEGITE)
  !
  ! Initial values
  !
  call gus_updunk(ITASK_BEGITE)
  !
  ! Properties on nodes
  !
  if( kfl_bendi_gus /= GUS_BEND_OFF ) then 
     call ker_proper('DENSI','NPOIN',dummi,dummi,densi_gus)
     call ker_proper('VISCO','NPOIN',dummi,dummi,visco_gus)
  end if
  
end subroutine gus_begite
