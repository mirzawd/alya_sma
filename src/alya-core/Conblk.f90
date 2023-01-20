!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Conblk()
  !-----------------------------------------------------------------------
  !****f* master/Conblk
  ! NAME
  !    Alya
  ! DESCRIPTION
  !    Increase block number
  ! USES
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,     only : ip
  use def_master,     only : ITASK_CONBLK
  use def_master,     only : iblok
  use def_master,     only : nblok
  use def_master,     only : itcou
  use def_master,     only : kfl_goblk
  use def_master,     only : kfl_gocou
  use def_coupli,     only : mcoup
  use def_coupli,     only : kfl_gozon
  use def_coupli,     only : coupling_driver_iteration
  use def_coupli,     only : coupling_driver_number_couplings
  use mod_messages,   only : livinf
  use mod_moduls,     only : moduls   
  use mod_ker_updpro, only : ker_updpro
  use mod_reset,      only : reset_check

  implicit none
  !
  ! Write message if block is in a coupling loop
  !
  if( mcoup > 0 ) then
     if( coupling_driver_number_couplings(iblok) /= 0 ) then
        call livinf(-13_ip,'END ZONAL COUPLING: ',coupling_driver_iteration(iblok))
     end if
  end if
  !
#ifdef COMMDOM
  call moduls(ITASK_CONBLK)
#endif 
  ! 
  ! Initialize
  !   
  coupling_driver_iteration(iblok) = 0
  kfl_gozon = 1
  iblok     = iblok+1 
  kfl_gocou = 1
  itcou     = 1
  if( iblok > nblok ) then
     kfl_goblk = 0
  end if
  if( nblok > 1 ) then
     call livinf(8_ip,' ',0_ip)
  end if

  call ker_updpro(ITASK_CONBLK) ! update property
  !
  ! Check reset
  !
  call reset_check()

end subroutine Conblk
