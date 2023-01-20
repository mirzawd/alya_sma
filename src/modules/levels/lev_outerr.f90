!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_outerr
  !------------------------------------------------------------------------
  !****f* Levels/lev_outerr
  ! NAME 
  !    lev_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !------------------------------------------------------------------------
  use      def_master
  use      def_levels
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  !
  ! Efficient redistancing is no longer working, should uncoment the 2 calls to elseap - I have chosen this manual solution because otherwise 
  ! elseap was beeing called in acses were it was not being needed because we do not have kfl_tyred_lev when it is called  
  !
  if( tyred_lev == 3 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'To use  tyred_lev == 3 uncoment !hhelseap & comment this')     
  end if

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,' ')

end subroutine lev_outerr
