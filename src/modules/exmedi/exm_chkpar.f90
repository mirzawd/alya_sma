!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_chkpar
  !-----------------------------------------------------------------------
  !****f* exmedi/exm_chkpar
  ! NAME 
  !    exm_chkpar
  ! DESCRIPTION
  !    This routine checks the master-slave-alone status.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_exmedi
  implicit none


  !checking...
  if (kfl_paral > 0) then           ! slaves
     imaster= .false.
     islave = .true.
  else if (kfl_paral == 0) then     ! master
     imaster= .true.
     islave = .false.
  else if (kfl_paral < 0) then      ! alone
     imaster= .false.
     islave = .false.
  end if
  
  weparal= .false.
  if (imaster .or. islave) weparal= .true.

end subroutine exm_chkpar
