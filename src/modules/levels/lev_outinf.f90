!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_outinf
  !-----------------------------------------------------------------------
  !****f* levels/lev_outinf
  ! NAME 
  !    lev_outinf
  ! DESCRIPTION
  !    This routine writes info in the wave equation files
  ! USES
  !    lev_inisol
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_levels
  implicit none

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then

     end if
  end if

end subroutine lev_outinf
      
