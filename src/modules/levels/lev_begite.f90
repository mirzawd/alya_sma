!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_begite
  !-----------------------------------------------------------------------
  !****f* Levels/lev_begite
  ! NAME 
  !    lev_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the level set equation
  !    equation
  ! USES
  !    lev_tittim
  !    lev_updbcs
  !    lev_inisol
  !    lev_updunk
  ! USED BY
  !    lev_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_levels
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_lev = 1 
  itinn(modul)  = 0
  if((kfl_timco_lev==1).and.(kfl_timco==0)) then
     if(dtinv>1.0d-10) dtime = 1.0_rp/dtinv
  endif
  call livinf(15_ip,' ',modul)
  if( itcou == 1 ) call lev_tistep()
  !
  ! Set up the solver parameters for the level set equation
  !
  call lev_inisol(1_ip)
  !
  ! Update boundary conditions
  !
  call lev_updbcs()
  !
  ! Obtain the initial guess for inner iterations
  !
  call lev_updunk(2_ip)

end subroutine lev_begite
