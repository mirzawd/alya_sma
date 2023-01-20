!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_doiter()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_doiter
  ! NAME 
  !    nsi_doiter
  ! DESCRIPTION
  !    This routine solves an iteration of the linearized incompressible NS
  !    equations.
  ! USES
  !    nsi_begite
  !    nsi_updbcs
  !    nsi_solite
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_nastin
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  implicit none
  
  if( kfl_stead_nsi == 0 ) then
     call timings_ini()
     call nsi_begite()
     call timings_end(ITASK_BEGITE)

     do while( kfl_goite_nsi == 1 )
        call nsi_solite()
        call nsi_endite(ITASK_ENDINN)
     end do

     call timings_ini()
     call nsi_endite(ITASK_ENDITE)
     call timings_end(ITASK_ENDITE)

  end if

end subroutine nsi_doiter

