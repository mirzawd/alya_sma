!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_doiter
  !-----------------------------------------------------------------------
  !****f* Levels/lev_doiter
  ! NAME 
  !    lev_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the level set equation.
  ! USES
  !    lev_begite
  !    lev_solite
  !    lev_endite
  ! USED BY
  !    LEVELS
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_levels
  use def_master
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  implicit none

  if(kfl_stead_lev==0) then
     call timings_ini()
     call lev_begite
     call timings_end(ITASK_BEGITE)
     do while(kfl_goite_lev==1)
        call lev_solite
        call lev_endite(one)
     end do
     call timings_ini()
     call lev_endite(two)
     call timings_end(ITASK_ENDITE)
  end if
  !
  ! Move mesh
  !
  !call lev_move_mesh_to_free_surface()

end subroutine lev_doiter
