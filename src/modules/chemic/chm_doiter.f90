!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_doiter()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_doiter
  ! NAME
  !    chm_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the transport equation.
  ! USES
  !    chm_begite
  !    chm_solite
  !    chm_endite
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_master,          only : ITASK_BEGITE, ITASK_ENDITE, momod, modul
  use def_chemic,          only : kfl_goite_chm, kfl_spray_chm, kfl_tisch_chm
  use mod_chm_explicit,    only : chm_explicit_solution
  use mod_chm_rk_explicit, only : chm_rk_explicit_solution
  use mod_chm_levSet,      only : chm_pseudo_time_reinit_levSet
  use mod_timings,         only : timings_ini
  use mod_timings,         only : timings_end
  implicit none

  external :: chm_begite
  external :: chm_endite
  external :: soldef

  if(momod(modul)%kfl_stead == 0) then

     call timings_ini()
     call chm_begite()
     call timings_end(ITASK_BEGITE)

     do while( kfl_goite_chm==1 )

        if (kfl_tisch_chm == 3) then
           call chm_explicit_solution()         ! Solve AB scheme
        else
           if(kfl_spray_chm == 2) then
              call chm_pseudo_time_reinit_levSet()

              call chm_rk_explicit_solution()   ! Solve surface density PDE

           else
              call chm_rk_explicit_solution()   ! Solve RK scheme
           end if
        endif

     end do

     call timings_ini()
     call chm_endite(ITASK_ENDITE)
     call timings_end(ITASK_ENDITE)

  end if


end subroutine chm_doiter
