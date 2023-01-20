!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_doiter()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_doiter
  ! NAME 
  !    tem_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the temperature equation.
  ! USES
  !    tem_begite
  !    tem_solite
  !    tem_endite
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_temper
  use mod_tem_explicit   , only : tem_explicit_solution 
  use mod_tem_rk_explicit, only : tem_rk_explicit_solution
  use mod_timings,         only : timings_ini
  use mod_timings,         only : timings_end
  use mod_ker_proper
  use def_kermod
  implicit none

  if( kfl_stead_tem == 0 ) then
     call timings_ini()
     call tem_begite()
     call timings_end(ITASK_BEGITE)
     do while( kfl_goite_tem == 1 )
        if (kfl_explicit_tem == 1) then
           if(kfl_tisch_tem==3) then 
              call tem_explicit_solution()
           else
              call tem_rk_explicit_solution()
           end if
        else 
           call tem_solite()
        end if
        call tem_endite(ITASK_ENDINN)

     end do
     call timings_ini()
     call tem_endite(ITASK_ENDITE)
     call timings_end(ITASK_ENDITE)
  end if
  
end subroutine tem_doiter
