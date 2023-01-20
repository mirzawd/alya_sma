!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_endste()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_endste
  ! NAME 
  !    tem_endste
  ! DESCRIPTION
  !    This routine ends a time step of the temperature equation.
  ! USES
  !    tem_cvgunk
  !    tem_updunk
  !    tem_output
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use mod_tem_therm_press, only : therm_press_update
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  !
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_tem==0.and.kfl_timei_tem==1) then
     call tem_cvgunk(ITASK_ENDSTE)
     call tem_updunk(ITASK_ENDSTE)
     call therm_press_update(2_ip)         ! PRTHE_NSI: Thermodynamic pressure
  end if
  !
  ! Compute averaged variables
  !  
  call tem_averag()
  !
  ! Coupling with dynamic solver
  !
  call tem_dyncou(2_ip)
  !
  ! If not steady, go on
  !
  if(kfl_stead_tem==0.and.kfl_timei_tem==1.and.kfl_conve(modul)==1) kfl_gotim = 1

end subroutine tem_endste
