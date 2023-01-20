!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_begite()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_begite
  ! NAME
  !    chm_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the transport
  !    equation
  ! USES
  !    chm_tittim
  !    chm_updbcs
  !    chm_inisol
  !    chm_updunk
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_master,        only : ITASK_BEGITE, itcou, modul, itinn
  use def_chemic,        only : kfl_dt_calc_CMC_chm, kfl_goite_chm, kfl_model_chm, kfl_overs_chm, kfl_solve_cond_CMC_chm,&
                                kfl_under_chm
  use def_kintyp,        only : ip
  use mod_messages,      only : livinf
  implicit none

  external :: chm_tistep
  external :: chm_updbcs
  external :: chm_updunk

  !
  ! Initializations
  !
  kfl_goite_chm = 1
  itinn(modul)  = 0
  kfl_under_chm = 0
  kfl_overs_chm = 0

  if(itcou==1) then
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        if (kfl_dt_calc_CMC_chm == 1_ip)  call chm_tistep
     else
        call chm_tistep()
     end if
  end if

  call livinf(15_ip,' ',modul)

  if (.not. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
     !
     ! Update boundary conditions
     !
     call chm_updbcs(ITASK_BEGITE)

     !
     ! Obtain the initial guess for inner iterations
     !
     call chm_updunk(ITASK_BEGITE)
  end if
end subroutine chm_begite

