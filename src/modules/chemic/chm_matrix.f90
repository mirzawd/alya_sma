!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_matrix()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_matrix
  ! NAME
  !    chm_matrix
  ! DESCRIPTION
  !    This routine assembles matrix and RHS
  ! USES
  !    chm_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    chm_solite
  !***
  !-----------------------------------------------------------------------
  use def_master,   only : INOTMASTER, pmatr
  use def_chemic,   only : dtinv_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, cputi_chm
  use def_kintyp,   only : ip, rp
  use def_solver,   only : solve_sol
  use mod_timings,  only : timings_assembly
  implicit none
  real(rp)    :: time1,time2,time3,time4,dtinv_tmp,time_elem,time_boun

  external    :: inisol
  external    :: cputim
  external    :: chm_massma
  external    :: chm_elmope_all
  external    :: chm_bouope

  !
  ! Solver initializations and preconditioners
  !
  call inisol()

  solve_sol(1) % xdiag = 1.0_rp/dtinv_chm
  if (.not. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
     call chm_massma(pmatr) !!!!!!!!!!!!! IS THIS ONLY FOR IMPLICIT? ¿PROVISIONAL?
  end if

  call cputim(time1)
  !
  ! Element assembly
  !
  if( INOTMASTER ) then

     call chm_elmope_all(1_ip)

  else
     dtinv_chm = 1.0e6_rp
  end if
  call cputim(time2)

  !
  ! Boundary assembly
  !
  call cputim(time3)
  if( INOTMASTER ) call chm_bouope()
  call cputim(time4)

  !
  ! Split assembly
  !
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_tmp = dtinv_chm
     dtinv_chm = 0.0_rp
  end if
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_chm = dtinv_tmp
  end if

  time_elem = time2-time1
  time_boun = time4-time3
  call timings_assembly(time_elem,time_boun,TYPE_OF_ASSEMBLY='ELEMENT, BOUNDARY')
  cputi_chm(1) = cputi_chm(1) + time_elem
  cputi_chm(2) = cputi_chm(2) + time_boun

end subroutine chm_matrix

