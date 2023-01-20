!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_timste
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_timste
  ! NAME 
  !    exm_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  !    exm_iniunk
  !    exm_updtss
  !    exm_updbcs
  !    exm_updunk
  ! USED BY
  !    Exmedi
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_exm_activation, only : exm_stim_compute

  implicit none


  !    Update applied currents appfi_exm (source terms) 
  !    using its definition function  and lapno_exm,  done in exm_comapp  

  if (kfl_paced_exm == 0) then
    !appfi_exm(1:npoin) = 0.0_rp
  else
     call exm_stim_compute() !calculate stimuli
  end if
  !
  ! Time step size.
  !
  call exm_updtss     

end subroutine exm_timste
