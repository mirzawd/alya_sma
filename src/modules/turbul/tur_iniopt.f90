!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_iniopt()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the adjoint equations of the turbulene equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_turbul
  use def_master
  use def_kermod, only : kfl_dvar_type,sens_mesh,kfl_cos_opt,costf
  use def_solver
  implicit none

  
  if (kfl_cos_opt == 1) costf = 0.0_rp
  if (kfl_dvar_type == 6) then
    sens_mesh = 0.0_rp
  endif
 
end subroutine tur_iniopt
