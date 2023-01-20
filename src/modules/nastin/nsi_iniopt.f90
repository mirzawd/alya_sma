!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_iniopt()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the adjoint equations of the NS equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_nastin
  use def_master
  use def_kermod, only : kfl_adj_prob,kfl_cos_opt,sens_mesh,costf,kfl_dvar_type
!  use def_kermod, only : sens
  use def_solver
  implicit none

    
  if (kfl_cos_opt == 1) costf = 0.0_rp
  
  if (kfl_adj_prob == 1) then
!     resdiff_nsi = 0.0_rp 
    dcost_dx_nsi = 0.0_rp
  endif
  
  if (kfl_dvar_type == 5) then
    sens_mesh = 0.0_rp
  endif
  
!   if ( kfl_coupl(ID_NASTIN,ID_TURBUL) == 0 ) then
!     if (kfl_dvar_type == 5) then
!       sens_mesh = 0.0_rp
!     else
!       sens = 0.0_rp
!     endif
!   endif
 
end subroutine nsi_iniopt
 
