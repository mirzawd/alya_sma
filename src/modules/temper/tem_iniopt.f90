!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_iniopt()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the adjoint equations of the temperature equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_temper
  use def_master
  use def_kermod, only : kfl_adj_prob,kfl_cos_opt,kfl_ndvars_opt,sens,costf
  use def_solver
  implicit none

  integer(ip)  :: idesvar
!  integer(ip)  :: izrhs
    
  if (kfl_cos_opt == 1) costf = 0.0_rp
  if (kfl_adj_prob == 1) then
  
      do idesvar = 1,kfl_ndvars_opt
	sens(idesvar) = 0.0_rp
	resdiff_tem(idesvar,:) = 0.0_rp
! 	do izrhs = 1,solve_sol(1) % nzrhs * solve_sol(1) % nseqn
! 	  resdiff_tem(idesvar,izrhs) = 0.0_rp
! 	end do
      end do
  
  endif
    
 
end subroutine tem_iniopt
 
