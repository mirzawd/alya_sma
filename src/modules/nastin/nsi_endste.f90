!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_endste()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_endste
  ! NAME 
  !    nsi_endste
  ! DESCRIPTION
  !    This routine ends a time step of the incompressible NS equations.
  ! USES
  !    nsi_cvgunk
  !    nsi_updunk
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_nastin
  use def_kermod, only : kfl_cos_opt,kfl_adj_prob
  use def_domain
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !    
  if( kfl_stead_nsi == 0 .and. kfl_timei_nsi == 1 ) then
     call nsi_cvgunk(ITASK_ENDSTE) ! Convergence
     call nsi_updunk(ITASK_ENDSTE) ! VELOC, PRESS
  end if
  !
  ! Compute averaged variables
  !
  call nsi_averag()
  !
  ! If not steady, go on
  !
  if(kfl_stead_nsi==0.and.kfl_timei_nsi==1.and.kfl_conve(modul)==1) kfl_gotim = 1
  !
  ! Calculate functional and sensitivities
  !
  if (kfl_cos_opt == 1) call nsi_costcal()
  if (kfl_adj_prob == 1) call nsi_senscal()
  
end subroutine nsi_endste
