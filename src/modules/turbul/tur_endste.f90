!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_endste()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_endste
  ! NAME 
  !    tur_endste
  ! DESCRIPTION
  !    This routine ends a time step of the turbulence equations.
  ! USES
  !    tur_cvgunk
  !    tur_updunk
  !    tur_output
  ! USED BY
  !    Turbul
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  use def_kermod, only : kfl_adj_prob,kfl_cos_opt

  implicit none  
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_tur==0.and.kfl_timei_tur==1) then
     call tur_cvgunk(ITASK_ENDSTE)
     call tur_updunk(ITASK_ENDSTE)
  end if
  !
  ! Compute averaged variables
  !
  call tur_averag()
  !
  ! Write restart file
  !
!  call tur_restar(two)
  !
  ! If not steady, go on
  !
  if(kfl_stead_tur==0.and.kfl_timei_tur==1.and.kfl_conve(modul)==1) kfl_gotim = 1
  !
  ! Calculate sensitivities
  !
  if (kfl_adj_prob == 1) call tur_senscal()
  if (kfl_cos_opt == 1) call tur_costcal()
  
end subroutine tur_endste
