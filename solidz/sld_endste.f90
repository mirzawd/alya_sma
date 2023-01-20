!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_endste.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine ends a time step of solidz
!> @details
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_endste()

  use def_kintyp,     only : ip, rp
  use def_master,     only : ITASK_ENDSTE
  use def_master,     only : modul, kfl_conve, kfl_gotim
  use def_solidz,     only : kfl_stead_sld, kfl_rigid_sld
  use def_solidz,     only : kfl_isoch_sld
  use mod_sld_energy, only : sld_updene
  use mod_sld_rbo,    only : sld_rbo_cvgunk, sld_rbo_updunk

  implicit none

  if( kfl_stead_sld == 0_ip ) then

     if ( kfl_rigid_sld == 0_ip ) then

        !-------------------------------------------------------------------
        !
        ! Deformable body
        !
        !-------------------------------------------------------------------
        !
        ! Compute convergence residual of the time evolution (that is,
        ! || d(n,*,*) - d(n-1,*,*)|| / ||d(n,*,*)||) and update unknowns
        !
        call sld_cvgunk(ITASK_ENDSTE) ! Convergence
        call sld_updunk(ITASK_ENDSTE) ! Update unknowns and state variables d(n-1,*,*) <-- d(n,*,*)
        call sld_updene(ITASK_ENDSTE) ! Update energies
        !
        ! Isochrones
        !
        if( kfl_isoch_sld == 1_ip ) call sld_isochr

     else if ( kfl_rigid_sld == 1_ip ) then

        !-------------------------------------------------------------------
        !
        ! Rigid body
        !
        !-------------------------------------------------------------------
        !
        ! Convergence
        !
        call sld_rbo_cvgunk(ITASK_ENDSTE)
        !
        ! Update RB unknowns: (:,TIME_N) <= (:,ITER_K)
        !
        call sld_rbo_updunk(ITASK_ENDSTE) ! Update unknowns: (:,TIME_N)   <= (:,ITER_K)

     end if

  end if
  !
  ! If not steady, go on
  !
  if( kfl_stead_sld == 0_ip .and. kfl_conve(modul) == 1_ip ) kfl_gotim = 1_ip

end subroutine sld_endste
