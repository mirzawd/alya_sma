!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    nsi_endste.f90
!> @author  Guillaume Houzeaux
!> @date    05/06/2013
!> @brief   End a time step
!> @details End a time step
!> @} 
!-----------------------------------------------------------------------
subroutine ker_endste()
  use def_kintyp,        only : ip 
  use def_master,        only : ITASK_ENDSTE
  use mod_ker_detection, only : ker_detection_boundaries
  use mod_eccoupling,    only : kfl_exmsld_ecc, eccou_manage_variables
  use mod_biofibers,     only : kfl_biofibers, biofibers
  implicit none
  !
  ! Feature detection
  !
  !call ker_detection_boundaries()
  !
  ! Velocity, temperature, concentration and displacement functions
  ! (:,3) <= (:,1)
  !
  call ker_velfun(ITASK_ENDSTE)
  call ker_temfun(ITASK_ENDSTE)
  call ker_confun(ITASK_ENDSTE)
  call ker_disfun(ITASK_ENDSTE)
  call ker_arefun(ITASK_ENDSTE)
  !
  ! Electro-mechanical coupling
  ! 
  if( kfl_exmsld_ecc )then
     call eccou_manage_variables(ITASK_ENDSTE)
  endif
  !
  ! Bio-fibers
  ! 
  if( kfl_biofibers )then
     call biofibers % update_fibers_at_nodes(ITASK_ENDSTE)
  endif

end subroutine ker_endste
