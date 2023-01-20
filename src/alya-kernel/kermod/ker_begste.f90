!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    nsi_begste.f90
!> @author  Guillaume Houzeaux
!> @date    05/06/2013
!> @brief   Begin a time step
!> @details Begin a time step
!> @} 
!-----------------------------------------------------------------------
subroutine ker_begste()
  
  use def_master,        only : ITASK_BEGSTE
  use mod_mass_matrix,   only : mass_matrix_consistent
  use mod_ker_subdomain, only : ker_subdomain_motion
  use mod_eccoupling,    only : eccou_save_log
  use mod_immersed,      only : cou_update_couplings  
  use mod_biofibers,     only : kfl_biofibers, biofibers
  use def_coupli,        only : kfl_dimbou

  implicit none
  !
  ! Transient fields
  !
  call calc_kx_tran_fiel()
  !
  ! Update fields: Velocity, temperature, concentration and displacement functions
  !
  call ker_velfun(ITASK_BEGSTE)
  call ker_temfun(ITASK_BEGSTE)
  call ker_confun(ITASK_BEGSTE)
  call ker_disfun(ITASK_BEGSTE)
  call ker_arefun(ITASK_BEGSTE)
  !
  ! Update subdomains
  !
  call ker_subdomain_motion()  
  !
  ! Update couplings 
  !
  if( kfl_dimbou ) call cou_update_couplings()
  !
  ! Weighted consistent mass... Kermod should have been read
  !
  call mass_matrix_consistent(CONSISTENT_WEIGHTED_MASS=.true.)
  !
  ! Bio-fibers
  !
  if( kfl_biofibers )then
     call biofibers % update_fibers_at_nodes(ITASK_BEGSTE)
  endif

end subroutine ker_begste
