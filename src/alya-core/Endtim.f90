!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Timste
!> @{
!> @file    Timste.f90
!> @author  Guillaume Houzeaux
!> @author  Adria Quintanas-Corominas
!> @brief
!> @details
!> @}
!-----------------------------------------------------------------------
subroutine Endtim()

  use def_master,       only : ITASK_ENDTIM
  use def_master,       only : iblok, nblok
  use mod_moduls,       only : moduls
  implicit none

  ! Compute the time step for each module
  call Kermod(ITASK_ENDTIM)
  do iblok = 1,nblok
     call moduls(ITASK_ENDTIM)
  end do

end subroutine Endtim