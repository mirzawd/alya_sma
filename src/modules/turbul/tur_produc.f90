!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Turbul
!> @{
!> @file    tur_produc.f90
!> @author  Guillaume Houzeaux
!> @date    12/06/2013
!> @brief   Production term
!> @details Project the production term using a L2 projection
!> @} 
!-----------------------------------------------------------------------
subroutine tur_produc()
  use def_master
  use def_turbul
  use mod_projec
  implicit none

  if( kfl_produ_tur == 1 ) then
     call projec_norm_symmetric_gradient_vector(advec,produ_tur)
  end if

end subroutine tur_produc
