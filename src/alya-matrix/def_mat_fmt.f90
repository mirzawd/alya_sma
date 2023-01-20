!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_csr.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   Matrix
!> @details Matrix formats
!-----------------------------------------------------------------------

module def_mat_fmt

  use def_mat
  use def_mat_csr, only : mat_csr
  use def_mat_coo, only : mat_coo
  use def_mat_ell, only : mat_ell
  use def_mat_den, only : mat_den
  use def_mat_blk, only : mat_blk
  use def_mat_dia, only : mat_dia
  use def_mat_tri, only : mat_tri
  use def_mat_sky, only : mat_sky
  use def_mat_bnd, only : mat_bnd

end module def_mat_fmt
!> @}
