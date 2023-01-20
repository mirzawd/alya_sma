!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    spcdec.f90
!> @author  Herbert Owen
!> @brief   Computes eigenval. and eigenvec. of real, symmetrix matrix A (3x3). Sorts them in accending order.
!> @details Output is a vector D
!!          containing the eigenvalues in ascending order, and a matrix V whose
!!          columns contain the corresponding eigenvectors.
!> @}
!------------------------------------------------------------------------
subroutine spcdec(A,D,V,NROT,kfl_wivec,callersub)
  
  use def_kintyp, only : ip,rp
  use mod_maths,  only : maths_eigen_3x3_symmetric_matrix

  implicit none
  real(rp),    intent(in)  :: A(3,3)
  real(rp),    intent(out) :: D(3), V(3,3)
  integer(ip), intent(out) :: NROT
  integer(ip), intent(in)  :: kfl_wivec   ! also obtain eigenvectors
  character(*)             :: callersub   ! who is calling spcdec

  call maths_eigen_3x3_symmetric_matrix(A,D,V,maxit=100_ip,toler=1.0e-06_rp)  

end subroutine spcdec

