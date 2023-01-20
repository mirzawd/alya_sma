!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine veloc_grad_tensor_lambda_max(gvelo,lambda_max)
  !--------------------------------------------------------------------------------
  !> @addtogroup Mathematics
  !> @{
  !> @file    veloc_grad_tensor_lambda_max.f90
  !> @author  Laura Nicolaou
  !> @brief   computes all the eigenvalues and eigenvectors of matrix A (3x3)
  !> @details returns max. modules of eigenvalues of veloc grad tensor to compute 
  !> 	      local fluid time scale for instantanteous particle Stokes number 
  !> @} 
  !--------------------------------------------------------------------------------
  use def_kintyp
  use def_domain, only    : ndime
  use mod_maths,  only    : maths_eigen_symmetric_matrix
  implicit none
  real(rp),    intent(in) :: gvelo(ndime,ndime)
  real(rp),    intent(out):: lambda_max

  real(rp)                :: mag(ndime)
  real(rp)                :: lambda(ndime)



  ! Eigenvalues (lambda) and Eigenvectors(V) of veloc grad tensor
  call  maths_eigen_symmetric_matrix(ndime,gvelo,lambda)
  
  ! Find largest modulus of the eigenvalues
  mag = abs(lambda)

  lambda_max = maxval(mag)

end subroutine veloc_grad_tensor_lambda_max

