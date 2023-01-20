!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine veloc_grad_tensor(gvelo,gpvol,Q_crit,lambda2)
  !------------------------------------------------------------------------
  !> @addtogroup Mathematics
  !> @{
  !> @file    veloc_grad_tensor.f90
  !> @author  hadrien calmet
  !> @brief   computes all the eigenvalues and eigenvectors of matrix A (3x3)
  !> @details blablabla
  !> @} 
  !------------------------------------------------------------------------
  use def_kintyp
  use def_domain,       only : ndime
  use mod_maths_solver, only : maths_eigen_symmetric_matrix
  
  implicit none
  real(rp),    intent(in) :: gvelo(ndime,ndime),gpvol
  real(rp),    intent(out):: Q_crit,lambda2
  integer(ip)             :: i,j,k

  real(rp)                :: S_ij(ndime,ndime), O_ij(ndime,ndime) 
  real(rp)                :: squareS_ij(ndime,ndime), squareO_ij(ndime,ndime)
  real(rp)                :: G__ij(ndime,ndime)
!  real(rp)                :: eigen(ndime)
  real(rp)                :: norm_S_ij, norm_O_ij
!  real(rp)                :: Ginv(ndime),angle(ndime)
  real(rp)                :: lambda(ndime),V(ndime,ndime)
!  real(rp)                :: pi
  !
  !Compute of rate-of-strain (Sij) and voticity tensor (\omega ij)
  !
  do i=1,ndime
     do j=1,ndime
        S_ij(i,j) = 0.0_rp 
        O_ij(i,j) = 0.0_rp
     enddo
  enddo

  do i=1,ndime
     do j=1,ndime
        S_ij(i,j) = S_ij(i,j) + 0.5_rp * gpvol* (gvelo(i,j) + gvelo(j,i))
        O_ij(i,j) = O_ij(i,j) + 0.5_rp * gpvol* (gvelo(i,j) - gvelo(j,i))
     enddo
  enddo
  !
  !Compute the norm L2 of rate-of-strain (Sij) and voticity tensor (\omega ij)
  !
  norm_S_ij=0.0_rp
  norm_O_ij=0.0_rp
  do i=1,ndime
     do j=1,ndime
        norm_S_ij = norm_S_ij + S_ij(i,j)*S_ij(i,j)
        norm_O_ij = norm_O_ij + O_ij(i,j)*O_ij(i,j)
     enddo
  enddo
  norm_S_ij = sqrt(norm_S_ij)
  norm_O_ij = sqrt(norm_O_ij)
  !
  ! Q criterion
  !
  Q_crit = 0.5_rp * ((norm_O_ij*norm_O_ij)-(norm_S_ij*norm_S_ij))
  !
  ! Lambda 2 criterion 
  !
  ! need S_ij^2 + O_ij^2 matrix
  do i=1,ndime
     do j=1,ndime
        squareS_ij(i,j)= 0.0_rp
        do k=1,ndime
           squareS_ij(i,j)= squareS_ij(i,j) + S_ij(i,k)*S_ij(k,j)
        enddo
     enddo
  enddo
  do i=1,ndime
     do j=1,ndime
        squareO_ij(i,j)= 0.0_rp
        do k=1,ndime
           squareO_ij(i,j)= squareO_ij(i,j) + O_ij(i,k)*O_ij(k,j)
        enddo
     enddo
  enddo
  do i=1,ndime
     do j=1,ndime
        G__ij(i,j)=0.0_rp
     enddo
  enddo
  do i=1,ndime
     do j=1,ndime
        G__ij(i,j)=squareS_ij(i,j) + squareO_ij(i,j)
     enddo
  enddo
  !
  !test yorgos 
  !

  ! Calculating the invariants of G
  ! I1 = tr(G)
  ! I2 = 1/2 * (tr(G)^2 - tr(G^2))
  ! I3 = det(G)

  !Ginv(1) = G__ij(1,1) + G__ij(2,2) + G__ij(3,3)
  !Ginv(2) =  G__ij(1,1)*G__ij(2,2) + G__ij(2,2)*G__ij(3,3) + G__ij(3,3)*G__ij(1,1) &
  !     - G__ij(1,2)*G__ij(1,2) - G__ij(2,3)*G__ij(2,3) - G__ij(1,3)*G__ij(1,3)
  !Ginv(3) =   G__ij(1,1) * (G__ij(2,2)*G__ij(3,3) - G__ij(2,3)*G__ij(2,3)) &
  !     - G__ij(1,2) * (G__ij(1,2)*G__ij(3,3) - G__ij(2,3)*G__ij(1,3)) &
  !     + G__ij(1,3) * (G__ij(1,2)*G__ij(2,3) - G__ij(2,2)*G__ij(1,3))

  ! Calculating angles from the invariants
  ! a1 = (I1^2)/9  - I2/3
  ! a2 = (I1^3)/27 - I1*I2/6 + I3/2
  ! a3 = 1/3 * arccos(a2/(a3^(3/2)))

  !angle(1) = Ginv(1)*Ginv(1)/9.0_rp - Ginv(2)/3.0_rp
  !angle(2) = Ginv(1)*Ginv(1)*Ginv(1)/27.0_rp - Ginv(1)*Ginv(2)/6.0_rp + Ginv(3)/2.0_rp
  !angle(3) = 0.3333333333_rp*(acos(angle(2)/(angle(1)**1.5_rp)))

  !pi = 4.0_rp*atan(1.0_rp)

  ! Calculating the eigenvalues and singular values
  ! lamda1= I1/3 + 2*sqrt(a1)*cos(a3)
  ! lamda2= I1/3 - 2*sqrt(a1)*cos(pi/3 + a3)
  ! lamda3= I1/3 - 2*sqrt(a1)*cos(pi/3 - a3)


  !eigen(1) = (Ginv(1)/3.0_rp + 2.0_rp*sqrt(angle(1))*cos(angle(3)))
  !eigen(2) = (Ginv(1)/3.0_rp - 2.0_rp*sqrt(angle(1))*cos(pi/3.0_rp + angle(3)))
  !eigen(3) = (Ginv(1)/3.0_rp - 2.0_rp*sqrt(angle(1))*cos(pi/3.0_rp - angle(3)))

  !eigen(3) = max(eigen(3), 0.0_rp)    ! The eigenvalues of a semi-definite positive matrix G are always non-negative
  ! Sorts the eigenvalues, lambda, into ascending order. 
  !call eigsrt(eigen,G__ij)
  !lambda2=eigen(2)

!----------------------------------------------
  ! Eigenvalues (lambda) and Eigenvectors(V) of the sum G__ij
  

  call maths_eigen_symmetric_matrix(ndime,G__ij,lambda,V,maxit=100_ip,toler=1.0e-06_rp)

  lambda2=lambda(2)



end subroutine veloc_grad_tensor

