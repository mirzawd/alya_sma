!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    eigsrt.f90
!> @author  Herbert Owen
!> @brief   Sorts the eigenvalues, D, into ascending order.
!> @details The eigenvectors, V, are rearranged accordingly.
!> @} 
!------------------------------------------------------------------------
subroutine eigsrt(D,V)

  use def_kintyp
  implicit none
  real(rp),    intent(inout) :: D(3), V(3,3)
  real(rp)                   :: P
  integer(ip)                :: idime, jdime, K
  
  do idime=1,2
     K = idime
     P = D(idime)
     do jdime = idime + 1,3
        if (D(jdime) >= P) then
           K = jdime
           P = D(jdime)
        end if
     end do

     if (K /= idime) then
        D(K) = D(idime)
        D(idime) = P
        do jdime=1,3
           P = V(jdime,idime)
           V(jdime,idime) = V(jdime,K)
           V(jdime,K) = P
        end do
     end if
  end do

end subroutine eigsrt


