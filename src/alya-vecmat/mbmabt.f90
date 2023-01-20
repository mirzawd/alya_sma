!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine mbmabt(a,b,c,n1,n2,n3)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix product A = B Ct, where
! A -> Mat(n1,n2), B -> Mat(n1,n3), C -> Mat(n2,n3)
!
!-----------------------------------------------------------------------
  use def_kintyp, only: ip,rp
  implicit none
  integer(ip), intent(in)  :: n1,n2,n3
  real(rp),    intent(in)  :: b(n1,n3), c(n2,n3)
  real(rp),    intent(out) :: a(n1,n2)
  integer(ip)              :: i,j,k
    
  do i=1,n1
     do j=1,n2
        a(i,j)=0.0_rp
        do k=1,n3
           a(i,j)=a(i,j)+b(i,k)*c(j,k)
        end do
     end do
  end do

end subroutine mbmabt
