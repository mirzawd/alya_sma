!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine mbmatb(a,b,c,n1,n2,n3)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix product A = Bt C, where
! A -> Mat(n1,n2), B -> Mat(n3,n1), C -> Mat(n3,n2)
!
!-----------------------------------------------------------------------
  use def_kintyp, only: ip,rp
  implicit none
  integer(ip), intent(in)  :: n1,n2,n3
  real(rp),    intent(out) :: a(n1,n2)
  integer(ip), intent(in)  :: b(n3,n1), c(n3,n2)    
  integer(ip)              :: i,j,k

  do i=1,n1
     do j=1,n2
        a(i,j)=0.0_rp
        do k=1,n3
           a(i,j)=a(i,j)+b(k,i)*c(k,j)
        end do
     end do
  end do

end subroutine mbmatb

