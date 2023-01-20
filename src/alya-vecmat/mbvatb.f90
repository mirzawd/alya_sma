!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine mbvatb(a,b,c,n1,n2)

!------------------------------------------------------------------------
!
! This routine evaluates the matrix-vector product A = Bt C, where
! A -> R(n1), B -> Mat(n2,n1), C -> R(n2)
!
!------------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: n1,n2
  real(rp),    intent(in)  :: b(n2,n1),c(n2)
  real(rp),    intent(out) :: a(n1)
  integer(ip)              :: i,j

  do i=1,n1
     a(i)=0.0_rp
     do j=1,n2
        a(i)=a(i)+b(j,i)*c(j)
     end do
  end do
  
end subroutine mbvatb
