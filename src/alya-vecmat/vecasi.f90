!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine vecasi(n,v1,v2)

!-----------------------------------------------------------------------
!
! Vector assign:    v2(i) = v1(i)   i=1..n
!
!-----------------------------------------------------------------------
  use      def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: n
  real(rp),    intent(in)  :: v1(n)
  real(rp),    intent(out) :: v2(n)
  integer(ip)              :: i

  do i=1,n
     v2(i)=v1(i)
  end do
  
end subroutine vecasi


