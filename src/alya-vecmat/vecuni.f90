!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine vecuni(n,v,rmod)

!-----------------------------------------------------------------------
!
! This routine computes the length of vector V and converts it to  
! a unit one 
!
!-----------------------------------------------------------------------
  use      def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)    :: n
  real(rp),    intent(inout) :: v(n)
  real(rp),    intent(out)   :: rmod
  integer(ip)                :: i

  rmod=0.0_rp
  do i=1,n
     rmod=rmod + v(i)*v(i)
  end do
  rmod=sqrt(rmod)
  if(rmod>epsilon(1.0_rp)) v=v/rmod
 
end subroutine vecuni
      
