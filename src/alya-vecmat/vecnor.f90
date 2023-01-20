!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine vecnor(v,n,vnor,inor)
!$acc routine seq
!-----------------------------------------------------------------------
!
! Compute the L1, L2 and L-inf norms of a vector V of length N
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: n,inor
  real(rp),    intent(in)  :: v(n)
  real(rp),    intent(out) :: vnor
  integer(ip)              :: i

  vnor=0.0_rp
  if(inor==0.or.inor==3) then
     do i=1,n
        if(abs(v(i))>vnor) vnor=abs(v(i))
     end do
  else if(inor==1) then
     do i=1,n
        vnor=vnor+abs(v(i))
     end do
  else if(inor==2) then
     if(n==2) then
        vnor=sqrt(v(1)*v(1)+v(2)*v(2))
     else if(n==3) then
        vnor=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
     else
        do i=1,n
           vnor=vnor+v(i)*v(i)
        end do
        vnor=sqrt(vnor)
     end if
  end if
  
end subroutine vecnor
