!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure subroutine vecres(norm,n,v1,v2,redif,zero)

!-----------------------------------------------------------------------
!
! Compute the relative difference between two vectors:
! 
! redif = ||v1 - v2|| / ||v1||      
!
!-----------------------------------------------------------------------
  use      def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: n,norm
  real(rp),    intent(in)  :: v1(n),v2(n),zero
  real(rp),    intent(out) :: redif
  integer(ip)              :: i
  real(rp)                 :: numer,denom,va,vo

  redif = 0.0_rp
  numer = 0.0_rp
  denom = 0.0_rp

  select case(norm)

  case(0)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = max(numer,abs(va-vo))
        denom = max(denom,abs(va))
     end do
     if(denom.gt.zero) redif = numer/denom

  case(1)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + abs(va-vo)
        denom = denom + abs(va)
     end do
     if(denom>zero) redif = numer/denom

  case(2)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + (va-vo)*(va-vo)
        denom = denom + va*va
     end do
     if(denom>zero) then
        redif = sqrt(numer/denom)
     else
        redif = 0.0_rp
     end if

  end select

end subroutine vecres
