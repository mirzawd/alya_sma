!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine invmtx(a,b,deter,nsize)

!-----------------------------------------------------------------------
!
! This routine inverts a square matrix A -> Mat(nsize,nsize). The
! inverse is stored in B. Its determinant is DETER
!
!
!-----------------------------------------------------------------------
  use def_kintyp, only: ip,rp
  implicit none
  integer(ip), intent(in)  :: nsize
  real(rp),    intent(in)  :: a(nsize,nsize)
  real(rp),    intent(out) :: b(nsize,nsize),deter
  real(rp)                 :: denom,t1,t2,t3,t4

  select case( nsize )

  case ( 1_ip )

     deter = a(1,1)
     if( abs(deter) == 0.0_rp ) return
     b(1,1) = 1.0_rp/a(1,1)

  case( 2_ip )

     deter = a(1,1)*a(2,2)-a(2,1)*a(1,2)
     if( abs(deter) == 0.0_rp ) return
     denom  = 1.0_rp/deter
     b(1,1) = a(2,2)*denom
     b(2,2) = a(1,1)*denom
     b(2,1) =-a(2,1)*denom
     b(1,2) =-a(1,2)*denom

  case(3)
     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
     deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
     if( abs(deter) == 0.0_rp ) return
     denom  = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
     b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
     b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
     b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
     b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
     b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

  case(4)
     t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
     t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
     t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
     t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
     deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
     if( abs(deter) == 0.0_rp ) return
     denom=1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(4,1) = t4*denom
     b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
          - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
          + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
     b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
          + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
          - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
     b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
          - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
          + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
     b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
          + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
          - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
     b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
          + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
          - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
     b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
          - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
          + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
     b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
          + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
          - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
     b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
          - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
          + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
     b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
          - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
          + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
     b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
          + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
          - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
     b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
          - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
          + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
     b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
          + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
          - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom


  case default
     b=a
     call invert(b,nsize,nsize)

  end select

end subroutine invmtx
