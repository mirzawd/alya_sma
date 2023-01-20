!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmlen(&
     ndime,pnode,deriv,tragl,elcod,hnatu,hleng)
  !-----------------------------------------------------------------------
  !****f* Domain/elmlen
  ! NAME
  !    elmlen
  ! DESCRIPTION
  !    Compute element lengths
  ! OUTPUT
  !    HLENG ... Element length with: 
  !              HLENG(1)     = Max length
  !              HLENG(NDIME) = Min length
  ! USED BY
  !    
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     : ip,rp
  implicit none 
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(out) :: tragl(ndime,ndime)
  real(rp),    intent(in)  :: hnatu,elcod(ndime,pnode)
  real(rp),    intent(in)  :: deriv(ndime,pnode)
  real(rp),    intent(out) :: hleng(ndime)
  integer(ip)              :: k
  real(rp)                 :: enor0,h_tem,gpdet,denom
  real(rp)                 :: xjacm(ndime,ndime),t1,t2,t3

  if( ndime == 1 ) then

     xjacm(1,1) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
     end do
     tragl(1,1) = 1.0_rp/xjacm(1,1)
     enor0      = tragl(1,1) * tragl(1,1)
     hleng(1)   = hnatu/sqrt(enor0)

  else if( ndime == 2 ) then

     xjacm(1,1) = 0.0_rp
     xjacm(1,2) = 0.0_rp
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
     end do

     gpdet      =  xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)          
     denom      =  1.0_rp/gpdet
     tragl(1,1) =  xjacm(2,2) * denom
     tragl(2,2) =  xjacm(1,1) * denom
     tragl(2,1) = -xjacm(2,1) * denom
     tragl(1,2) = -xjacm(1,2) * denom  

     enor0      =  tragl(1,1) * tragl(1,1) + tragl(1,2) * tragl(1,2) 
     hleng(1)   =  hnatu/sqrt(enor0)
     enor0      =  tragl(2,1) * tragl(2,1) + tragl(2,2) * tragl(2,2) 
     hleng(2)   =  hnatu/sqrt(enor0)

     if( hleng(2) > hleng(1) ) then
        h_tem    = hleng(2)
        hleng(2) = hleng(1)
        hleng(1) = h_tem
     end if

  else if( ndime == 3 ) then

     xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
     xjacm(1,2) = 0.0_rp ! tragl = xjacm^-1
     xjacm(1,3) = 0.0_rp 
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     xjacm(2,3) = 0.0_rp
     xjacm(3,1) = 0.0_rp
     xjacm(3,2) = 0.0_rp
     xjacm(3,3) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(1,3) = xjacm(1,3) + elcod(1,k) * deriv(3,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
        xjacm(2,3) = xjacm(2,3) + elcod(2,k) * deriv(3,k)
        xjacm(3,1) = xjacm(3,1) + elcod(3,k) * deriv(1,k)
        xjacm(3,2) = xjacm(3,2) + elcod(3,k) * deriv(2,k)
        xjacm(3,3) = xjacm(3,3) + elcod(3,k) * deriv(3,k)
     end do

     t1         =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2         = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3         =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet      =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
     denom      =  1.0_rp / gpdet
     tragl(1,1) =  t1 * denom
     tragl(2,1) =  t2 * denom
     tragl(3,1) =  t3 * denom
     tragl(2,2) =  ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
     tragl(3,2) =  (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
     tragl(3,3) =  ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
     tragl(1,2) =  (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
     tragl(1,3) =  ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
     tragl(2,3) =  (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom
     !
     ! Element length HLENG
     !
     enor0    =   tragl(1,1) * tragl(1,1) &
          &     + tragl(1,2) * tragl(1,2) &
          &     + tragl(1,3) * tragl(1,3)
     hleng(1) = hnatu/sqrt(enor0)
     enor0    =   tragl(2,1) * tragl(2,1) &
          &     + tragl(2,2) * tragl(2,2) &
          &     + tragl(2,3) * tragl(2,3)
     hleng(2) = hnatu/sqrt(enor0)
     enor0    =   tragl(3,1) * tragl(3,1) &
          &     + tragl(3,2) * tragl(3,2) &
          &     + tragl(3,3) * tragl(3,3) 
     hleng(3) = hnatu/sqrt(enor0) 
     !
     ! Sort hleng: hleng(1)=max; hleng(ndime)=min
     !     
     if( hleng(2) > hleng(1) ) then
        h_tem    = hleng(2)
        hleng(2) = hleng(1)
        hleng(1) = h_tem
     end if
     if( hleng(3) > hleng(1) ) then
        h_tem    = hleng(3)
        hleng(3) = hleng(1)
        hleng(1) = h_tem
     end if
     if( hleng(3) > hleng(2) ) then
        h_tem    = hleng(3)
        hleng(3) = hleng(2)
        hleng(2) = h_tem
     end if

  end if

end subroutine elmlen
