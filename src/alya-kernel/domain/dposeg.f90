!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine dposeg(xcoor,coor1,coor2,dista,proje,ifoun)
  !-----------------------------------------------------------------------
  ! NAME
  !    segdis
  ! DESCRIPTION
  !    Minimun distance between a point and a segment
  !    xcoor : point xcoorinates 
  !    coor1,coor2 : defines the segment
  !    ndime: dimension
  !    dista: distance
  !    proje: projection of the point on the segment
  !    ifoun = 1 the projection point is inside the segment
  !    ifoun = 0 the projection point is outside the segment
  ! USED BY
  !    pofadi
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only           :  ip,rp 
  use def_domain, only           :  ndime
  implicit none
  integer(ip),   intent(out)     :: ifoun
  real(rp),      intent(in)      :: xcoor(ndime),coor1(ndime),coor2(ndime)
  real(rp),      intent(out)     :: dista,proje(ndime)
  integer(ip)                    :: idime  
  real(rp)                       :: numer,denom,dsegm

  numer = 0.0_rp
  denom = 0.0_rp
  do idime = 1,ndime
      numer = numer + (coor2(idime) - coor1(idime)) * (xcoor(idime) - coor1(idime))
      denom = denom + (coor2(idime) - coor1(idime)) * (coor2(idime) - coor1(idime))
  end do  

  dsegm = numer / denom
  
  if( dsegm < -0.01_rp ) then     
     ifoun        = 0_ip
     proje(1)     = coor1(1)
     proje(2)     = coor1(2)
     proje(ndime) = coor1(ndime)

  else if ( dsegm >= -0.01_rp .and. dsegm <= 1.01_rp ) then

     ifoun        = 1_ip
     proje(1)     = coor1(1)     + dsegm * (coor2(1)     - coor1(1))
     proje(2)     = coor1(2)     + dsegm * (coor2(2)     - coor1(2))
     proje(ndime) = coor1(ndime) + dsegm * (coor2(ndime) - coor1(ndime))

  else

     ifoun        = 0_ip
     proje(1)     = coor2(1)
     proje(2)     = coor2(2)
     proje(ndime) = coor2(ndime)

  end if

  dista = 0.0_rp
  do idime = 1,ndime
     dista = dista + (xcoor(idime) - proje(idime))*(xcoor(idime) - proje(idime))
  end do

end subroutine dposeg
