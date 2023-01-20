!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine instri(pnodb,point,facoo,ifoun)
  !-----------------------------------------------------------------------
  ! NAME
  !    insiib
  ! DESCRIPTION
  !    Determine if a point is inside a triangle using the same side technique
  !    point: point
  !    facoo(dimension,vertices): triangle coordinates
  !    ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  ! USED BY
  !    pofadi,faliib
  !***
  !----------------------------------------------------------------------- 

  use def_kintyp, only           :  ip,rp
  use def_domain, only           :  ndime

  implicit   none
  integer(ip),   intent(out)    :: ifoun
  integer(ip),   intent(in)     :: pnodb
  integer(ip)                   :: test1,test2,test3
  real(rp),      intent(in)     :: point(ndime)
  real(rp),      intent(in)     :: facoo(ndime,*)
  !
  ! The same side technique is applied at each side of the triangle
  !
  call samsid(facoo(1,1),facoo(1,2),point,facoo(1,3),ndime,test1)
  call samsid(facoo(1,1),facoo(1,3),point,facoo(1,2),ndime,test2)
  call samsid(facoo(1,2),facoo(1,3),point,facoo(1,1),ndime,test3)

  if ( test1 == 1 .and. test2 == 1 .and. test3 == 1) then
     ifoun = 1
  else
     ifoun = 0
  end if

  if( pnodb == 4 .and. ifoun == 0 ) then
     !
     ! Check other triangle
     !
     call samsid(facoo(1,1),facoo(1,3),point,facoo(1,4),ndime,test1)
     call samsid(facoo(1,1),facoo(1,4),point,facoo(1,3),ndime,test2)
     call samsid(facoo(1,3),facoo(1,4),point,facoo(1,1),ndime,test3)

     if ( test1 == 1 .and. test2 == 1 .and. test3 == 1 ) then
        ifoun = 1
     else
        ifoun = 0
     end if
  end if


end subroutine instri

subroutine instr2(pnodb,point,facoo,ifoun,bari1,bari2,ntria)
  !-----------------------------------------------------------------------
  ! NAME
  !    insiib
  ! DESCRIPTION
  !    Determine if a point is inside a triangle using the same side technique
  !    point: point
  !    facoo(dimension,vertices): triangle coordinates
  !    ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  ! USED BY
  !    pofadi,faliib
  !***
  !----------------------------------------------------------------------- 

  use def_kintyp, only           :  ip,rp
  use def_domain, only           :  ndime
  use def_master, only           :  zeror

  implicit   none
  integer(ip),   intent(in)     :: pnodb
  integer(ip),   intent(out)    :: ifoun
  integer(ip),   intent(out)    :: ntria
  real(rp),      intent(in)     :: point(ndime)
  real(rp),      intent(in)     :: facoo(ndime,*)
  real(rp),      intent(out)    :: bari1,bari2
  real(rp)                      :: v0(3),v1(3),v2(3)
  real(rp)                      :: dot00,dot01,dot02,dot11,dot12
  real(rp)                      :: invDenom
  !
  ! 3D
  !
  v0(1) = facoo(1,3) - facoo(1,1)
  v0(2) = facoo(2,3) - facoo(2,1)
  v0(3) = facoo(3,3) - facoo(3,1)

  v1(1) = facoo(1,2) - facoo(1,1)
  v1(2) = facoo(2,2) - facoo(2,1)
  v1(3) = facoo(3,2) - facoo(3,1)

  v2(1) = point(1)   - facoo(1,1)
  v2(2) = point(2)   - facoo(2,1)
  v2(3) = point(3)   - facoo(3,1)

  dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
  dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
  dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
  dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
  dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
  !
  ! Compute barycentric coordinates
  !
  ntria = 0
  if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
     invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
     bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
     bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
     !
     ! Check
     !
     if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
        ntria = 1
        ifoun = 1
     else
        ifoun = 0
     end if
  else
     ifoun = 0
  end if
  if( pnodb == 4 .and. ifoun == 0 ) then
  v0(1) = facoo(1,4) - facoo(1,1)
  v0(2) = facoo(2,4) - facoo(2,1)
  v0(3) = facoo(3,4) - facoo(3,1)

  v1(1) = facoo(1,3) - facoo(1,1)
  v1(2) = facoo(2,3) - facoo(2,1)
  v1(3) = facoo(3,3) - facoo(3,1)

  v2(1) = point(1)   - facoo(1,1)
  v2(2) = point(2)   - facoo(2,1)
  v2(3) = point(3)   - facoo(3,1)

   
     dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
     dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
     dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
     dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
     dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  

     if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
        invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
        bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
        bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
        !
        ! Check
        !
        if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
           ifoun = 1
           ntria = 2
        else
           ifoun = 0
        end if
     else
        ifoun = 0
     end if
  end if

end subroutine instr2
