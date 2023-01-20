!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_funcre(param,npara,ifuge,timev,coord,coord_ini,displ)

  !------------------------------------------------------------------------
  !
  ! This function yields a parabolic or periodic evolution
  !
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_parame, only     :  pi
  implicit none
  integer(ip), intent(in)  :: npara
  integer(ip), intent(in)  :: ifuge
  real(rp),    intent(in)  :: param(npara)
  real(rp),    intent(in)  :: timev
  real(rp),    intent(in)  :: coord(ndime)
  real(rp),    intent(in)  :: coord_ini(ndime)
  real(rp),    intent(out) :: displ(ndime)
  integer(ip)              :: idime
  real(rp)                 :: deltx,deltz
  real(rp)                 :: coor0(3),betaa,radiu,alpha,w,d,A


  if( ifuge == 0 ) then
     return

  else if( ifuge == 1 ) then

  else if( ifuge == 2 ) then

  else if( ifuge == 3 ) then

  else if( ifuge == 4 ) then
 

  else if( ifuge == 5 ) then


  else if( ifuge == 6 ) then
     !
     ! Translation
     !
     do idime = 1,ndime
        displ(idime) = param(idime) 
     end do

  else if( ifuge == 7 ) then
     !
     ! Rotation
     !
     coor0(1) = param(1)
     coor0(2) = param(2)
     coor0(3) = param(3)
     betaa    = param(4)*pi/180.0_rp
     deltx    = coord(1)-coor0(1)
     deltz    = coord(3)-coor0(3)
     radiu    = sqrt(  (deltx*deltx) + (deltz*deltz) )     
     alpha    = atan2( deltz , deltx )
     
     displ(1) =  ( radiu * cos(alpha+betaa) ) + coor0(1)
     displ(2) =  coord(2)
     displ(3) =  ( radiu * sin(alpha+betaa) ) + coor0(3) 
     do idime = 1,ndime
        displ(idime) = displ(idime) - coord(idime)
     end do

  else if( ifuge == 8 ) then
     !
     ! Sin
     !
     w          = 1.0_rp
     A          = param(1)
     d          = 0.5_rp*A*(1.0_rp-cos(2.0_rp*pi*w*timev))
     displ(1)   = coord_ini(1)-coord(1)+d
     displ(2)   = 0.0_rp

  else if( ifuge == 9 ) then
     !
     ! Linear
     !
     do idime = 1,ndime
        displ(idime) = param(idime) 
     end do

  end if

end subroutine ale_funcre

