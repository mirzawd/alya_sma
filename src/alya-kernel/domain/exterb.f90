!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exterb(iimbo,pnodb,iboib,bouno)
  !-----------------------------------------------------------------------
  !****f* exterb/exterb
  ! NAME
  !    exterb
  ! DESCRIPTION
  !    This routines computes the exterior normal to the RB - very similar to exteib 
  ! USED BY
  !    domain
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  rbbou
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: iimbo,pnodb,iboib
  real(rp),    intent(out) :: bouno(ndime)
  integer(ip)              :: p1,p2,p3
  real(rp)                 :: xfact,vec(3,3)

  if( pnodb == 2 ) then
     p1       = rbbou(iimbo) % lnoib(1,iboib)
     p2       = rbbou(iimbo) % lnoib(2,iboib)
     vec(1,1) = rbbou(iimbo) % cooib(1,p2) - rbbou(iimbo) % cooib(1,p1)
     vec(2,1) = rbbou(iimbo) % cooib(2,p2) - rbbou(iimbo) % cooib(2,p1)
     bouno(1) =  vec(2,1)
     bouno(2) = -vec(1,1)
  else if( pnodb == 3 .or. pnodb==4 ) then
     p1       = rbbou(iimbo) % lnoib(1,iboib)
     p2       = rbbou(iimbo) % lnoib(2,iboib)
     p3       = rbbou(iimbo) % lnoib(3,iboib)
     call nortri(p1,p2,p3,rbbou(iimbo) % cooib,vec,ndime)
     bouno(1) = 0.5_rp*vec(1,3)
     bouno(2) = 0.5_rp*vec(2,3)
     bouno(3) = 0.5_rp*vec(3,3)
  else
     call runend('EXTERB: COULD NOT COMPUTE EXTERIOR NORMAL')
  end if
  call vecuni(ndime,bouno,xfact)


end subroutine exterb
