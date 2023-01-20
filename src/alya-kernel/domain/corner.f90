!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine corner(exwor,ipoin,ibott)

!------------------------------------------------------------------------
!
! This check if a corner point is of step of bottom type
!   
!          jpoin         jpoin
!            O             O<--------O ipoin
!           /|\                     /|
!  (domain)  |                     / |  (domain)
!            |                   \/_ |     
!            |                      \|/
!  O<--------O ipoin                 O 
!             \
!             _\|
!
!     Bottom type            Step type
!     angle>90               0<angle<90
!
!------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  r_bou,c_bou,coord,ndime
  implicit none
  integer(ip), intent(in)  :: ipoin
  real(rp),    intent(in)  :: exwor(ndime)
  integer(ip), intent(out) :: ibott
  integer(ip)              :: ista0,ista1,icolu,jpoin
  real(rp)                 :: vnorm,vecau(3)
      
  ibott=1                                                 ! By default, 
     
  ista0=r_bou(ipoin)
  ista1=r_bou(ipoin+1)-1
  icolu=ista0
  do while(icolu/=ista1)
     icolu=icolu+1
     jpoin=c_bou(icolu)
     if(ipoin/=jpoin) then
        vecau(    1) = coord(    1,jpoin) - coord(    1,ipoin)
        vecau(    2) = coord(    2,jpoin) - coord(    2,ipoin)
        vecau(ndime) = coord(ndime,jpoin) - coord(ndime,ipoin)
        call vecuni(ndime,vecau,vnorm)            
        vnorm=dot_product(vecau(1:ndime),exwor(1:ndime))
        if(vnorm>0.1_rp) then                             ! ipoin is of step type            
           icolu=ista1
           ibott=0
           return
        end if
     end if
  end do
  
end subroutine corner
