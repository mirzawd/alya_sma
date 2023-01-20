!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_rotunk(irose,istpo,ifipo,vetar)
!-----------------------------------------------------------------------
!****f* Nastal/nsm_rotunk
! NAME 
!    nsa_rotunk
! DESCRIPTION
!    This routine rotates ahead or rotates back a nodal vector 
!    using the appropiate rotation matrix. 
!    The sense of rotation is given by input integer irose:
!            irose = -1 --> rotate ahead, from GCB to LCB
!            irose =  1 --> rotate back, from LCB to GCB 
!    
!    Recall that:
!    GCB : Global Cartesian Basis
!    LCB : Local Curvilinear Basis
!    
! USES
!    mbvab0
! USED BY
!    
!***
!-----------------------------------------------------------------------
  use      def_domain
  use      def_master
  implicit none
  integer(ip) :: ipoin,istpo,ifipo,ibopo
  integer(ip) :: irose
  real(rp)    :: venew(ndime),veold(ndime),vetar(ndime)

  if (istpo < ifipo) return 
!
! 
!
  do ipoin=istpo,ifipo

     ibopo=lpoty(ipoin)
     
     if(ibopo>0) then

        veold(1:ndime)=vetar(1:ndime)
        !
        ! Target Vector boundary conditions linked to exnor
        !           
        call sld_rotvec(irose,veold,exnor(1,1,ibopo),venew,ndime)

        vetar(1:ndime)=venew(1:ndime)
     end if

  end do

end subroutine sld_rotunk

