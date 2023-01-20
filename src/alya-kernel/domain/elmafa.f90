!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup ???
!> @{
!> @file    elmafa.f90
!> @author  gguillam
!> @date    2022-11-15
!> @brief   Computes the area for each face of hexahedron
!> @details Computes the area for each face of hexahedron
!>          
!> @} 
!-----------------------------------------------------------------------

subroutine elmafa(pelty,pnode,elcod,aface)

  use def_master
  use def_domain
  use def_elmtyp
  use def_elmgeo, only : element_type
  
  implicit none

  external                 :: runend
  
  integer(ip), intent(in)  :: pelty
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(out) :: aface(ndime)
  
  integer(ip)              :: iface,inodf,pnodf
  real(rp)                 :: xloc(4)
  real(rp)                 :: yloc(4)
  real(rp)                 :: zloc(4)
  real(rp)                 :: area_face, a
  real(rp)                 :: lface(6)
  integer(ip)              :: k
  real(rp)                 :: dshxi(pnode,ndime)
  real(rp)                 :: dXdXi,dXdEta,dXdZeta
  real(rp)                 :: dYdXi,dYdEta,dYdZeta
  real(rp)                 :: dZdXi,dZdEta,dZdZeta

  real(rp)                 :: enor0,h_tem,gpdet,denom
  real(rp)                 :: xjacm(ndime,ndime),t1,t2,t3
    
  if( pelty == HEX08 ) then

     !-------------------------------------------------------------------
     !
     ! HEX08
     !
     !-------------------------------------------------------------------

     pnodf = 4_ip

     do iface = 1,element_type(pelty) % number_faces

        area_face = 0.0_rp
        !
        ! Gauss point locations in each face
        !
        if( iface == 1_ip ) then
           xloc(1) = -sqrt(1.0_rp/3.0_rp)
           yloc(1) = -sqrt(1.0_rp/3.0_rp)
           zloc(1) = -1.0_rp
           xloc(2) =  sqrt(1.0_rp/3.0_rp)
           yloc(2) = -sqrt(1.0_rp/3.0_rp)
           zloc(2) = -1.0_rp
           xloc(3) =  sqrt(1.0_rp/3.0_rp)
           yloc(3) =  sqrt(1.0_rp/3.0_rp)
           zloc(3) = -1.0_rp
           xloc(4) = -sqrt(1.0_rp/3.0_rp)
           yloc(4) =  sqrt(1.0_rp/3.0_rp)
           zloc(4) = -1.0_rp
        else if( iface == 2_ip ) then
           xloc(1) = -sqrt(1.0_rp/3.0_rp)
           yloc(1) = -sqrt(1.0_rp/3.0_rp)
           zloc(1) =  1.0_rp
           xloc(2) =  sqrt(1.0_rp/3.0_rp)
           yloc(2) = -sqrt(1.0_rp/3.0_rp)
           zloc(2) =  1.0_rp
           xloc(3) =  sqrt(1.0_rp/3.0_rp)
           yloc(3) =  sqrt(1.0_rp/3.0_rp)
           zloc(3) =  1.0_rp
           xloc(4) = -sqrt(1.0_rp/3.0_rp)
           yloc(4) =  sqrt(1.0_rp/3.0_rp)
           zloc(4) =  1.0_rp
        else if( iface == 3_ip ) then
           xloc(1) = -sqrt(1.0_rp/3.0_rp)
           yloc(1) = -1.0_rp
           zloc(1) = -sqrt(1.0_rp/3.0_rp)
           xloc(2) =  sqrt(1.0_rp/3.0_rp)
           yloc(2) = -1.0_rp
           zloc(2) = -sqrt(1.0_rp/3.0_rp)
           xloc(3) =  sqrt(1.0_rp/3.0_rp)
           yloc(3) = -1.0_rp
           zloc(3) =  sqrt(1.0_rp/3.0_rp)
           xloc(4) = -sqrt(1.0_rp/3.0_rp)
           yloc(4) = -1.0_rp
           zloc(4) = sqrt(1.0_rp/3.0_rp)
        else if( iface == 4 ) then
           xloc(1) =  1.0_rp
           yloc(1) = -sqrt(1.0_rp/3.0_rp)
           zloc(1) = -sqrt(1.0_rp/3.0_rp)
           xloc(2) =  1.0_rp
           yloc(2) =  sqrt(1.0_rp/3.0_rp)
           zloc(2) = -sqrt(1.0_rp/3.0_rp)
           xloc(3) =  1.0_rp
           yloc(3) =  sqrt(1.0_rp/3.0_rp)
           zloc(3) =  sqrt(1.0_rp/3.0_rp)
           xloc(4) =  1.0_rp
           yloc(4) = -sqrt(1.0_rp/3.0_rp)
           zloc(4) =  sqrt(1.0_rp/3.0_rp)
        else if( iface == 5_ip ) then
           xloc(1) = -sqrt(1.0_rp/3.0_rp)
           yloc(1) =  1.0_rp
           zloc(1) = -sqrt(1.0_rp/3.0_rp)
           xloc(2) =  sqrt(1.0_rp/3.0_rp)
           yloc(2) =  1.0_rp
           zloc(2) = -sqrt(1.0_rp/3.0_rp)
           xloc(3) =  sqrt(1.0_rp/3.0_rp)
           yloc(3) = 1.0_rp
           zloc(3) =  sqrt(1.0_rp/3.0_rp)
           xloc(4) = -sqrt(1.0_rp/3.0_rp)
           yloc(4) =  1.0_rp
           zloc(4) =  sqrt(1.0_rp/3.0_rp)
        else if( iface == 6_ip ) then
           xloc(1) = -1.0_rp
           yloc(1) = -sqrt(1.0_rp/3.0_rp)
           zloc(1) = -sqrt(1.0_rp/3.0_rp)
           xloc(2) = -1.0_rp
           yloc(2) =  sqrt(1.0_rp/3.0_rp)
           zloc(2) = -sqrt(1.0_rp/3.0_rp)
           xloc(3) = -1.0_rp
           yloc(3) =  sqrt(1.0_rp/3.0_rp)
           zloc(3) =  sqrt(1.0_rp/3.0_rp)
           xloc(4) = -1.0_rp
           yloc(4) = -sqrt(1.0_rp/3.0_rp)
           zloc(4) =  sqrt(1.0_rp/3.0_rp)
        end if
        !
        ! Calculate derivatives
        !
        do inodf = 1,pnodf
           ! Hex shape function derivatives
           ! Replace derivatives by the ones from alya!!!
           dshxi(1,1) = -1.0_rp/8.0_rp*(1.0_rp - yloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(1,2) = -1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(1,3) = -1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp - yloc(inodf))
           dshxi(2,1) =  1.0_rp/8.0_rp*(1.0_rp - yloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(2,2) = -1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(2,3) = -1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp - yloc(inodf))
           dshxi(3,1) =  1.0_rp/8.0_rp*(1.0_rp + yloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(3,2) =  1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(3,3) = -1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp + yloc(inodf))
           dshxi(4,1) = -1.0_rp/8.0_rp*(1.0_rp + yloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(4,2) =  1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp - zloc(inodf))
           dshxi(4,3) = -1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp + yloc(inodf))
           dshxi(5,1) = -1.0_rp/8.0_rp*(1.0_rp - yloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(5,2) = -1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(5,3) =  1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp - yloc(inodf))
           dshxi(6,1) =  1.0_rp/8.0_rp*(1.0_rp - yloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(6,2) = -1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(6,3) =  1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp - yloc(inodf))
           dshxi(7,1) =  1.0_rp/8.0_rp*(1.0_rp + yloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(7,2) =  1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(7,3) =  1.0_rp/8.0_rp*(1.0_rp + xloc(inodf))*(1.0_rp + yloc(inodf))
           dshxi(8,1) = -1.0_rp/8.0_rp*(1.0_rp + yloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(8,2) =  1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp + zloc(inodf))
           dshxi(8,3) =  1.0_rp/8.0_rp*(1.0_rp - xloc(inodf))*(1.0_rp + yloc(inodf))

           dXdXi   = 0.0_rp
           dXdEta  = 0.0_rp
           dXdZeta = 0.0_rp
           dYdXi   = 0.0_rp
           dYdEta  = 0.0_rp
           dYdZeta = 0.0_rp
           dZdXi   = 0.0_rp
           dZdEta  = 0.0_rp
           dZdZeta = 0.0_rp
           do k = 1,pnode
              dXdXi   = dXdXi   + dshxi(k,1)*elcod(1,k)
              dXdEta  = dXdEta  + dshxi(k,2)*elcod(1,k)
              dXdZeta = dXdZeta + dshxi(k,3)*elcod(1,k)

              dYdXi   = dYdXi   + dshxi(k,1)*elcod(2,k)
              dYdEta  = dYdEta  + dshxi(k,2)*elcod(2,k)
              dYdZeta = dYdZeta + dshxi(k,3)*elcod(2,k)

              dZdXi   = dZdXi   + dshxi(k,1)*elcod(3,k)
              dZdEta  = dZdEta  + dshxi(k,2)*elcod(3,k)
              dZdZeta = dZdZeta + dshxi(k,3)*elcod(3,k)
           end do
           !
           ! Jacobian of the mapping
           !
           if(      iface == 1 .or. iface == 2 ) then
              ! zeta = constant on this face
              a = sqrt((dYdXi*dZdEta - dYdEta*dZdXi)**2+(dXdXi*dZdEta - dXdEta*dZdXi)**2+(dXdXi*dYdEta - dXdEta*dYdXi)**2)
           else if( iface == 3 .or. iface == 4 ) then
              ! eta = constant on this face
              a = sqrt((dYdXi*dZdZeta - dYdZeta*dZdXi)**2+(dXdXi*dZdZeta - dXdZeta*dZdXi)**2+(dXdXi*dYdZeta - dXdZeta*dYdXi)**2)
           else if( iface == 5 .or. iface == 6 ) then
              ! xi = constant on this face
              a = sqrt((dYdEta*dZdZeta - dYdZeta*dZdEta)**2+(dXdEta*dZdZeta - dXdZeta*dZdEta)**2+(dXdEta*dYdZeta - dXdZeta*dYdEta)**2)
           end if

           area_face = area_face + a

        end do
        lface(iface) = area_face
     end do

     aface(1)     = maxval(lface)
     aface(ndime) = minval(lface)

  else

     !-------------------------------------------------------------------
     !
     ! UNKNOWN ELEMENT
     !
     !-------------------------------------------------------------------

     call runend('ELMAFA: ELEMENT TYPE NOT PROGRAMMED')

  end if


end subroutine elmafa
