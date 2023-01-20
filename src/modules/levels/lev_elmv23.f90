!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv23(&
     pnode,pgaus,pelty,elcod,ellev,voloc)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv23
  ! NAME 
  !    lev_elmv23
  ! DESCRIPTION
  !    Compute in a 2D triangular element the volume of 
  !    the phase with positive level set (the liquid normally)
  !                                    
  ! USES
  !    
  ! USED BY
  !    lev_calvol
  !    lev_l1norm
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp  
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: pelty
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)              :: compl,compt,type
  integer(ip)              :: inode,jnode,igaus
  real(rp)                 :: l1,lp,volel
  real(rp)                 :: x1,y1,x2,y2
  real(rp)                 :: xa,ya,xb,yb,xc,yc
  real(rp)                 :: xjaci(9),xjacm(9) 
  real(rp)                 :: gpcar(ndime,mnode,pgaus)                ! dNk/dxj
  real(rp)                 :: gpvol,gpdet                             ! |J|*w,|J|

  compt = 0_ip
  volel = 0.0_rp
  !
  ! Determine type of element
  ! type = 1 empty element:        (phi<0) in each node
  ! type = 2 full element:         (phi>=0) in each node
  ! type = 3 partly full element:  (phi<0) in some nodes and (phi>=0) in some nodes
  !
  do inode = 1,pnode
     if( ellev(inode) >= 0.0_rp ) then
        compt = compt + 1
     end if
  end do

  compl = 0_ip

  if( compt == 0_ip ) then
     !
     ! All are air
     !
     type  = 1_ip
     voloc = 0.0_rp

  else if( compt == pnode ) then
     !
     ! All are water
     !
     type = 2_ip
     do igaus = 1,pgaus
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
        gpvol = elmar(pelty) % weigp(igaus)*gpdet
        volel = volel + gpvol
     end do
     voloc = volel

  else 
     !
     ! Some are air and some are water
     !
     type = 3_ip
     do igaus = 1,pgaus
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
        gpvol = elmar(pelty) % weigp(igaus) * gpdet
        volel = volel + gpvol
     end do

     do inode = 1,pnode

        if( inode == pnode ) then
           jnode = 1 
        else
           jnode = inode + 1
        endif

        if( compt == 1 ) then
           if( ellev(inode) >= 0.0_rp ) then
              xc = elcod(1,inode)
              yc = elcod(2,inode)
           end if
        else if( compt == 2 ) then 
           if( ellev(inode) < 0.0_rp ) then
              xc = elcod(1,inode)
              yc = elcod(2,inode)
           end if
        else
           print*,'IMPOSSIBLE!'
        end if

        if( ellev(inode)*ellev(jnode) < 0.0_rp ) then
           compl = compl+1
           !
           ! Compute the intersection of the elements with the surface 
           !
           l1 = abs(ellev(inode))
           lp = abs(ellev(inode)-ellev(jnode))
           x1 = elcod(1,inode)
           x2 = elcod(1,jnode)
           y1 = elcod(2,inode)
           y2 = elcod(2,jnode)

           if( compl == 1 ) then
              xa = x1 * (1-l1/lp) + x2 * l1/lp  
              ya = y1 * (1-l1/lp) + y2 * l1/lp 
           else if( compl == 2 ) then
              xb = x1 * (1-l1/lp) + x2 * l1/lp  
              yb = y1 * (1-l1/lp) + y2 * l1/lp 
           end if
        end if
     end do

     if( compt == 1 ) then
        voloc = 0.5_rp*abs((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
     else if( compt == 2 ) then
        voloc = volel-0.5_rp*abs((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
     end if

  end if

end subroutine lev_elmv23



