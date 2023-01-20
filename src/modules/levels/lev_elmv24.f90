!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv24(&
     pnode,pgaus,pelty,elcod,ellev,voloc)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv24
  ! NAME 
  !    lev_elmv24
  ! DESCRIPTION
  !    Compute in a 2D quadrilateral element the volume of 
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
  use      def_parame
  use      def_master
  use      def_domain

  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,pelty
  real(rp),    intent(in)  :: elcod(ndime,pnode),ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)             :: compl,compt,type, iprem,iseco
  integer(ip)             :: inode,jnode,igaus
  real(rp)                :: l1,l2,lp,volel
  real(rp)                :: x1,y1,x2,y2
  real(rp)                :: xa,ya,xb,yb
  real(rp)                :: xjaci(ndime,ndime) 
  real(rp)                :: xjacm(ndime,ndime) 
  real(rp)                :: gpcar(ndime,mnode,pgaus)                ! dNk/dxj
  real(rp)                :: gpvol,gpdet                             ! |J|*w,|J|
  real(rp)                :: tragl(ndime,ndime)     ! Stabilization
  real(rp)                :: hleng(ndime)  

  compt=0_ip
  volel=0_rp

  !
  ! Determine type of element
  ! type= 1 empty element (phi<0) in each node
  ! type= 2 full element (phi>0) in each node
  ! type= 3 partly full element (phi<0) in some nodes and (phi>0) in some nodes
  !
  do inode=1,pnode
     if(ellev(inode)>=0.0_rp) then
        compt=compt+1
     endif
  end do

  compl=0_ip

  if(compt==0_ip) then

     type=1_ip
     voloc=0_rp

  else if (compt==pnode) then

     type=2_ip
     !
     ! hleng and tragl at center of gravity
     !
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)

     do igaus=1,pgaus
        !
        ! Cartesian derivatives, and volume: GPCAR, PGVOL
        !
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
        gpvol=elmar(pelty)%weigp(igaus)*gpdet
        volel=volel+gpvol
     end do
     voloc=volel

  else 

     type=3_ip

     !! attention specific case for element cartesian quadraliteral 
     !
     ! hleng and tragl at center of gravity
     !
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)

     do igaus=1,pgaus
        !
        ! Cartesian derivatives, and volume: GPCAR, PGVOL
        !
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
        gpvol=elmar(pelty)%weigp(igaus)*gpdet
        volel=volel+gpvol
     end do

     do inode=1,pnode

        if(inode==pnode) then
           jnode=1 
        else
           jnode=inode+1
        endif

        if(ellev(inode)*ellev(jnode)<=0.0_rp) then

           compl=compl+1
           if(compl==1) then
              iprem=jnode   
           else if(compl==2) then
              iseco=inode   
           endif
           !
           ! Compute the intersection of the elements with the surface 
           !
           l1=abs(ellev(inode))
           lp=abs(ellev(inode)-ellev(jnode))
           x1=elcod(1,inode)
           x2=elcod(1,jnode)
           y1=elcod(2,inode)
           y2=elcod(2,jnode)

           if(compl==1) then
              xa = x1*(1-l1/lp)+x2*l1/lp  
              ya = y1*(1-l1/lp)+y2*l1/lp 
           else if(compl==2) then
              xb = x1*(1-l1/lp)+x2*l1/lp  
              yb = y1*(1-l1/lp)+y2*l1/lp 
           endif

        endif

     end do

     !
     ! Treatement of different geometrical cases
     !
     if(iprem==iseco) then 
        !
        ! Case interface separation forms one triangle and one
        !
        x1=elcod(1,iprem)
        y1=elcod(2,iprem)
        ! attention 
        l1=sqrt((x1-xa)*(x1-xa)+(y1-ya)*(y1-ya))
        lp=sqrt((x1-xb)*(x1-xb)+(y1-yb)*(y1-yb))

        if(ellev(iprem)<=0.0_rp) then
           !
           ! Triangle empty
           !
           voloc=volel-l1*lp*0.5_rp
        else
           !
           ! Triangle full
           !
           voloc=l1*lp*0.5_rp

        endif

     else
        !
        ! Case interface separation forms two quadrilaterals
        !
        x1=elcod(1,iprem)
        y1=elcod(2,iprem)
        x2=elcod(1,iseco)
        y2=elcod(2,iseco)
        ! attention 
        l1=sqrt((x1-xa)*(x1-xa)+(y1-ya)*(y1-ya))
        l2=sqrt((x2-xb)*(x2-xb)+(y2-yb)*(y2-yb))
        lp=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

        if(ellev(iprem)<=0.0_rp) then
           !
           ! Trapeze empty
           !
           voloc=volel-0.5_rp*(l1+l2)*lp
        else
           !
           ! Trapeze full
           !
           voloc=0.5_rp*(l1+l2)*lp

        endif


     endif

  endif



end subroutine lev_elmv24



