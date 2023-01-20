!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv38(&
     pnode,pgaus,pelty,elcod,ellev,tetra,voloc)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv38
  ! NAME 
  !    lev_elmv38
  ! DESCRIPTION
  !    Compute in a 3D hexaedra element the volume of 
  !    the phase with positive level set (the liquid normally)
  !                                    
  ! USES
  !    
  ! USED BY
  !    lev_calvol
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp  
  use      def_parame
  use      def_master
  use      def_domain

  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,pelty,tetra(4_ip,6_ip)
  real(rp),    intent(in)  :: elcod(ndime,pnode),ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)             :: compl,compt,type,typ1,typ2,typ3,typa
  integer(ip)             :: inode,jnode,igaus,idime,itetr,pnodt
  integer(ip)             :: ino1l,ino2l,ino1g,ino2g
  real(rp)                :: elcot(ndime,4_ip)
  real(rp)                :: ellet(4_ip)
  real(rp)                :: l1,lp
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  real(rp)                :: x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8
  real(rp)                :: phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8
  real(rp)                :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
  real(rp)                :: volel,voaux 
  real(rp)                :: xjacm(ndime,ndime) 
  real(rp)                :: gpvol,gpdet                             ! |J|*w,|J|
  real(rp)                :: tragl(ndime,ndime)     ! Stabilization
  real(rp)                :: hleng(ndime)  

  compt=0_ip
  voloc=0_rp        
  volel=0_rp
  pnodt=4_ip

  typ1=0_ip
  typ2=0_ip
  typ3=0_ip
  typa=0_ip

  x1 = elcod(1,1)
  y1 = elcod(2,1)
  z1 = elcod(3,1)

  x2 = elcod(1,2)
  y2 = elcod(2,2)
  z2 = elcod(3,2)

  x3 = elcod(1,3)
  y3 = elcod(2,3)
  z3 = elcod(3,3)

  x4 = elcod(1,4)
  y4 = elcod(2,4)
  z4 = elcod(3,4)

  x5 = elcod(1,5)
  y5 = elcod(2,5)
  z5 = elcod(3,5)

  x6 = elcod(1,6)
  y6 = elcod(2,6)
  z6 = elcod(3,6)

  x7 = elcod(1,7)
  y7 = elcod(2,7)
  z7 = elcod(3,7)

  x8 = elcod(1,8)
  y8 = elcod(2,8)
  z8 = elcod(3,8)

  phi1 = ellev(1)
  phi2 = ellev(2)
  phi3 = ellev(3)
  phi4 = ellev(4)
  phi5 = ellev(5)
  phi6 = ellev(6)
  phi7 = ellev(7)
  phi8 = ellev(8)

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
        call jacdet(&
             ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
             xjacm,gpdet)
        gpvol=elmar(pelty)%weigp(igaus)*gpdet
        volel=volel+gpvol
     end do
     voloc=volel

  else 

     type=3_ip

     ! loop on the tetraedra
     do itetr=1,6_ip

        do inode=1,pnodt
           ellet(inode)=ellev(tetra(inode,itetr))
           do idime=1,ndime
              elcot(idime,inode)=elcod(idime,tetra(inode,itetr))
           end do
        end do

        x1 = elcot(1,1)
        y1 = elcot(2,1)
        z1 = elcot(3,1)

        x2 = elcot(1,2)
        y2 = elcot(2,2)
        z2 = elcot(3,2)

        x3 = elcot(1,3)
        y3 = elcot(2,3)
        z3 = elcot(3,3)

        x4 = elcot(1,4)
        y4 = elcot(2,4)
        z4 = elcot(3,4)

        phi1 = ellet(1)
        phi2 = ellet(2)
        phi3 = ellet(3)
        phi4 = ellet(4)

        ! We count the number of point in the liquid
        compl = 0_ip

        do inode=1,pnodt
           if(ellet(inode)>=0_rp) then
              compl = compl + 1_ip
           endif
        end do

        if((compl==1_ip).or.(compl==3_ip)) then

           if(compl==1_ip) then

              do inode=1,pnodt
                 if(ellet(inode)>=0_rp) then
                    jnode = inode
                 endif
              end do

           else if(compl==3_ip) then

              do inode=1,pnodt
                 if(ellet(inode)<0_rp) then
                    jnode = inode
                 endif
              end do

              ! Compute volume of tetrahedra
              call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,volel)

           endif


           x1=elcot(1,jnode)
           y1=elcot(2,jnode)
           z1=elcot(3,jnode)

           compt=0_ip

           do inode=1,pnodt

              if(inode/=jnode) then

                 compt=compt+1

                 !
                 ! Compute the intersection of the elements with the surface 
                 !
                 l1=abs(ellet(jnode))
                 lp=abs(ellet(jnode)-ellet(inode))
                 x2=elcot(1,inode)
                 y2=elcot(2,inode)
                 z2=elcot(3,inode)

                 if(compt==1) then
                    xa = x1*(1-l1/lp)+x2*l1/lp  
                    ya = y1*(1-l1/lp)+y2*l1/lp 
                    za = z1*(1-l1/lp)+z2*l1/lp 
                 else if(compt==2) then
                    xb = x1*(1-l1/lp)+x2*l1/lp  
                    yb = y1*(1-l1/lp)+y2*l1/lp 
                    zb = z1*(1-l1/lp)+z2*l1/lp 
                 else if(compt==3) then
                    xc = x1*(1-l1/lp)+x2*l1/lp  
                    yc = y1*(1-l1/lp)+y2*l1/lp 
                    zc = z1*(1-l1/lp)+z2*l1/lp 
                 endif

              endif

           enddo


           ! Compute volume of tetraedra of gaz or liquid
           call voltet(x1,y1,z1,xa,ya,za,xb,yb,zb,xc,yc,zc,voaux)

           if(compl==1_ip) then
              voloc = voloc + voaux
           else if(compl==3_ip) then
              voloc = voloc + volel - voaux
           endif


        else if(compl==2_ip) then

           compt=0_ip
           do inode=1,pnodt
              if(ellet(inode)>=0_rp) then
                 if(compt==0_ip) then 
                    ino1l = inode
                 else if(compt==1_ip) then 
                    ino2l = inode
                 endif
                 compt=compt+1_ip 
              endif
           end do

           compt=0_ip
           do inode=1,pnodt
              if(ellet(inode)<0_rp) then
                 if(compt==0_ip) then 
                    ino1g = inode
                 else if(compt==1_ip) then 
                    ino2g = inode
                 endif
                 compt=compt+1_ip 
              endif
           end do

           l1=abs(ellet(ino1l))
           lp=abs(ellet(ino1l)-ellet(ino1g))
           x1=elcot(1,ino1l)
           y1=elcot(2,ino1l)
           z1=elcot(3,ino1l)
           x2=elcot(1,ino1g)
           y2=elcot(2,ino1g)
           z2=elcot(3,ino1g)

           xa = x1*(1-l1/lp)+x2*l1/lp  
           ya = y1*(1-l1/lp)+y2*l1/lp 
           za = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellet(ino1l))
           lp=abs(ellet(ino1l)-ellet(ino2g))
           x1=elcot(1,ino1l)
           y1=elcot(2,ino1l)
           z1=elcot(3,ino1l)
           x2=elcot(1,ino2g)
           y2=elcot(2,ino2g)
           z2=elcot(3,ino2g)

           xb = x1*(1-l1/lp)+x2*l1/lp  
           yb = y1*(1-l1/lp)+y2*l1/lp 
           zb = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellet(ino2l))
           lp=abs(ellet(ino2l)-ellet(ino1g))
           x1=elcot(1,ino2l)
           y1=elcot(2,ino2l)
           z1=elcot(3,ino2l)
           x2=elcot(1,ino1g)
           y2=elcot(2,ino1g)
           z2=elcot(3,ino1g)

           xc = x1*(1-l1/lp)+x2*l1/lp  
           yc = y1*(1-l1/lp)+y2*l1/lp 
           zc = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellet(ino2l))
           lp=abs(ellet(ino2l)-ellet(ino2g))
           x1=elcot(1,ino2l)
           y1=elcot(2,ino2l)
           z1=elcot(3,ino2l)
           x2=elcot(1,ino2g)
           y2=elcot(2,ino2g)
           z2=elcot(3,ino2g)

           xd = x1*(1-l1/lp)+x2*l1/lp  
           yd = y1*(1-l1/lp)+y2*l1/lp 
           zd = z1*(1-l1/lp)+z2*l1/lp 

           ! volume ino1l a b inod2l d c divided in three tetraedra

           x1=elcot(1,ino1l)
           y1=elcot(2,ino1l)
           z1=elcot(3,ino1l)
           x2=xa              
           y2=ya            
           z2=za            
           x3=xb              
           y3=yb            
           z3=zb            
           x4=xd              
           y4=yd            
           z4=zd 

           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
           voloc = voloc + voaux

           x1=elcot(1,ino1l)
           y1=elcot(2,ino1l)
           z1=elcot(3,ino1l)
           x2=elcot(1,ino2l)
           y2=elcot(2,ino2l)
           z2=elcot(3,ino2l)
           x3=xc              
           y3=yc            
           z3=zc            
           x4=xd              
           y4=yd            
           z4=zd 

           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
           voloc = voloc + voaux

           x1=elcot(1,ino1l)
           y1=elcot(2,ino1l)
           z1=elcot(3,ino1l)
           x2=xa              
           y2=ya            
           z2=za            
           x3=xc              
           y3=yc            
           z3=zc            
           x4=xd              
           y4=yd            
           z4=zd 

           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
           voloc = voloc + voaux

        endif

     end do

  endif

  if(type==1_ip) then
     typ1 = 1_ip
  else if(type==2_ip) then
     typ2 = 1_ip
  else if(type==3_ip) then
     typ3 = 1_ip
  endif

  typa = typ1 + typ2 + typ3


end subroutine lev_elmv38
