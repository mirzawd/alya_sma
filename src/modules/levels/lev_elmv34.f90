!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv34(&
     pnode,elcod,ellev,voloc)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv34
  ! NAME 
  !    lev_elmv34
  ! DESCRIPTION
  !    Compute in a 3D tetraedral element the volume of 
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
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: elcod(ndime,pnode),ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)             :: compl,compt
  integer(ip)             :: inode,jnode
  integer(ip)             :: ino1l,ino2l,ino1g,ino2g
  real(rp)                :: l1,lp
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  real(rp)                :: phi1,phi2,phi3,phi4
  real(rp)                :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
  real(rp)                :: volel,voaux   

  compt=0_ip
  voloc=0_rp        
  volel=0_rp

  do inode=1,pnode
     if(ellev(inode)>=0.0_rp) then
        compt=compt+1
     endif
  end do

  compl=0_ip


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

  phi1 = ellev(1)
  phi2 = ellev(2)
  phi3 = ellev(3)
  phi4 = ellev(4)

  if(compt==0_ip) then

     voloc=0_rp

  else if (compt==pnode) then

     ! Compute volume of tetrahedra
     call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voloc)

  else 

     ! We count the number of point in the liquid
     compl = 0_ip

     do inode=1,pnode
        if(ellev(inode)>=0_rp) then
           compl = compl + 1_ip
        endif
     end do


     if((compl==1_ip).or.(compl==3_ip)) then

        if(compl==1_ip) then

           do inode=1,pnode
              if(ellev(inode)>=0_rp) then
                 jnode = inode
              endif
           end do

        else if(compl==3_ip) then

           do inode=1,pnode
              if(ellev(inode)<0_rp) then
                 jnode = inode
              endif
           end do

           ! Compute volume of tetrahedra
           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,volel)

        endif


        x1=elcod(1,jnode)
        y1=elcod(2,jnode)
        z1=elcod(3,jnode)

        compt=0_ip

        do inode=1,pnode

           if(inode/=jnode) then

              compt=compt+1

              !
              ! Compute the intersection of the elements with the surface 
              !
              l1=abs(ellev(jnode))
              lp=abs(ellev(jnode)-ellev(inode))
              x2=elcod(1,inode)
              y2=elcod(2,inode)
              z2=elcod(3,inode)


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
        do inode=1,pnode
           if(ellev(inode)>=0_rp) then
              if(compt==0_ip) then 
                 ino1l = inode
              else if(compt==1_ip) then 
                 ino2l = inode
              endif
              compt=compt+1_ip 
           endif
        end do

        compt=0_ip
        do inode=1,pnode
           if(ellev(inode)<0_rp) then
              if(compt==0_ip) then 
                 ino1g = inode
              else if(compt==1_ip) then 
                 ino2g = inode
              endif
              compt=compt+1_ip 
           endif
        end do

        l1=abs(ellev(ino1l))
        lp=abs(ellev(ino1l)-ellev(ino1g))
        x1=elcod(1,ino1l)
        y1=elcod(2,ino1l)
        z1=elcod(3,ino1l)
        x2=elcod(1,ino1g)
        y2=elcod(2,ino1g)
        z2=elcod(3,ino1g)

        xa = x1*(1-l1/lp)+x2*l1/lp  
        ya = y1*(1-l1/lp)+y2*l1/lp 
        za = z1*(1-l1/lp)+z2*l1/lp 

        l1=abs(ellev(ino1l))
        lp=abs(ellev(ino1l)-ellev(ino2g))
        x1=elcod(1,ino1l)
        y1=elcod(2,ino1l)
        z1=elcod(3,ino1l)
        x2=elcod(1,ino2g)
        y2=elcod(2,ino2g)
        z2=elcod(3,ino2g)

        xb = x1*(1-l1/lp)+x2*l1/lp  
        yb = y1*(1-l1/lp)+y2*l1/lp 
        zb = z1*(1-l1/lp)+z2*l1/lp 

        l1=abs(ellev(ino2l))
        lp=abs(ellev(ino2l)-ellev(ino1g))
        x1=elcod(1,ino2l)
        y1=elcod(2,ino2l)
        z1=elcod(3,ino2l)
        x2=elcod(1,ino1g)
        y2=elcod(2,ino1g)
        z2=elcod(3,ino1g)

        xc = x1*(1-l1/lp)+x2*l1/lp  
        yc = y1*(1-l1/lp)+y2*l1/lp 
        zc = z1*(1-l1/lp)+z2*l1/lp 

        l1=abs(ellev(ino2l))
        lp=abs(ellev(ino2l)-ellev(ino2g))
        x1=elcod(1,ino2l)
        y1=elcod(2,ino2l)
        z1=elcod(3,ino2l)
        x2=elcod(1,ino2g)
        y2=elcod(2,ino2g)
        z2=elcod(3,ino2g)

        xd = x1*(1-l1/lp)+x2*l1/lp  
        yd = y1*(1-l1/lp)+y2*l1/lp 
        zd = z1*(1-l1/lp)+z2*l1/lp 

        ! volume ino1l a b inod2l d c divided in three tetraedra

        x1=elcod(1,ino1l)
        y1=elcod(2,ino1l)
        z1=elcod(3,ino1l)
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

        x1=elcod(1,ino1l)
        y1=elcod(2,ino1l)
        z1=elcod(3,ino1l)
        x2=elcod(1,ino2l)
        y2=elcod(2,ino2l)
        z2=elcod(3,ino2l)
        x3=xc              
        y3=yc            
        z3=zc            
        x4=xd              
        y4=yd            
        z4=zd 

        call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
        voloc = voloc + voaux

        x1=elcod(1,ino1l)
        y1=elcod(2,ino1l)
        z1=elcod(3,ino1l)
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

  endif

end subroutine lev_elmv34
