!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv36(&
     pnode,elcod,ellev,voloc,tetra)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv36
  ! NAME 
  !    lev_elmv36
  ! DESCRIPTION
  !    Compute in a 3D prism element the volume of 
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
  integer(ip), intent(in)  :: pnode,tetra(4_ip,3_ip)
  real(rp),    intent(in)  :: elcod(ndime,pnode),ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)             :: compl,compt,type
  integer(ip)             :: inode,jnode, pnodt, itetr
  integer(ip)             :: ino1l,ino2l,ino1g,ino2g
  real(rp)                :: l1,lp
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  real(rp)                :: phi1,phi2,phi3,phi4
  real(rp)                :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
  real(rp)                :: volel,voaux



  voloc = 0_rp        
  volel = 0_rp
  pnodt = 4_ip

  do itetr = 1, 3_ip

     compt=0_ip

     do inode=1,pnodt
        if(ellev(tetra(inode,itetr))>=0.0_rp) then
           compt=compt+1
        endif
     end do

     compl=0_ip


     x1 = elcod(1,tetra(1,itetr))
     y1 = elcod(2,tetra(1,itetr))
     z1 = elcod(3,tetra(1,itetr))

     x2 = elcod(1,tetra(2,itetr))
     y2 = elcod(2,tetra(2,itetr))
     z2 = elcod(3,tetra(2,itetr))

     x3 = elcod(1,tetra(3,itetr))
     y3 = elcod(2,tetra(3,itetr))
     z3 = elcod(3,tetra(3,itetr))

     x4 = elcod(1,tetra(4,itetr))
     y4 = elcod(2,tetra(4,itetr))
     z4 = elcod(3,tetra(4,itetr))

     phi1 = ellev(tetra(1,itetr))
     phi2 = ellev(tetra(2,itetr))
     phi3 = ellev(tetra(3,itetr))
     phi4 = ellev(tetra(4,itetr))

     if(compt==0_ip) then

        type=1_ip

     else if (compt==pnodt) then

        type=2_ip

        ! Compute volume of tetrahedra
        call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
        voloc = voloc + voaux

     else 

        type=3_ip


        ! We count the number of point in the liquid
        compl = 0_ip

        do inode=1,pnodt
           if(ellev(tetra(inode,itetr))>=0_rp) then
              compl = compl + 1_ip
           endif
        end do


        if((compl==1_ip).or.(compl==3_ip)) then

           if(compl==1_ip) then

              do inode=1,pnodt
                 if(ellev(tetra(inode,itetr))>=0_rp) then
                    jnode = inode
                 endif
              end do

           else if(compl==3_ip) then

              do inode=1,pnodt
                 if(ellev(tetra(inode,itetr))<0_rp) then
                    jnode = inode
                 endif
              end do

              ! Compute volume of tetrahedra
              call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,volel)

           endif


           x1=elcod(1,tetra(jnode,itetr))
           y1=elcod(2,tetra(jnode,itetr))
           z1=elcod(3,tetra(jnode,itetr))

           compt=0_ip

           do inode=1,pnodt

              if(inode/=jnode) then

                 compt=compt+1

                 !
                 ! Compute the intersection of the elements with the surface 
                 !
                 l1=abs(ellev(tetra(jnode,itetr)))
                 lp=abs(ellev(tetra(jnode,itetr))-ellev(tetra(inode,itetr)))
                 x2=elcod(1,tetra(inode,itetr))
                 y2=elcod(2,tetra(inode,itetr))
                 z2=elcod(3,tetra(inode,itetr))

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
              if(ellev(tetra(inode,itetr))>=0_rp) then
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
              if(ellev(tetra(inode,itetr))<0_rp) then
                 if(compt==0_ip) then 
                    ino1g = inode
                 else if(compt==1_ip) then 
                    ino2g = inode
                 endif
                 compt=compt+1_ip 
              endif
           end do

           l1=abs(ellev(tetra(ino1l,itetr)))
           lp=abs(ellev(tetra(ino1l,itetr))-ellev(tetra(ino1g,itetr)))
           x1=elcod(1,tetra(ino1l,itetr))
           y1=elcod(2,tetra(ino1l,itetr))
           z1=elcod(3,tetra(ino1l,itetr))
           x2=elcod(1,tetra(ino1g,itetr))
           y2=elcod(2,tetra(ino1g,itetr))
           z2=elcod(3,tetra(ino1g,itetr))

           xa = x1*(1-l1/lp)+x2*l1/lp  
           ya = y1*(1-l1/lp)+y2*l1/lp 
           za = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino1l,itetr)))
           lp=abs(ellev(tetra(ino1l,itetr))-ellev(tetra(ino2g,itetr)))
           x1=elcod(1,tetra(ino1l,itetr))
           y1=elcod(2,tetra(ino1l,itetr))
           z1=elcod(3,tetra(ino1l,itetr))
           x2=elcod(1,tetra(ino2g,itetr))
           y2=elcod(2,tetra(ino2g,itetr))
           z2=elcod(3,tetra(ino2g,itetr))

           xb = x1*(1-l1/lp)+x2*l1/lp  
           yb = y1*(1-l1/lp)+y2*l1/lp 
           zb = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino2l,itetr)))
           lp=abs(ellev(tetra(ino2l,itetr))-ellev(tetra(ino1g,itetr)))
           x1=elcod(1,tetra(ino2l,itetr))
           y1=elcod(2,tetra(ino2l,itetr))
           z1=elcod(3,tetra(ino2l,itetr))
           x2=elcod(1,tetra(ino1g,itetr))
           y2=elcod(2,tetra(ino1g,itetr))
           z2=elcod(3,tetra(ino1g,itetr))

           xc = x1*(1-l1/lp)+x2*l1/lp  
           yc = y1*(1-l1/lp)+y2*l1/lp 
           zc = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino2l,itetr)))
           lp=abs(ellev(tetra(ino2l,itetr))-ellev(tetra(ino2g,itetr)))
           x1=elcod(1,tetra(ino2l,itetr))
           y1=elcod(2,tetra(ino2l,itetr))
           z1=elcod(3,tetra(ino2l,itetr))
           x2=elcod(1,tetra(ino2g,itetr))
           y2=elcod(2,tetra(ino2g,itetr))
           z2=elcod(3,tetra(ino2g,itetr))

           xd = x1*(1-l1/lp)+x2*l1/lp  
           yd = y1*(1-l1/lp)+y2*l1/lp 
           zd = z1*(1-l1/lp)+z2*l1/lp 

           ! volume ino1l a b inod2l d c divided in three tetraedra

           x1=elcod(1,tetra(ino1l,itetr))
           y1=elcod(2,tetra(ino1l,itetr))
           z1=elcod(3,tetra(ino1l,itetr))
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

           x1=elcod(1,tetra(ino1l,itetr))
           y1=elcod(2,tetra(ino1l,itetr))
           z1=elcod(3,tetra(ino1l,itetr))
           x2=elcod(1,tetra(ino2l,itetr))
           y2=elcod(2,tetra(ino2l,itetr))
           z2=elcod(3,tetra(ino2l,itetr))
           x3=xc              
           y3=yc            
           z3=zc            
           x4=xd              
           y4=yd            
           z4=zd 

           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
           voloc = voloc + voaux

           x1=elcod(1,tetra(ino1l,itetr))
           y1=elcod(2,tetra(ino1l,itetr))
           z1=elcod(3,tetra(ino1l,itetr))
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

  enddo

!  typa = typ1 + typ2 + typ3

end subroutine lev_elmv36
