!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmv35(&
     pnode,elcod,ellev,voloc,tetra)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_elmv35
  ! NAME 
  !    lev_elmv35
  ! DESCRIPTION
  !    Compute in a 3D pyramid element the volume of 
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
  integer(ip), intent(in)  :: pnode,tetra(4_ip,2_ip)
  real(rp),    intent(in)  :: elcod(ndime,pnode),ellev(pnode)
  real(rp),    intent(out) :: voloc
  integer(ip)             :: compl,compt,type
  integer(ip)             :: inode,jnode, pnodt, ipyra
  integer(ip)             :: ino1l,ino2l,ino1g,ino2g
  real(rp)                :: l1,lp
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  real(rp)                :: phi1,phi2,phi3,phi4
  real(rp)                :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
  real(rp)                :: volel,voaux



  voloc = 0_rp        
  volel = 0_rp
  pnodt = 4_ip

  do ipyra = 1, 2_ip

     compt=0_ip

     do inode=1,pnodt
        if(ellev(tetra(inode,ipyra))>=0.0_rp) then
           compt=compt+1
        endif
     end do

     compl=0_ip


     x1 = elcod(1,tetra(1,ipyra))
     y1 = elcod(2,tetra(1,ipyra))
     z1 = elcod(3,tetra(1,ipyra))

     x2 = elcod(1,tetra(2,ipyra))
     y2 = elcod(2,tetra(2,ipyra))
     z2 = elcod(3,tetra(2,ipyra))

     x3 = elcod(1,tetra(3,ipyra))
     y3 = elcod(2,tetra(3,ipyra))
     z3 = elcod(3,tetra(3,ipyra))

     x4 = elcod(1,tetra(4,ipyra))
     y4 = elcod(2,tetra(4,ipyra))
     z4 = elcod(3,tetra(4,ipyra))

     phi1 = ellev(tetra(1,ipyra))
     phi2 = ellev(tetra(2,ipyra))
     phi3 = ellev(tetra(3,ipyra))
     phi4 = ellev(tetra(4,ipyra))

     if(compt==0_ip) then

        type=1_ip

     else if (compt==pnodt) then

        type=2_ip

        ! Compute volume of tetrahedra
        call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
        voloc = voloc + voaux

     else 

        type=3_ip

        ! We count the number of points in the liquid
        compl = 0_ip

        do inode=1,pnodt
           if(ellev(tetra(inode,ipyra))>=0_rp) then
              compl = compl + 1_ip
           endif
        end do


        if((compl==1_ip).or.(compl==3_ip)) then

           if(compl==1_ip) then

              do inode=1,pnodt
                 if(ellev(tetra(inode,ipyra))>=0_rp) then
                    jnode = inode
                 endif
              end do

           else if(compl==3_ip) then

              do inode=1,pnodt
                 if(ellev(tetra(inode,ipyra))<0_rp) then
                    jnode = inode
                 endif
              end do

              ! Compute volume of tetrahedra
              call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,volel)
           endif


           x1=elcod(1,tetra(jnode,ipyra))
           y1=elcod(2,tetra(jnode,ipyra))
           z1=elcod(3,tetra(jnode,ipyra))

           compt=0_ip

           do inode=1,pnodt

              if(inode/=jnode) then

                 compt=compt+1

                 !
                 ! Compute the intersection of the elements with the surface 
                 !
                 l1=abs(ellev(tetra(jnode,ipyra)))
                 lp=abs(ellev(tetra(jnode,ipyra))-ellev(tetra(inode,ipyra)))
                 x2=elcod(1,tetra(inode,ipyra))
                 y2=elcod(2,tetra(inode,ipyra))
                 z2=elcod(3,tetra(inode,ipyra))

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
              if(ellev(tetra(inode,ipyra))>=0_rp) then
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
              if(ellev(tetra(inode,ipyra))<0_rp) then
                 if(compt==0_ip) then 
                    ino1g = inode
                 else if(compt==1_ip) then 
                    ino2g = inode
                 endif
                 compt=compt+1_ip 
              endif
           end do

           l1=abs(ellev(tetra(ino1l,ipyra)))
           lp=abs(ellev(tetra(ino1l,ipyra))-ellev(tetra(ino1g,ipyra)))
           x1=elcod(1,tetra(ino1l,ipyra))
           y1=elcod(2,tetra(ino1l,ipyra))
           z1=elcod(3,tetra(ino1l,ipyra))
           x2=elcod(1,tetra(ino1g,ipyra))
           y2=elcod(2,tetra(ino1g,ipyra))
           z2=elcod(3,tetra(ino1g,ipyra))

           xa = x1*(1-l1/lp)+x2*l1/lp  
           ya = y1*(1-l1/lp)+y2*l1/lp 
           za = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino1l,ipyra)))
           lp=abs(ellev(tetra(ino1l,ipyra))-ellev(tetra(ino2g,ipyra)))
           x1=elcod(1,tetra(ino1l,ipyra))
           y1=elcod(2,tetra(ino1l,ipyra))
           z1=elcod(3,tetra(ino1l,ipyra))
           x2=elcod(1,tetra(ino2g,ipyra))
           y2=elcod(2,tetra(ino2g,ipyra))
           z2=elcod(3,tetra(ino2g,ipyra))

           xb = x1*(1-l1/lp)+x2*l1/lp  
           yb = y1*(1-l1/lp)+y2*l1/lp 
           zb = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino2l,ipyra)))
           lp=abs(ellev(tetra(ino2l,ipyra))-ellev(tetra(ino1g,ipyra)))
           x1=elcod(1,tetra(ino2l,ipyra))
           y1=elcod(2,tetra(ino2l,ipyra))
           z1=elcod(3,tetra(ino2l,ipyra))
           x2=elcod(1,tetra(ino1g,ipyra))
           y2=elcod(2,tetra(ino1g,ipyra))
           z2=elcod(3,tetra(ino1g,ipyra))

           xc = x1*(1-l1/lp)+x2*l1/lp  
           yc = y1*(1-l1/lp)+y2*l1/lp 
           zc = z1*(1-l1/lp)+z2*l1/lp 

           l1=abs(ellev(tetra(ino2l,ipyra)))
           lp=abs(ellev(tetra(ino2l,ipyra))-ellev(tetra(ino2g,ipyra)))
           x1=elcod(1,tetra(ino2l,ipyra))
           y1=elcod(2,tetra(ino2l,ipyra))
           z1=elcod(3,tetra(ino2l,ipyra))
           x2=elcod(1,tetra(ino2g,ipyra))
           y2=elcod(2,tetra(ino2g,ipyra))
           z2=elcod(3,tetra(ino2g,ipyra))

           xd = x1*(1-l1/lp)+x2*l1/lp  
           yd = y1*(1-l1/lp)+y2*l1/lp 
           zd = z1*(1-l1/lp)+z2*l1/lp 

           ! volume ino1l a b inod2l d c divided in three tetraedra

           x1=elcod(1,tetra(ino1l,ipyra))
           y1=elcod(2,tetra(ino1l,ipyra))
           z1=elcod(3,tetra(ino1l,ipyra))
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

           x1=elcod(1,tetra(ino1l,ipyra))
           y1=elcod(2,tetra(ino1l,ipyra))
           z1=elcod(3,tetra(ino1l,ipyra))
           x2=elcod(1,tetra(ino2l,ipyra))
           y2=elcod(2,tetra(ino2l,ipyra))
           z2=elcod(3,tetra(ino2l,ipyra))
           x3=xc              
           y3=yc            
           z3=zc            
           x4=xd              
           y4=yd            
           z4=zd 

           call voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,voaux)
           voloc = voloc + voaux

           x1=elcod(1,tetra(ino1l,ipyra))
           y1=elcod(2,tetra(ino1l,ipyra))
           z1=elcod(3,tetra(ino1l,ipyra))
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

end subroutine lev_elmv35
