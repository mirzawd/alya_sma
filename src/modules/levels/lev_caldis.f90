!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_caldis
  !-----------------------------------------------------------------------
  !****f* Levels/lev_caldis
  ! NAME 
  !    lev_caldis
  ! DESCRIPTION
  !    Compute the level set function redistanciation with Sussman equation
  ! USES
  !    
  ! USED BY
  !    lev_redieq
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  implicit none
  integer(ip)             :: tetra(4,6),tetrp(4_ip,3_ip),tetpy(4_ip,2_ip),ipasq,ipasp,ipapy
  integer(ip)             :: ielcr,ij,ii,idime,inode,jnode,ielem,ipoin
  integer(ip)             :: pelty,pnode,compt,compg, flag
  integer(ip)             :: inod1,inod2,signn,signp,sigtn,sigtp
  real(rp)                :: intlo(ndime,6) ! local interface piece
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode),inter(3,4)
  real(rp)                :: d,dx,dy,ratio,px,py,pz,l1,lp,dista
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,xa,ya

  if( INOTMASTER ) then

     if(ndime==2) then

        ! compute distance to interface on crossed element 
        do ielcr=1,nelcr_lev

           ielem=elcro_lev(ielcr)
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           !
           ! Gather operations
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode)=fleve(ipoin,1)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           compt=0_ip

           do inode=1,pnode

              if(inode==pnode) then
                 jnode=1 
              else
                 jnode=inode+1
              endif

              if(ellev(inode)*ellev(jnode)<=0.0_rp) then
                 compt=compt+1

                 !
                 ! Compute the intersection of the elements with the surface 
                 !
                 l1=abs(ellev(inode))
                 lp=abs(ellev(inode)-ellev(jnode))
                 x1=elcod(1,inode)
                 x2=elcod(1,jnode)
                 y1=elcod(2,inode)
                 y2=elcod(2,jnode)

                 intlo(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                 intlo(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 

              endif

           end do

           !
           ! Compute the distance 2d
           !
           do inode=1,pnode

              ipoin = lnods(inode,ielem)
              px    = coord(1,ipoin)
              py    = coord(2,ipoin)
              x1    = intlo(1,1)
              y1    = intlo(2,1)
              x2    = intlo(1,2)
              y2    = intlo(2,2)
              dx    = x2-x1
              dy    = y2-y1

              ratio = ((px-x1)*dx+(py-y1)*dy)/(dx*dx+dy*dy)

              if( ratio < 0.0_rp ) then 
                 d = sqrt((px-x1)*(px-x1)+(py-y1)*(py-y1))
              else if( ratio > 1.0_rp ) then 
                 d = sqrt((px-x2)*(px-x2)+(py-y2)*(py-y2))
              else
                 xa = (1.0_rp-ratio) * x1 + ratio * x2
                 ya = (1.0_rp-ratio) * y1 + ratio * y2
                 d  = sqrt((px-xa)*(px-xa)+(py-ya)*(py-ya))
              endif

              dista            = dista_lev(ipoin)
              dista            = min(dista,d)
              dista_lev(ipoin) = dista

           end do

        enddo

     else if(ndime==3) then

        ipasq = 0_ip
        ipasp = 0_ip
        ipapy = 0_ip

        ! compute distance to interface on crossed element 
        do ielcr=1,nelcr_lev

           ielem=elcro_lev(ielcr)
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           flag = 0_ip

           if((pnode==8_ip).and.(ipasq==0_ip)) then

              tetra(1,1)=1
              tetra(2,1)=3
              tetra(3,1)=4
              tetra(4,1)=5

              tetra(1,2)=3
              tetra(2,2)=4
              tetra(3,2)=5
              tetra(4,2)=8

              tetra(1,3)=1
              tetra(2,3)=2
              tetra(3,3)=3
              tetra(4,3)=5

              tetra(1,4)=2
              tetra(2,4)=3
              tetra(3,4)=5
              tetra(4,4)=6

              tetra(1,5)=3
              tetra(2,5)=5
              tetra(3,5)=6
              tetra(4,5)=7

              tetra(1,6)=3
              tetra(2,6)=5
              tetra(3,6)=7
              tetra(4,6)=8

              ipasq = 1

           else if((pnode==6_ip).and.(ipasp==0_ip)) then

              ipasp = 1

              tetrp(1,1)=1
              tetrp(2,1)=2
              tetrp(3,1)=3
              tetrp(4,1)=4

              tetrp(1,2)=2
              tetrp(2,2)=3
              tetrp(3,2)=4
              tetrp(4,2)=5

              tetrp(1,3)=3
              tetrp(2,3)=4
              tetrp(3,3)=5
              tetrp(4,3)=6

           else if((pnode==5_ip).and.(ipapy==0_ip)) then

              ipapy = 1

              tetpy(1,1)=1
              tetpy(2,1)=2
              tetpy(3,1)=3
              tetpy(4,1)=5

              tetpy(1,2)=1
              tetpy(2,2)=3
              tetpy(3,2)=4
              tetpy(4,2)=5

           endif
           !
           ! Gather operations
           !
           if(pnode==8_ip) then

              signn=0_ip
              signp=0_ip

              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 ellev(inode)=fleve(ipoin,1)
                 if(ellev(inode)>=0_rp) then
                    signp=signp+1
                 else 
                    signn=signn+1
                 endif

                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do

              if(signp/=8.and.signn/=8) then

                 do ii=1,6

                    sigtn=0_ip
                    sigtp=0_ip
                    compt=0_ip
                    inod1=0_ip

                    do ij=1,4
                       if(ellev(tetra(ij,ii))>=0_rp) then
                          sigtp=sigtp+1
                       else 
                          sigtn=sigtn+1
                       endif
                    end do

                    if((sigtn==1_ip).or.(sigtp==1_ip)) then
                       !  research of one interface triangle

                       if(sigtn==1_ip) then
                          do ij=1,4
                             if(ellev(tetra(ij,ii))<0_rp) then
                                inod1=ij
                             endif
                          end do
                       else
                          do ij=1,4
                             if(ellev(tetra(ij,ii))>=0_rp) then
                                inod1=ij
                             endif
                          end do
                       endif

                       do ij=1,4

                          if(ij/=inod1) then

                             compt=compt+1
                             !
                             ! Compute the intersection of the elements with the surface 
                             !
                             l1=abs(ellev(tetra(inod1,ii)))
                             lp=abs(ellev(tetra(inod1,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod1,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod1,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod1,ii))
                             z2=elcod(3,tetra(ij,ii))

                             intlo(1,compt) =  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             intlo(2,compt) =  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             intlo(3,compt) =  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                    else if(sigtn==2) then

                       !  research of two interface triangles
                       do ij=1,4
                          if(ellev(tetra(ij,ii))<0_rp) then
                             if(inod1==0) then
                                inod1=ij
                             else 
                                inod2=ij
                             endif
                          endif
                       end do

                       do ij=1,4
                          if(ij/=inod1.and.ij/=inod2) then
                             !
                             ! Compute the intersection of the elements with the surface 
                             !

                             compt=compt+1
                             l1=abs(ellev(tetra(inod1,ii)))
                             lp=abs(ellev(tetra(inod1,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod1,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod1,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod1,ii))
                             z2=elcod(3,tetra(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 


                             compt=compt+1
                             l1=abs(ellev(tetra(inod2,ii)))
                             lp=abs(ellev(tetra(inod2,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod2,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod2,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod2,ii))
                             z2=elcod(3,tetra(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                       compg=0_ip
                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,1)
                       intlo(2,compg)=  inter(2,1)
                       intlo(3,compg)=  inter(3,1) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,4)
                       intlo(2,compg)=  inter(2,4)
                       intlo(3,compg)=  inter(3,4) 

                    endif

                    if((sigtn>0_ip).and.(sigtp>0_ip)) then

                       do inode=1,pnode

                          ipoin=lnods(inode,ielem)
                          px=coord(1,ipoin)
                          py=coord(2,ipoin)
                          pz=coord(3,ipoin)

                          x1=intlo(1,1)
                          y1=intlo(2,1)
                          z1=intlo(3,1)

                          x2=intlo(1,2)
                          y2=intlo(2,2)
                          z2=intlo(3,2)

                          x3=intlo(1,3)
                          y3=intlo(2,3)
                          z3=intlo(3,3)

                          call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                          dista=dista_lev(ipoin)
                          dista=min(dista,d)
                          dista_lev(ipoin)=dista

                          if(sigtn==2) then

                             x1=intlo(1,4)
                             y1=intlo(2,4)
                             z1=intlo(3,4)

                             x2=intlo(1,5)
                             y2=intlo(2,5)
                             z2=intlo(3,5)

                             x3=intlo(1,6)
                             y3=intlo(2,6)
                             z3=intlo(3,6)

                             call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                             dista=dista_lev(ipoin)
                             dista=min(dista,d)
                             dista_lev(ipoin)=dista

                          endif

                       end do

                    endif


                 end do

              endif

           else if(pnode==6_ip) then

              signn=0_ip
              signp=0_ip

              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 ellev(inode)=fleve(ipoin,1)
                 if(ellev(inode)>=0_rp) then
                    signp=signp+1
                 else 
                    signn=signn+1
                 endif

                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do

              if(signp/=6_ip.and.signn/=6_ip) then

                 do ii=1,3

                    sigtn=0_ip
                    sigtp=0_ip
                    compt=0_ip
                    inod1=0_ip

                    do ij=1,4
                       if(ellev(tetrp(ij,ii))>=0_rp) then
                          sigtp=sigtp+1
                       else 
                          sigtn=sigtn+1
                       endif
                    end do

                    if((sigtn==1_ip).or.(sigtp==1_ip)) then
                       !  research of one interface triangle

                       if(sigtn==1_ip) then
                          do ij=1,4
                             if(ellev(tetrp(ij,ii))<0_rp) then
                                inod1=ij
                             endif
                          end do
                       else
                          do ij=1,4
                             if(ellev(tetrp(ij,ii))>=0_rp) then
                                inod1=ij
                             endif
                          end do
                       endif

                       do ij=1,4

                          if(ij/=inod1) then

                             compt=compt+1
                             !
                             ! Compute the intersection of the elements with the surface 
                             !
                             l1=abs(ellev(tetrp(inod1,ii)))
                             lp=abs(ellev(tetrp(inod1,ii))-ellev(tetrp(ij,ii)))
                             x1=elcod(1,tetrp(inod1,ii))
                             x2=elcod(1,tetrp(ij,ii))
                             y1=elcod(2,tetrp(inod1,ii))
                             y2=elcod(2,tetrp(ij,ii))
                             z1=elcod(3,tetrp(inod1,ii))
                             z2=elcod(3,tetrp(ij,ii))

                             intlo(1,compt) =  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             intlo(2,compt) =  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             intlo(3,compt) =  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                    else if(sigtn==2) then

                       !  research of two interface triangles
                       do ij=1,4
                          if(ellev(tetrp(ij,ii))<0_rp) then
                             if(inod1==0) then
                                inod1=ij
                             else 
                                inod2=ij
                             endif
                          endif
                       end do

                       do ij=1,4
                          if(ij/=inod1.and.ij/=inod2) then
                             !
                             ! Compute the intersection of the elements with the surface 
                             !

                             compt=compt+1
                             l1=abs(ellev(tetrp(inod1,ii)))
                             lp=abs(ellev(tetrp(inod1,ii))-ellev(tetrp(ij,ii)))
                             x1=elcod(1,tetrp(inod1,ii))
                             x2=elcod(1,tetrp(ij,ii))
                             y1=elcod(2,tetrp(inod1,ii))
                             y2=elcod(2,tetrp(ij,ii))
                             z1=elcod(3,tetrp(inod1,ii))
                             z2=elcod(3,tetrp(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 


                             compt=compt+1
                             l1=abs(ellev(tetrp(inod2,ii)))
                             lp=abs(ellev(tetrp(inod2,ii))-ellev(tetrp(ij,ii)))
                             x1=elcod(1,tetrp(inod2,ii))
                             x2=elcod(1,tetrp(ij,ii))
                             y1=elcod(2,tetrp(inod2,ii))
                             y2=elcod(2,tetrp(ij,ii))
                             z1=elcod(3,tetrp(inod2,ii))
                             z2=elcod(3,tetrp(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                       compg=0_ip
                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,1)
                       intlo(2,compg)=  inter(2,1)
                       intlo(3,compg)=  inter(3,1) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,4)
                       intlo(2,compg)=  inter(2,4)
                       intlo(3,compg)=  inter(3,4) 

                    endif

                    if((sigtn>0_ip).and.(sigtp>0_ip)) then

                       do inode=1,pnode

                          ipoin=lnods(inode,ielem)
                          px=coord(1,ipoin)
                          py=coord(2,ipoin)
                          pz=coord(3,ipoin)

                          x1=intlo(1,1)
                          y1=intlo(2,1)
                          z1=intlo(3,1)

                          x2=intlo(1,2)
                          y2=intlo(2,2)
                          z2=intlo(3,2)

                          x3=intlo(1,3)
                          y3=intlo(2,3)
                          z3=intlo(3,3)

                          call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                          dista=dista_lev(ipoin)
                          dista=min(dista,d)
                          dista_lev(ipoin)=dista

                          if(sigtn==2) then

                             x1=intlo(1,4)
                             y1=intlo(2,4)
                             z1=intlo(3,4)

                             x2=intlo(1,5)
                             y2=intlo(2,5)
                             z2=intlo(3,5)

                             x3=intlo(1,6)
                             y3=intlo(2,6)
                             z3=intlo(3,6)

                             call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                             dista=dista_lev(ipoin)
                             dista=min(dista,d)
                             dista_lev(ipoin)=dista

                          endif

                       end do

                    endif


                 end do

              endif

           else if(pnode==5_ip) then

              signn=0_ip
              signp=0_ip

              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 ellev(inode)=fleve(ipoin,1)
                 if(ellev(inode)>=0_rp) then
                    signp=signp+1
                 else 
                    signn=signn+1
                 endif

                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do

              if(signp/=5_ip.and.signn/=5_ip) then

                 do ii=1,2

                    sigtn=0_ip
                    sigtp=0_ip
                    compt=0_ip
                    inod1=0_ip

                    do ij=1,4
                       if(ellev(tetpy(ij,ii))>=0_rp) then
                          sigtp=sigtp+1
                       else 
                          sigtn=sigtn+1
                       endif
                    end do

                    if((sigtn==1_ip).or.(sigtp==1_ip)) then
                       !  research of one interface triangle

                       if(sigtn==1_ip) then
                          do ij=1,4
                             if(ellev(tetpy(ij,ii))<0_rp) then
                                inod1=ij
                             endif
                          end do
                       else
                          do ij=1,4
                             if(ellev(tetpy(ij,ii))>=0_rp) then
                                inod1=ij
                             endif
                          end do
                       endif

                       do ij=1,4

                          if(ij/=inod1) then

                             compt=compt+1
                             !
                             ! Compute the intersection of the elements with the surface 
                             !
                             l1=abs(ellev(tetpy(inod1,ii)))
                             lp=abs(ellev(tetpy(inod1,ii))-ellev(tetpy(ij,ii)))
                             x1=elcod(1,tetpy(inod1,ii))
                             x2=elcod(1,tetpy(ij,ii))
                             y1=elcod(2,tetpy(inod1,ii))
                             y2=elcod(2,tetpy(ij,ii))
                             z1=elcod(3,tetpy(inod1,ii))
                             z2=elcod(3,tetpy(ij,ii))

                             intlo(1,compt) =  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             intlo(2,compt) =  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             intlo(3,compt) =  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                    else if(sigtn==2) then

                       !  research of two interface triangles
                       do ij=1,4
                          if(ellev(tetpy(ij,ii))<0_rp) then
                             if(inod1==0) then
                                inod1=ij
                             else 
                                inod2=ij
                             endif
                          endif
                       end do

                       do ij=1,4
                          if(ij/=inod1.and.ij/=inod2) then
                             !
                             ! Compute the intersection of the elements with the surface 
                             !

                             compt=compt+1
                             l1=abs(ellev(tetpy(inod1,ii)))
                             lp=abs(ellev(tetpy(inod1,ii))-ellev(tetpy(ij,ii)))
                             x1=elcod(1,tetpy(inod1,ii))
                             x2=elcod(1,tetpy(ij,ii))
                             y1=elcod(2,tetpy(inod1,ii))
                             y2=elcod(2,tetpy(ij,ii))
                             z1=elcod(3,tetpy(inod1,ii))
                             z2=elcod(3,tetpy(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 


                             compt=compt+1
                             l1=abs(ellev(tetpy(inod2,ii)))
                             lp=abs(ellev(tetpy(inod2,ii))-ellev(tetpy(ij,ii)))
                             x1=elcod(1,tetpy(inod2,ii))
                             x2=elcod(1,tetpy(ij,ii))
                             y1=elcod(2,tetpy(inod2,ii))
                             y2=elcod(2,tetpy(ij,ii))
                             z1=elcod(3,tetpy(inod2,ii))
                             z2=elcod(3,tetpy(ij,ii))

                             inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                             inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                             inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                          endif

                       end do

                       compg=0_ip
                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,1)
                       intlo(2,compg)=  inter(2,1)
                       intlo(3,compg)=  inter(3,1) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,2)
                       intlo(2,compg)=  inter(2,2)
                       intlo(3,compg)=  inter(3,2) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,3)
                       intlo(2,compg)=  inter(2,3)
                       intlo(3,compg)=  inter(3,3) 

                       compg=compg+1_ip
                       intlo(1,compg)=  inter(1,4)
                       intlo(2,compg)=  inter(2,4)
                       intlo(3,compg)=  inter(3,4) 

                    endif

                    if((sigtn>0_ip).and.(sigtp>0_ip)) then

                       do inode=1,pnode

                          ipoin=lnods(inode,ielem)
                          px=coord(1,ipoin)
                          py=coord(2,ipoin)
                          pz=coord(3,ipoin)

                          x1=intlo(1,1)
                          y1=intlo(2,1)
                          z1=intlo(3,1)

                          x2=intlo(1,2)
                          y2=intlo(2,2)
                          z2=intlo(3,2)

                          x3=intlo(1,3)
                          y3=intlo(2,3)
                          z3=intlo(3,3)

                          call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                          dista=dista_lev(ipoin)
                          dista=min(dista,d)
                          dista_lev(ipoin)=dista

                          if(sigtn==2) then

                             x1=intlo(1,4)
                             y1=intlo(2,4)
                             z1=intlo(3,4)

                             x2=intlo(1,5)
                             y2=intlo(2,5)
                             z2=intlo(3,5)

                             x3=intlo(1,6)
                             y3=intlo(2,6)
                             z3=intlo(3,6)

                             call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                             dista=dista_lev(ipoin)
                             dista=min(dista,d)
                             dista_lev(ipoin)=dista

                          endif

                       end do

                    endif


                 end do

              endif

           else if(pnode==4_ip) then

              signn=0_ip
              signp=0_ip
              compt=0_ip
              inod1=0_ip

              do inode=1,pnode

                 ipoin=lnods(inode,ielem)
                 ellev(inode)=fleve(ipoin,1)

                 if(ellev(inode)>=0.0_rp) then
                    signp=signp+1
                 else 
                    signn=signn+1
                 endif

                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do

              end do

              if((signn==1_ip).or.(signp==1_ip)) then

                 !  research of one interface triangle
                 if(signn==1_ip) then
                    do ij=1,4
                       if(ellev(ij)<0.0_rp) then
                          inod1=ij
                       endif
                    end do
                 else if(signp==1_ip) then
                    do ij=1,4
                       if(ellev(ij)>=0.0_rp) then
                          inod1=ij
                       endif
                    end do
                 endif

                 do ij=1,4
                    if(ij/=inod1) then

                       compt=compt+1
                       !
                       ! Compute the intersection of the elements with the surface 
                       !
                       l1=abs(ellev(inod1))
                       lp=abs(ellev(inod1)-ellev(ij))
                       x1=elcod(1,inod1)
                       x2=elcod(1,ij)
                       y1=elcod(2,inod1)
                       y2=elcod(2,ij)
                       z1=elcod(3,inod1)
                       z2=elcod(3,ij)

                       intlo(1,compt) =  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                       intlo(2,compt) =  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                       intlo(3,compt) =  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                    endif

                 end do


              else if(signp==2) then

                 !  research of two interface triangles
                 do ij=1,4
                    if(ellev(ij)<0.0_rp) then
                       if(inod1==0) then
                          inod1=ij
                       else 
                          inod2=ij
                       endif
                    endif
                 end do

                 do ij=1,4

                    if(ij/=inod1.and.ij/=inod2) then
                       !
                       ! Compute the intersection of the elements with the surface 
                       !

                       compt=compt+1
                       l1=abs(ellev(inod1))
                       lp=abs(ellev(inod1)-ellev(ij))
                       x1=elcod(1,inod1)
                       x2=elcod(1,ij)
                       y1=elcod(2,inod1)
                       y2=elcod(2,ij)
                       z1=elcod(3,inod1)
                       z2=elcod(3,ij)

                       inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                       inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                       inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 


                       compt=compt+1
                       l1=abs(ellev(inod2))
                       lp=abs(ellev(inod2)-ellev(ij))
                       x1=elcod(1,inod2)
                       x2=elcod(1,ij)
                       y1=elcod(2,inod2)
                       y2=elcod(2,ij)
                       z1=elcod(3,inod2)
                       z2=elcod(3,ij)

                       inter(1,compt)=  x1*(1.0_rp-l1/lp)+x2*l1/lp  
                       inter(2,compt)=  y1*(1.0_rp-l1/lp)+y2*l1/lp 
                       inter(3,compt)=  z1*(1.0_rp-l1/lp)+z2*l1/lp 

                    endif

                 end do

                 compg=0_ip
                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,1)
                 intlo(2,compg)=  inter(2,1)
                 intlo(3,compg)=  inter(3,1) 

                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,2)
                 intlo(2,compg)=  inter(2,2)
                 intlo(3,compg)=  inter(3,2) 

                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,3)
                 intlo(2,compg)=  inter(2,3)
                 intlo(3,compg)=  inter(3,3) 

                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,2)
                 intlo(2,compg)=  inter(2,2)
                 intlo(3,compg)=  inter(3,2) 

                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,3)
                 intlo(2,compg)=  inter(2,3)
                 intlo(3,compg)=  inter(3,3) 

                 compg=compg+1_ip
                 intlo(1,compg)=  inter(1,4)
                 intlo(2,compg)=  inter(2,4)
                 intlo(3,compg)=  inter(3,4) 

              endif

              if((signn>0_ip).and.(signp>0_ip))then
                 do inode=1,pnode

                    ipoin=lnods(inode,ielem)
                    px=coord(1,ipoin)
                    py=coord(2,ipoin)
                    pz=coord(3,ipoin)

                    x1=intlo(1,1)
                    y1=intlo(2,1)
                    z1=intlo(3,1)

                    x2=intlo(1,2)
                    y2=intlo(2,2)
                    z2=intlo(3,2)

                    x3=intlo(1,3)
                    y3=intlo(2,3)
                    z3=intlo(3,3)

                    call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                    dista=dista_lev(ipoin)
                    dista=min(dista,d)
                    dista_lev(ipoin)=dista

                    if(signn==2) then

                       x1=intlo(1,1)
                       y1=intlo(2,1)
                       z1=intlo(3,1)

                       x2=intlo(1,2)
                       y2=intlo(2,2)
                       z2=intlo(3,2)

                       x3=intlo(1,3)
                       y3=intlo(2,3)
                       z3=intlo(3,3)

                       call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                       dista=dista_lev(ipoin)
                       dista=min(dista,d)
                       dista_lev(ipoin)=dista

                    endif

                 enddo

              endif

           endif

        end do

     endif

  endif

end subroutine lev_caldis
