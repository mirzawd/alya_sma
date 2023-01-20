!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_outgau
  !****f* Levels/lev_outgau
  ! NAME 
  !    lev_outgau
  ! DESCRIPTION
  !    Postprocess Interface Gauges
  ! USES
  !    
  ! USED BY
  !    lev_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  implicit none

  integer(ip)             :: tetra(4,6)
  integer(ip)             :: ielem, idime, inode, ipoin, ii, ij
  integer(ip)             :: pelty, pnode
  integer(ip)             :: signn, signp, sigtn, sigtp           ! sign inside the element 
  integer(ip)             :: compl, inod1, inod2
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode), inter(3,4)
  real(rp)                :: l1, lp, x1, y1, z1, x2, y2, z2
  real(rp)                :: trian(ndime,ndime)    ! interface triangle

  if( INOTMASTER ) then
     if(ndime==2) then

        print*,' lev_outgau 2d case not implemented '

     else if(ndime==3) then

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

        !
        ! Loop over the elements to construct discrete interface 
        ! collection of triangle (3d case)
        !
        do ielem=1, nelem

           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           pnode=nnode(pelty)

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
                    compl=0_ip
                    inod1=0_ip

                    do ij=1,4

                       if(ellev(tetra(ij,ii))>=0_rp) then
                          sigtp=sigtp+1
                       else 
                          sigtn=sigtn+1
                       endif

                    end do

                    if(sigtn==1) then

                       !  research of one interface triangle
                       do ij=1,4
                          if(ellev(tetra(ij,ii))<0_rp) then
                             inod1=ij
                          endif
                       end do

                       do ij=1,4

                          if(ij/=inod1) then
                             compl = compl + 1

                             !
                             ! Compute the intersection of the elements 
                             ! with the surface 
                             !
                             l1=abs(ellev(tetra(inod1,ii)))
                             lp=abs(ellev(tetra(inod1,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod1,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod1,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod1,ii))
                             z2=elcod(3,tetra(ij,ii))
                             trian(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                             trian(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                             trian(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 
                          endif

                       enddo

                       call lev_calgau(trian)

                    else if(sigtp==1) then

                       !  research of one interface triangle
                       do ij=1,4
                          if(ellev(tetra(ij,ii))>=0_rp) then
                             inod1=ij
                          endif
                       end do

                       do ij=1,4

                          if(ij/=inod1) then

                             compl=compl+1
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
                             trian(1,compl) =  x1*(1-l1/lp)+x2*l1/lp  
                             trian(2,compl) =  y1*(1-l1/lp)+y2*l1/lp 
                             trian(3,compl) =  z1*(1-l1/lp)+z2*l1/lp 

                          endif

                       end do

                       call lev_calgau(trian)

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

                             compl=compl+1
                             l1=abs(ellev(tetra(inod1,ii)))
                             lp=abs(ellev(tetra(inod1,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod1,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod1,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod1,ii))
                             z2=elcod(3,tetra(ij,ii))

                             inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                             inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                             inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                             compl=compl+1
                             l1=abs(ellev(tetra(inod2,ii)))
                             lp=abs(ellev(tetra(inod2,ii))-ellev(tetra(ij,ii)))
                             x1=elcod(1,tetra(inod2,ii))
                             x2=elcod(1,tetra(ij,ii))
                             y1=elcod(2,tetra(inod2,ii))
                             y2=elcod(2,tetra(ij,ii))
                             z1=elcod(3,tetra(inod2,ii))
                             z2=elcod(3,tetra(ij,ii))

                             inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                             inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                             inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                          endif

                       end do

                       trian(1,1)=  inter(1,1)
                       trian(2,1)=  inter(2,1)
                       trian(3,1)=  inter(3,1) 
                       trian(1,2)=  inter(1,2)
                       trian(2,2)=  inter(2,2)
                       trian(3,2)=  inter(3,2) 
                       trian(1,3)=  inter(1,3)
                       trian(2,3)=  inter(2,3)
                       trian(3,3)=  inter(3,3) 

                       call lev_calgau(trian)

                       trian(1,1)=  inter(1,2)
                       trian(2,1)=  inter(2,2)
                       trian(3,1)=  inter(3,2) 
                       trian(1,2)=  inter(1,3)
                       trian(2,2)=  inter(2,3)
                       trian(3,2)=  inter(3,3) 
                       trian(1,3)=  inter(1,4)
                       trian(2,3)=  inter(2,4)
                       trian(3,3)=  inter(3,4) 

                       call lev_calgau(trian)

                    endif

                 enddo

              endif

           else if (pnode==4_ip) then

              print*,' lev_outgau tetraedral element case not implemented '

           endif

        enddo

     endif
  endif

end subroutine lev_outgau
