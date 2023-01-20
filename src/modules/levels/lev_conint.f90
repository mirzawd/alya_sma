!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_conint()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_conint
  ! NAME 
  !    lev_conint
  ! DESCRIPTION
  !    Compute the dscrete interface : vector of segments in 2D
  !                                    vector of triangles in 3D
  ! USES
  !    
  ! USED BY
  !    lev_redist
  !
  ! TASK 1 : compute segment/triangle of interface
  ! TASK 2 : compute segment/triangle of interface orientated
  ! TASK 3 : deallocate memory for interface
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_levels
  use def_domain
  use def_solver
  use mod_memchk

  implicit none
  integer(ip)             :: ielem, idime, icomp     ! Indices and dimensions
  integer(ip)             :: pelty, pnode , p1 , p2 , p3
  integer(ip)             :: inode, jnode, ipoin, compt, compl, compg , compi
  integer(ip)             :: ii, ij, inod1, inod2         ! element crossed by interface = 1 
  integer(ip)             :: signn, signp, sigtn, sigtp   ! sign inside the element 
  integer(ip)             :: tetra(4_ip,6_ip), tetrp(4_ip,3_ip), tetpy(4_ip,2_ip)
  integer(ip)             :: ipasq, ipasp, ipapy

  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode),inter(3,4)
  real(rp)                :: l1,lp,x1,y1,x2,y2,z1,z2
  real(rp)                :: gpgrl(3),gpcar(ndime,mnode),xjaci(9),xjacm(9),gpdet,vec(3,3)
  real(rp)                :: xx,nx,ny,nz,coord_xxx(3)
  real(rp),    pointer    :: coord_tmp(:,:) 
  integer(ip), pointer    :: lnodb_tmp(:,:) 
  integer(ip), pointer    :: lebsu_tmp(:)    ! element (whole mesh numeration) to which the surface belongs
  integer(ip), pointer    :: lauxi(:)        ! auxiliary vector used to build lsbel_lev 
  integer(4)              :: istat

  nullify(coord_tmp)
  nullify(lnodb_tmp)
  nullify(lebsu_tmp)
  nullify(lauxi)

  !
  ! Does not work with divisor. maybe do a recursive fringe geometry
  !
  if( ndivi > 0 .and. tyred_lev > 2 ) call runend('LEV_CONINT: NOT READIY FOR AUTOMATIC DIVISION')

  compt=0_ip

  if( INOTMASTER ) then

     ! 
     ! Create array lnuew in the serial case so that // & serial cases can be writen in a unified way
     !
     if( ISEQUEN ) then
        if( tyred_lev > 2 ) then
           allocate(lnuew(nelem), stat=istat)
           call memchk( zero, istat, memor_dom, 'lnuew', 'lev_conint', lnuew ) 
           do ielem=1,nelem
              lnuew(ielem) = ielem
           end do
        end if
     end if

     if(ndime==2) then
        !
        ! Loop over the elements to construct discrete interface 
        ! collection of segments (2d case)
        !
        do ielem=1,nelem
           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           if( lelch(ielem) == ELFEM .or. lelch(ielem) == ELCUT) then
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

              compl=0_ip

              do inode=1,pnode
                 if(inode==pnode) then
                    jnode=1 
                 else
                    jnode=inode+1
                 endif
                 if(ellev(inode)*ellev(jnode)<=0.0_rp) then
                    ! count interface points
                    compl=compl+1
                    if(compl==3) then
                       compl=1
                    endif
                    if(compl==1) then
                       ! count interface segments
                       compt=compt+1
                    endif
                 endif
              end do

           end if
        end do

        !
        ! registering of discrete interface vectors size
        !
        if(kfl_paral==-1) then
           nboun_lev=compt
           npoin_lev=ndime*nboun_lev
        else if(kfl_paral>0) then
           nredi_lev(1)=compt
           nredi_lev(2)=ndime*nredi_lev(1)
           nredt_lev(1)=compt
           nredt_lev(2)=ndime*nredt_lev(1)
        endif

        compt=0_ip
        compg=0_ip
        !
        ! Memory allocation for discrete interface vectors
        ! connectivity (lnodb_lev) and points coordinates (coord_lev)
        !
        ! Slave:  LNODP_LEV, COORP_LEV
        ! Sequen: LNODB_LEV, COORD_LEV
        !
        call lev_memall(2_ip)
        if( ISEQUEN ) then
           coord_tmp => coord_lev
           lnodb_tmp => lnodb_lev
           lebsu_tmp => lebsu_lev
        else if( ISLAVE ) then
           coord_tmp => coorp_lev
           lnodb_tmp => lnodp_lev 
           lebsu_tmp => lebsp_lev             
        end if
        !
        ! filling of discrete interface vectors 
        ! connectivity (lnodb_lev) and points coordinates (coord_lev)
        !

        do ielem=1,nelem
           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           if( lelch(ielem) == ELFEM .or. lelch(ielem) == ELCUT) then
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

              compl=0_ip

              do inode=1,pnode
                 if(inode==pnode) then
                    jnode=1 
                 else
                    jnode=inode+1
                 endif
                 if(ellev(inode)*ellev(jnode)<=0.0_rp) then
                    compg=compg+1
                    compl=compl+1

                    if(compl==3) then
                       compl=1
                    endif
                    if(compl==1) then
                       compt=compt+1
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
                    lnodb_tmp(compl,compt)= compg
                    coord_tmp(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                    coord_tmp(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp
                    if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                 endif
              end do

           end if
        end do

     else if(ndime==3) then

        ipasq = 0_ip
        ipasp = 0_ip
        ipapy = 0_ip

        !
        ! Loop over the elements to construct discrete interface 
        ! collection of triangle (3d case)
        !
        do ielem=1, nelem

           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           if( lelch(ielem) == ELFEM .or. lelch(ielem) == ELCUT) then
              pnode=nnode(pelty)

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

              ! Case Q1
              if(pnode==8_ip) then

                 signn=0_ip
                 signp=0_ip
                 !
                 ! Gather operations
                 !
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


                 if(signp/=8.and.signn/=8) then


                    do ii=1,6
                       sigtn=0_ip
                       sigtp=0_ip

                       do ij=1,4
                          if(ellev(tetra(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       !! We count the number of interface triangle
                       if(sigtp==1.or.sigtn==1) then
                          compt=compt+1 
                       else if(sigtp==2.or.sigtn==2) then
                          compt=compt+2 
                       endif
                    end do
                 endif

                 ! Case Prism
              else if(pnode==6_ip) then

                 signn=0_ip
                 signp=0_ip
                 !
                 ! Gather operations
                 !
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

                 if(signp/=6_ip.and.signn/=6_ip) then

                    do ii=1,3
                       sigtn=0_ip
                       sigtp=0_ip

                       do ij=1,4
                          if(ellev(tetrp(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       !! We count the number of interface triangle
                       if(sigtp==1.or.sigtn==1) then
                          compt=compt+1 
                       else if(sigtp==2.or.sigtn==2) then
                          compt=compt+2 
                       endif
                    end do

                 endif

                 ! Case Pyramid
              else if(pnode==5_ip) then

                 signn=0_ip
                 signp=0_ip
                 !
                 ! Gather operations
                 !
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

                 if(signp/=5_ip.and.signn/=5_ip) then

                    do ii=1,2
                       sigtn=0_ip
                       sigtp=0_ip

                       do ij=1,4
                          if(ellev(tetpy(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       !! We count the number of interface triangle
                       if(sigtp==1.or.sigtn==1) then
                          compt=compt+1 
                       else if(sigtp==2.or.sigtn==2) then
                          compt=compt+2 
                       endif
                    end do

                 endif

              else if(pnode==4_ip) then

                 signn=0_ip
                 signp=0_ip
                 !
                 ! Gather operations
                 !
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

                 !! We count the number od interface triangle
                 if(signp==1.or.signn==1) then
                    compt=compt+1 
                 else if(signp==2.or.signn==2) then
                    compt=compt+2 
                 endif

              endif
           end if

        end do

        !
        ! registering of discrete interface vectors size
        !
        if( ISEQUEN ) then
           nboun_lev=compt
           npoin_lev=ndime*nboun_lev
        else if( ISLAVE ) then
           nredi_lev(1)=compt
           nredi_lev(2)=ndime*nredi_lev(1)
           nredt_lev(1)=compt
           nredt_lev(2)=ndime*nredt_lev(1)
        endif

        compt=0_ip
        compg=0_ip
        gpgrl(1)=0.0_rp
        gpgrl(2)=0.0_rp
        gpgrl(3)=0.0_rp
        
        !
        ! Memory allocation for discrete interface vectors
        ! connectivity (lnodb_lev) and points coordinates (coord_lev)
        !
        ! Slave:  LNODP_LEV, COORP_LEV
        ! Sequen: LNODB_LEV, COORD_LEV
        !
        call lev_memall(2_ip)
        if( ISEQUEN ) then
           coord_tmp => coord_lev
           lnodb_tmp => lnodb_lev
           lebsu_tmp => lebsu_lev
        else if( ISLAVE ) then
           coord_tmp => coorp_lev
           lnodb_tmp => lnodp_lev  
           lebsu_tmp => lebsp_lev            
        end if
        !
        ! filling of discrete interface vectors 
        ! connectivity (lnodb_lev) and points coordinates (coord_lev)
        !
        do ielem=1,nelem

           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           if( lelch(ielem) == ELFEM .or. lelch(ielem) == ELCUT) then
              pnode=nnode(pelty)
              compi=compt+1

              ! Case Q1
              if(pnode==8_ip) then

                 signn=0_ip
                 signp=0_ip
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

                 if(signp/=8.and.signn/=8) then
                    do ii=1,6
                       sigtn=0_ip
                       sigtp=0_ip
                       compl=0_ip
                       inod1=0_ip

                       do ij=1,4
                          if(ellev(tetra(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       if(sigtn==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetra(ij,ii))<0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg) =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg) =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg) =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)
                             endif

                          end do

                       else if(sigtp==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetra(ij,ii))>=0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg)=  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)    ! it will repeat this asigment once for each point of the surface but this creates no problem

                             endif
                          end do

                       else if(sigtn==2) then

                          !  research of two interface triangles
                          do ij=1,4
                             if(ellev(tetra(ij,ii))<0.0_rp) then
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

                          compg=compg+1
                          compt=compt+1

                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,1)
                          coord_tmp(2,compg)=  inter(2,1)
                          coord_tmp(3,compg)=  inter(3,1) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3)

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem) 

                          compg=compg+1
                          compt=compt+1
                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,4)
                          coord_tmp(2,compg)=  inter(2,4)
                          coord_tmp(3,compg)=  inter(3,4) 

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                       endif


                    end do

                 endif

                 ! Case Prism
              else if(pnode==6_ip) then

                 signn=0_ip
                 signp=0_ip
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

                 if(signp/=6_ip.and.signn/=6_ip) then

                    do ii=1,3
                       sigtn=0_ip
                       sigtp=0_ip
                       compl=0_ip
                       inod1=0_ip

                       do ij=1,4
                          if(ellev(tetrp(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       if(sigtn==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetrp(ij,ii))<0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg) =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg) =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg) =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                             endif

                          end do

                       else if(sigtp==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetrp(ij,ii))>=0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg)=  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                             endif
                          end do

                       else if(sigtn==2) then

                          !  research of two interface triangles
                          do ij=1,4
                             if(ellev(tetrp(ij,ii))<0.0_rp) then
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
                                l1=abs(ellev(tetrp(inod1,ii)))
                                lp=abs(ellev(tetrp(inod1,ii))-ellev(tetrp(ij,ii)))
                                x1=elcod(1,tetrp(inod1,ii))
                                x2=elcod(1,tetrp(ij,ii))
                                y1=elcod(2,tetrp(inod1,ii))
                                y2=elcod(2,tetrp(ij,ii))
                                z1=elcod(3,tetrp(inod1,ii))
                                z2=elcod(3,tetrp(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                                compl=compl+1
                                l1=abs(ellev(tetrp(inod2,ii)))
                                lp=abs(ellev(tetrp(inod2,ii))-ellev(tetrp(ij,ii)))
                                x1=elcod(1,tetrp(inod2,ii))
                                x2=elcod(1,tetrp(ij,ii))
                                y1=elcod(2,tetrp(inod2,ii))
                                y2=elcod(2,tetrp(ij,ii))
                                z1=elcod(3,tetrp(inod2,ii))
                                z2=elcod(3,tetrp(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                             endif

                          end do

                          compg=compg+1
                          compt=compt+1

                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,1)
                          coord_tmp(2,compg)=  inter(2,1)
                          coord_tmp(3,compg)=  inter(3,1) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3) 

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                          compg=compg+1
                          compt=compt+1
                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,4)
                          coord_tmp(2,compg)=  inter(2,4)
                          coord_tmp(3,compg)=  inter(3,4) 

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                       endif


                    end do

                 endif

                 ! Case Pyramid
              else if(pnode==5_ip) then

                 signn=0_ip
                 signp=0_ip
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

                 if(signp/=5_ip.and.signn/=5_ip) then

                    do ii=1,2
                       sigtn=0_ip
                       sigtp=0_ip
                       compl=0_ip
                       inod1=0_ip

                       do ij=1,4
                          if(ellev(tetpy(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       if(sigtn==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetpy(ij,ii))<0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg) =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg) =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg) =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                             endif

                          end do

                       else if(sigtp==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(tetpy(ij,ii))>=0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
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

                                lnodb_tmp(compl,compt)= compg
                                coord_tmp(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                coord_tmp(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                coord_tmp(3,compg)=  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                                if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                             endif
                          end do

                       else if(sigtn==2) then

                          !  research of two interface triangles
                          do ij=1,4
                             if(ellev(tetpy(ij,ii))<0.0_rp) then
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
                                l1=abs(ellev(tetpy(inod1,ii)))
                                lp=abs(ellev(tetpy(inod1,ii))-ellev(tetpy(ij,ii)))
                                x1=elcod(1,tetpy(inod1,ii))
                                x2=elcod(1,tetpy(ij,ii))
                                y1=elcod(2,tetpy(inod1,ii))
                                y2=elcod(2,tetpy(ij,ii))
                                z1=elcod(3,tetpy(inod1,ii))
                                z2=elcod(3,tetpy(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                                compl=compl+1
                                l1=abs(ellev(tetpy(inod2,ii)))
                                lp=abs(ellev(tetpy(inod2,ii))-ellev(tetpy(ij,ii)))
                                x1=elcod(1,tetpy(inod2,ii))
                                x2=elcod(1,tetpy(ij,ii))
                                y1=elcod(2,tetpy(inod2,ii))
                                y2=elcod(2,tetpy(ij,ii))
                                z1=elcod(3,tetpy(inod2,ii))
                                z2=elcod(3,tetpy(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                             endif

                          end do

                          compg=compg+1
                          compt=compt+1

                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,1)
                          coord_tmp(2,compg)=  inter(2,1)
                          coord_tmp(3,compg)=  inter(3,1) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3) 

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                          compg=compg+1
                          compt=compt+1
                          lnodb_tmp(1,compt)= compg
                          coord_tmp(1,compg)=  inter(1,2)
                          coord_tmp(2,compg)=  inter(2,2)
                          coord_tmp(3,compg)=  inter(3,2) 

                          compg=compg+1
                          lnodb_tmp(2,compt)= compg
                          coord_tmp(1,compg)=  inter(1,3)
                          coord_tmp(2,compg)=  inter(2,3)
                          coord_tmp(3,compg)=  inter(3,3) 

                          compg=compg+1
                          lnodb_tmp(3,compt)= compg
                          coord_tmp(1,compg)=  inter(1,4)
                          coord_tmp(2,compg)=  inter(2,4)
                          coord_tmp(3,compg)=  inter(3,4) 

                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                       endif


                    end do

                 endif


                 ! Case P1   
              else if(pnode==4_ip) then

                 signn=0_ip
                 signp=0_ip
                 compl=0_ip
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

                 if(signn==1) then

                    !  research of one interface triangle
                    do ij=1,4
                       if(ellev(ij)<0.0_rp) then
                          inod1=ij
                       endif
                    end do

                    do ij=1,4
                       if(ij/=inod1) then
                          compg=compg+1
                          compl=compl+1
                          if(compl==1) then
                             compt=compt+1
                          endif
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

                          lnodb_tmp(compl,compt)= compg
                          coord_tmp(1,compg)=  x1*(1-l1/lp)+x2*l1/lp  
                          coord_tmp(2,compg)=  y1*(1-l1/lp)+y2*l1/lp 
                          coord_tmp(3,compg)=  z1*(1-l1/lp)+z2*l1/lp 
                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                       endif

                    end do

                 else if(signp==1) then
                    !  research of one interface triangle
                    do ij=1,4
                       if(ellev(ij)>=0.0_rp) then
                          inod1=ij
                       endif
                    end do

                    do ij=1,4
                       if(ij/=inod1) then
                          compg=compg+1
                          compl=compl+1
                          if(compl==1) then
                             compt=compt+1
                          endif
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

                          lnodb_tmp(compl,compt)= compg
                          coord_tmp(1,compg)=  x1*(1-l1/lp)+x2*l1/lp  
                          coord_tmp(2,compg)=  y1*(1-l1/lp)+y2*l1/lp 
                          coord_tmp(3,compg)=  z1*(1-l1/lp)+z2*l1/lp 
                          if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

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

                          compl=compl+1
                          l1=abs(ellev(inod1))
                          lp=abs(ellev(inod1)-ellev(ij))
                          x1=elcod(1,inod1)
                          x2=elcod(1,ij)
                          y1=elcod(2,inod1)
                          y2=elcod(2,ij)
                          z1=elcod(3,inod1)
                          z2=elcod(3,ij)

                          inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                          inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                          inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                          compl=compl+1
                          l1=abs(ellev(inod2))
                          lp=abs(ellev(inod2)-ellev(ij))
                          x1=elcod(1,inod2)
                          x2=elcod(1,ij)
                          y1=elcod(2,inod2)
                          y2=elcod(2,ij)
                          z1=elcod(3,inod2)
                          z2=elcod(3,ij)

                          inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                          inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                          inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                       endif

                    end do

                    compg=compg+1
                    compt=compt+1

                    lnodb_tmp(1,compt)= compg
                    coord_tmp(1,compg)=  inter(1,1)
                    coord_tmp(2,compg)=  inter(2,1)
                    coord_tmp(3,compg)=  inter(3,1) 

                    compg=compg+1
                    lnodb_tmp(2,compt)= compg
                    coord_tmp(1,compg)=  inter(1,2)
                    coord_tmp(2,compg)=  inter(2,2)
                    coord_tmp(3,compg)=  inter(3,2) 

                    compg=compg+1
                    lnodb_tmp(3,compt)= compg
                    coord_tmp(1,compg)=  inter(1,3)
                    coord_tmp(2,compg)=  inter(2,3)
                    coord_tmp(3,compg)=  inter(3,3)

                    if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                    compg=compg+1
                    compt=compt+1
                    lnodb_tmp(1,compt)= compg
                    coord_tmp(1,compg)=  inter(1,2)
                    coord_tmp(2,compg)=  inter(2,2)
                    coord_tmp(3,compg)=  inter(3,2) 

                    compg=compg+1
                    lnodb_tmp(2,compt)= compg
                    coord_tmp(1,compg)=  inter(1,3)
                    coord_tmp(2,compg)=  inter(2,3)
                    coord_tmp(3,compg)=  inter(3,3) 

                    compg=compg+1
                    lnodb_tmp(3,compt)= compg
                    coord_tmp(1,compg)=  inter(1,4)
                    coord_tmp(2,compg)=  inter(2,4)
                    coord_tmp(3,compg)=  inter(3,4) 

                    if( tyred_lev > 2 ) lebsu_tmp(compt) = lnuew(ielem)

                 endif

              endif
              !
              ! Reorder interface to point always in the same direction
              !
              if( compt >= compi ) then
                 call elmder(&
                         pnode,ndime,elmar(pelty)%dercg,elcod,&
                         gpcar,gpdet,xjacm,xjaci)
                 gpgrl(1) = 0.0_rp
                 gpgrl(2) = 0.0_rp
                 gpgrl(3) = 0.0_rp
                 do inode = 1,pnode
                    gpgrl(1) = gpgrl(1) + gpcar(1,inode) * ellev(inode)
                    gpgrl(2) = gpgrl(2) + gpcar(2,inode) * ellev(inode)
                    gpgrl(3) = gpgrl(3) + gpcar(3,inode) * ellev(inode)                              
                 end do
              end if
              do ii = compi,compt
                 p1 = lnodb_tmp(1,ii)
                 p2 = lnodb_tmp(2,ii)
                 p3 = lnodb_tmp(3,ii)
                 call nortri(p1,p2,p3,coord_tmp,vec,ndime)
                 nx = vec(1,3)
                 ny = vec(2,3)
                 nz = vec(3,3)
                 xx = gpgrl(1) * nx + gpgrl(2) * ny + gpgrl(3) * nz
                 if( xx < 0.0_rp ) then
                    coord_xxx(1)    = coord_tmp(1,p2)
                    coord_xxx(2)    = coord_tmp(2,p2)
                    coord_xxx(3)    = coord_tmp(3,p2)
                    coord_tmp(1,p2) = coord_tmp(1,p3) 
                    coord_tmp(2,p2) = coord_tmp(2,p3) 
                    coord_tmp(3,p2) = coord_tmp(3,p3) 
                    coord_tmp(1,p3) = coord_xxx(1)    
                    coord_tmp(2,p3) = coord_xxx(2)    
                    coord_tmp(3,p3) = coord_xxx(3)    
                    !lnodb_tmp(1,ii) = p1
                    !lnodb_tmp(2,ii) = p3
                    !lnodb_tmp(3,ii) = p2
                 end if
              end do


           end if

        end do

     endif

  endif

  if( IPARALL ) then

     ! First
     ! master needs to know the discrete interface data 
     ! of each slave (connectivity and points coordinates vectors size) 
     !
     ! Master: NREDM_LEV
     !
     call lev_parall(6_ip)
     ! Second 
     ! sizes for the discrete interface on the whole domain
     ! has to be computed
     call lev_parall(7_ip)

     !   parallelization case
     !   each slave learns total discrete interface vectors size
     nboun_lev=nredt_lev(1)
     npoin_lev=nredt_lev(2)  

     !   parallelization case : allocation of total discrete interface vectors 
     !   lnodb_lev (connectivity) and coord_lev ( interface points coordinates)
     !
     ! Slave & Master:  LNODB_LEV, COORD_LEV
     !
     call lev_memall(3_ip)

     ! master gathers total interface points coordinates vector 
     call lev_parall(8_ip)

     ! master sends total interface points coordinates vector 
     ! to each slave
     call lev_parall(9_ip)

     ! master gathers vector that indicates each surface to which element belongs (lebsu)
     call lev_parall(10_ip)
     ! master sends vector that indicates each surface to which element belongs
     ! to each slave
     call lev_parall(11_ip)

  endif

!!$     if( ISLAVE ) then
  if( IPARALL ) then
     !   parallelization case
     !   lnodb_lev has to be written new
     compt = 1_ip
     do icomp=1,nboun_lev
        do idime=1,ndime
           lnodb_lev(idime,icomp) = compt
           compt = compt+1
        end do
     end do

  endif
  !
  ! Deallocate
  ! Master: LNODB_LEV, COORD_LEV, NREDM_LEV
  ! Slave:  LNODP_LEV, COORP_LEV
  !
  call lev_memall(7_ip)
  !
  ! Build psbel_lev & lsbel_lev
  ! psbel_lev has already been allocated in lev_memall(1) of fixed size nelwh
  !
  if(( tyred_lev == 3 ).or.( tyred_lev == 5 ).or.( tyred_lev == 6 )) then
     if( INOTMASTER ) then
        allocate(lauxi(nboun_lev),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LAUXI ','lev_conint',lauxi)
        allocate(lsbel_lev(nboun_lev),stat=istat)    ! deallocated in lev_redisc
        call memchk(zero,istat,mem_modul(1:2,modul),'LSBEL_LEV ','lev_conint',lsbel_lev)
        !
        ! count number of surfaces that belong to each element
        !
        do ielem=1,nelwh+1
           psbel_lev(ielem) = 0
        end do
        do ii = 1,nboun_lev
           ielem = lebsu_lev(ii)
           lauxi(ii) = psbel_lev(ielem+1)    ! in lauxi we have 0 for the first surface that corresponds to ielem, 1 for the second ...
           psbel_lev(ielem+1) = psbel_lev(ielem+1) + 1
        end do
        !
        ! obtain psbel_lev
        !
        psbel_lev(1) = 1
        do ii = 2,nelwh+1
           psbel_lev(ii) = psbel_lev(ii-1) + psbel_lev(ii)
        end do
        !
        ! obtain lsbel_lev
        !
        do ii = 1,nboun_lev
           ielem = lebsu_lev(ii)
           lsbel_lev( psbel_lev(ielem) + lauxi(ii) ) = ii
        end do

        call memchk(two,istat,mem_modul(1:2,modul),'LAUXI','lev_conint',lauxi)
        deallocate(lauxi, stat=istat)
        if(istat/=0) call memerr(two,'LAUXI','lev_conint',0_ip)

     end if
  end if

  if( INOTMASTER ) then
     ! 
     ! Deallocate lnuew 
     !
     if( ISEQUEN ) then
        if( tyred_lev > 2 ) then
           call memchk(two,istat,memor_dom,'lnuew','lev_conint',lnuew)
           deallocate(lnuew,stat=istat)
           if(istat/=0) call memerr(two,'lnuew','lev_conint',0_ip)
        end if
     end if
  end if

end subroutine lev_conint
