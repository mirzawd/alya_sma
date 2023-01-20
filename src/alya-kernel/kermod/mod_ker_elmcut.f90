!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ker_elmcut

  use def_kintyp, only : ip,rp,lg
  use def_elmtyp, only : ELCUT,ELFEM, TET04, HEX08, PYR05
  use def_master, only : fleve,INOTMASTER
  use def_domain, only : nelem,lnnod,lnods
  use def_domain, only : coord,lelch,ndime
  use def_domain, only : ltype,lgaus,ngaus
  use mod_cutele, only : onecut
  use mod_cutele, only : cutel
  use mod_cutele, only : pgaus_sub
  
  private

  public :: ker_elmcut_free_surface

contains

  subroutine ker_elmcut_free_surface()
    implicit none
    integer(ip)  :: ielem,inode,ipoin,pnode,jnode
    integer(ip)  :: jpoin,izdom,knode,idime,iedge,jedge
    integer(ip)  :: pair_of_edges(2,12),touch(12)
    real(rp)     :: regresion_plane(3,3),inv_regre(3,3),b_regresion(3)
    real(rp)     :: deter_regre,A_coef,B_coef,C_coef
    real(rp)     :: deltaf(12),deltaf_max
    real(rp)     :: fi,fj,t,norma(3),norma1(3),norma2(3),zeror
    real(rp)     :: xx(3,12),xx_avg(3),pcoor(3),max_fleve,min_fleve
    logical(lg)  :: lplus,lminu

    if( INOTMASTER ) then
       !
       ! Detect cut elements
       !
       zeror = epsilon(1.0_rp)
       norma = 0.0_rp
       pair_of_edges = 0_ip
       do ielem = 1,nelem
          lelch(ielem) = ELFEM
          lgaus(ielem) = ngaus(abs(ltype(ielem)))
       end do
       do ielem = 1,nelem
          lplus = .false.
          lminu = .false.
          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( fleve(ipoin,1) > 0.0_rp ) then
                lplus = .true.
             else
                lminu= .true.
             end if
          end do
          if( lplus .and. lminu ) lelch(ielem) = ELCUT
          max_fleve = -100000000000.0_rp
          min_fleve = 100000000000.0_rp
          pnode = lnnod(ielem)
          do inode = 1,pnode
             max_fleve= max(max_fleve,fleve(lnods(inode,ielem),1))
             min_fleve= min(min_fleve,fleve(lnods(inode,ielem),1))
          end do
          if(max_fleve * min_fleve == 0.0_rp)lelch(ielem)= ELFEM 
       end do
       !
       ! Cut elements
       !
       do ielem = 1,nelem
          if( lelch(ielem) == ELCUT ) then
             if(ndime==2)then
                ! print*,'ielemmmmmmmmmmmmmm',ielem
                knode = 0
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   if( inode == pnode ) then
                      jnode = 1
                   else
                      jnode = inode + 1
                   end if
                   ipoin = lnods(inode,ielem)
                   jpoin = lnods(jnode,ielem)
                   fi    = fleve(ipoin,1)
                   fj    = fleve(jpoin,1)
                   t     = -fi / ( fj - fi )
                   !if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                   !if( t >= 0.0_rp .and. t < 1.0_rp ) then
                   !if( t > 0.0_rp .and. t < 1.0_rp ) then
                   do while ( knode < 2 .and. ( t >= 0.0_rp .and. t <= 1.0_rp ) )
                      knode = knode + 1
                      do idime = 1,ndime
                         xx(idime,knode) = coord(idime,ipoin) + t * ( coord(idime,jpoin) - coord(idime,ipoin) )
                      end do
                   end do
                   !end if
                end do
                !print*,'mod_ker_elmcut',ielem,knode,fleve(lnods(1,ielem),1),fleve(lnods(2,ielem),1),fleve(lnods(3,ielem),1),fleve(lnods(4,ielem),1),t
                !if( knode /= 2 ) then
                !   call runend('PROBLEM CUTTING ELEMENT')
                !else
                !!if( knode == 0 .or. knode == 1 ) then
                !!   lelch(ielem) = ELFEM
                !!else
                norma(1) = - ( xx(2,2) - xx(2,1) ) 
                norma(2) =   ( xx(1,2) - xx(1,1) ) 
                norma    = norma / (sqrt(dot_product(norma,norma))+zeror)
                
                pcoor    = 0.0_rp
                do inode = 1,knode
                   pcoor(1:ndime) = pcoor(1:ndime) + xx(1:ndime,inode)
                end do
                pcoor(1:ndime) = pcoor(1:ndime) / real(knode,rp)
                call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub
                !!end if
             else
                if(ltype(ielem)==TET04)then
                   pair_of_edges(1,1) = 1
                   pair_of_edges(2,1) = 2
                   pair_of_edges(1,2) = 1
                   pair_of_edges(2,2) = 3
                   pair_of_edges(1,3) = 1
                   pair_of_edges(2,3) = 4
                   pair_of_edges(1,4) = 2
                   pair_of_edges(2,4) = 3
                   pair_of_edges(1,5) = 2
                   pair_of_edges(2,5) = 4
                   pair_of_edges(1,6) = 3
                   pair_of_edges(2,6) = 4
                   knode = 0
                   do iedge=1,6
                      ipoin = lnods(pair_of_edges(1,iedge),ielem)
                      jpoin = lnods(pair_of_edges(2,iedge),ielem)
                      fi    = fleve(ipoin,1)
                      fj    = fleve(jpoin,1)
                      t     = -fi / ( fj - fi )
                      if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                         knode = knode + 1
                         do idime = 1,ndime
                            xx(idime,knode) = coord(idime,ipoin) + t * ( coord(idime,jpoin) - coord(idime,ipoin) )
                         end do
                      end if
                   end do
                   pcoor    = 0.0_rp
                   do inode = 1,knode
                      pcoor(1:ndime) = pcoor(1:ndime) + xx(1:ndime,inode)
                   end do
                   pcoor(1:ndime) = pcoor(1:ndime) / real(knode,rp)

                   if (knode==3)then
                      norma(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma    = norma / (sqrt(dot_product(norma,norma))+zeror)
                   else if (knode==4)then
                      norma1(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma1(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma1(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma1    = norma1 / (sqrt(dot_product(norma1,norma1))+zeror)
                      norma2(1)= (xx(2,3)-xx(2,1))*(xx(3,4)-xx(3,1))- (xx(2,4)-xx(2,1))*(xx(3,3)-xx(3,1))
                      norma2(2)= (xx(3,3)-xx(3,1))*(xx(1,4)-xx(1,1))- (xx(1,3)-xx(1,1))*(xx(3,4)-xx(3,1))
                      norma2(3)= (xx(1,3)-xx(1,1))*(xx(2,4)-xx(2,1))- (xx(2,3)-xx(2,1))*(xx(1,4)-xx(1,1))
                      norma1    = norma2 / (sqrt(dot_product(norma2,norma2))+zeror)
                      norma = (norma1 + norma2)/2.0_rp
                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')

                      !call runend('PROGRAMARLO-MOD_KER_ELMCUT')
                      !print*,'programarlo'
                   else
                      call runend('PROBLEM CUTTING ELEMENT')
                   end if
                   call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                   lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub
                else if (ltype(ielem)==HEX08)then
                   !print*,'fleve',fleve(lnods(1:8,ielem),1)
                   touch = 0_ip
                   pair_of_edges(1,1) = 1
                   pair_of_edges(2,1) = 2
                   pair_of_edges(1,2) = 1
                   pair_of_edges(2,2) = 4
                   pair_of_edges(1,3) = 1
                   pair_of_edges(2,3) = 5
                   pair_of_edges(1,4) = 2
                   pair_of_edges(2,4) = 3
                   pair_of_edges(1,5) = 2
                   pair_of_edges(2,5) = 6
                   pair_of_edges(1,6) = 3
                   pair_of_edges(2,6) = 4
                   pair_of_edges(1,7) = 3
                   pair_of_edges(2,7) = 7
                   pair_of_edges(1,8) = 4
                   pair_of_edges(2,8) = 8
                   pair_of_edges(1,9) = 5
                   pair_of_edges(2,9) = 6
                   pair_of_edges(1,10) = 5
                   pair_of_edges(2,10) = 8
                   pair_of_edges(1,11) = 6
                   pair_of_edges(2,11) = 7
                   pair_of_edges(1,12) = 7
                   pair_of_edges(2,12) = 8
                   knode = 0
                   do iedge=1,12
                      ipoin = lnods(pair_of_edges(1,iedge),ielem)
                      jpoin = lnods(pair_of_edges(2,iedge),ielem)
                      fi    = fleve(ipoin,1)
                      fj    = fleve(jpoin,1)
                      t     = -fi / ( fj - fi )
                      if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                         knode = knode + 1
                         do idime = 1,ndime
                            xx(idime,knode) = coord(idime,ipoin) + t * ( coord(idime,jpoin) - coord(idime,ipoin) )
                         end do
                         !print*,'iedge-t',iedge,t,knode,xx(:,knode)
                      end if
                      deltaf(iedge) = abs(fi-fj)
                      touch(iedge) = 1
                   end do
                   do iedge=1,11
                      do jedge = iedge+1,12 
                         if(deltaf(iedge)>deltaf(jedge))then
                            deltaf_max = deltaf(iedge)
                            deltaf(iedge) = deltaf(jedge)
                            deltaf(jedge) = deltaf_max
                         end if
                      end do
                   end do
                   !print*,'delta-ordenado',deltaf
                   pcoor    = 0.0_rp
                   do inode = 1,knode
                      pcoor(1:ndime) = pcoor(1:ndime) + xx(1:ndime,inode)
                   end do
                   pcoor(1:ndime) = pcoor(1:ndime) / real(knode,rp)

                   if (knode==3)then
                      norma(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma    = norma / (sqrt(dot_product(norma,norma))+zeror)

!!$                      print*,'knode 3',knode,ielem, 'norma-ok',norma
!!$                      !regresion_plane(1,1)= knode
!!$                      regresion_plane = 0.0_rp
!!$                      b_regresion = 0.0_rp
!!$                      A_coef = 0.0_rp
!!$                      B_coef = 0.0_rp
!!$                      C_coef = 0.0_rp
!!$                     do inode=1,knode
!!$                         regresion_plane(1,1) = regresion_plane(1,1) + xx(1,inode)* xx(1,inode)
!!$                         regresion_plane(1,2) = regresion_plane(1,2) + xx(1,inode)* xx(2,inode)
!!$                         regresion_plane(1,3) = regresion_plane(1,3) + xx(1,inode)
!!$                         regresion_plane(2,1) = regresion_plane(2,1) + xx(1,inode) * xx(2,inode) 
!!$                         regresion_plane(2,2) = regresion_plane(2,2) + xx(2,inode) * xx(2,inode)                
!!$                         regresion_plane(2,3) = regresion_plane(2,3) + xx(2,inode) 
!!$                         regresion_plane(3,1) = regresion_plane(3,1) + xx(1,inode)
!!$                         regresion_plane(3,2) = regresion_plane(3,2) + xx(2,inode) 
!!$                         regresion_plane(3,3) = regresion_plane(3,3) + 1
!!$                         b_regresion(1) = b_regresion(1) + xx(1,inode) * xx(3,inode)
!!$                         b_regresion(2) = b_regresion(2) + xx(2,inode) * xx(3,inode)
!!$                         b_regresion(3) = b_regresion(3) + xx(3,inode)
!!$                      end do
!!$                      call invmtx(regresion_plane,inv_regre,deter_regre,3)
!!$                      print*,inv_regre
!!$                      !
!!$                      !regresion plane : Z = AX + BY + CZ
!!$                      !
!!$                      do izdom = 1,3
!!$                         A_coef = A_coef + inv_regre(1,izdom) * b_regresion(izdom)
!!$                         B_coef = B_coef + inv_regre(2,izdom) * b_regresion(izdom)
!!$                         C_coef = C_coef + inv_regre(3,izdom) * b_regresion(izdom)
!!$                      end do
!!$                      norma(1) = A_coef
!!$                      norma(2) = B_coef
!!$                      norma(3) = -1.0_rp
!!$                      pcoor = 0.0_rp
!!$                      do inode = 1,knode
!!$                         pcoor(1:ndime-1) = pcoor(1:ndime-1) + xx(1:ndime-1,inode)
!!$                      end do
!!$                      pcoor(1:ndime-1) = pcoor(1:ndime-1) / real(knode,rp)                         
!!$                      pcoor(3) = C_coef+A_coef*pcoor(1)+B_coef*pcoor(2)
!!$
!!$                      print*,'normal-punto', norma, pcoor
!!$                      stop
                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub


                   else if (knode==4)then
                      norma1(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma1(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma1(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma1    = norma1 / (sqrt(dot_product(norma1,norma1))+zeror)
                      norma2(1)= (xx(2,3)-xx(2,1))*(xx(3,4)-xx(3,1))- (xx(2,4)-xx(2,1))*(xx(3,3)-xx(3,1))
                      norma2(2)= (xx(3,3)-xx(3,1))*(xx(1,4)-xx(1,1))- (xx(1,3)-xx(1,1))*(xx(3,4)-xx(3,1))
                      norma2(3)= (xx(1,3)-xx(1,1))*(xx(2,4)-xx(2,1))- (xx(2,3)-xx(2,1))*(xx(1,4)-xx(1,1))
                      norma1    = norma2 / (sqrt(dot_product(norma2,norma2))+zeror)
                      norma = (norma1 + norma2)/2.0_rp
                      
!!$                      print*,'knode 4',knode,ielem
!!$                      print*,'xx-1',xx(:,1)
!!$                      print*,'xx-2',xx(:,2)
!!$                      print*,'xx-3',xx(:,3)
!!$                      print*,'xx-4',xx(:,4)
!!$                      print*,'norma-ok',norma
!!$                      print*,'pcoor-ok',pcoor
!!$                      !regresion_plane(1,1)= knode
!!$                      regresion_plane = 0.0_rp
!!$                      b_regresion = 0.0_rp
!!$                      A_coef = 0.0_rp
!!$                      B_coef = 0.0_rp
!!$                      C_coef = 0.0_rp
!!$                    
!!$                      do inode=1,knode
!!$                         regresion_plane(1,1) = regresion_plane(1,1) + xx(1,inode)* xx(1,inode)
!!$                         regresion_plane(1,2) = regresion_plane(1,2) + xx(1,inode)* xx(2,inode)
!!$                         regresion_plane(1,3) = regresion_plane(1,3) + xx(1,inode)
!!$                         regresion_plane(2,1) = regresion_plane(2,1) + xx(1,inode) * xx(2,inode) 
!!$                         regresion_plane(2,2) = regresion_plane(2,2) + xx(2,inode) * xx(2,inode)                
!!$                         regresion_plane(2,3) = regresion_plane(2,3) + xx(2,inode) 
!!$                         regresion_plane(3,1) = regresion_plane(3,1) + xx(1,inode)
!!$                         regresion_plane(3,2) = regresion_plane(3,2) + xx(2,inode) 
!!$                         regresion_plane(3,3) = regresion_plane(3,3) + 1
!!$                         b_regresion(1) = b_regresion(1) + xx(1,inode) * xx(3,inode)
!!$                         b_regresion(2) = b_regresion(2) + xx(2,inode) * xx(3,inode)
!!$                         b_regresion(3) = b_regresion(3) + xx(3,inode)
!!$                      end do
!!$                      print*,'R',regresion_plane
!!$                      print*,'b',b_regresion
!!$                     call invmtx(regresion_plane,inv_regre,deter_regre,3)
!!$                     print*,inv_regre
!!$                      !
!!$                      !regresion plane : Z = AX + BY + CZ
!!$                      !
!!$                      do izdom = 1,3
!!$                         A_coef = A_coef + inv_regre(1,izdom) * b_regresion(izdom)
!!$                         B_coef = B_coef + inv_regre(2,izdom) * b_regresion(izdom)
!!$                         C_coef = C_coef + inv_regre(3,izdom) * b_regresion(izdom)
!!$                      end do
!!$                      norma(1) = A_coef
!!$                      norma(2) = B_coef
!!$                      norma(3) = -1.0_rp
!!$                      pcoor = 0.0_rp
!!$                      do inode = 1,knode
!!$                         pcoor(1:ndime-1) = pcoor(1:ndime-1) + xx(1:ndime-1,inode)
!!$                      end do
!!$                      pcoor(1:ndime-1) = pcoor(1:ndime-1) / real(knode,rp)                         
!!$                      pcoor(3) = C_coef+A_coef*pcoor(1)+B_coef*pcoor(2)
!!$
!!$                      print*,'PLANO', A_coef,B_coef,C_coef
!!$                      print*,'normal-punto', norma, pcoor
!!$                      stop

                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub

                   else
                      regresion_plane = 0.0_rp
                      b_regresion = 0.0_rp
                      A_coef = 0.0_rp
                      B_coef = 0.0_rp
                      C_coef = 0.0_rp
                      xx_avg = 0.0_rp

                      do inode=1,knode
                         regresion_plane(1,1) = regresion_plane(1,1) + (xx(1,inode)-xx_avg(1))* (xx(1,inode)-xx_avg(1))
                         regresion_plane(1,2) = regresion_plane(1,2) + (xx(1,inode)-xx_avg(1))* (xx(2,inode)-xx_avg(2))
                         regresion_plane(1,3) = regresion_plane(1,3) + (xx(1,inode)-xx_avg(1))
                         regresion_plane(2,1) = regresion_plane(2,1) + (xx(1,inode)-xx_avg(1)) * (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(2,2) = regresion_plane(2,2) + (xx(2,inode)-xx_avg(2)) * (xx(2,inode)-xx_avg(2))                
                         regresion_plane(2,3) = regresion_plane(2,3) + (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(3,1) = regresion_plane(3,1) + (xx(1,inode)-xx_avg(1))
                         regresion_plane(3,2) = regresion_plane(3,2) + (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(3,3) = regresion_plane(3,3) + 1
                         b_regresion(1) = b_regresion(1) + (xx(1,inode)-xx_avg(1)) * (xx(3,inode)-xx_avg(3))
                         b_regresion(2) = b_regresion(2) + (xx(2,inode)-xx_avg(2)) * (xx(3,inode)-xx_avg(3))
                         b_regresion(3) = b_regresion(3) + (xx(3,inode)-xx_avg(3))
                      end do
                      call invmtx(regresion_plane,inv_regre,deter_regre,3)
                      !print*,inv_regre
                      !
                      !regresion plane : Z = AX + BY + CZ
                      !
                      do izdom = 1,3
                         A_coef = A_coef + inv_regre(1,izdom) * b_regresion(izdom)
                         B_coef = B_coef + inv_regre(2,izdom) * b_regresion(izdom)
                         C_coef = C_coef + inv_regre(3,izdom) * b_regresion(izdom)
                      end do
                      norma(1) = A_coef
                      norma(2) = B_coef
                      norma(3) = -1.0_rp
                      pcoor = 0.0_rp
                      do inode = 1,knode
                         pcoor(1:ndime-1) = pcoor(1:ndime-1) + xx(1:ndime-1,inode)
                      end do
                      pcoor(1:ndime-1) = pcoor(1:ndime-1) / real(knode,rp)                         
                      pcoor(3) = C_coef+A_coef*pcoor(1)+B_coef*pcoor(2)

                      !print*,'normal-punto', norma, pcoor
                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub
                      !call runend('PROBLEM CUTTING ELEMENT')
                   end if
                else if (ltype(ielem)==PYR05)then
                   touch = 0_ip
                   pair_of_edges(1,1) = 1
                   pair_of_edges(2,1) = 2
                   pair_of_edges(1,2) = 2
                   pair_of_edges(2,2) = 3
                   pair_of_edges(1,3) = 3
                   pair_of_edges(2,3) = 4
                   pair_of_edges(1,4) = 4
                   pair_of_edges(2,4) = 1
                   pair_of_edges(1,5) = 1
                   pair_of_edges(2,5) = 5
                   pair_of_edges(1,6) = 2
                   pair_of_edges(2,6) = 5
                   pair_of_edges(1,7) = 3
                   pair_of_edges(2,7) = 5
                   pair_of_edges(1,8) = 4
                   pair_of_edges(2,8) = 5
                    knode = 0
                   do iedge=1,8
                      ipoin = lnods(pair_of_edges(1,iedge),ielem)
                      jpoin = lnods(pair_of_edges(2,iedge),ielem)
                      fi    = fleve(ipoin,1)
                      fj    = fleve(jpoin,1)
                      t     = -fi / ( fj - fi )
                      if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                         knode = knode + 1
                         do idime = 1,ndime
                            xx(idime,knode) = coord(idime,ipoin) + t * ( coord(idime,jpoin) - coord(idime,ipoin) )
                         end do
                         !print*,'iedge-t',iedge,t,knode,xx(:,knode)
                      end if
                      deltaf(iedge) = abs(fi-fj)
                      touch(iedge) = 1
                   end do
                   do iedge=1,7
                      do jedge = iedge+1,8 
                         if(deltaf(iedge)>deltaf(jedge))then
                            deltaf_max = deltaf(iedge)
                            deltaf(iedge) = deltaf(jedge)
                            deltaf(jedge) = deltaf_max
                         end if
                      end do
                   end do
                   pcoor    = 0.0_rp
                   do inode = 1,knode
                      pcoor(1:ndime) = pcoor(1:ndime) + xx(1:ndime,inode)
                   end do
                   pcoor(1:ndime) = pcoor(1:ndime) / real(knode,rp)

                   if (knode==3)then
                      norma(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma    = norma / (sqrt(dot_product(norma,norma))+zeror)


                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub


                   else if (knode==4)then
                      norma1(1)= (xx(2,2)-xx(2,1))*(xx(3,3)-xx(3,1))- (xx(2,3)-xx(2,1))*(xx(3,2)-xx(3,1))
                      norma1(2)= (xx(3,2)-xx(3,1))*(xx(1,3)-xx(1,1))- (xx(1,2)-xx(1,1))*(xx(3,3)-xx(3,1))
                      norma1(3)= (xx(1,2)-xx(1,1))*(xx(2,3)-xx(2,1))- (xx(2,2)-xx(2,1))*(xx(1,3)-xx(1,1))
                      norma1    = norma1 / (sqrt(dot_product(norma1,norma1))+zeror)
                      norma2(1)= (xx(2,3)-xx(2,1))*(xx(3,4)-xx(3,1))- (xx(2,4)-xx(2,1))*(xx(3,3)-xx(3,1))
                      norma2(2)= (xx(3,3)-xx(3,1))*(xx(1,4)-xx(1,1))- (xx(1,3)-xx(1,1))*(xx(3,4)-xx(3,1))
                      norma2(3)= (xx(1,3)-xx(1,1))*(xx(2,4)-xx(2,1))- (xx(2,3)-xx(2,1))*(xx(1,4)-xx(1,1))
                      norma1    = norma2 / (sqrt(dot_product(norma2,norma2))+zeror)
                      norma = (norma1 + norma2)/2.0_rp
                      

                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub

                   else
                      regresion_plane = 0.0_rp
                      b_regresion = 0.0_rp
                      A_coef = 0.0_rp
                      B_coef = 0.0_rp
                      C_coef = 0.0_rp
                      xx_avg = 0.0_rp

                      do inode=1,knode
                         regresion_plane(1,1) = regresion_plane(1,1) + (xx(1,inode)-xx_avg(1))* (xx(1,inode)-xx_avg(1))
                         regresion_plane(1,2) = regresion_plane(1,2) + (xx(1,inode)-xx_avg(1))* (xx(2,inode)-xx_avg(2))
                         regresion_plane(1,3) = regresion_plane(1,3) + (xx(1,inode)-xx_avg(1))
                         regresion_plane(2,1) = regresion_plane(2,1) + (xx(1,inode)-xx_avg(1)) * (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(2,2) = regresion_plane(2,2) + (xx(2,inode)-xx_avg(2)) * (xx(2,inode)-xx_avg(2))                
                         regresion_plane(2,3) = regresion_plane(2,3) + (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(3,1) = regresion_plane(3,1) + (xx(1,inode)-xx_avg(1))
                         regresion_plane(3,2) = regresion_plane(3,2) + (xx(2,inode)-xx_avg(2)) 
                         regresion_plane(3,3) = regresion_plane(3,3) + 1
                         b_regresion(1) = b_regresion(1) + (xx(1,inode)-xx_avg(1)) * (xx(3,inode)-xx_avg(3))
                         b_regresion(2) = b_regresion(2) + (xx(2,inode)-xx_avg(2)) * (xx(3,inode)-xx_avg(3))
                         b_regresion(3) = b_regresion(3) + (xx(3,inode)-xx_avg(3))
                      end do
                      call invmtx(regresion_plane,inv_regre,deter_regre,3)
                      !print*,inv_regre
                      !
                      !regresion plane : Z = AX + BY + CZ
                      !
                      do izdom = 1,3
                         A_coef = A_coef + inv_regre(1,izdom) * b_regresion(izdom)
                         B_coef = B_coef + inv_regre(2,izdom) * b_regresion(izdom)
                         C_coef = C_coef + inv_regre(3,izdom) * b_regresion(izdom)
                      end do
                      norma(1) = A_coef
                      norma(2) = B_coef
                      norma(3) = -1.0_rp
                      pcoor = 0.0_rp
                      do inode = 1,knode
                         pcoor(1:ndime-1) = pcoor(1:ndime-1) + xx(1:ndime-1,inode)
                      end do
                      pcoor(1:ndime-1) = pcoor(1:ndime-1) / real(knode,rp)                         
                      pcoor(3) = C_coef+A_coef*pcoor(1)+B_coef*pcoor(2)

                      !print*,'normal-punto', norma, pcoor
                      call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
                      lgaus(ielem) = size(cutel(ielem) % l) * pgaus_sub
                      !call runend('PROBLEM CUTTING ELEMENT')
                   end if
                end if
             end if
          end if
       end do
       !print*,'SALGO DE BUCLE DE ELEMENTOS-MOD_KER_ELMCUT'
    end if
  end subroutine ker_elmcut_free_surface

end module mod_ker_elmcut

