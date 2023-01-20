!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_projec()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_projec
  ! NAME 
  !    tur_projec
  ! DESCRIPTION
  !
  ! A value 'f/y' is wanted on the wall (where y=0)
  ! It is computed by setting f=0 on the wall and projecting its gradients:
  ! f/y = grad(f).n
  !
  ! 1. f = sqrt(k)       => |grad(f).n| = sqrt(k)/y
  ! 2. f = 1/y           => |grad(f).n| = 1/y^2
  ! 3. f = y^2*eps^{1/4} => |grad(f).n| = y*eps^{1/4}
  !
  ! USED BY
  !    tur_boucon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_turbul
  use mod_memory
  implicit none
  integer(ip)          :: ipoin,idime,inode,ielem,ibopo,jnode,knode
  integer(ip)          :: pnode,pelty
  real(rp)             :: detjm,gpvol,cartc(ndime,mnode) 
  real(rp)             :: elcod(ndime,mnode)
  real(rp)             :: elk12(mnode),elono(mnode),eleps(mnode)
  real(rp)             :: xjaci(9),xjacm(9),eps,rmass
  integer(ip), pointer :: lelex(:)
  real(rp),    pointer :: grk12_tmp(:)
  real(rp),    pointer :: grono_tmp(:)
  real(rp),    pointer :: greps_tmp(:)

  nullify(lelex)
  nullify(grk12_tmp)
  nullify(grono_tmp)
  nullify(greps_tmp)

  if ( (kfl_paral/=0) .and. ( kfl_grk12_tur/=0 .or. kfl_grono_tur/=0 .or. kfl_greps_tur/=0 ) ) then
     !
     ! LELEX: =1 if element has a node IPOIN such that LPOTY(IPOIN)>0
     !
     call memory_alloca(mem_modul(1:2,modul),'LELEX','tur_projec',lelex,nelem)

     do ielem=1,nelem
        inode=0
        pnode=nnode(ltype(ielem)) 
        do while(inode<pnode)
           inode=inode+1
           ipoin=lnods(inode,ielem)
           if(lpoty(ipoin)>0) then
              lelex(ielem)=1
              inode=pnode
           end if
        end do
     end do
     !
     ! Initialization
     !
     if(kfl_grk12_tur/=0) then
        call memory_alloca(mem_modul(1:2,modul),'GRK12_TMP','tur_projec',grk12_tmp,npoin)
        do ibopo=1,nbopo
           grk12_tur(ibopo)=0.0_rp
        end do
     end if
     if(kfl_grono_tur/=0) then
        call memory_alloca(mem_modul(1:2,modul),'GRONO_TMP','tur_projec',grono_tmp,npoin)
        do ibopo=1,nbopo
           grono_tur(ibopo)=0.0_rp
        end do
     end if
     if(kfl_greps_tur/=0) then
        call memory_alloca(mem_modul(1:2,modul),'GREPS_TMP','tur_projec',greps_tmp,npoin)
        do ibopo=1,nbopo
           greps_tur(ibopo)=0.0_rp
        end do
     end if
     !
     ! Project gradients
     !
     elements: do ielem=1,nelem

        if(lelex(ielem)==1) then
           pelty=ltype(ielem) 
           pnode=nnode(pelty)
           !           if( pelty == PYR05 ) call runend('TUR_PROJEC: NOT CODED FOR PYR05')
           if( pelty .ne. PYR05 ) then
              !
              ! Gather vectors
              !
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
                 if(lpoty(ipoin)==0) then
                    if(kfl_grk12_tur/=0) elk12(inode) = sqrt(untur(1,ipoin,1)) ! sqrt(k)/y
                    if(kfl_grono_tur/=0) then
                       if(walld(ipoin)/=0.0_rp) then
                          elono(inode) = 1.0_rp/walld(ipoin)               ! 1/y
                       else
                          elono(inode) = 0.0_rp
                       end if
                    end if
                    if(kfl_greps_tur/=0) then
                       eps = 0.09_rp*untur(1,ipoin,1)*untur(2,ipoin,1)         ! eps=Cmu*k*w
                       eleps(inode) = walld(ipoin)**2*eps**0.25_rp           ! y^2*eps^{1/4}
                    end if
                 else
                    elk12(inode) = 0.0_rp
                    elono(inode) = 0.0_rp
                    eleps(inode) = 0.0_rp
                 end if
              end do
              !
              ! Loop over Gauss points (which are nodes)
              !
              knode = pnode
              if( lelch(ielem) == ELEXT ) knode = 1
              gauss_points: do inode=1,knode
                 ipoin=lnods(inode,ielem) 
                 ibopo=lpoty(ipoin)
                 call elmder(&
                      pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                      elcod,cartc,detjm,xjacm,xjaci)
                 gpvol=elmar(pelty)%weigc(inode)*detjm
                 if(kfl_naxis==1) then
                    if(elcod(1,inode)==0.0_rp) then
                       gpvol=gpvol*twopi*1.0e-12_rp
                    else
                       gpvol=gpvol*twopi*elcod(1,inode)
                    end if
                 end if
                 !
                 ! Gradients
                 !
                 if(ibopo>0) then
                    if(kfl_grk12_tur/=0) then
                       do idime=1,ndime
                          do jnode=1,pnode
                             grk12_tmp(ipoin)=grk12_tmp(ipoin)+gpvol&
                                  *cartc(idime,jnode)*elk12(jnode)&
                                  *exnor(idime,1,ibopo)
                          end do
                       end do
                    end if
                    if(kfl_grono_tur/=0) then
                       do idime=1,ndime
                          do jnode=1,pnode
                             grono_tmp(ipoin)=grono_tmp(ipoin)+gpvol&
                                  *cartc(idime,jnode)*elono(jnode)&
                                  *exnor(idime,1,ibopo)
                          end do
                       end do
                    end if
                    if(kfl_greps_tur/=0) then
                       do idime=1,ndime
                          do jnode=1,pnode
                             greps_tmp(ipoin)=greps_tmp(ipoin)+gpvol&
                                  *cartc(idime,jnode)*eleps(jnode)&
                                  *exnor(idime,1,ibopo)
                          end do
                       end do
                    end if
                 end if
              end do gauss_points
           end if
        end if

     end do elements
     !
     ! Periodicity and Parall service
     !
     if(kfl_grk12_tur/=0) call rhsmod(1_ip,grk12_tmp)
     if(kfl_grono_tur/=0) call rhsmod(1_ip,grono_tmp)  
     if(kfl_greps_tur/=0) call rhsmod(1_ip,greps_tmp)  
     !
     ! Solve diagonal system
     !
     do ipoin=1,npoin
        ibopo=lpoty(ipoin)
        if(ibopo>0) then
           rmass=1.0_rp/vmass(ipoin)
           if(kfl_grk12_tur/=0) grk12_tur(ibopo)=abs(grk12_tmp(ipoin))*rmass
           if(kfl_grono_tur/=0) grono_tur(ibopo)=abs(grono_tmp(ipoin))*rmass
           if(kfl_greps_tur/=0) greps_tur(ibopo)=abs(greps_tmp(ipoin))*rmass
        end if
     end do

     call memory_deallo(mem_modul(1:2,modul),'LELEX'    ,'tur_projec',lelex    )
     call memory_deallo(mem_modul(1:2,modul),'GRK12_TMP','tur_projec',grk12_tmp)
     call memory_deallo(mem_modul(1:2,modul),'GRONO_TMP','tur_projec',grono_tmp)
     call memory_deallo(mem_modul(1:2,modul),'GREPS_TMP','tur_projec',greps_tmp)
 
  end if

end subroutine tur_projec

