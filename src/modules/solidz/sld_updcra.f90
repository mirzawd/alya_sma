!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updcra.f90
!> @author  Denny 
!> @date    16/11/1966
!> @brief   Updating the crack propagation
!> @details Updating the crack propagation
!> @} 
!-----------------------------------------------------------------------
subroutine sld_updcra()
  use def_parame
  use def_master
  use def_domain
  use def_solidz
  use mod_cutele
  use mod_communications,        only : PAR_MAX
  implicit none
  integer(ip)       :: nctip,ifacg,ielem,pmate,idime,ielty,ipoin,inode,jdime,pnode
!  integer(ip)       :: iface,pface
  logical(lg)       :: debugging
  real(rp)          :: tcrit,weigh,towei,sgsur
  real(rp)          :: sigma(ndime,ndime)
  
  debugging = .false.

  if (kfl_xfeme_sld > 0) then
     !
     !  Project stress from GP to nodes and exchange in parallel for nodes at the boundary of subdomains
     !
     call sld_pstres()
     !
     !  In case of propagation from initial/existing crack
     !
     nctip = 1_ip
     do while (nctip > 0_ip) 
        !
        !  Wind forward the crack (propagation) at the beginning of time step
        !
        if (INOTMASTER) then
           do ifacg = 1,nfacg
              do idime = 1,ndime
                 crtip_sld(idime,ifacg) = 0.0_rp
              end do
           end do
           call sld_elmope(6_ip)                               ! get the position and normal of new element crack
           call cutele(1_ip,cranx_sld,crapx_sld,lcrkf_sld)     ! define and perform element cut
           call sld_cratip(nctip)                              ! update the new crack tip position
        else
            nctip = 0_ip
        end if
        !
        !  Exchange crack tip information across subdomains in parallel
        !
        call PAR_MAX(nctip)
        ncrak_sld = ncrak_sld + nctip
        if (INOTMASTER .and. nctip > 0_ip) then
           call parari('SLX',NFACE_TYPE,nfacg,lcrkf_sld)
           do ifacg = 1,nfacg
              lcrkf_sld(ifacg) = min(1_ip,lcrkf_sld(ifacg))
           end do
!           call pararr('SLX',NFACE_TYPE,ndime*nfacg,cockf_sld)
           call pararr('SLX',NFACE_TYPE,ndime*nfacg,crtip_sld)
           do ifacg = 1,nfacg
              do idime = 1,ndime
                 cockf_sld(idime,ifacg) = cockf_sld(idime,ifacg) + crtip_sld(idime,ifacg)
              end do
           end do
!           do ielem = 1,nelem
!              ielty = ltype(ielem)
!              pface = nface(ielty)
!              do iface=1,pface
!                 ifacg = lelfa(ielem) % l(iface)
!                 if (lcrkf_sld(ifacg) == 1) then
!                    write(*,*)'subdomain :',kfl_paral
!                    write(*,*)'element :',ielem,leenr_sld(ielem),crapx_sld(:,ielem)
!                    write(*,*)'face :',iface,ifacg,cockf_sld(:,ifacg)
!                    write(*,*)'face :',iface,ifacg,crtip_sld(:,ifacg)
!                 end if
!              end do
!           end do
        end if
     end do
     !
     !  In case of crack nucleation due to maximum tensile stress reaching critical value
     !
     call PAR_MAX(ncrak_sld)
     if (ncrak_sld == 0) then
        if (INOTMASTER) then
           call sld_elmope(6_ip)                     ! get the position and normal of new element crack
        else
           sgult_sld = 0.0_rp
        end if
        call PAR_MAX(sgult_sld)
        if (INOTMASTER) then


    
           do ielem = 1,nelem
              pmate = 1
              if( nmate > 1 ) then
                 pmate = lmate_sld(ielem)
              end if
              tcrit = parch_sld(1,pmate)
              ielty = ltype(ielem)
              pnode = nnode(ielty)
              !
              ! Calculate sigma (weighted averaged Cauchy stress tensor)
              !
              sigma = 0.0_rp
              towei = 0.0_rp
              do inode = 1,pnode
                 weigh = 1.0_rp / real(pnode,rp) 
                 towei = towei + weigh
                 ipoin = lnods(inode,ielem)           
                 do idime = 1,ndime
                    do jdime = 1,ndime
                       sigma(idime,jdime) = sigma(idime,jdime) +  weigh * nopio_sld(jdime+(idime-1)*ndime,ipoin)
                    end do
                 end do
              end do

              !Von Misses surface
              if (ndime==2) then
                 sgsur = sqrt(0.5_rp*((sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))))
              else
                 sgsur = sqrt(0.5_rp*((sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))+(sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))+(sigma(3,3)-sigma(1,1))*(sigma(3,3)-sigma(1,1))))
              end if

              if (sgult_sld > 1.0_rp*tcrit .and. sgmax_sld(ielem) >= sgult_sld .and. sgsur > 0 .and. leenr_sld(ielem) == 0_ip) then
                 leenr_sld(ielem) = -1
                 if( kfl_cohes_sld /= 0 .and. kfl_elcoh > 0) lecoh_sld(ielem) = 1
                 lelch(ielem) = ELCUT
                 
                 if (debugging) write(*,'(a29)')             '--| ALYA     INITIATION CRACK'
                 if (debugging) write(*,'(a20,1x,i8,1x,a20)')'--| ALYA     ELEMENT',ielem,'EXCEEDS CRIT. STRESS'
                 if (debugging) write(*,'(a27,1x,e12.5)')    '--| ALYA     ELEMENT STRESS',sgmax_sld(ielem)
                 if (debugging) write(*,'(a30,3(1x,f12.7))') '--| ALYA     CRACK ORIENTATION',(cranx_sld(idime,ielem),idime=1,ndime)
                 if (debugging) write(*,'(a36,1x,f12.5)')    '--| ALYA     CRACK ORIENTATION ANGLE',&
                                atan(-cranx_sld(1,ielem)/cranx_sld(2,ielem))*45.0_rp/atan(1.0_rp)
                 if (debugging) write(*,'(a30,3(1x,f12.7))') '--| ALYA     CRACK POSITION   ',(crapx_sld(idime,ielem),idime=1,ndime)

              end if
           end do
           do ifacg = 1,nfacg
              do idime = 1,ndime
                 crtip_sld(idime,ifacg) = 0.0_rp
              end do
           end do
           call cutele(1_ip,cranx_sld,crapx_sld,lcrkf_sld)     ! define and perform element cut
           call sld_cratip(nctip)                              ! update the new crack tip position
           !
           !  Exchange crack tip information across subdomains in parallel
           !
           call parari('SLX',NFACE_TYPE,nfacg,lcrkf_sld)
           do ifacg = 1,nfacg
              lcrkf_sld(ifacg) = min(1_ip,lcrkf_sld(ifacg))
           end do
           call pararr('SLX',NFACE_TYPE,ndime*nfacg,crtip_sld)
           do ifacg = 1,nfacg
              do idime = 1,ndime
                 cockf_sld(idime,ifacg) = cockf_sld(idime,ifacg) + crtip_sld(idime,ifacg)
              end do
           end do
        end if
     end if
     
  end if

end subroutine sld_updcra

