!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_uptcoh.f90
!> @author  Eva Casoni
!> @date    11/06/2013
!> @brief   Updating traction in cohesive elements
!> @details Updating traction in cohesive elements
!> @} 
!-----------------------------------------------------------------------
subroutine sld_uptcoh()
  use def_parame
  use def_master
  use def_domain
  use def_solidz
  use def_elmtyp
  use def_domain  !EVA
  implicit none
  integer(ip)       :: ielem, ipoin, jpoin, kpoin,  ipoin1, ipoin2, idime, jdime, inode, iperi
  integer(ip)       :: pmate, pnode, pelty
  logical(lg)       :: debugging
  real(rp)          :: midcoord(ndime,mnode), norm(ndime), vec1(ndime), vec2(ndime)
  real(rp)          :: tract(ndime)
  real(rp)          :: tracnx(mnode), tractx(mnode)
  real(rp)          :: beta, eta, traceff, tcrit, modno
  
  debugging = .false.

  if (kfl_elcoh > 0) then

  pmate = 1
  if( nmate > 1 ) then
    pmate = lmate_sld(ielem)
  end if
  tcrit = parch_sld(1,pmate)

     if (INOTMASTER) then
     ! 
     ! Compute normal in the midplane of the cohesive element (same for all nodes of the face) 
     !
     do ielem=1,nelem
        if (lelch(ielem)==ELCOH) then
           pelty = ltype(ielem) 
           pnode = nnode(pelty)
           do inode =1, pnode/2
              ipoin1 = lnods(inode,ielem)
              ipoin2 = lnods(pnode/2+inode,ielem)
              midcoord(:,inode) = 0.5_rp*(coord(:,ipoin1)+coord(:,ipoin2))
           end do
           if (ndime==3) then
              vec1 = midcoord(:,1)-midcoord(:,2)
              vec2 = midcoord(:,2)-midcoord(:,3)
              norm(1) = abs(vec1(2)*vec2(3) - vec1(3)*vec2(2)) 
              norm(2) = -abs(vec1(1)*vec2(3)-vec1(3)*vec2(1))
              norm(3) = abs(vec1(1)*vec2(2)-vec1(2)*vec2(1))
              modno = sqrt(norm(1)**2+norm(2)**2+norm(3)**2)
              norm(:) = norm(:)/modno
           else
              call runend('sld_uptcoh: Cohesive elements not implemented for 2D')
           end if
           cohnx_sld(:,ielem) = norm
        end if
     end do

     !
     !  Project stress from GP to nodes and exchange in parallel for nodes at the boundary of subdomains (return frist Piola-Kirchoff)
     !
     call sld_pstres()

     !
     !  Compute surface traction t=sigma*n and the tangencial and normal components
     !
     do ielem=1,nelem
        if (lelch(ielem)==ELCOH) then
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           do inode=1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 tract(idime) = 0.0_rp
                 do jdime = 1,ndime
                    tract(idime) = tract(idime) + nopio_sld(jdime+(idime-1)*ndime,ipoin)*cohnx_sld(jdime,ielem) 
                 end do
              end do
              !print*, ipoin, nopio_sld(:,ipoin)
              ! Compute normal and tangetial components of the surface traction
              tracnx(inode) = 0.0_rp
              tractx(inode) = 0.0_rp
              do idime=1,ndime
                 tracnx(inode) = tracnx(inode) + cohnx_sld(idime,ielem)*tract(idime)
              end do
              do idime=1,ndime
                 tractx(inode) = tractx(inode) + tract(idime)*tract(idime)
              end do
              tractx(inode) = sqrt(tractx(inode) - tracnx(inode)*tracnx(inode))   
           end do
           ! 
           ! Check criteria of cohesive fracture            !
            pmate = 1
            if( nmate > 1 ) then
               pmate = lmate_sld(ielem)
            end if
           if (lawch_sld(pmate)==900) then
              beta = parch_sld(3,pmate)
           else if (lawch_sld(pmate)==901) then
              beta = parch_sld(4,pmate)
           else if (lawch_sld(pmate)==902) then
              beta = parch_sld(4,pmate)
           else if (lawch_sld(pmate)==903) then
              beta = parch_sld(5,pmate)
           end if
           do inode=1,pnode
              ipoin = lnods(inode,ielem)   
              if (tracnx(inode) >= 0.0_rp) then
                 traceff = sqrt(tracnx(inode)*tracnx(inode) + (beta**(-2))*(tractx(inode)*tractx(inode)))
              else 
                 if (kfl_cohft_sld == 0) then
                    traceff = 0.0_rp
                 else
                    call runend('sld_uptcoh: Not implemented for friction')
                 !
                 ! In case of considering friction 
                 !
                    if (abs(tractx(inode))-eta*abs(tracnx(inode)) < 0) then
                       traceff = 0
                    else
                       traceff = (abs(tractx(inode))-eta*abs(tracnx(inode)))/beta
                    end if
                 end if
              end if   
              treff_sld(ipoin) = traceff
              ! 
              ! Check fracture cohesive criteria
              !
              if (traceff > tcrit) then
                 nocoh_sld(ipoin) = 1
                 ! Break periodicity between the nodes of the cohesive element
                 do iperi = 1,nperi
                    jpoin = lperi(1,iperi)    !master
                    kpoin = lperi(2,iperi)    !slave
                    if (jpoin > 0 .and. kpoin > 0) then 
                       if (ipoin == jpoin .or. ipoin == kpoin) then
                          lperi(1,iperi) = -lperi(1,iperi)
                          lperi(2,iperi) = -lperi(2,iperi)
                          kfl_domar = 1
                       end if
                    end if 
                 end do
                 if( kfl_cohes_sld /= 0 .and. kfl_elcoh > 0 ) lecoh_sld(ielem) = 1
              end if
           end do
        end if
     end do

  end if
end if

end subroutine sld_uptcoh

