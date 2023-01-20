!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_updtss()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_updtss
  ! NAME 
  !    tur_updtss
  ! DESCRIPTION
  !    This routine computes the time step size for the turbulence
  !    equations.
  ! USED BY
  !    tur_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_turbul
  use mod_communications, only : PAR_MIN
  use def_kermod,          only     :  kfl_adj_prob
  implicit none 
  integer(ip) :: ielem,idime,inode,ipoin,iunkn
  integer(ip) :: pnode,pelty,porde 
  real(rp)    :: dtcri,dtmin
  real(rp)    :: chale(2),hleng(3),tragl(9),chave(6)
  real(rp)    :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)    :: eltur(nturb_tur,mnode),eledd(mnode)
  real(rp)    :: elfle(mnode)

  if( kfl_timei_tur /= 0 ) then

     dtmin = 1e6_rp

     if( INOTMASTER ) then

        do ielem = 1,nelem
           pelty = ltype(ielem)

           if( lelch(ielem) /= ELHOL ) then
              pnode = nnode(pelty)
              porde = lorde(pelty)
              !
              ! ELCOD, ELVEL, ELTUR, ELEDD: Gather
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do iunkn = 1,nturb_tur
                    if (kfl_adj_prob == 0) then
		      eltur(iunkn,inode) = untur(iunkn,ipoin,1)
                    else
		      eltur(iunkn,inode) = untur_forw(iunkn,ipoin,1)
                    endif
                 end do
                 eledd(inode) = turmu(ipoin)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                    if (kfl_adj_prob == 0) then
                      elvel(idime,inode) = advec(idime,ipoin,1)
                    else
                      elvel(idime,inode) = veloc_forw(idime,ipoin,1)                    
                    endif
                 end do
              end do
              !
              ! ELFLE
              !
              if( kfl_colev_tur /= 0 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    elfle(inode) = fleve(ipoin,kfl_colev_tur)
                 end do
              end if
              !
              ! HLENG and TRAGL at center of gravity
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              !
              ! CHALE: Compute the characteristic length
              ! 
              call elmchl(&
                   tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                   porde,hnatu(pelty),kfl_advec_tur,kfl_ellen_tur)
              do iunkn_tur=1,min(nturb_tur,3_ip)
                 call tur_elmtss(&
                      pelty,pnode,ielem,elvel,eledd,dtcri,eltur,elfle,&
                      lnods(:,ielem),chale,hleng)
                 dtmin=min(dtmin,dtcri)
              end do
           end if
        end do
     end if
     !
     ! Look for minimum over whole mesh
     !
     call PAR_MIN(dtmin,'IN MY CODE')
    
     dtcri_tur = dtmin
     if(dtcri_tur/=0.0_rp) then
        dtinv_tur = 1.0_rp/(dtcri_tur*safet_tur)
     else
        dtinv_tur = 0.0_rp
     end if
     if(kfl_timco==1) dtinv=max(dtinv,dtinv_tur)
   
  end if

end subroutine tur_updtss
