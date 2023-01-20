!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_heatRelease_lookup()
   !------------------------------------------------------------------------
   ! lookup and integration of heat release rate
   !------------------------------------------------------------------------
   use def_chemic,             only : kfl_hrr_fw_chm, kfl_hrr_col_chm, hrr_int_chm
   use def_master,             only : INOTEMPTY
   use def_master,             only : mem_modul,modul
   use def_kermod,             only : lookup_fw
   use def_domain,             only : ltype,lnods,ndime,nelem,ngaus,nnode,elmar
   use def_domain,             only : mnode,mgaus
   use def_domain,             only : coord
   use def_kintyp,             only : ip,rp,r2p
   use mod_memory_basic,       only : memory_alloca
   use mod_memory_basic,       only : memory_deallo
   use mod_communications,     only : PAR_SUM
   type(r2p),pointer         :: aux_r2p(:)
   integer(ip)               :: ielem,igaus,ipoin,inode,idime
   integer(ip)               :: pelty,pgaus,pnode
   real(rp)                  :: elcod(ndime,mnode)
   real(rp)                  :: gpdet, gpvol
   real(rp)                  :: gpcar(ndime,mnode,mgaus)
   real(rp)                  :: xjaci(ndime,ndime),xjacm(ndime,ndime)

   external                  :: chm_post_gp_lookup
   external                  :: elmder
   external                  :: pararr

   !
   ! Initialize
   !
   hrr_int_chm = 0.0_rp

   if ( INOTEMPTY ) then
      !
      ! Allocate auxiliary variable on elements
      !
      nullify(aux_r2p)
      call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_heatRelease_lookup',aux_r2p,nelem)
      do ielem = 1,nelem
         pelty = ltype(ielem)
         pgaus = ngaus(pelty)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_heatRelease_lookup',aux_r2p(ielem)%a,pgaus,&
             lookup_fw( kfl_hrr_fw_chm ) % main_table % nvar)
      end do

      !
      ! Lookup from HRR table
      !
      call chm_post_gp_lookup(aux_r2p, lookup_fw( kfl_hrr_fw_chm ))

      !
      ! Integrate Heat Release Rate
      !
      do ielem = 1,nelem
         pelty = ltype(ielem)
         pnode = nnode(pelty)
         pgaus = ngaus(pelty)

         do inode = 1,pnode
            ipoin = lnods(inode,ielem)
            do idime = 1,ndime
               elcod(idime,inode) = coord(idime,ipoin)
            end do
         end do


         do igaus = 1,pgaus
            !
            ! Get volume associated to Gaussian integration point
            !
            call elmder(&
                 pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                 elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
            gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg

            !
            ! Volumetric integral of heat release rate
            !
            hrr_int_chm = hrr_int_chm + gpvol * aux_r2p(ielem) % a(igaus, kfl_hrr_col_chm)

         enddo
      end do
      !
      ! Deallocate
      !
      call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_heatRelease_lookup',aux_r2p)
   endif

   !
   ! Sum heat release rates of subdomains
   !
   call PAR_SUM(hrr_int_chm)



end subroutine chm_heatRelease_lookup
