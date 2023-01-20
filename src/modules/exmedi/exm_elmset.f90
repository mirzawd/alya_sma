!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_elmset(iesec,ieset)
   !-----------------------------------------------------------------------
   !****f* Exmedi/exm_elmset
   ! NAME
   !    exm_elmset
   ! DESCRIPTION
   !    This routine computes variables on an element set W.
   ! USES
   ! USED BY
   !-----------------------------------------------------------------------
   use def_kintyp_basic,   only : ip,rp
   use def_master,         only : postp
   use def_master,         only : elmag,vconc
   use def_domain,         only : nelem, nnode, ngaus, elmar, ndime, coord, mnode, lnods, mgaus, ltype, leset
   implicit none

   integer(ip), intent(in)  :: iesec !< element of the set
   integer(ip), intent(in)  :: ieset !< tag of the set

   real(rp), pointer        :: s_intra(:), s_calci(:), s_sodiu(:), s_potas(:), setvo(:)
   integer(ip)              :: ielem, inode, pnode, igaus, pgaus, ipoin, pelty, i_current, nvabi
   real(rp)                 :: gpdet, gpvol
   real(rp)                 :: gpcar(ndime,mnode,mgaus)
   real(rp)                 :: xjacm(ndime,ndime), xjaci(ndime,ndime)
   real(rp)                 :: elcod(ndime,mnode)
   real(rp)                 :: e_intra(mnode), e_vconc(11,mnode)
   real(rp)                 :: gp_intra, gp_vconc(11)

   !----------------------------------------------------------------------
   ! Initialization
   !----------------------------------------------------------------------
   nvabi =  postp(1) % nvaes+1
   setvo => postp(1) % veset( nvabi:nvabi ,ieset)
   setvo =  0.0_rp  ! Set volume

   if( postp(1) % npp_setse(1) /= 0 ) s_intra => postp(1) % veset(1:1,ieset)
   if( postp(1) % npp_setse(2) /= 0 ) s_calci => postp(1) % veset(2:2,ieset)
   if( postp(1) % npp_setse(3) /= 0 ) s_sodiu => postp(1) % veset(3:3,ieset)
   if( postp(1) % npp_setse(4) /= 0 ) s_potas => postp(1) % veset(4:4,ieset)

   if( postp(1) % npp_setse(1) /= 0 ) s_intra = 0.0_rp
   if( postp(1) % npp_setse(2) /= 0 ) s_calci = 0.0_rp
   if( postp(1) % npp_setse(3) /= 0 ) s_sodiu = 0.0_rp 
   if( postp(1) % npp_setse(4) /= 0 ) s_potas = 0.0_rp

   elements: do ielem = 1,nelem

      if( leset(ielem) == iesec ) then
         !
         ! Element properties and dimensions
         !
         pelty = ltype(ielem)
         if( pelty > 0 ) then
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            !
            ! Gather operations
            !
            do inode = 1,pnode
               ipoin = lnods(inode,ielem)
               elcod(1:ndime,inode) = coord(1:ndime,ipoin)
               e_intra(inode) = elmag(ipoin,1)
               do i_current=1,11
                 e_vconc(i_current,inode)= vconc(i_current,ipoin,1)
               enddo
            end do

           gaus: do igaus=1, pgaus
               call elmder(&
                    pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                    elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
               gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
               setvo = setvo + gpvol

               gp_intra = 0.0_rp
               gp_vconc(:)=0.0_rp
               do inode = 1,pnode
                    gp_intra = gp_intra + elmar(pelty) % shape(inode,igaus) * e_intra(inode)
                    do i_current=1,11
                        gp_vconc(i_current) = gp_vconc(i_current) + elmar(pelty) % shape(inode,igaus) * e_vconc(i_current,inode)
                    enddo
               enddo

               !
               ! Intra
               !
               if( postp(1) % npp_setse(1) /= 0 ) then
                    s_intra = s_intra + gpvol*gp_intra 
               endif
               !
               ! Calci
               !
               if( postp(1) % npp_setse(2) /= 0 ) then
                    s_calci = s_calci + gpvol*gp_vconc(1)
               endif
               !
               ! Sodiu
               !
               if( postp(1) % npp_setse(3) /= 0 ) then
                    s_sodiu = s_sodiu + gpvol*gp_vconc(3)
               endif
               !
               ! Potas
               !
               if( postp(1) % npp_setse(4) /= 0 ) then
                    s_potas = s_potas + gpvol*gp_vconc(4)
               endif

           enddo gaus


         endif ! from if (pelty>0)



      endif ! from if(leset(ielem) == iesec)
   enddo elements



end subroutine exm_elmset
