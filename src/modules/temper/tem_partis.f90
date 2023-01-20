!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_partis.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling with Partis module
!> @details Exchange heat source in energy equation 
!> @}
!------------------------------------------------------------------------

subroutine tem_partis()
  use def_master
  use def_domain
  use mod_parall,  only : commd
  implicit none
  integer(ip) :: ipoin,jj
  real(rp)    :: multiplicity_loc

  if( associated(heat_sink) ) then
     jj    = 0 
     do ipoin = 1,npoin
        if (ipoin <= npoi1) then
           multiplicity_loc = 1.0_rp
        else
           jj = jj + 1
           multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
        endif
        !
        ! Add heat source term
        !
        rhsid(ipoin) = rhsid(ipoin) + heat_sink(ipoin)/multiplicity_loc
     end do
  end if
  if( associated(mass_sink) ) then
     jj    = 0 
     do ipoin = 1,npoin
        if (ipoin <= npoi1) then
           multiplicity_loc = 1.0_rp
        else
           jj = jj + 1
           multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
        endif
        !
        ! Subtract dilution term:
        !
        rhsid(ipoin)  = rhsid(ipoin) - mass_sink(ipoin)/multiplicity_loc * unkno(ipoin)
     end do
  end if

  !real(rp)    :: elheask(mnode)
  !real(rp)    :: elmassk(mnode)
  !real(rp)    :: gpvol(mgaus)
  !real(rp)    :: gpheask(mgaus)
  !real(rp)    :: gpmassk(mgaus)
  !real(rp)    :: gpcar(ndime,mnode,mgaus) 
  !real(rp)    :: gphes(ntens,mnode) 
  !real(rp)    :: elcod(ndime,mnode) 
  !integer(ip) :: pnode,pgaus,pelty
  !integer(ip) :: ielem,igaus,inode
  !!
  !! Loop over elements
  !!
  !elements: do ielem = 1,nelem
  !   ! 
  !   ! Element properties and dimensions
  !   !
  !   pelty = ltype(ielem)
  !   if( pelty > 0 ) then
  !      pnode = nnode(pelty)
  !      pgaus = ngaus(pelty)
  !      !
  !      ! Gather operations
  !      !
  !      do inode = 1,pnode
  !         ipoin = lnods(inode,ielem)
  !         elcod(1:ndime,inode) = coord(1:ndime,ipoin)
  !      enddo
  !      elheask(:) = 0.0_rp
  !      if( associated(heat_sink) ) then
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            elheask(inode) = heat_sink(ipoin)
  !         enddo
  !      endif
  !      elmassk(:) = 0.0_rp
  !      if( associated(mass_sink) ) then
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            elmassk(inode) = mass_sink(ipoin)
  !         enddo
  !      endif
  !      call elmcar(&
  !           pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
  !           elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
  !           gphes,ielem)
  !         
  !      gpheask(:) = 0.0_rp
  !      gpmassk(:) = 0.0_rp
  !      do igaus = 1,pgaus     
  !         do inode = 1,pnode
  !            gpheask(igaus) = gpheask(igaus)                   &
  !                          + elmar(pelty) % shape(inode,igaus) &
  !                          * elheask(inode)              
  !            gpmassk(igaus) = gpmassk(igaus)                   &
  !                          + elmar(pelty) % shape(inode,igaus) &
  !                          * elmassk(inode)              
  !         enddo
  !      enddo

  !      !
  !      ! Add to right hand side
  !      !
  !      do igaus = 1,pgaus     
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            if( kfl_fixno_tem(1,ipoin) <= 0 ) then
  !               !
  !               ! Add source term:
  !               !
  !               rhsid(ipoin) = rhsid(ipoin) +  elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpheask(igaus)  
  !               !
  !               ! Subtract dilution term:
  !               !
  !               rhsid(ipoin) = rhsid(ipoin) -  elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpmassk(igaus) * unkno(ipoin)
  !            end if
  !         enddo
  !      end do
  !   endif

  !end do elements

end subroutine tem_partis
