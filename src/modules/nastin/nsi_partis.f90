!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_partis.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling with Partis module
!> @details Remove momentum from the Navier-Stokes equations
!>          \verbatim
!>          x_p ......... particle position
!>          Ff .......... force fluid on paticle
!>          Navier-Stokes = - \int_W (Ff/Vp).v.delta(x-x_p) dw,
!>                        = - (Ff(xp)/Vp).v(xp)*Vp
!>          So that the nodeal force Fi is:
!>                     Fi = - Fp * Ni
!>          \endverbatim
!> @}
!------------------------------------------------------------------------

subroutine nsi_partis()
  use def_master
  use def_domain
  use def_nastin
  use mod_parall,  only : commd
  implicit none
  integer(ip) :: idime,idofn,ipoin,jj
  real(rp)    :: multiplicity_loc

  if( associated(momentum_sink) ) then
     jj    = 0 
     do ipoin = 1,npoin
        if (ipoin <= npoi1) then
           multiplicity_loc = 1.0_rp
        else
           jj = jj + 1
           multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
        endif
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime+idime
           if( kfl_fixno_nsi(idime,ipoin) <= 0 ) then
              !
              ! Add source term:
              !
              rhsid(idofn) = rhsid(idofn) +  momentum_sink(idime,ipoin) / multiplicity_loc
           end if
        end do
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
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime+idime
           if( kfl_fixno_nsi(idime,ipoin) <= 0 ) then
              !
              ! Subtract dilution term:
              !
              rhsid(idofn)  = rhsid(idofn) -   mass_sink(ipoin) / multiplicity_loc * unkno(idofn)
           end if
        end do
     end do
  end if


  !real(rp)    :: elmomsk(ndime,mnode)
  !real(rp)    :: elmassk(mnode)
  !real(rp)    :: gpvol(mgaus)
  !real(rp)    :: gpmomsk(ndime,mgaus)
  !real(rp)    :: gpmassk(mgaus)
  !real(rp)    :: gpcar(ndime,mnode,mgaus) 
  !real(rp)    :: gphes(ntens,mnode) 
  !real(rp)    :: elcod(ndime,mnode) 
  !integer(ip) :: pnode,pgaus,pelty
  !integer(ip) :: ielem,igaus,inode

  !block
  !  real(rp) :: fostalicska,csereplec,sipcsont,dummr
  !  fostalicska = 0.0_rp
  !  csereplec = 0.0_rp
  !  sipcsont = 0.0_rp
  !  dummr = 0.0_rp
  !  if ( associated(gesca) )call memgen(2_ip,(ndime+1)*npoin,zero)
  !  nullify(gesca)
  !  call memgen(zero,(ndime+1)*npoin,zero)
  !  do ipoin = 1, npoin
  !     idofn = (ipoin-1)*ndime + 1
  !     gesca(ipoin) = rhsid(idofn)
  !  enddo
  !  call rhsmod(1_ip,gesca)
  !  !call rhsmod(1_ip,rhsid)
  !  if( INOTEMPTY  ) then
  !     call norm2x(1_ip,gesca,fostalicska)
  !  else
  !     call norm2x(1_ip,dummr,fostalicska)
  !  end if
  !  if( INOTSLAVE ) print*,'before rhsid=',fostalicska
  !end block

  !
  ! Loop over elements
  !
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
  !      elmomsk(:,:) = 0.0_rp
  !      if( associated(momentum_sink) ) then
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            elmomsk(1:ndime,inode) = momentum_sink(1:ndime,ipoin)
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
  !      gpmomsk(1:ndime,:) = 0.0_rp
  !      gpmassk(:)         = 0.0_rp
  !      do igaus = 1,pgaus     
  !         do inode = 1,pnode
  !            gpmomsk(1:ndime,igaus) = gpmomsk(1:ndime,igaus)   &
  !                          + elmar(pelty) % shape(inode,igaus) &
  !                          * elmomsk(1:ndime,inode)              
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
  !            do idime = 1,ndime
  !               idofn = (ipoin-1)*ndime + idime
  !               if( kfl_fixno_nsi(idime,ipoin) <= 0 ) then
  !                  !
  !                  ! Add source term:
  !                  !
  !                  rhsid(idofn) = rhsid(idofn) +  elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpmomsk(idime,igaus)  
  !                  !
  !                  ! Subtract dilution term:
  !                  !
  !                  rhsid(idofn) = rhsid(idofn) -  elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpmassk(igaus) * unkno(idofn)
  !               end if
  !            enddo
  !         enddo
  !      end do
  !   endif

  !end do elements



  !block
  !  real(rp) :: fostalicska,csereplec,sipcsont,dummr
  !  fostalicska = 0.0_rp
  !  csereplec = 0.0_rp
  !  sipcsont = 0.0_rp
  !  dummr = 0.0_rp
  !  if ( associated(gesca) )call memgen(2_ip,(ndime+1)*npoin,zero)
  !  nullify(gesca)
  !  call memgen(zero,(ndime+1)*npoin,zero)
  !  do ipoin = 1, npoin
  !     idofn = (ipoin-1)*ndime + 1
  !     gesca(ipoin) = rhsid(idofn)
  !  enddo
  !  call rhsmod(1_ip,gesca)
  !  !call rhsmod(1_ip,rhsid)
  !  if( INOTEMPTY  ) then
  !     call norm2x(1_ip,gesca,fostalicska)
  !  else
  !     call norm2x(1_ip,dummr,fostalicska)
  !  end if
  !  if( INOTSLAVE ) print*,'rhsid=',fostalicska
  !end block


end subroutine nsi_partis
