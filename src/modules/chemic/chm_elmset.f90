!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_elmset(iesec,ieset)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_elmset
  ! NAME
  !    chm_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are:
  !                                      where w=dv/dx-du/dy
  !    SETVE and SETVR are normailzed further on in chm_outset
  ! USES
  !    chm_elmgat
  !    elmder
  ! USED BY
  !    chm_outset
  !***
  !-----------------------------------------------------------------------
  use def_master,     only : postp, conce, veset
  use def_domain,     only : ndime, mnode, mgaus, elmar, nelem, leset, ltype, nnode, ngaus, lnods, coord
  use def_kintyp,     only : ip, rp
  use def_chemic,     only : mass_gp, kfl_model_chm, nclas_chm, hrr_chm
  use mod_ker_proper, only : ker_proper
  implicit none
  integer(ip), intent(in)  :: iesec,ieset
  real(rp),    pointer     :: setvo(:)
  real(rp),    pointer     :: setmass(:)
  real(rp),    pointer     :: sethr(:)
  real(rp),    pointer     :: setconc(:)
  integer(ip)              :: pnode,pgaus,pelty,nvabi,ispec
  integer(ip)              :: ielem,igaus,idime,inode,ipoin,iclas
  integer(ip)              :: dummi
  real(rp)                 :: gpcar(ndime,mnode,mgaus)
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: elcon(mnode,nclas_chm)
  real(rp)                 :: elhrr(mnode)
  real(rp)                 :: gpcon(mgaus,nclas_chm)
  real(rp)                 :: gpden(mgaus)
  real(rp)                 :: gphr(mgaus)
  real(rp)                 :: gpvol,gpdet

  external                 :: elmder

  !
  ! Initialization
  !
  nvabi =  postp(1) % nvaes+1
  setvo => veset( nvabi:nvabi ,ieset)
  setvo =  0.0_rp  ! Set volume

  if( postp(1) % npp_setse(1) /= 0 ) setmass => postp(1) % veset(1:1,ieset)
  if( postp(1) % npp_setse(2) /= 0 ) sethr   => postp(1) % veset(2:2,ieset)
  if( postp(1) % npp_setse(3) /= 0 ) setconc => postp(1) % veset(3:6,ieset)

  if( postp(1) % npp_setse(1) /= 0 ) setmass = 0.0_rp
  if( postp(1) % npp_setse(2) /= 0 ) sethr   = 0.0_rp
  if( postp(1) % npp_setse(3) /= 0 ) setconc = 0.0_rp
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( leset(ielem) == iesec ) then
        !
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Gather operations
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do ispec = 1,nclas_chm
              elcon(inode,ispec) = conce(ipoin,ispec,1)
           end do
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        if( kfl_model_chm == 3 ) then
           do inode = 1,pnode
              elhrr(inode) = hrr_chm(ipoin)
           enddo
        endif
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
        !
        ! Gauss point values
        !
        do igaus = 1,pgaus
           do ispec = 1,nclas_chm
              gpcon(igaus,ispec) = 0.0_rp
           end do
        end do
        do ispec = 1,nclas_chm
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gpcon(igaus,ispec) = gpcon(igaus,ispec) &
                      +elmar(pelty)%shape(inode,igaus)   &
                      *elcon(inode,ispec)
              end do
           end do
        end do

        if( kfl_model_chm == 1 ) then
           !
           ! Flamelet model: hr is substituted by omega_Yc
           !
           do igaus = 1,pgaus
              gphr(igaus) = mass_gp(ielem) % a(igaus,1,1)
           end do
        elseif( kfl_model_chm == 3 ) then
           !
           ! Finite rate chemistry model
           !
           do igaus = 1,pgaus
              gphr(igaus) = 0.0_rp
              do inode = 1,pnode
                 gphr(igaus) = gphr(igaus) +elmar(pelty)%shape(inode,igaus) * elhrr(inode)
              enddo
           end do
        else
           do igaus = 1,pgaus
              gphr(igaus) = 0.0_rp
           end do
        endif

        !
        ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
        !
        do igaus = 1,pgaus
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
           gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
           setvo = setvo + gpvol

           if( postp(1) % npp_setse(1) /= 0 ) then
              !
              ! Mass
              !
              setmass = setmass + gpvol * gpden(igaus)
           end if

           if( postp(1) % npp_setse(2) /= 0 ) then
              !
              ! Heat realease
              !
              sethr = sethr + gpvol * gphr(igaus)
           end if

           if( postp(1) % npp_setse(3) /= 0 ) then
              !
              ! Mass fraction
              !
              do iclas = 1,min(nclas_chm,size(setconc, KIND=ip))
                 setconc(iclas) = setconc(iclas) + gpvol * gpden(igaus) * gpcon(igaus,iclas)
              enddo

           end if


        end do
     end if

  end do elements

end subroutine chm_elmset

