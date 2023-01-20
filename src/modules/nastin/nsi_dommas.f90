!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_dommas(xmass)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_dommas
  ! NAME 
  !    nsi_dommas
  ! DESCRIPTION
  !    This routine computes the mass of the domain:
  !     +-
  !     |   div(u) dw
  !    -+ W
  ! USES
  !    elmder
  ! USED BY
  !    nsi_cvgunk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin 
  implicit none
  real(rp),    intent(out) :: xmass
  integer(ip)              :: pnode,pgaus,pelty
  integer(ip)              :: ielem,igaus,idime,inode,ipoin
  real(rp)                 :: gpcar(ndime,mnode,mgaus) 
  real(rp)                 :: xjaci(9),xjacm(9) 
  real(rp)                 :: elvel(ndime,mnode),elcod(ndime,mnode)
  real(rp)                 :: gpvol,gpdet,xvolu
  real(rp),    target      :: xvalu(2)

  xmass = 0.0_rp
  xvolu = 0.0_rp

  if( INOTMASTER ) then
     !
     ! Loop over elements
     !
     elements: do ielem = 1,nelem
        ! 
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! ELCOD and ELVEL: Gather operations
           !
           if( NSI_MONOLITHIC ) then
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elvel(idime,inode) = unkno((ipoin-1)*(ndime+1)+idime)
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
           else
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elvel(idime,inode) = unkno((ipoin-1)*ndime+idime)
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
           end if
           !
           ! XMASS and XVOLU: mass and volume
           !
           do igaus = 1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
              gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
              xvolu = xvolu + gpvol
              do inode = 1,pnode
                 do idime = 1,ndime
                    xmass = xmass + gpvol*elvel(idime,inode)*gpcar(idime,inode,igaus)
                 end do
              end do
           end do
        end if

     end do elements

  end if

  if( IPARALL ) then
     nparr    =  2
     xvalu(1) =  xmass
     xvalu(2) =  xvolu
     parre    => xvalu
     call par_operat(3_ip)
     xmass    =  xvalu(1)
     xvolu    =  xvalu(2)
  end if
  xmass = abs(xmass / (xvolu**(1.0_rp/real(ndime,rp))**(real(ndime,rp)-1.0_rp)))

end subroutine nsi_dommas

