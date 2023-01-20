!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_asspro(&
     itask,pnode,pevat,pgaus,lnods,gpden,gpvis,gpvol,gpsha,&
     elrhs,rhsid)
  !-----------------------------------------------------------------------
  !****f* nastin/tur_asspro
  ! NAME 
  !    tur_asspro
  ! DESCRIPTION
  !    Assembly of the RHS
  ! USES
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp
  implicit none
  integer(ip),  intent(in)    :: itask,pnode,pevat,pgaus
  integer(ip),  intent(in)    :: lnods(pnode)
  real(rp),     intent(in)    :: gpden(pgaus)
  real(rp),     intent(in)    :: gpvis(pgaus)
  real(rp),     intent(in)    :: gpvol(pgaus)
  real(rp),     intent(in)    :: gpsha(pnode,pgaus)
  real(rp),     intent(out)   :: elrhs(pevat)
  real(rp),     intent(inout) :: rhsid(*)
  integer(ip)                 :: inode,ipoin,idofl,idofg,igaus
  real(rp)                    :: fact1,fact2

  do inode = 1,pevat
     elrhs(inode) = 0.0_rp
  end do
  
  if( itask == 10 ) then
     !
     ! Density
     !
     do igaus = 1,pgaus
        fact1 = gpden(igaus) * gpvol(igaus)
        do inode = 1,pnode
           elrhs(inode) = elrhs(inode) + fact1 * gpsha(inode,igaus) 
        end do
     end do

  else if( itask == 11 ) then
     !
     ! Viscosity
     !
     do igaus = 1,pgaus
        fact1 = gpvis(igaus) * gpvol(igaus)
        do inode = 1,pnode
           elrhs(inode) = elrhs(inode) + fact1 * gpsha(inode,igaus) 
        end do
     end do
     
  else if( itask == 12 ) then
     !
     ! Density and viscosity
     !
     do igaus = 1,pgaus
        fact1 = gpden(igaus) * gpvol(igaus)
        fact2 = gpvis(igaus) * gpvol(igaus)
        do inode = 1,pnode
           idofl = (inode-1)*2 + 1
           elrhs(idofl) = elrhs(idofl) + fact1 * gpsha(inode,igaus) 
           idofl = idofl + 1
           elrhs(idofl) = elrhs(idofl) + fact2 * gpsha(inode,igaus) 
        end do
     end do
     
  end if

  !----------------------------------------------------------------------
  !
  ! Assembly
  !
  !----------------------------------------------------------------------

  if( itask == 10 .or. itask == 11 ) then
     do inode = 1,pnode
        ipoin        = lnods(inode)
        rhsid(ipoin) = rhsid(ipoin) + elrhs(inode)
     end do

  else if( itask == 12 ) then
     do inode = 1,pnode
        ipoin        = lnods(inode)
        idofg        = (ipoin-1)*2
        idofl        = (inode-1)*2
        idofg        = idofg + 1
        idofl        = idofl + 1
        rhsid(idofg) = rhsid(idofg) + elrhs(idofl)
        idofg        = idofg + 1
        idofl        = idofl + 1
        rhsid(idofg) = rhsid(idofg) + elrhs(idofl)
     end do

  end if

end subroutine tur_asspro

