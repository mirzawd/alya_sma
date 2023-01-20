!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_pprdes(lnods,pnode,fddes,gddes,pgaus,gpvol,gpsha)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_pprdes
  ! NAME
  !   tur_pprdes
  ! DESCRIPTION
  !    Smoothes the FDDES and GDDES values to make it usable for the postprocessing
  ! OUTPUT 
  !
  ! USES
  ! USED BY
  !    tur_elmop2
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  fddes_tur, gddes_tur
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: fddes(pgaus)
  real(rp),    intent(in)    :: gddes(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  integer(ip)                :: inode, igaus
  real(rp)                   :: elfde(1,pnode)  !element fddes
  real(rp)                   :: elgde(1,pnode)  !element gddes

  do inode = 1,pnode
    elfde(1,inode) = 0.0_rp
    elgde(1,inode) = 0.0_rp
  end do

  do igaus = 1,pgaus
     do inode = 1,pnode
         elfde(1,inode) =  elfde(1,inode) + fddes(igaus) * gpvol(igaus) * gpsha(inode,igaus)
         elgde(1,inode) =  elgde(1,inode) + gddes(igaus) * gpvol(igaus) * gpsha(inode,igaus)
     end do
  end do

  !
  ! Assemble
  !

  call assrhs(1_ip,pnode,lnods,elfde,fddes_tur)
  call assrhs(1_ip,pnode,lnods,elgde,gddes_tur)

end subroutine tur_pprdes
