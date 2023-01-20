!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_pprsst(itask,lnods,pnode,pprva,pgaus,gpvol,gpsha)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_pprdes
  ! NAME
  !   tur_pprdes
  ! DESCRIPTION
  !    Smoothes values to make it usable for the postprocessing (SST-Model)
  ! OUTPUT 
  !
  ! USES
  ! USED BY
  !    tur_elmop2
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  sstf1_tur, sasso_tur
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: pprva(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  integer(ip)                :: inode, igaus
  real(rp)                   :: elppr(1,pnode)  !element postprocess value

  do inode = 1,pnode
    elppr(1,inode) = 0.0_rp
  end do

  do igaus = 1,pgaus
     do inode = 1,pnode
         elppr(1,inode) =  elppr(1,inode) + pprva(igaus) * gpvol(igaus) * gpsha(inode,igaus)
     end do
  end do
  !
  ! Assemble
  !
  select case (itask)
  case (1_ip)
     call assrhs(1_ip,pnode,lnods,elppr,sstf1_tur)
  case (2_ip)
     call assrhs(1_ip,pnode,lnods,elppr,sasso_tur)
  end select

end subroutine tur_pprsst
