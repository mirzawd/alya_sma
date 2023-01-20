!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_radiat(&
     ielem,pnode,pgaus,lnods,gpsha,gprhs)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_radiat
  ! NAME
  !   tem_radiat
  ! DESCRIPTION
  !    Couple to the heat source from radiation
  ! USES
  ! USED BY
  !    tem_elmop2
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  radso,kfl_coupl,ID_TEMPER,ID_RADIAT
  implicit none 
  integer(ip), intent(in)    :: pnode,pgaus,lnods(pnode),ielem
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(inout) :: gprhs(pgaus)
  real(rp)                   :: gprad(pgaus)
  integer(ip)                :: inode,igaus,ipoin

  !
  ! Coupling with RADIATion
  !
  if (kfl_coupl(ID_TEMPER,ID_RADIAT) >= 1 ) then
     ! Initialize
     do igaus=1,pgaus
        gprad(igaus)=0.0_rp
     end do
     ! Gather
     do inode=1,pnode
        ipoin=lnods(inode) 
        do igaus=1,pgaus
           gprad(igaus)=gprad(igaus)&
                +gpsha(inode,igaus)*radso(ipoin,1)
        end do
     end do
     ! Add to RHS
     do igaus=1,pgaus
        gprhs(igaus)= gprhs(igaus)+gprad(igaus)
     end do
  endif

end subroutine tem_radiat
