!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_chemic(&
     ielem,pgaus,gprhs)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_chemic
  ! NAME
  !   tem_radiat
  ! DESCRIPTION
  !    Couple to the heat source from chemical processes
  ! USES
  ! USED BY
  !    tem_elmop2 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  div_enthalpy_transport,kfl_coupl, &
                                ID_TEMPER,ID_CHEMIC
  implicit none 
  integer(ip), intent(in)    :: ielem,pgaus
  real(rp),    intent(inout) :: gprhs(pgaus)
  integer(ip)                :: igaus

  !
  ! Coupling with CHEMIC
  !
  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .and. associated(div_enthalpy_transport)) then
     !
     ! RHS terms: Divergence of enthalpy transport by diffusion
     !
     do igaus=1,pgaus
        gprhs(igaus) = gprhs(igaus) + div_enthalpy_transport(ielem)%a(igaus,1,1)
     end do

  endif

end subroutine tem_chemic
