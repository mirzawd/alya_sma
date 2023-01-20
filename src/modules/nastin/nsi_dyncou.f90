!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_dyncou(itask)
  !------------------------------------------------------------------------
  !****f* Temper/nsi_dyncou
  ! NAME 
  !    nsi_dyncou
  ! DESCRIPTION
  !    This routine manages the dynamic coupling of nastin module
  ! USES
  ! USED BY
  !    nsi_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask

  select case(kfl_dynco_nsi)

  case(1_ip)
     !
     ! Marek: Cooling system + typical room
     !
     call nsi_dync01(itask)

  end select

end subroutine nsi_dyncou
