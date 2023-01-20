!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_dyncou(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_dyncou
  ! NAME 
  !    tem_dyncou
  ! DESCRIPTION
  !    This routine manages the dynamic coupling of temper module
  ! USES
  ! USED BY
  !    tem_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_temper
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask

  select case(kfl_dynco_tem)

  case(1_ip)
     !
     ! Marek: Cooling system + typical room
     !
     call tem_dync01(itask)

  end select

end subroutine tem_dyncou
