!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_ifconf(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_ifconf
  ! NAME 
  !    nsi_ifconf
  ! DESCRIPTION
  !    This routine checks if the flow is confined or not and look for a node
  !    to impose pressure 
  ! USED BY
  !    nsi_solite
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iperi

  return

  select case(itask)

  case(2_ip)

     if( INOTMASTER .and. kfl_confi_nsi == 1 .and. nodpr_nsi /= 0 .and. nperi > 0 ) then
        do iperi = 1,nperi
           if( lperi(2,iperi) == nodpr_nsi ) then
              nodpr_nsi = lperi(1,iperi)
           end if
        end do

     end if

  end select

end subroutine nsi_ifconf
