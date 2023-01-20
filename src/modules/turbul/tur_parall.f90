!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_parall
  ! NAME
  !    tur_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask

  if( IPARALL ) then

     select case(itask)

     case(1_ip)
        !
        ! Exchange data read in tur_reaphy, tur_reanut and tur_reaous
        ! always using MPI, even if this is a partition restart
        !
        call tur_sendat(1_ip)

     case(2_ip)
        !
        ! Exchange data read in tur_reabcs
        !
        call tur_sendat(2_ip)

     case(3_ip)

     end select

  end if

end subroutine tur_parall
