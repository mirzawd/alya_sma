!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmini(pevat,pgaus,wmatr,wrhsi)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmini
  ! NAME 
  !    nsi_elmini
  ! DESCRIPTION
  !    Initialize matrix and RHS
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pevat,pgaus
  real(rp),    intent(out) :: wmatr(pevat,pevat,pgaus)
  real(rp),    intent(out) :: wrhsi(pevat,pgaus)
  integer(ip)              :: ievat,jevat,igaus
  !
  ! Initialization
  !
  do igaus=1,pgaus
     do ievat=1,pevat
        wrhsi(ievat,igaus)=0.0_rp
        do jevat=1,pevat
           wmatr(jevat,ievat,igaus)=0.0_rp
        end do
     end do
  end do

end subroutine nsi_elmini
