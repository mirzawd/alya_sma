!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mixing(itask,wemix,slmix,whmix,value)
  !-----------------------------------------------------------------------
  !****f* mathru/mixing
  ! NAME
  !   mixing
  ! DESCRIPTION
  !   Computes a mixing function
  ! OUTPUT 
  !   WEIGH
  ! USES
  ! USED BY
  !    tur_tworod
  !    tur_updedd
  !***
  !-----------------------------------------------------------------------
  use      def_kintyp 
  implicit none
  integer(ip), intent(in)  :: itask
  real(rp),    intent(in)  :: slmix,value,whmix
  real(rp),    intent(out) :: wemix

  select case(itask)

  case(1)
     !
     ! From 0 to 1. Equal to 0.5 at 1.
     !
     wemix  = 0.5_rp*(tanh(slmix*(value-whmix))+1.0_rp)

  end select

end subroutine mixing
