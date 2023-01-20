!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine codfix(ndime,kfl_fixno)
!-----------------------------------------------------------------------
!****f* wtool/codifx
! NAME 
!    codfix
! DESCRIPTION
!    This routine codes fixno
! USES
! USED BY
!    
!***
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)    :: ndime
  integer(ip), intent(inout) :: kfl_fixno(ndime)
  integer(ip)                :: icode,idime
  real(rp)                   :: rcode

  icode = kfl_fixno(1)
  do idime = 1,ndime
     rcode = icode/10**(ndime-idime)
     kfl_fixno(idime) = int(rcode)
     icode = icode - kfl_fixno(idime)*10**(ndime-idime)
  end do

end subroutine codfix

