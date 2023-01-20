!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_qp2rp
  use def_kintyp_basic, only: rp,qp
  implicit none
    
  public :: qp2rp

contains

  function qp2rp(value_in) result(value_out)

    implicit none
    real(qp), intent(in) :: value_in
    real(rp)             :: value_out
    
#if !defined Q16 && !defined R4
    value_out = value_in
#else
    value_out = real(value_in, rp)
#endif
  end function qp2rp

end module mod_qp2rp
