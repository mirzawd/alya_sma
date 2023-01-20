!-----------------------------------------------------------------!
! Copyright 2005 - 2021 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



pure elemental function safe_exp_cond(x, enabled) result(r)
    ! enabled -- if 0, normal exponent is called, if 1 - safe exponent
    use def_kintyp, only : rp, lg, ip
    implicit none
    real(rp), intent(in)              :: x
    integer(ip), intent(in)           :: enabled
    real(rp)                          :: r
    real(rp), parameter               :: tol = 100.0_rp
    real(rp), parameter               :: expptol = exp(tol)
    real(rp), parameter               :: expntol = exp(-tol)
          
    r = exp( real(enabled,rp)*sign(min(abs(x),tol),x) + real(1-enabled,rp)*x )
end function safe_exp_cond



pure elemental function safe_exp(x) result(r)
    use def_kintyp, only : rp
    implicit none
    real(rp), intent(in)              :: x
    real(rp)                          :: r
    real(rp), parameter               :: tol = 100.0_rp
          
    r = exp(sign(min(abs(x),tol),x))
end function safe_exp