!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_comtem
  !-----------------------------------------------------------------------
  !****f* Temper/solve_tem
  ! NAME 
  !    solve_tem
  ! DESCRIPTION
  !    This routine computes the temperature from the enthalpy using the 
  !    polynomial coefficients for the specific heat
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_ker_proper
  use mod_physics, only : physics_H_2_TCp
  implicit none
  integer(ip)             :: ipoin,ivalu
  real(rp)                :: cploc(6,2),ent_loc,acval

  !
  ! total enthalpy equation, source term implicit in cp coefficients
  !
  do ipoin=1,npoin
    !
    ! Check if adiabatic calculation
    !
    if (kfl_adiab_tem /= 0) then           ! Enthalpy is given by the mixture fraction
       acval   = min(1.0_rp,max(0.0_rp,(conce(ipoin,3,1))))
       ent_loc = acval * cfi_hmax_tem + (1.0_rp - acval) * cfi_hmin_tem
    else
       ent_loc = therm(ipoin,1)            ! Enthalpy is recomputed each time-step
    endif

    do ivalu = 1,6
      cploc(ivalu,1) = sphec(ipoin,ivalu,1)
      cploc(ivalu,2) = sphec(ipoin,ivalu,2)
    end do
               
    call physics_H_2_TCp(ent_loc, cploc, tempe(ipoin,1), sphek(ipoin,1)) 

  end do

end subroutine tem_comtem
