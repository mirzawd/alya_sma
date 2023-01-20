!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_calcEnthalpyBC()
!-----------------------------------------------------------------------
!****f* Temper/tem_calcEnthalpyBC
! NAME 
!    tem_endite
! DESCRIPTION
!    This routine re-computes the enthalpy after the properties were updated in chemic to 
!    allow obtaining correct temperature BC for isothermal cases.
! USES
!    physics_T_2_HCp
! USED BY
!    tem_endite
!    tem_endste
!    tem_updbcs
!***
!-----------------------------------------------------------------------
  use def_kintyp,   only : ip,rp
  use def_master,   only : sphec
  use def_domain,   only : npoin
  use def_temper,   only : kfl_fixno_tem,bvess_tem
  use def_temper,   only : kfl_regim_tem,kfl_plepp_tem
  use def_temper,   only : posit_tem,negat_tem,bvtem_tem
  use def_temper,   only : kfl_posit_tem,kfl_negat_tem
  use mod_physics,  only : physics_T_2_HCp
  use mod_communications, only : PAR_MIN,PAR_MAX
  implicit none
  integer(ip) :: ipoin,ivalu
  real(rp)    :: dummr,hnew,cploc(6,2)

  if( kfl_regim_tem == 4 .and. kfl_plepp_tem /= 4 ) then
     
     posit_tem = -huge(1.0_rp)
     negat_tem =  huge(1.0_rp)
     
     do ipoin = 1,npoin
        if( kfl_fixno_tem(1,ipoin) == 1 ) then
           do ivalu = 1,6
              cploc(ivalu,1) = sphec(ipoin,ivalu,1)
              cploc(ivalu,2) = sphec(ipoin,ivalu,2)
           end do 
           dummr = 0.0_rp
           call physics_T_2_HCp(bvtem_tem(1,ipoin,1),cploc,hnew,dummr)
           bvess_tem(1,ipoin,1) = hnew
           posit_tem            = max(posit_tem,hnew)
           negat_tem            = min(negat_tem,hnew)
        end if
     end do
     !
     ! Clipping
     !
     if( kfl_posit_tem == 1 ) call PAR_MAX(posit_tem,'IN MY CODE')
     if( kfl_negat_tem == 1 ) call PAR_MIN(negat_tem,'IN MY CODE')

  end if

end subroutine tem_calcEnthalpyBC
