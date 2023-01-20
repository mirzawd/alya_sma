!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> tur_costcal.f90
!> @file tur_costcal.f90 
!> @fn tur_costcal 
!> This subroutine calculates the costf by summing between subdomain and demostrating
!>

subroutine tur_costcal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : costf,kfl_cost_type
  use def_domain
  use mod_communications, only : PAR_SUM
  use def_turbul
  implicit none

  if (kfl_cost_type == 6) then
    call PAR_SUM(costf,  'IN MY CODE' )
    if( INOTMASTER) then
        
    else
      print*, "costf tur paral", costf, (costf-36.4591404223040_rp)/0.000001_rp
    endif

  endif
                                     

end subroutine tur_costcal
