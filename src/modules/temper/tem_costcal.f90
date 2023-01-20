!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> tem_costcal.f90
!> @file tem_costcal.f90 
!> @fn tem_costcal 
!> This subroutine calculates the costf by summing between subdomain and demostrating
!>

subroutine tem_costcal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : costf,kfl_cost_type
  use def_domain
  use mod_communications, only : PAR_SUM
  use def_temper
  implicit none

  
    
  call PAR_SUM(costf,  'IN MY CODE' )
!   print*, "costf temp", costf, (costf-340.102760851646_rp)/0.00001_rp
!   print*, "costf temp", costf, (costf-1578.56301580740_rp)/0.00001_rp
  if (kfl_cost_type == 1) print*, "costf temp", costf, (costf-775190.831160077_rp)/0.00001_rp
    

end subroutine tem_costcal
