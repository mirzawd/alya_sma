!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_temperature_from_enthalpy.f90
!> @author  houzeaux
!> @date    2019-11-27
!> @brief   Compute temperature from enthalpy
!> @details Compute temperature from enthalpy
!> @} 
!-----------------------------------------------------------------------

subroutine tem_temperature_from_enthalpy()

  use def_master
  use def_temper
  use def_domain
  implicit none
    
  if (kfl_regim_tem == 4_ip) then
     !
     ! Compute temperature from enthalpy and chemistry
     !
     if (kfl_lookg_tem > 0 ) then
        !
        ! Lookup was done on Gauss points. 
        ! Temperature is calculated on Gauss points from the enthalpy and NASA polynomials
        ! Nodal temperature is projected, and boundary temperatures are imposed.
        !
        call tem_gp_comtem()
        
     else
        !
        ! Lookup was done on nodes.
        ! Temperature is calculated on nodes from enthalpy and NASA polynomials.
        !
        call tem_comtem()

     end if
     !
     ! Impose boundary conditions to overwrite what0s been computed
     !
     call tem_calcEnthalpyBC()
     
  end if

end subroutine tem_temperature_from_enthalpy
