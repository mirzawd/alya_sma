!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_plugin.f90
!> @date    20/04/2017
!> @author  Vishal Mehta
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine exm_plugin(icoup)

   use def_kintyp,         only :  ip,rp, lg
   use def_coupli,         only :  coupling_type
   use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
   use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
   use mod_communications, only :  PAR_SUM,PAR_BARRIER
   use mod_memory,         only :  memory_alloca
   use mod_memory,         only :  memory_deallo
   use mod_eccoupling

   implicit none

   integer(ip), intent(in) :: icoup
   character(5)            :: variable

   variable = coupling_type(icoup) % variable 

   if( variable == 'ECCOU' )then
      ! Electro-mechanical coupling
      call eccou_do_plugin(icoup)   

   else  if( variable == 'DISPL' ) then
      !Lagrangian mechano-electric coupling
      call eccou_do_lagrangian(icoup)   

   end if

end subroutine exm_plugin
