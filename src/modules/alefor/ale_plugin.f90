!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    ale_plugin.f90
!> @author  J.C. Cajas
!> @date    04/11/2014
!> @brief   Plugin for zonal coupling with other modules
!> @details Plugin for zonal coupling with other modules. 
!> @        Used in the FSI coupling with Nastin and Solidz.
!> @} 
!----------------------------------------------------------------------

subroutine ale_plugin(icoup)

  use def_master,        only :  INOTMASTER
  use def_master,        only :  kfl_fixno_ale
  use def_master,        only :  bvess_ale
  use def_coupli,        only :  coupling_type
  use def_domain,        only :  npoin
  use def_domain,        only :  ndime
  use def_master,        only :  modul,mem_modul
  use def_kintyp,        only :  ip,rp
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,     only :  COU_SET_FIXITY_ON_TARGET
  use mod_memory,        only :  memory_deallo
  use mod_memory,        only :  memory_alloca
  use mod_matrix,        only :  matrix_initialize

  implicit none

  real(rp),    pointer    :: svalu(:,:)

  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  integer(ip), save       :: ipass = 0

  nullify(svalu)
  variable = coupling_type(icoup) % variable

  if( variable == 'ALEFO' ) then
     !
     ! Coupling with solidz
     !
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'SVALU','ale_plugin',svalu,ndime,npoin)
        svalu = 0.0_rp
     else
        call memory_alloca(mem_modul(1:2,modul),'SVALU','ale_plugin',svalu,1_ip,1_ip)
        svalu = 0.0_rp
     end if
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('ALEFO',modul,kfl_fixno_ale)
     end if
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_ale,svalu,kfl_fixno_ale)

  end if

  if( associated(svalu) ) call memory_deallo(mem_modul(1:2,modul),'SVALU','ale_plugin',svalu)
  
end subroutine ale_plugin
!> @} 
!-----------------------------------------------------------------------
