!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis 
!> @{
!> @file    pts_plugin.f90
!> @date    29/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine pts_plugin(icoup)

  use def_kintyp,    only :  ip,rp
  use def_domain,    only :  ndime
  use def_kermod,    only :  kfl_vefun,kfl_tefun
  use def_coupli,    only :  coupling_type
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use def_partis,    only :  kfl_momentum_sink_pts
  use def_partis,    only :  kfl_heat_sink_pts
  use def_partis,    only :  kfl_mass_sink_pts
  use def_master,    only :  INOTMASTER
  !
  ! Possible variables => 
  ! 
  use def_master,    only :  advec
  use def_master,    only :  therm
  use def_master,    only :  conce
  use def_master,    only :  momentum_sink
  use def_master,    only :  heat_sink
  use def_master,    only :  mass_sink
  use def_partis,    only :  nclas_pts
  use def_domain,    only :  npoin
  use mod_memory,    only :  memory_deallo
  use mod_memory,    only :  memory_alloca
  use def_master,    only :  mem_modul,modul
  implicit none
  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp), pointer       :: conce_communicate(:,:,:)
  integer(ip)             :: iclas,ipoin
  real(rp), pointer       :: dumm1(:)
  real(rp), pointer       :: dumm2(:,:)

  nullify(dumm1,dumm2)
  
  variable = coupling_type(icoup) % variable 
   
  if( variable == 'VELOC' .or. variable == 'ADVEC' ) then  
     !
     ! Velocity/Advection
     ! 
     if( kfl_vefun /= 99 ) call runend('PTS_PLUGIN: VELOCITY FUNCITON SHOULD BE 99 TO USE WITH PLUGIN')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,advec)

  else if( variable == 'TEMPE' .or. variable == 'ENTHA') then
     !
     ! Temperature or Enthalpy depending on the equation solved in Temper
     ! 
     if( kfl_tefun /= 99 ) call runend('PTS_PLUGIN: TEMPERATURE FUNCTiON SHOULD BE 99 TO USE WITH PLUGIN')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,therm)
     
  else if( variable == 'CONCE' ) then
     !
     ! Concentration
     ! 
     if( kfl_tefun /= 99 ) call runend('PTS_PLUGIN: concentration FUNCTiON SHOULD BE 99 TO USE WITH PLUGIN')

     nullify(conce_communicate)
     if (INOTMASTER) then
        call memory_alloca(mem_modul(1:2,modul),'CONCE_COMMUNICATE','pts_plugin',conce_communicate,nclas_pts,npoin,1_ip)
     else
        call memory_alloca(mem_modul(1:2,modul),'CONCE_COMMUNICATE','pts_plugin',conce_communicate,nclas_pts,1_ip,1_ip)
     endif
     conce_communicate = -1.0_rp

     call COU_INTERPOLATE_NODAL_VALUES(icoup,nclas_pts,conce_communicate)

     do iclas = 1,nclas_pts
        do ipoin = 1,npoin
            conce(ipoin,iclas,1) = conce_communicate(iclas,ipoin,1)
        end do
     end do

     call memory_deallo(mem_modul(1:2,modul),'CONCE_COMMUNICATE','pts_plugin',conce_communicate)

     
  else if( variable == 'MOMSK' ) then
     !
     ! Momentum sink
     ! 
     if( kfl_momentum_sink_pts == 0 ) call runend('PTS_PLUGIN: ACTIVATE MOMSK OPTION')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,dumm2,momentum_sink)
     
  else if( variable == 'HEASK' ) then
     !
     ! Heat sink
     ! 
     if( kfl_heat_sink_pts == 0 ) call runend('PTS_PLUGIN: ACTIVATE HEASK OPTION')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,dumm1,heat_sink)
     
  else if( variable == 'MASSK' ) then
     !
     ! Mass sink
     ! 
     if( kfl_mass_sink_pts == 0 ) call runend('PTS_PLUGIN: ACTIVATE MASSK OPTION')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,dumm1,mass_sink)
     
  end if

end subroutine pts_plugin
