!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_plugin.f90
!> @author  Guillaume Houzeaux
!> @date    13/04/2014
!> @brief   Plugin with coupling
!> @details Plugin for the zonal coupling
!> @}
!-----------------------------------------------------------------------

subroutine chm_plugin(icoup)
  !
  ! Obligatory variables
  !
  use def_kintyp,    only :  ip,rp
  use def_coupli,    only :  coupling_type
  use def_master,    only :  INOTMASTER
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_memory,    only :  memory_deallo
  use mod_memory,    only :  memory_alloca
  use def_master,    only :  mem_modul,modul
  use def_domain,    only :  npoin
  use def_chemic,    only :  nclas_chm
  !
  ! Possible variables
  !
  use def_master,    only :  mass_sink, conce

  use mod_chm_operations_CMC, only : chm_data_CMC_to_CFD_domain_CMC, &
                                     chm_data_CFD_to_CMC_domain_CMC
  implicit none

  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp), pointer       :: dummr(:,:,:), dummr2(:,:)
  real(rp), pointer       :: conce_communicate(:,:,:)

  integer(ip)             :: ipoin,iclas

  nullify(conce_communicate)
  nullify(dummr)
  nullify(dummr2)

  variable = coupling_type(icoup) % variable

  if( variable == 'CONCE' ) then
     !
     ! Concentration
     !
     !!!call COU_INTERPOLATE_NODAL_VALUES(icoup,nclas_chm,dummr,conce)
     if (INOTMASTER) then
        call memory_alloca(mem_modul(1:2,modul),'CONCE_COMMUNICATE','chm_plugin',conce_communicate,nclas_chm,npoin,1_ip)
     else
        call memory_alloca(mem_modul(1:2,modul),'CONCE_COMMUNICATE','chm_plugin',conce_communicate,nclas_chm,1_ip,1_ip)
     endif

     do iclas = 1,nclas_chm
        do ipoin = 1,npoin
            conce_communicate(iclas,ipoin,1) = conce(ipoin,iclas,1)
        end do
     end do

     call COU_INTERPOLATE_NODAL_VALUES(icoup,nclas_chm,dummr,conce_communicate)

     call memory_deallo(mem_modul(1:2,modul),'CONCE_COMMUNICATE','chm_plugin',conce_communicate)

  else if( variable == 'MASSK' ) then
     !
     ! Mass sink
     !
     if( .not. associated(mass_sink) ) then
         if (INOTMASTER) then
             call memory_alloca(mem_modul(1:2,modul),'MASS_SINK','chm_plugin',mass_sink,npoin)
         else
             call memory_alloca(mem_modul(1:2,modul),'MASS_SINK','chm_plugin',mass_sink,1_ip)
         endif
     endif
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,mass_sink)


  !
  ! Communications from CMC to CFD
  !

  else if( variable == 'CM2CF' ) then
     !
     ! Send information from CMC to CFD: density, laminar viscosity,
     ! conductivity, cp and temperature
     !
     call chm_data_CMC_to_CFD_domain_CMC(icoup)


  !
  ! Communications from CFD to CMC
  !
  else if( variable == 'CF2CM' ) then
     !
     ! Send information from CFD to CMC: mixing variables and velocity
     !
     call chm_data_CFD_to_CMC_domain_CMC(icoup)

  end if

end subroutine chm_plugin
