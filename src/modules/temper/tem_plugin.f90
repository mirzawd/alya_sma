!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_plugin.f90
!> @author  Guillaume Houzeaux
!> @date    13/04/2014
!> @brief   Plugin with coupling
!> @details Plugin for the zonal coupling
!> @} 
!-----------------------------------------------------------------------

subroutine tem_plugin(icoup)
  !
  ! Obligatory variables 
  !
  use def_kintyp,         only :  ip,rp
  use def_coupli,         only :  coupling_type
  use def_master,         only :  solve_sol,modul, INOTMASTER
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
  use mod_memory,         only :  memory_deallo
  use mod_memory,         only :  memory_alloca
  use mod_matrix,         only :  matrix_initialize
  use def_master,         only :  mem_modul,modul
  use def_domain,         only :  npoin
  !
  ! Possible variables 
  !
  use def_master,         only :  therm
  use def_master,         only :  advec
  use def_master,         only :  heat_sink
  use def_temper,         only :  bvess_tem

  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_master,         only : inotmaster
  use mod_parall
  use def_temper,         only : kfl_fixno_tem
  use def_domain
  implicit none 
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp),    pointer    :: dummr(:)
  real(rp),    pointer    :: dumm3(:,:,:)
  integer(ip), save       :: ipass=0

  nullify(dummr)
  
  variable = coupling_type(icoup) % variable 

  if( variable == 'TEMPE' .or. variable == 'ENTHA'  ) then   
     !
     ! Temperature 
     !
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('TEMPE',modul,kfl_fixno_tem)
     end if
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,bvess_tem,therm,kfl_fixno_tem)      
     
  else if( variable == 'RESID' ) then
     !
     ! Residual
     !
     call matrix_initialize(solve_sol(1) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,solve_sol(1) % bvnat,solve_sol(1) % reaction,kfl_fixno_tem)
    
  else if( variable == 'VELOC' ) then
     !
     ! Velocity
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,advec,dumm3)
    
  else if( variable == 'FAKE ' ) then
     !
     ! This variable is defined to allow for the use of the coupling tools
     !
  else if( variable == 'HEASK' ) then
     !
     ! Heat sink
     ! 
     if( .not. associated(heat_sink) ) then 
         if (INOTMASTER) then 
             call memory_alloca(mem_modul(1:2,modul),'HEAT_SINK','tem_plugin',heat_sink,npoin)
         else
             call memory_alloca(mem_modul(1:2,modul),'HEAT_SINK','tem_plugin',heat_sink,1_ip)
         endif
     endif
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,heat_sink)
  end if

end subroutine tem_plugin

