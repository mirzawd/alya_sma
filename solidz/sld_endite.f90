!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_endite.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Convergence checks and unknowns update
!> @details
!>          \verbatim
!>          ITASK = 1 The end of an internal iteration
!>                = 2 The end of an external iteration
!>          \endverbatim
!> @todo
!>           <GGU> Revision sld_coupli (obsolete for nastin and solidz)
!> @}
!------------------------------------------------------------------------

subroutine sld_endite(itask)

  use def_kintyp,            only : ip, rp
  use def_master,            only : itinn, modul, intost, INOTSLAVE, INOTMASTER
  use def_master,            only : ITASK_ENDINN, ITASK_ENDITE, TIME_N
  use def_master,            only : solve_sol
  use def_master,            only : ittim
  use def_domain,            only : kfl_elcoh
  use mod_messages,          only : livinf
  use def_solidz,            only : miinn_sld
  use def_solidz,            only : kfl_xfeme_sld, last_iters_sld
  use def_solidz,            only : kfl_psmat_sld
  use mod_sld_energy,        only : sld_updene
  use mod_sld_fe2,           only : fe2_update_vars, fe2_write_profiling
  use mod_sld_cardiac_cycle, only : kfl_cardiac_cycle
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_compute_cavity_volume
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_manage_cavity_pressure
  use def_solidz,            only : lun_sysnet_heart_res_sld, lun_sysnet_system_res_sld
  use mod_communications,    only : PAR_BROADCAST
#if defined COMMDOM && COMMDOM==2
  use def_solidz,                only : kfl_rigid_sld,kfl_conta_sld, kfl_contf_sld
  use mod_sld_pdn_contact_plepp, only : kfl_pdnco_sld,kfl_pdncf_sld
  use mod_sld_pdn_contact_plepp, only : sld_pdn_contact_cvgunk
#endif
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet,        only : kfl_sysnet, sysnet_compute_volumes, sysnet_compute_pvloop, sysnet_write_alya_res
  use mod_ker_sysnet,        only : sysnet_alya_couple, sysnet_update_phase,cavities_sysnet, iphase_sysnet
#endif

  implicit none

  integer(ip), intent(in) :: itask   !< 1, inner iterations; 2, outer iterations
  integer(ip)             :: maxiter
  character(300)          :: messa

  !
  ! Cardiac Cycle ON: Calculate the volume inside the cavities and windkessel pressure if needed
  !
  if( kfl_cardiac_cycle )then
     call sld_cardiac_cycle_compute_cavity_volume(itask)
     call sld_cardiac_cycle_manage_cavity_pressure(itask)
  end if

  select case(itask)

  case( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the internal iteration
     !
     !-------------------------------------------------------------------

     call sld_cvgunk(ITASK_ENDINN)   ! Residual
     call sld_updunk(ITASK_ENDINN)   ! Update: u(n,i,j-1) <-- u(n,i,j)
     !
     ! Live info iterations
     !
     maxiter = solve_sol(1) % miter
     messa = &
          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_sld))//&
          ' IT: '//trim(intost(last_iters_sld))//'/'//trim(intost(maxiter))//')'
     call livinf(-3_ip,messa,1_ip)
     call livinf(56_ip,' ',modul)
     !
     ! Convergence and Timings
     !
     call sld_cvgunk(0_ip)
     !
     ! Convergence of PDN-contact algorithm
     !
#if defined COMMDOM && COMMDOM==2
     if( (kfl_pdnco_sld > 0 .or. kfl_conta_sld == 3_ip) .and. &
         (kfl_pdncf_sld .or. kfl_contf_sld == 1) .and. kfl_rigid_sld == 0_ip ) then
        call sld_pdn_contact_cvgunk()
     end if
#endif
     !
     ! Matrix output
     !
     if( kfl_psmat_sld == 1_ip ) call sld_outmat()

  case( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the external iteration
     !
     !-------------------------------------------------------------------

     call livinf(16_ip,' ',itinn(modul))

     call sld_cvgunk(ITASK_ENDITE)  ! Residual || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||)
     call sld_updunk(ITASK_ENDITE)  ! Update unknowns: u(n,i-1,*) <-- u(n,i,*)
     call sld_updene(ITASK_ENDITE)  ! Update energies
     !
     ! Coupling
     !
     call sld_coupli(ITASK_ENDITE)
     !
     ! Cohesive elements and X-FEM
     !
     if (kfl_elcoh > 0 .or. kfl_xfeme_sld > 0) call sld_updcoh(2_ip)
     !
     ! Write MicroPP profiling
     !
     call fe2_write_profiling(ittim)
     !
     ! Update MicroPP internal variables
     !
     call fe2_update_vars()

     !
     ! Sysnet ON: Calculate the volume inside the cavities and windkessel pressure if needed (only at the final iteration)
     !
#ifndef PROPER_ELEM_PRIVATE_OFF
     if( kfl_sysnet)then
        ! compute volumes
        call sysnet_compute_volumes(ITASK_ENDITE)
        ! update sysnet volumes coming from alya
        call sysnet_alya_couple(1_ip)
        ! update phase
        call sysnet_update_phase(ITASK_ENDITE)
        !
        ! During PRELOAD phase, sysnet coupling is bridged, storing the funcar and unksys values coming
        ! from the startup phase
        !
        if (adjustl(trim(iphase_sysnet(1))) == "pre") then    
           !
           ! Preload phase, write alya_res files
           !
        else 
           !
           ! Regular coupled phase
           !
           if (INOTSLAVE) then
              ! update pvloop (only the master)
              call sysnet_compute_pvloop(ITASK_ENDITE)
              ! update alya pressures coming from sysnet, converting from mm/hg to CGS
              call sysnet_alya_couple(2_ip)
           end if
           ! broadcast new pressures to the slaves
           if (INOTMASTER) then
              cavities_sysnet(1) % pres(TIME_N) = 0.0_rp
              cavities_sysnet(2) % pres(TIME_N) = 0.0_rp
           end if
           call PAR_BROADCAST(cavities_sysnet(1) % pres(TIME_N))
           call PAR_BROADCAST(cavities_sysnet(2) % pres(TIME_N))
           messa = &
                '    SOLIDZ CALL SYSNET: COMPUTE CYCLE'
           call livinf(0_ip,messa,1_ip)
        end if
     end if
#endif

  end select

end subroutine sld_endite
