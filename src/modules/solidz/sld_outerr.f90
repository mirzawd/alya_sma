!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outerr.f90
!> @author  Guillaume Houzeaux
!> @date    
!> @brief   This routine checks if there are errors and warnings
!> @details This routine checks if there are errors and warnings
!>    USED BY
!>       sld_turnon
!> @}
!-----------------------------------------------------------------------

subroutine sld_outerr()
  
  use def_master,                only : INOTMASTER
  use def_master,                only : postp
  use def_master,                only : zeror
  use def_kermod,                only : kfl_cutel
  use def_domain,                only : xfiel
  use def_domain,                only : nbset,neset,nnset
  use def_domain,                only : nmate
  use mod_messages,              only : messages_live
  use def_solidz
  use mod_sld_interface_element, only : ELINT_LAW_TURON
  use mod_sld_interface_element, only : ELINT_LAW_TURON_2018
  use mod_communications,        only : PAR_MAX
  
  implicit none

  integer(ip)    :: ierro=0,iwarn=0,jerro
  integer(ip)    :: imate
  
  iwarn = 0
  ierro = 0
  !
  ! Total number of materials
  !
  if( nmate == 0 ) then
     ierro = ierro + 1
     call messages_live('MATERIALS MUST BE DECLARED EXPLICITLY IN THE DOMAIN (*.dom.dat FILE)','ERROR')     
  end if
  !
  ! Plane stress/strain conditions
  !
  if( kfl_plane_sld == 1_ip ) then
     do imate = 1,nmate_sld
        if( lawst_sld(imate) .ne. 100_ip ) then
           ierro = ierro + 1
           call messages_live("PLANE STRESS CONDITION IS ONLY PROGRAMMED FOR ISOLIN MODEL")
        end if
     end do
  end if
  !
  ! Material density 
  !
  if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then ! Dynamic problem
     if(      kfl_rigid_sld == 0_ip ) then
        do imate = 1,nmate_sld
           if( densi_sld(1,imate) < zeror ) then
              ierro = ierro + 1
              call messages_live(&
                   'DENSITY VALUE IS NEGATIVE, NULL OR IT IS NOT DEFINED. IT IS MANDATORY FOR DYNAMIC PROBLEMS','ERROR')
           end if
        end do
     else if( kfl_rigid_sld == 1_ip ) then
        do imate = 1,nmate_sld
           if( densi_sld(1,imate) > zeror ) then
              iwarn = iwarn + 1
              call messages_live('DENSITY IS NOT USED FOR RIGID BODIES','WARNING')
           end if
        end do
     end if
  else                                            ! Static problem (density not required)
     do imate = 1,nmate_sld
        if( densi_sld(1,imate) > zeror ) then
           iwarn = iwarn + 1
           if(      kfl_rigid_sld == 0_ip) then
              call messages_live('DENSITY IS NOT USED FOR STATIC PROBLEMS USING DEFORMABLE BODIES','WARNING')
           else if( kfl_rigid_sld == 1_ip ) then
              call messages_live('DENSITY IS NOT USED FOR RIGID BODIES','WARNING')
           end if
        end if
     end do
  end if
  !
  ! Material properties
  !
  do imate = 1,nmate_sld
     if(      lawst_sld(imate) == 100_ip ) then
        if( parco_sld(1,imate) <= zeror ) then
           ierro = ierro + 1
           call messages_live('ISOTROPIC ELASTIC MATERIAL MODEL CANNOT HAVE NEGATIVE OR NULL YOUNG MODULUS','ERROR')
        end if
     else if( lawst_sld(imate) == 151_ip ) then
        if( any(parco_sld(1:3,imate) <= zeror) ) then
           ierro = ierro + 1
           call messages_live('ORTHOTROPIC ELASTIC MATERIAL MODEL CANNOT HAVE NEGATIVE OR NULL YOUNG MODULI','ERROR')
        end if
     else if( lawst_sld(imate) == 154_ip ) then
        if( any(parco_sld(1:2,imate) <= zeror) ) then
           ierro = ierro + 1
           call messages_live('BESSA MATERIAL MODEL CANNOT HAVE NEGATIVE OR NULL YOUNG MODULI','ERROR')
        end if
     else if( lawch_sld(imate) == ELINT_LAW_TURON ) then
        if( any(parch_sld(1:5,imate) <= zeror) ) then
           ierro = ierro + 1
           call messages_live(&
                'GIC,GIIC,TAUI,KP and ETA VALUES FROM COHESIVE MODEL MUST BE GREATER OR EQUAL THAN ZERO','ERROR')
        end if
     else if( lawch_sld(imate) == ELINT_LAW_TURON_2018 ) then
        if( any(parch_sld(1:6,imate) <= zeror) ) then
           ierro = ierro + 1
           call messages_live(&
                'GIC,GIIC,TAUI,TAUII,ETA and KP VALUES FROM COHESIVE MODEL MUST BE GREATER OR EQUAL THAN ZERO','ERROR')
        end if
     end if
  end do
  !
  ! Material orientation at element level
  !
  if( kfl_fiber_sld < 4_ip ) then
     do imate=1,nmate_sld
        if( lawst_sld(imate) == 151_ip .or. &
            lawst_sld(imate) == 152_ip .or. &
            lawst_sld(imate) == 154_ip .or. &
           (lawst_sld(imate) == 200_ip .and. lawco_sld(imate) == 2_ip) ) then
           ierro = ierro + 1
           call messages_live('ORTHOTROPIC MATERIAL MODELS REQUIRE MATERIAL ORIENTATIONS','ERROR')
        end if
     end do
  end if
  !
  ! Initial condition for velocity
  !
  if( kfl_timei_sld == SLD_STATIC_PROBLEM .and. kfl_invel_sld > 0_ip ) then
     iwarn = iwarn + 1
     call messages_live('INITIAL CONDITION FOR VELOCITY IS NOT USED FOR STATIC PROBLEMS','WARNING')
  end if
  !
  ! Check initial conditions for state dependent variables 
  !
  if( kfl_sdvar_sld == 0_ip .and. kfl_insdv_sld > 0_ip) then
     iwarn = iwarn + 1
     call messages_live('INITIAL CONDITION FOR STATE DEPENDENT VARIABLES IS ONLY FOR MATERIALS WITH SDVAR','WARNING')
  end if
  !
  ! Set postprocess
  !
  if( neset == 0 .and. maxval(postp(1) % npp_setse) > 0 ) then
     iwarn = iwarn + 1
     call messages_live('ELEMENT SETS HAVE NOT BEEN DEFINED: CANNOT OUTPUT ELEMENT SET RESULTS','WARNING')
  end if
  if( nbset == 0 .and. maxval(postp(1) % npp_setsb) > 0 ) then
     iwarn = iwarn + 1
     call messages_live('BOUNDARY SETS HAVE NOT BEEN DEFINED: CANNOT OUTPUT BOUNDARY SET RESULTS','WARNING')
  end if
  if( nnset == 0 .and. maxval(postp(1) % npp_setsn) > 0 ) then
     iwarn = iwarn + 1
     call messages_live('NODE SETS HAVE NOT BEEN DEFINED: CANNOT OUTPUT NODE SET RESULTS','WARNING')
  end if
  !
  ! Time step strategy
  !
  if( kfl_savdt_sld == 1 .and. kfl_celen_sld == 1 ) then
     ierro = ierro + 1
     call messages_live('SAVE TIME STEP OPTION IS INCOMPATIBLE WITH THE MIN_LENGHT USING UPDATED','ERROR')
  end if
  !
  ! X-FEM
  !
  if( kfl_cutel == 0 .and. kfl_xfeme_sld == 1 ) then
     ierro = ierro+1
     call messages_live('PUT CUT ELEMENTS ON IN THE KERNAL DATA FILE WHEN USING X-FEM','ERROR')
  end if
  !
  ! Volume force
  !
  if( kfl_vofor_sld > 0 ) then
     jerro = 0
     if( INOTMASTER ) then
        if( size(xfiel) < kfl_vofor_sld ) then
           jerro = 1
        else if( .not. associated(xfiel(kfl_vofor_sld) % a)) then
           jerro = 1
        end if
     end if
     call PAR_MAX(jerro)
     if( jerro == 1 ) then
        ierro = ierro + 1
        call messages_live('EXTERNAL VOLUME FORCE FIELD DOES NOT EXIST','ERROR')
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,'NULL')

end subroutine sld_outerr
