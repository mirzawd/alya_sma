!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outerr()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_outerr
  ! NAME
  !    nsi_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use def_coupli,              only : mcoup
  use mod_communications,      only : PAR_MAX
  use mod_outfor,              only : outfor
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_output_postprocess,  only : output_postprocess_cancel_variable_postprocess
  use mod_maths_arrays,        only : maths_findloc
  use mod_messages,            only : messages_live
  use mod_arrays,              only : arrays_number
  implicit none
  integer(ip)    :: ierro=0,iwarn=0,iboun,ipoin,idime,ibopo,ipara
!  integer(ip)    :: iaux,ielty
  integer(ip)    :: dumm1,dumm2,ii
  character(200) :: wmess
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
  ! Check the transient evolution
  !
  if( kfl_timei /= 0 ) then
     if( kfl_timei_nsi == 0 ) then
        iwarn = iwarn + 1
        call messages_live('STEADY NAVIER STOKES EQUATIONS IN A TRANSIENT CALCULATION','WARNING')
     end if
  end if
  !
  ! Coupling with levels
  !
  if( kfl_colev_nsi /= 0 .and. kfl_modul(ID_LEVELS) == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'LEVELS MODULE SHOULD BE PUT TO ON')
  end if
  !
  ! Dirichlet Matrix or Algorithm needed when there are periodic nodes
  !
  if( (kfl_matdi_nsi /= NSI_DIRICHLET_MATRIX .and. kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM) .and. nperi > 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'DIRICLET MATRIX OR ALGORITHM NEEDED WHEN THERE ARE PERIODIC NODES')
  end if
  !
  ! Boussinesq without temperature
  !
  if( kfl_cotem_nsi==1 .and. kfl_modul(ID_TEMPER) == 0 .and. mcoup == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'BOUSSINESQ COUPLING IS IMPOSSIBLE IF TEMPER MODULE IS NOT SOLVED')
  end if
  !
  ! Boussinesq without gravity vector
  !
  if( kfl_cotem_nsi==1 ) then
     if( &
          abs(gravb_nsi(    1))<zensi .and. &
          abs(gravb_nsi(    2))<zensi .and. &
          abs(gravb_nsi(ndime))<zensi ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'BOUSSINESQ COUPLING REQUIRES TO DEFINE THE GRAVITY VECTOR')
     end if
  end if
  !
  ! Low Mach without temperature
  
!  if( kfl_regim_nsi==3 .and. kfl_modul(ID_TEMPER) == 0 ) then
!     ierro = ierro + 1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'LOW MACH REGIME IS IMPOSSIBLE IF TEMPER MODULE IS NOT SOLVED')
!  end if
! not any more true with spray
  
  ! Low Mach with only one Gauss point
  !
  if( kfl_regim_nsi==3 .and. mgaus == 1 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'LOW MACH REGIME WITH 1 GP: NOT CODED')
  end if
  !
  ! Low Mach only developed with KERMOD
  !
  if( kfl_regim_nsi==3 .and. kfl_prope == 0 ) then
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'WARNING: LOW MACH REGIME ONLY TESTED WITH PROPERTIES FROM KERNEL (KERMOD)')
  end if
!  !
!  ! Turbulence without solving TURBUL
!  !
!  if( kfl_cotur_nsi==1 .and. kfl_modul(4) == 0 ) then
!     ierro = ierro + 1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'TURBULENCE COUPLING IS IMPOSSIBLE IF TURBUL MODULE IS NOT SOLVED')
!  end if
  !
  ! Turbulence without convective term
  !
  if( turmu_ker % kfl_exist /= 0_ip .and. kfl_advec_nsi == 0 ) then
     iwarn = iwarn + 1
     kfl_advec_nsi=1
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'TURBULENCE MODELING REQUIRES A CONVECTIVE TERM: IT WAS TURNED ON AUTOMATICALLY')
  end if
  !
  ! Turbulence without viscous term
  !
  if( turmu_ker % kfl_exist /= 0_ip .and. kfl_visco_nsi == 0 ) then
     iwarn = iwarn + 1
     kfl_visco_nsi=1
     fvins_nsi=1.0_rp
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'TURBULENCE MODELING REQUIRES A VISCOUS TERM: IT WAS TURNED ON AUTOMATICALLY')
  end if
  !
  ! Turbulence with Laplacian viscous term
  !
  if( turmu_ker % kfl_exist /= 0_ip .and. fvins_nsi == 0.0_rp ) then
     iwarn = iwarn + 1
     kfl_visco_nsi=1
     !fvins_nsi=1.0_rp
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'TURBULENCE MODELING REQUIRES THE VISCOUS TERM IN DIVERGENCE FORM. SHOULD BE CHANGED IN NSI-DAT FILE')
  end if
  if( turmu_ker % kfl_exist /= 0_ip .and. kfl_grtur_nsi /= 0 ) then
     iwarn = iwarn + 1
     kfl_grtur_nsi=0
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'CANNOT INCLUDE -2/3*rho*K TERM IF TURBULENCE MODELING IS OFF')
  end if
  !
  ! Wall law conditions with zero wall distance
  !
  if( delta_nsi<=zensi .and. INOTMASTER  .and. kfl_delta /= 1 ) then
     boundaries: do iboun=1,nboun
        if( kfl_fixbo_nsi(iboun)==3 ) then
           iwarn = iwarn + 1
           call outfor(2_ip,momod(modul)%lun_outpu,&
                'WALL DISTANCE IS ZERO: WALL LAW CONDITION IS REPLACED BY A SLIP CONDITION')
           exit boundaries
        end if
     end do boundaries
  end if
  !
  ! Time integration scheme and accuracy
  !
  if( kfl_timei_nsi==1 .and. kfl_tisch_nsi==1 .and. kfl_tiacc_nsi>2 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'WRONG TIME INTEGRATION ORDER USING TRAPEZOIDAL RULE. MUST BE 1 OR 2.')
  end if
  if( kfl_timei_nsi==1 .and. kfl_tisch_nsi==2 ) then   ! BDF
     if ( kfl_timco > 1 ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'BDF MAKES NO SENSE WITH LOCAL TIME STEPS')
     end if
     if ( kfl_timco /= 0 .and. kfl_tiacc_nsi > 2 ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'BDF OF ORDER HIGHER THAN 2 ONLY READY FOR PRESCRIBED TIME STEP')
     end if
  end if
  !
  ! Time tracking and integration scheme
  !
  if( kfl_timei_nsi==1 .and. kfl_sgsti_nsi /= 0 .and. kfl_tisch_nsi==2 .and. kfl_sgsac_nsi/=1) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'CANNOT TRACK THE SUBGRID SCALES WITH A BDF SCHEME AND SECOND ORDER ACCURACY')
  end if
  !
  ! Time tracking of the subscales
  !
  if( kfl_timei_nsi == 0 .and. kfl_sgsti_nsi /= 0 ) then
     iwarn = iwarn + 1
     kfl_sgsti_nsi=0
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'CANNOT TRACK THE SUBGRID SCALES IN TIME FOR STATIONARY PROBLEM')
  end if
  !
  ! Convection tracking of the subscales
  !
  if( kfl_advec_nsi == 0 .and. kfl_sgsco_nsi /= 0 ) then
     iwarn = iwarn + 1
     kfl_sgsco_nsi=0
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'CANNOT TRACK THE SUBGRID SCALES IN CONVECTION FOR STOKES PROBLEM')
  end if

  if( kfl_sgsli_nsi == 2 .and.  kfl_stabi_nsi /= 0  ) then

     !call outfor(1_ip,momod(modul)%lun_outpu,&
     !     'NEWTON RAPHSON LINEARIZATION OF SUBGRID SCALES ONLY READY WITH ASGS')
  end if
  if( kfl_sgsli_nsi == 2 .and.  kfl_taust_nsi /= 1  ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'NEWTON RAPHSON LINEARIZATION OF SUBGRID SCALES ONLY READY WITH CODINA TAU')
  end if
  !
  ! Non-inertial boundary conditions
  !
  if( INOTMASTER ) then
     if( kfl_local_nsi==1 ) then
        nodes1: do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if( ibopo /= 0 ) then
              if( kfl_fixno_nsi(1,ipoin)==9 .and. kfl_fixrs_nsi(ipoin) /= 0 ) then
                 ierro = ierro + 1
                 call outfor(1_ip,momod(modul)%lun_outpu,&
                      'CANNOT PRESCRIBE NON-INERTIAL DIRICHLET BOUNDARY CONDITION IN LOCAL AXES.'&
                      //'CHECK NODE '//intost(ipoin))
                 exit nodes1
              end if
           end if
        end do nodes1
     end if
     nodes2: do ipoin=1,npoin
        dumm1=0
        dumm2=0
        do idime=1,ndime
           if( kfl_fixno_nsi(idime,ipoin)==9 ) then
              dumm1=1
           else
              dumm2=1
           end if
        end do
        if( dumm1 /= 0 .and. dumm2 /= 0 ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,&
                'WHEN PRESCRIBING A NON-INERTIAL DIRICHLET BOUNDARY, ALL VELOCITY D.O.F'&
                //char(13)//' MUST HAVE THIS CONDITION. CHECK NODE '//intost(ipoin))
           exit nodes2
        else if( dumm1==1 .and. kfl_conbc_nsi == 0 ) then
           if( kfl_funno_nsi(ipoin) /= 0 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,&
                   'CANNOT APPLY A TIME FUNCTION TO A NON-INERTIAL DIRICHLET BOUNDARY,'&
                   //' CHECK NODE '//intost(ipoin))
              exit nodes2
           end if
        else if( dumm1==1 .and. kfl_fvfua_nsi /= 0 .and. kfl_conbc_nsi==1 ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,&
                'NON-INERTIAL BOUNDARY CONDITIONS WITH A TIME DEPENDENT ROTATION CANNOT BE CONSTANT')
           exit nodes2
        end if
     end do nodes2
  end if
  !
  ! Flow rate
  !
  if( associated(kfl_flow_rate_codes_nsi) ) then
     if( maxval(kfl_flow_rate_codes_nsi(1:mflow_nsi)) > 0 .and. kfl_conbc_nsi == 1 ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'CONSTRAINING FLOW RATE REQUIRES NON-CONSTANT BOUNDARY CONDITIONS')
     end if
  end if
  !
  ! Solver
  !
  if(      (solve(1)%nkryd == 0 .and. solve(1)%kfl_algso==8)&
       .or.(solve(1)%nkryd == 0 .and. solve(1)%kfl_algso==8) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'KRYLOV DIOMENSION MUST BE  > 0 WHEN USING GMRES SOLVER')
  end if
  !
  ! Laplacian preconditioning
  !
  if( (kfl_predi_nsi==1.or.kfl_predi_nsi==2) .and. kfl_timei_nsi == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'LAPLACIAN PRECONDITIONING ONLY FOR TRANSIENT PROBLEM')
  end if
  !
  ! Fast assembly
  !
  if( kfl_assem_nsi == 5 ) then
     if( corio_nsi > zeror ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'FAST ASSEMBLY INCOMPATIBLE WITH ROTATION')
     end if
     if( fvins_nsi==2.0_rp ) then
        iwarn = iwarn + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'FAST ASSEMBLY INCOMPATIBLE WITH COMPLETE FORM OF VISCOUS TERM.'//&
             ' SWITCHED TO DIVERGENCE FORM AUTOMATICALLY')
     end if
  end if
  !
  ! Stabilization parameters
  !
  do ipara = 1,size(staco_nsi,KIND=ip)
     if( staco_nsi(ipara) /= 0.0_rp .and. staco_nsi(ipara) /= 1.0_rp  ) then
        iwarn = iwarn+1
        wmess = 'CHECK NAVIER-STOKES STABILIZATION PARAMETERS'
        call outfor(2_ip,momod(modul)%lun_outpu,trim(wmess))
     end if
  end do
  !
  ! Comfort postprocess
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('PMV  ')) .or. output_postprocess_check_variable_postprocess(arrays_number('PPD  ')) ) then
     if(&
          cloth_nsi == -500.0_rp .or. &
          metab_nsi == -500.0_rp .or. &
          wetme_nsi == -500.0_rp .or. &
          ambie_nsi == -500.0_rp .or. &
          radia_nsi == -500.0_rp .or. &
          relat_nsi == -500.0_rp ) then
        ierro = ierro+1
        wmess = 'COMFORT CONSTANT WERE NOT DEFINED: CANNOT POSTPROCESS PMV NOR PPD'
        call outfor(1_ip,momod(modul)%lun_outpu,trim(wmess))
     end if
  end if
  !
  ! Schur complement solver
  !
  if( NSI_SCHUR_COMPLEMENT .and. kfl_sosch_nsi == 1  ) then
     iwarn = iwarn + 1
     kfl_sosch_nsi = 2
     call outfor(2_ip,momod(modul)%lun_outpu,'ORTHOMIN(1) IS EITHER MOMENTUM OR CONTINUITY PRESERVING')
  end if
  !
  ! Boundary conditions
  !
  if( kfl_matdi_nsi == 1 .and. NSI_MONOLITHIC ) then
     iwarn = iwarn + 1
     kfl_matdi_nsi = 0
     call outfor(2_ip,momod(modul)%lun_outpu,'DIRICHLET ONLY ON ELEMENT WHEN DING MONOLITHIC')
  end if

  !if (nodpr_nsi > npoin .and. INOTMASTER ) then
  !   ierro = ierro + 1
  !   call outfor(2_ip,momod(modul)%lun_outpu,'PRESSURE FIXED ON A NON-EXISTING NODE')
  !endif
  !
  ! Fast assembly+monolithic
  !
  !if(kfl_assem_nsi==1.and.NSI_MONOLITHIC ) then
  !   ierro = ierro + 1
  !   wmess='FAST ASSEMBLY NOT COMPATIBLE WITH MONOLITHIC APPROACH'
  !   call outfor(1_ip,momod(modul)%lun_outpu,trim(wmess))
  !end if
  if( kfl_normc_nsi == 4 .and. kfl_refer_nsi == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'REFERENCE RESIDUAL NORM ONLY POSSIBLE IF REFERENCE SOLUTION IS OUTPUT')
  end if
  !
  ! Internal force
  !
  if( kfl_intfo_nsi == 1 .and. kfl_matdi_nsi == 0 ) then
     ierro = ierro + 1
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'REAIDUAL INTERNAL FORCE ONLY IF DIRICHLET BC ARE IMPOSED DIRECTLY IN MATRICES')
  end if
  !
  ! Decoupled SGS and OSS-like stabilizations
  !
  if( kfl_sgscp_nsi == 1 .and. kfl_stabi_nsi /= 0 ) then
     ierro = ierro + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'OSS-LIKE STABILIZATION INCOMPATBILE WITH COUPLED SGS')
  end if
  !
  ! Properties
  !
  if( kfl_prope == 0 ) then
     ierro = ierro + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'PROPERTIES SHOULD BE DECALRED IN KERMOD')
  end if
  !
  ! Hydrostatic pressure
  !
  if( kfl_inipr_nsi == 2 .and. kfl_hydro_nsi == 0 ) then
     ierro = ierro + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'HYDROSTATIC PRESSURE MUST BE COMPUTED TO IMPOSE IT AS INITIAL SOLUTION')
  end if
  !
  ! Element cut + tracking
  !
  if( kfl_cutel == 1 .and. ( kfl_sgsco_nsi /= 0 .or. kfl_sgsti_nsi /= 0 ) ) then
     ierro = ierro + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'ELEMENT CUT NOT COMPATIBLE WITH SUBGRID SCALE TRACKING')
  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess
  !
  !----------------------------------------------------------------------
  !
  ! DCG groups
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('GROUP')) .and. solve(2)%kfl_algso /= 2 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('GROUP'))
     iwarn = iwarn+1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS GROUPS IF DEFLATED CG SOLVER WAS NOT CHOSEN')
  end if
  if( output_postprocess_check_variable_postprocess(arrays_number('LINEL')) .and. solve(2)%kfl_preco /= 4 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('LINEL'))
     iwarn = iwarn+1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS LINELET IF LINELET PRECONDITIONER WAS NOT CHOSEN')
  end if
  !
  ! Streamlines
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('STREA')) .and.ndime /= 2 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('STREA'))
     iwarn = iwarn+1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS STREAMLINES IN 3D')
  end if
  !
  ! Projection
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('VEPRO')) .and. kfl_stabi_nsi == 0 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('VEPRO'))
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS PROJECTION')
  end if
  if( output_postprocess_check_variable_postprocess(arrays_number('PRPRO')) .and. kfl_stabi_nsi == 0 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('PRPRO'))
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS PROJECTION')
  end if
  !
  ! For the cases without no slip wall what should be outputed here has not been thought nor programmed
  !
  if( kfl_noslw_ker == 0 ) then
     if( output_postprocess_check_variable_postprocess(arrays_number('AVVAF')) .or. &
         output_postprocess_check_variable_postprocess(arrays_number('AVNTR')) .or. &
         output_postprocess_check_variable_postprocess(arrays_number('AVGTR')) ) then
!        output_postprocess_check_variable_postprocess(arrays_number('FANSW')) ) then
        ierro = ierro + 1
        call outfor(2_ip,momod(modul)%lun_outpu,'AVNTR, AVGTR, VAFOR, AND AVVAF CAN ONLY BE POSTPROCESSED WITH NO SLIP WALL LAW')
     end if
  end if
  !
  ! For the cases without exchange location, or two layer, or no slip wall what should be outputed here has not been thought nor programmed
  !
  if( kfl_waexl_ker == 0 .and. kfl_twola_ker == 0 .and. kfl_noslw_ker == 0 ) then
     if( output_postprocess_check_variable_postprocess(arrays_number('NOTRA'))) then
        ierro = ierro + 1
        call outfor(2_ip,momod(modul)%lun_outpu,'NOTRA CAN ONLY BE POSTPROCESSED WITH EXCHANGE LOCATION, OR TWO LAYER, OR NO SLIP WALL LAW')
     end if
  end if
  !
  ! For exchange location (no slip walllaw included) elsest must be ON
  !
  if( kfl_waexl_ker /= 0_ip .and. kfl_elses == 0_ip ) then
     if( output_postprocess_check_variable_postprocess(arrays_number('NOTRA'))) then
        ierro = ierro + 1
        call outfor(2_ip,momod(modul)%lun_outpu,'For exchange location (no slip walllaw included) elsest must be ON')
     end if
  end if
  !
  ! Internal force
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('INTFO')) .and. kfl_intfo_nsi == 0 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('INTFO'))
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS INTERNAL FORCE')
  end if

  !
  ! Consistent mass matrix for end of step mass correction when there are periodic nodes
  ! To do: add call solver_preprocess in nsi matrix  & probably solver_postprocess in nsi_solite
  !
  if( kfl_corre_nsi == 3 .and. nperi > 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'PERIODIC NODES NOT READY WITH CONSISTENT MASS MATRIX FOR MASS CORRECTION')
  end if
  !
  ! Consistent mass matrix & DIRICHLET MATRIX
  ! To do: should no be difficult to add
  !
  if( kfl_corre_nsi == 3 .and. kfl_matdi_nsi /= 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'DIRICHLET MATRIX NOT READY WITH CONSISTENT MASS MATRIX FOR MASS CORRECTION')
  end if
  !
  ! Consistent mass matrix & extension elements
  ! MIISING ADD RUNEND . I COULD NOT FIND THE FLAG FOR EXTENSION ELEMENTS
  !
  !
  ! kfl_ini_ts_guess_order_nsi for the moment only valid up to 2
  !
  if( kfl_ini_ts_guess_order_nsi > 2 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'kfl_ini_ts_guess_order_nsi > 2 not ready')
  end if

  !
  ! INCLUDE DT IN TAU Does not make sense with tracking
  !
  if( (kfl_taust_nsi == 5.or.kfl_taust_nsi == 6) .and. ( kfl_sgsti_nsi/= 0 .or. kfl_sgscp_nsi /= 0) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'INCLUDE DT IN TAU Does not make sense with tracking' )
  end if

  !
  ! INCLUDE DT IN TAU not sure if it makes sense with Local time step
  !
  if( (kfl_taust_nsi == 5.or.kfl_taust_nsi == 6) .and. kfl_timco == 2 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'INCLUDE DT IN TAU - not sure if it makes sense with Local time step' )
  end if

  !
  ! INCLUDE DT IN TAU not sure if it makes sense ORTOMIN TAU
  !
  if( (kfl_taust_nsi == 5.or.kfl_taust_nsi == 6) .and. kfl_predi_nsi==3 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'INCLUDE DT IN TAU not sure if it makes sense ORTOMIN TAU' )
  end if
  !
  ! Pressure should be prescribed to 0
  !
  if( NSI_FRACTIONAL_STEP .and. gamma_nsi /= 1.0_rp .and. kfl_confi_nsi >= 0 .and. &
       nodpr_nsi > 0 .and. valpr_nsi /= 0.0_rp .AND. kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
 !    ierro = ierro + 1
 !    call outfor(1_ip,momod(modul)%lun_outpu,&
 !         'PRESSURE SHOULD BE PRESCRIBED TO ZERO' )
  end if
  !
  ! Impossible combination
  !
  if(  NSI_FRACTIONAL_STEP .and. kfl_stabi_nsi /= NSI_GALERKIN .and. &
       gamma_nsi /= 1.0_rp .and. kfl_nota1_nsi == 0 .and. staco_nsi(4) < 1.0e-8_rp ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'WHEN STABILIZING FRACTIONAL STEP, GAMMA SHOULD BE 1' )
  end if
  !
  ! Fraction step stabilized with TAU
  !
  if( kfl_press_stab_nsi >= 1 .and. gamma_nsi /= 0.0_rp ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'TAU STABILIZATION ONLY POSSIBLE WITH GAMMA=0' )
  end if
  !
  ! Galerking and SGS tracking impossible
  !
  if(  ( kfl_sgsti_nsi /= 0 .or. kfl_sgsco_nsi /= 0 ) .and. &
       ( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'ALGEBRAIC SPLIT OSS IMPOSSIBLE WITH SUBGRID SCALE TRACKING' )
  end if
  !
  ! ALgebraic split and pressure Schur complement preconditioner
  !
  if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS .and. &
       ( kfl_predi_nsi /= 7 .and. kfl_predi_nsi /= 8 .and. kfl_predi_nsi /= 9 ) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'ALGEBRAIC SPLIT OSS IMPOSSIBLE WITH SELECT PRESSURE SCHUR COMPLEMENT PRECONDITIONER' )
  end if
  !
  ! Galerking and limiter impossible
  !
  if( kfl_stabi_nsi == -1 .and. kfl_limit_nsi /= 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'GALERKIN IMPOSSIBLE WITH LIMITER' )
  end if
  !
  ! Enrichement not possible with Galerkin and without penalization
  !
  if( kfl_bubbl_nsi /= 0 .and. kfl_stabi_nsi == NSI_GALERKIN .and. kfl_penal_nsi == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'ENRICHEMENT IMPOSSIBLE WITH GALERKIN WITHOUT PENALIZATION' )
  end if
  !
  ! Fractional step
  !
  if( kfl_algor_nsi == 2 .and. (kfl_tisch_nsi /= 3 .and. kfl_tisch_nsi /= 4)) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'FRACTIONAL STEP SHOULD USE ADAMS-BASHOFRTH OR RUNGE-KUTTA TIME INTEGRATION SCHEME' )
  end if
  !
  ! Fractional step and Dirichlet conditions
  !
  if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM .and. NSI_FRACTIONAL_STEP ) then
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'MAYBE YOU SHOULD USE THE OPTION DIRICHLET: ALGORITHM FOR EFFICIENCY REASONS')
  end if
  !if( kfl_grad_div_nsi /= 0 .and. kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
  !   ierro = ierro + 1
  !   call outfor(1_ip,momod(modul)%lun_outpu,&
  !        'GRAD-DIV FORMULATION MUST BE USED WITH ALGORITHM DIRICHLET OPTION' )
  !end if
  if( .not. NSI_FRACTIONAL_STEP .and. kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'ALGORITHMIC DIRICHLET CONDITIONS DOES NOT WORK WITHOUT FRACTIONAL STEP')
  end if
  if( NSI_FRACTIONAL_STEP .and. kfl_massm_nsi == NSI_CONSISTENT_MASS .and. kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
    ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'YOU SHOULD USE DIRICHLET ALGORITHM FOR CONSISTENT MASS MATRIX IN NSI_FRACTIONAL_STEP METHOD')
  end if
  if( NSI_FRACTIONAL_STEP .and. kfl_massm_nsi == NSI_CONSISTENT_MASS .and. kfl_conma_weighted == 0 ) then
    ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'YOU SHOULD ACTIVE THE WEIGHTED MASS MATRIX OPTION IN KER.DAT FILE FOR CONSISTENT MASS MATRIX IN NSI_FRACTIONAL_STEP METHOD')
  end if
! if( kfl_grad_div_nsi /= 0 ) then
!    dumm1 = 0
!    if( INOTMASTER ) then
!       !if( kfl_regim_nsi == 3 .or. bemol_nsi > 0.0_rp ) dumm1 = 1    ! LOWMA(3) -- bemol_nsi - Integration of convective term by parts
!       do iboun = 1,nboun                                            ! CODES BOUNDARIES
!          if( kfl_fixbo_nsi(iboun) > 0 .and. (kfl_fixbo_nsi(iboun) /= 3 ) ) dumm1 = 1
!       end do
!    end if
!    call PAR_MAX(dumm1)
!    if( dumm1 > 0 ) then
!       ierro = ierro + 1
!       call outfor(1_ip,momod(modul)%lun_outpu,&
!            'GRAD-DIV FORMULATION IS INCOMPATIBLE WITH IMPLICIT NATURAL BOUDNARY CONDITION' )
!    end if
! end if
  !
  ! GPU assembly -- nsi_element_operations_fast
  !
!!$  if( kfl_assem_nsi == 5 ) then    
!!$     if( kfl_grad_div_nsi == 0 ) then
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly only ready with GRAD-DIV FORMULATION' )
!!$     end if
!!$     if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly only ready with NSI_DIRICHLET_ALGORITHM' )
!!$     end if
!!$     if( kfl_convection_type_nsi /= NSI_CONVECTION_EMAC ) then
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly only ready with EMAC scheme' )
!!$     end if
!!$     if( abs(fvins_nsi - 1.0_rp) > 1e-7_rp ) then 
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly only ready DIVERGENCE FORM' )
!!$     end if
!!$     if( kfl_timco == 2 ) then
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly not ready for LOCAL time step' )
!!$     end if
!!$     !
!!$     ! Only linear elements for the moment
!!$     ! For the moment I have only inlined the part corresponding to porde = 1 from nsi_rhodt_rhotau_nu_vector in nsi...fast
!!$     ! The idea is to minimize ifs but this one perhaps could be allowed
!!$     !
!!$     iaux = 0
!!$     do ielty = 1,nelty
!!$        if( lexis(ielty) == 1 ) then
!!$           if ( lorde(ielty) /= 1 ) iaux = iaux+1
!!$        end if
!!$     end do
!!$     if( iaux /=0 ) then
!!$        ierro = ierro + 1
!!$        call outfor(1_ip,momod(modul)%lun_outpu,&
!!$             'GPU assembly not ready for non linear elements - could be easilly addapted with an if' )
!!$     end if
!!$  end if
  !
  ! kfl_grad_div/=0 is only ready for NSI_FRACTIONAL_STEP.
  ! In some special cases it could be used (weighted laplacian... Q = - tau * L  & Split OSS)
  ! but several files need to be adapted, for example:   in mod_ndi_element_operations L909
  ! we would need to add an if kfl_grad_div as is done on fractional step case so that lapla_nsi is not calculated.
  !  
  !if ( kfl_grad_div_nsi/=0 .and. ( .not. NSI_FRACTIONAL_STEP) ) then
  !   ierro = ierro + 1
  !   call outfor(1_ip,momod(modul)%lun_outpu,&
  !        'kfl_grad_div/=0 is only ready for NSI_FRACTIONAL_STEP ')
  !end if
  !
  ! GPU and GPU2 assemblies only work with GRAD_DIV
  !
  if ( kfl_grad_div_nsi == 0_ip .and. ( kfl_assem_nsi == 5_ip ) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'GPU and GPU2 assemblies only work with GRAD_DIV ')
  end if  
  !
  ! Stop by witness not ready with restart - extra values would ned to be saved
  !
!  if ( kfl_stop_by_wit_nsi == 1 .and. kfl_rstar /= 0 ) then
!    ierro = ierro + 1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'Stop by witness not ready with restart')
!  end if


  !
  ! Exchange location and no slip wall can not be used at the same time
  ! Now I have actually set wall exchange when no slip wall 
  !
!  if ((kfl_noslw_ker/=0) .and.(kfl_waexl_ker/=0)) then
!    ierro = ierro + 1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'Exchange location and no slip wall can not be used at the same time')   ! later we might do nsw that uses velocity further away from the wall
!     ! or have some boundaries with exchage location an others with NSWL -- nsi_memall L458 (if(kfl_waexl_ker==0))  would need to be changed
  !  end if

  !
  ! No slip wall needs INTERNAL_FORCES:      RESIDUAL
  !
  if ( kfl_noslw_ker /= 0 .and. kfl_intfo_nsi /= 1) then
    ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'No slip wall needs INTERNAL_FORCES:      RESIDUAL')   
  end if
  !
  ! No slip wall needs INTERNAL_FORCES:      RESIDUAL
  !
  if ( kfl_asbou_nsi == 5_ip .and. kfl_ustar /= 0_ip ) then
    ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'BOUNDARY ASSEMBLY GPU2 is only ready for smooth wall law')
     ! else rougness would need to be defined of size (DEF_VECT,pgaub) and calculated accordingly in mod_nsi_boundary_operations_fast
     ! for the moment I have just eliminated the need of using an include and called frivel directly - much cleaner
  end if
  !
  ! Correction to fractional step when gravity forces are present is only ready with fractional step or semi-implicit scheme  -- see waht happens in the case viscous is treated implicitly
  !
  if ( ( kfl_algor_nsi /= 2_ip ) .and.  kfl_fsgrb_nsi == 1_ip ) then
    ierro = ierro + 1
    call outfor(1_ip,momod(modul)%lun_outpu,&
          'Correction to fractional step when gravity forces are present is only ready with explicit fractional step - semi implicit needs to be implemented')
  end if
  !
  ! Pressure level
  !
  if( kfl_press_lev_nsi /= 0 ) then
     ii = maths_findloc(lbsec,kfl_press_lev_nsi)
     if( ii == 0 ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'NSI_UPDBCS: WRONG SET TO COMPUTE PRESSURE LEVEL')
     end if
     if( postp(1) % npp_setsb(1) == 0 ) then
        iwarn = iwarn + 1
        call messages_live('ACTIVATING AUTOMATICALLY MEAN PRESSURE ON BOUNDARY SETS TO LEVEL PRESSURE','WARNING')
     end if
  end if
  
  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,'NULL')

end subroutine nsi_outerr
