!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outerr()
!------------------------------------------------------------------------
!****f* Chemic/chm_outerr
! NAME
!    chm_outerr
! DESCRIPTION
!    This routine checks if there are erros and warnings
! USES
! USED BY
!    chm_turnon
!***
!------------------------------------------------------------------------
  use def_master,       only : momod, kfl_timei, intost, modul
  use def_kintyp,       only : ip, rp
  use def_chemic,       only : mixedEq_eqs_chm, mixedEq_groups_chm, kfl_model_chm, kfl_premix_chm, kfl_tab_fw_chm, kfl_timei_chm,&
                               kfl_ufpv_chm, kfl_varYc_chm, kfl_varZ_chm, nclas_chm, ngrou_chm, diffu_chm, table_fw, &
                               kfl_tab_fw_chm_diff, kfl_multimod_chm
  use mod_outfor,       only : outfor
  use def_kermod,       only : turmu_ker
  use def_kermod,       only : kfl_chemic_vect
  use def_kermod,       only : lookup_fw
  use mod_messages,     only : messages_live
  use mod_chm_mixedEq,  only : CHM_EQ_Z
  use mod_chm_mixedEq,  only : CHM_EQ_ZVAR
  use mod_chm_mixedEq,  only : CHM_EQ_ZZ
  use mod_chm_mixedEq,  only : CHM_EQ_YC
  use mod_chm_mixedEq,  only : CHM_EQ_YCVAR
  use mod_chm_mixedEq,  only : CHM_EQ_YCYC
  use mod_chm_mixedEq,  only : CHM_EQ_PASSIVE
  use mod_chm_mixedEq,  only : CHM_EQ_REACTIVE
  use mod_chm_mixedEq,  only : CHM_EQ_SECTIONAL
  use mod_chm_mixedEq,  only : CHM_GR_PASSIVE
  use mod_chm_mixedEq,  only : CHM_GR_CONTROL
  use mod_chm_mixedEq,  only : CHM_GR_REACTIVE
  use mod_chm_mixedEq,  only : CHM_GR_SECTIONAL
  use mod_chm_mixedEq,  only : CHM_SRC_OFF
  use mod_chm_mixedEq,  only : CHM_SRC_TABLE
  use mod_chm_mixedEq,  only : CHM_SRC_DSM
  use mod_chm_mixedEq,  only : chm_mixedEq_getEqType
  use mod_chm_mixedEq,  only : chm_mixedEq_getSrcType
  implicit none
  integer(ip)   :: ierro,iwarn
  integer(ip)   :: idimt
  integer(ip)   :: cmean_present
  integer(ip)   :: cvar_present
  integer(ip)   :: chist_present
  integer(ip)   :: zmean_present
  integer(ip)   :: zvar_present
  integer(ip)   :: imean_present
  integer(ip)   :: iequa, igrou, ncontrol, nsectional

  external      :: errors

  ierro = 0_ip
  iwarn = 0_ip
  !
  ! TRANSIENT PROBLEM
  !
  if( kfl_timei /= 0 .and. kfl_timei_chm == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'STEADY CHEMIC IN A TRANSIENT CALCULATION')
  end if
  !
  ! Turbulent diffusion
  !
  if(turmu_ker % kfl_exist /= 0_ip .and. abs(diffu_chm(1,1)) < 1e-16_rp ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE DEFINE A TURBULENT SCHMIDT NUMBER TO RUN WITH A TURBULENT VISCOSITY MODEL')
  elseif (turmu_ker % kfl_exist == 0_ip .and. abs(diffu_chm(1,1)) > 1e-16_rp) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE DO NOT DEFINE A TURBULENT SCHMIDT NUMBER IF RUNNING WITHOUT A TURBULENT&
         & VISCOSITY MODEL')
  endif

  !
  ! MIXED EQUATION MODEL
  ! No vectorization:
  !
  if (kfl_chemic_vect == 0 .and. kfl_model_chm == 2) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE TURN ON CHEMIC VECTORIZATION IN .ker.dat: CHMVECTORIZATION: ON')
  endif

  !
  ! Variance main equation is wrong:
  !
  do iequa = 1, nclas_chm
     if (CHM_EQ_ZVAR  == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
        &CHM_EQ_ZZ    == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
        &CHM_EQ_YCVAR == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
        &CHM_EQ_YCYC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
        if ( mixedEq_eqs_chm(iequa) % kfl_ieq_mean <= 0_ip .or. &
            &mixedEq_eqs_chm(iequa) % kfl_ieq_mean >  nclas_chm) then
            !
            ! Out of range
            !
            ierro = ierro + 1
            call outfor(1_ip,momod(modul)%lun_outpu,'FOR VARIANCE EQUATION '//trim(intost(iequa))//' THE INDEX OF MEAN EQUTION: '//&
                trim(intost(mixedEq_eqs_chm(iequa) % kfl_ieq_mean))//' IS OUT OF RANGE: 0<IEQ_MEAN<='//trim(intost(nclas_chm)))
        else
            !
            ! Mean equation is wrong type:
            !
            if ((CHM_EQ_ZVAR  == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
               &CHM_EQ_ZZ    == mixedEq_eqs_chm(iequa) % kfl_eqtype) .and.&
               &CHM_EQ_Z     /= mixedEq_eqs_chm( mixedEq_eqs_chm(iequa) % kfl_ieq_mean) % kfl_eqtype) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'MEAN EQUATION: '//trim(intost(mixedEq_eqs_chm(iequa) % kfl_ieq_mean))//&
                   ' IS NOT A MIXTURE FRACTION EQUATION FOR VARIANCE EQUATION '//trim(intost(iequa)))

            elseif((CHM_EQ_YCVAR == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
                  &CHM_EQ_YCYC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) .and.&
                  &CHM_EQ_YC    /= mixedEq_eqs_chm( mixedEq_eqs_chm(iequa) % kfl_ieq_mean) % kfl_eqtype) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'MEAN EQUATION: '//trim(intost(mixedEq_eqs_chm(iequa) % kfl_ieq_mean))//&
                   ' IS NOT A PROGRESS VARIABLE EQUATION FOR VARIANCE EQUATION '//trim(intost(iequa)))
            endif

        endif
     endif
  enddo

  !
  ! Eqtype is YCVAR OR YCYC but the source term is undefined
  !
  do iequa = 1, nclas_chm
     if (CHM_EQ_YCVAR == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
        &CHM_EQ_YCYC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
        !
        ! Source type has to be tabulated
        !
        if ( CHM_SRC_TABLE /= mixedEq_eqs_chm(iequa) % kfl_source_type ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' IS FOR PROGRESS VARIABLE VARIACE OR&
              & PROGRESS VARIABLE SQUARE BUT THE SOURCE TYPE IS NOT TABULATED SOURCE: '//&
              trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa))))
        endif
     endif
  enddo


  !
  ! Eqtype is YC but the source term is undefined
  ! (ONLY WARNING IN CASE SOMEBODY WANTS TO RUN FOR A WHILE TO SMOOTH THE FIELDS)
  !
  do iequa = 1, nclas_chm
     if (CHM_EQ_YC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
        !
        ! Source type has to be defined
        !
        if ( CHM_SRC_OFF == mixedEq_eqs_chm(iequa) % kfl_source_type ) then
           call messages_live('EQUATION: '//trim(intost(iequa))//' IS FOR PROGRESS VARIABLE BUT THE SOURCE TERM IS OFF' ,'WARNING')
           call messages_live('EQUATION: '//trim(intost(iequa))//' IS FOR PROGRESS VARIABLE BUT THE SOURCE TERM IS OFF' ,'WARNING')
           call messages_live('EQUATION: '//trim(intost(iequa))//' IS FOR PROGRESS VARIABLE BUT THE SOURCE TERM IS OFF' ,'WARNING')
           call messages_live('EQUATION: '//trim(intost(iequa))//' IS FOR PROGRESS VARIABLE BUT THE SOURCE TERM IS OFF, ARE YOU&
              & REALLY SURE ABOUT THIS?' ,'WARNING')
        endif
     endif
  enddo


  !
  ! Source type is tabulated by framework is undefined
  !
  do iequa = 1, nclas_chm
     if (CHM_SRC_TABLE == mixedEq_eqs_chm(iequa) % kfl_source_type .or. &
        &CHM_EQ_YCVAR == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
        &CHM_EQ_YCYC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
        !
        ! Source framework has to be defined
        !
        if ( mixedEq_eqs_chm(iequa) % kfl_source_fw == 0_ip .and.      &
             mixedEq_eqs_chm(iequa) % kfl_premsource_fw == 0_ip .and.  &
             mixedEq_eqs_chm(iequa) % kfl_diffsource_fw == 0_ip) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' THE SOURCE FRAMEWORK IS NOT DEFINED')
        endif
        !
        ! Source column has to be defined
        !
        if ( mixedEq_eqs_chm(iequa) % kfl_source_col == 0_ip     .and. &
             mixedEq_eqs_chm(iequa) % kfl_diffsource_col == 0_ip .and. &
             mixedEq_eqs_chm(iequa) % kfl_premsource_col == 0_ip) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' THE SOURCE COLUMN IS NOT DEFINED')
        endif
     endif
  enddo


  !
  ! Disallow certain uses of the mixed equations model
  !
  ncontrol = 0
  nsectional = 0
  do igrou = 1,ngrou_chm
     if (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_PASSIVE) then
        !
        ! Passive groups
        !
        do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
           !
           ! Only allow passive equations
           !
           if (mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_PASSIVE) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A PASSIVE GROUP BUT IS OF&
                  & TYPE: '//trim(chm_mixedEq_getEqType(mixedEq_eqs_chm(iequa))))
           endif

           !
           ! Only allow if source terms are off
           !
           if (mixedEq_eqs_chm(iequa) % kfl_source_type /= CHM_SRC_OFF) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A PASSIVE GROUP BUT HAS&
                  & SOURCE TERM TYPE: '//trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa))))
           endif
        enddo


     elseif (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_CONTROL) then
        !
        ! Control group
        !
        ncontrol = ncontrol + 1
        do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
           !
           ! Only allow control equations
           !
           if (mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_Z     .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_ZVAR  .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_ZZ    .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_YC    .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_YCVAR .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_YCYC  .and. &
               mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_PASSIVE ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A CONTROL GROUP BUT IS OF&
                  & TYPE: '//trim(chm_mixedEq_getEqType(mixedEq_eqs_chm(iequa))))
           endif
           if (mixedEq_eqs_chm(iequa) % kfl_eqtype == CHM_EQ_PASSIVE) then
               call messages_live('EQAUTION: '//trim(intost(iequa))//' BELONGS TO A CONTROL GROUP BUT IS OF TYPE: '//&
                  trim(chm_mixedEq_getEqType(mixedEq_eqs_chm(iequa)))//', ARE YOU USING THE OLD FLAMELET METHOD?','WARNING')
           endif
        enddo



     elseif (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_REACTIVE) then
        !
        ! Reactive group
        !
        do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
           !
           ! Only allow reactive equations
           !
           if (mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_REACTIVE) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A REACTIVE GROUP BUT IS OF&
                  & TYPE: '//trim(chm_mixedEq_getEqType(mixedEq_eqs_chm(iequa))))
           endif
        enddo




     elseif (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_SECTIONAL) then
        !
        ! Sectional group
        !
        nsectional = nsectional + 1
        do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
           !
           ! Only allow sectional equations
           !
           if (mixedEq_eqs_chm(iequa) % kfl_eqtype /= CHM_EQ_SECTIONAL) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A SECTIONAL GROUP BUT IS OF&
                  & TYPE: '//trim(chm_mixedEq_getEqType(mixedEq_eqs_chm(iequa))))
           endif

           !
           ! Only allow if source terms are of type DSM
           !
           if (mixedEq_eqs_chm(iequa) % kfl_source_type /= CHM_SRC_DSM) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'EQUATION: '//trim(intost(iequa))//' BELONGS TO A SECTIONAL GROUP BUT HAS&
                  & SOURCE TERM TYPE: '//trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa))))
           endif
        enddo
     endif
  enddo

  if (ncontrol > 1) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'THERE IS MORE THAN ONE CONTROL GROUPS: '//trim(intost(ncontrol)))
  endif

  if (nsectional > 1) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'THERE IS MORE THAN ONE SECTIONAL GROUPS: '//trim(intost(nsectional)))
  endif




  !
  ! Table framework
  !
  if (kfl_tab_fw_chm >= 0) then
     cmean_present = 0_ip
     cvar_present  = 0_ip
     chist_present = 0_ip
     zmean_present = 0_ip
     zvar_present  = 0_ip
     imean_present = 0_ip
     if (kfl_multimod_chm == 1_ip) then
        do idimt = 1, lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
           select case (lookup_fw(kfl_tab_fw_chm_diff) % main_table % coords(idimt) % name)
           case ('CMEAN','C    ')
               cmean_present  = idimt
           case ('CVAR ')
               cvar_present   = idimt
           case ('CHIST')
               chist_present  = idimt
           case ('ZMEAN','Z    ')
               zmean_present  = idimt
           case ('ZVAR ')
               zvar_present   = idimt
           case ('IMEAN','I    ')
               imean_present  = idimt
           end select
        enddo
     else
        do idimt = 1, table_fw % main_table % ndim
           select case (table_fw % main_table % coords(idimt) % name)
           case ('CMEAN','C    ')
               cmean_present  = idimt
           case ('CVAR ')
               cvar_present   = idimt
           case ('CHIST')
               chist_present  = idimt
           case ('ZMEAN','Z    ')
               zmean_present  = idimt
           case ('ZVAR ')
               zvar_present   = idimt
           case ('IMEAN','I    ')
               imean_present  = idimt
           end select
        enddo
     end if 

     if(turmu_ker % kfl_exist == 0_ip) then
        !
        ! Turbulent models
        !
        if (cvar_present > 0) then
           if ( table_fw % kfl_scale(cvar_present) /= -1 .and. table_fw % kfl_scale(cvar_present) /= 3 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE TURN OFF CVAR SCALING IN TABLE FRAMEWORK')
           endif
        endif
        if (zvar_present > 0) then
           if ( table_fw % kfl_scale(zvar_present) /= -1  .and. table_fw % kfl_scale(zvar_present) /= 3 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE TURN OFF ZVAR SCALING IN TABLE FRAMEWORK')
           endif
        endif
     end if

     if( kfl_ufpv_chm /= 0 .and. chist_present==0) then
        !
        ! Scalar dissipation rate is in table
        !
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'CHIST SHOULD BE IN TABLE FOR UFPV CACULATION')
     end if

     if( kfl_varZ_chm /= 0 ) then
        !
        ! Z variance model
        !
        if (zvar_present==0) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SHOULD BE IN TABLE')
        else
            if ( abs(kfl_varZ_chm) == 1 .and.  table_fw % kfl_scale(zvar_present) /= (100+zmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SCALING SHOULD BE VARIANCE IN TABLE FRAMEWORK')
            endif
            if ( abs(kfl_varZ_chm) == 2 .and.  table_fw % kfl_scale(zvar_present) /= -1*(100+zmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SCALING SHOULD BE SQUARE IN TABLE FRAMEWORK')
            endif
        end if
     end if

     if( kfl_varYc_chm /= 0 ) then
        !
        ! Yc variance model
        !
        if (cvar_present==0) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SHOULD BE IN TABLE')
        else
            if (kfl_varYc_chm==1 .and. table_fw % kfl_scale(cvar_present) /= (100+cmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SCALING SHOULD BE VARIANCE IN TABLE FRAMEWORK')
            endif
            if (kfl_varYc_chm==2 .and. table_fw % kfl_scale(cvar_present) /= -1*(100+cmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SCALING SHOULD BE SQUARE IN TABLE FRAMEWORK')
            endif
        end if
     end if

     if( kfl_premix_chm /= 0 .and. zmean_present > 0) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'ZMEAN SHOULD NOT BE PRESENT FOR PREMIXED CALCULATION')
     end if

     if( kfl_premix_chm == 0 .and. zmean_present == 0 .and. kfl_model_chm == 1) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'Z OR ZMEAN SHOULD BE PRESENT FOR NON-PREMIXED CALCULATION')
     end if

     !
     ! Chack if enthalpy is needed for lookup
     !
     if( imean_present > 0 ) then
        if ( table_fw % kfl_scale(imean_present) /= -1 .and. table_fw % kfl_scale(imean_present) /= 3&
            .and. table_fw % kfl_needs_enthalpy == 0 ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'ENTHALPY IS A DIMENSION OF THE TABLE BUT IT WILL NOT BE GATHERED FOR LOOKUP')
        endif
     endif

  endif



  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------
  call errors(3_ip,ierro,iwarn,'NULL')

end subroutine chm_outerr
