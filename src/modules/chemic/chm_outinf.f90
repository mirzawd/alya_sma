!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outinf()
  !-----------------------------------------------------------------------
  !****f* partis/chm_outinf
  ! NAME
  !    chm_outinf
  ! DESCRIPTION
  !    This routine writes informtation
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master,      only : INOTSLAVE, kfl_rstar, kfl_htran, momod, modul, intost, retost
  use def_kintyp,      only : ip, rp
  use def_chemic,      only : table_fw, mixedEq_groups_chm, mixedEq_eqs_chm, kfl_annfw_has_sources, kfl_cpCoefHT_end_tab_index_chm,&
                              kfl_cpCoefHT_tab_index_chm, kfl_cpCoefLT_end_tab_index_chm, kfl_cpCoefLT_tab_index_chm,&
                              kfl_DtRho_tab_index_chm, kfl_icmean_chm, kfl_icvar_chm, kfl_izmean_chm, kfl_izvar_chm,&
                              kfl_k_tab_index_chm, kfl_max_src_annfw_chm, kfl_max_srcfw_chm, kfl_min_src_annfw_chm,&
                              kfl_min_srcfw_chm, kfl_model_chm, kfl_mu_tab_index_chm, kfl_T_tab_index_chm, kfl_tab_fw_chm,&
                              kfl_W_tab_index_chm, nclas_chm, ngrou_chm, kfl_fw_src_equa_list, kfl_annfw_src_equa_list,&
                              kfl_fw_has_sources
  use mod_chm_mixedEq, only : chm_mixedEq_getGroupType
  use mod_chm_mixedEq, only : chm_mixedEq_getEqType
  use mod_chm_mixedEq, only : chm_mixedEq_getSrcType
  use mod_chm_mixedEq, only : chm_mixedEq_getIniType
  use mod_chm_mixedEq, only : CHM_SRC_OFF
  use mod_chm_mixedEq, only : CHM_SRC_TABLE
  use mod_chm_mixedEq, only : CHM_SRC_DSM
  use mod_chm_mixedEq, only : CHM_SRC_ANN
  use mod_chm_mixedEq, only : CHM_INI_OFF
  use mod_chm_mixedEq, only : CHM_INI_FIELD
  use mod_chm_mixedEq, only : CHM_INI_CONST
  implicit none
  integer(ip)   :: idimt,iequa,igrou,ifw,ii

  if( INOTSLAVE ) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then


        !
        ! Information about tabulated properties
        !
        if (kfl_tab_fw_chm > 0_ip) then
           !
           ! Print table dimensions
           !
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'---------------------------------------------'
           write(momod(modul) % lun_outpu,*)'TABULATED CHEMISTRY: THERMOCHEMICAL DATABASE '
           write(momod(modul) % lun_outpu,*)'---------------------------------------------'
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'TABLE DIMENSIONS:'

           do idimt = 1,table_fw % main_table % ndim
              write(momod(modul) % lun_outpu,'(5X,A5,1X,I5,1X,A4,1X,I5)') table_fw % main_table % coords(idimt) % name ,&
                table_fw % main_table % coords(idimt) % leng, 'IEQ:',  table_fw % kfl_chm_control(idimt)
           enddo
        endif

        !
        ! Information about all equations
        !
        if (kfl_model_chm == 2) then
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'----------------------'
           write(momod(modul) % lun_outpu,*)'MIXED EQUATIONS MODEL '
           write(momod(modul) % lun_outpu,*)'----------------------'
        endif
        write(momod(modul) % lun_outpu,*)''
        write(momod(modul) % lun_outpu,*)'EQUATION GROUPS:'
        do igrou = 1,ngrou_chm
           !
           ! Go through groups
           !
           write(momod(modul) % lun_outpu,*)'----------------------'
           !
           ! Group name and type
           !
           write(momod(modul) % lun_outpu,'(5X,A7,I5,2X,A5,1X,A)') 'IGROUP=',igrou,                                           &
                                                                   & mixedEq_groups_chm(igrou) % name,                        &
                                                                   & trim(chm_mixedEq_getGroupType(mixedEq_groups_chm(igrou)))
           write(momod(modul) % lun_outpu,'(19X,"RANGE",1X,A,"--",A)') trim(intost(mixedEq_groups_chm(igrou) % i_start)), &
                                                                   &   trim(intost(mixedEq_groups_chm(igrou) % i_end))

           if (mixedEq_groups_chm(igrou) % kfl_therm_phor /= 0) then
              write(momod(modul) % lun_outpu,'(19X,A)') 'THERMOPHORETIC DIFFUSION ON'
           endif

           do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
              !
              ! Go through equations of one group
              !
              if (iequa>0 .and. iequa<= nclas_chm) then
                 !
                 ! Name and type
                 !
                 write(momod(modul) % lun_outpu,'(7X,A6,I5,1X,A5,1X,A)') 'IEQUA=', iequa,                                          &
                                                                          & mixedEq_eqs_chm (iequa) % name,                        &
                                                                          & trim(chm_mixedEq_getEqType(mixedEq_eqs_chm (iequa)))

                 !
                 ! If variance equation, print which is the mean equation
                 !
                 if (mixedEq_eqs_chm (iequa) % kfl_ieq_mean /= 0_ip) then
                    write(momod(modul) % lun_outpu,'(19X,A,1X,A)')                                                                 &
                                                                          & 'INDEX OF MEAN EQUATION:',                             &
                                                                          & trim(intost(mixedEq_eqs_chm(iequa) % kfl_ieq_mean))
                 endif

                 !
                 ! Source terms
                 !
                 if (mixedEq_eqs_chm (iequa) % kfl_source_type == CHM_SRC_TABLE) then
                    write(momod(modul) % lun_outpu,'(19X,A5,5(1X,A))')&
                                                             &                                                                     &
                                                                  & 'SRC:', trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa))),  &
                                                                  & 'FW:',  trim(intost(mixedEq_eqs_chm (iequa) % kfl_source_fw)), &
                                                                  & 'COL:', trim(intost(mixedEq_eqs_chm (iequa) % kfl_source_col))
                 elseif (mixedEq_eqs_chm (iequa) % kfl_source_type == CHM_SRC_ANN) then
                    write(momod(modul) % lun_outpu,'(19X,A5,5(1X,A))')&
                                                                  &                                                                &
                                                               & 'SRC:', trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa))),     &
                                                               & 'ANNFW:',  trim(intost(mixedEq_eqs_chm (iequa) % kfl_source_ann)),&
                                                               & 'COL:', trim(intost(mixedEq_eqs_chm (iequa) % kfl_source_col))
                 elseif (mixedEq_eqs_chm (iequa) % kfl_source_type == CHM_SRC_OFF .or. &
                         mixedEq_eqs_chm (iequa) % kfl_source_type == CHM_SRC_DSM) then
                    write(momod(modul) % lun_outpu,'(19X,A5,1X,A)')                                                                &
                                                                  & 'SRC:', trim(chm_mixedEq_getSrcType(mixedEq_eqs_chm(iequa)))
                 endif

                 !
                 ! Initialization
                 !
                 if (mixedEq_eqs_chm (iequa) % kfl_ini_type == CHM_INI_OFF) then
                    write(momod(modul) % lun_outpu,'(19X,A5,2(1X,A))')                                                             &
                                                                 & 'INI:', trim(chm_mixedEq_getIniType(mixedEq_eqs_chm (iequa))),  &
                                                                 & '(INITIALIZED TO ZERO)'

                 elseif (mixedEq_eqs_chm (iequa) % kfl_ini_type == CHM_INI_FIELD) then
                    write(momod(modul) % lun_outpu,'(19X,A5,2(1X,A))')                                                             &
                                                                 & 'INI:', trim(chm_mixedEq_getIniType(mixedEq_eqs_chm (iequa))),  &
                                                                 & trim(intost(mixedEq_eqs_chm (iequa) % kfl_ini_field))
                 elseif (mixedEq_eqs_chm (iequa) % kfl_ini_type == CHM_INI_CONST) then
                    write(momod(modul) % lun_outpu,'(19X,A5,2(1X,A))')                                                             &
                                                                 & 'INI:', trim(chm_mixedEq_getIniType(mixedEq_eqs_chm (iequa))),  &
                                                                 & trim(retost(mixedEq_eqs_chm (iequa) % ini_value))
                 endif

                 !
                 ! Lewis number
                 !
                 if (mixedEq_eqs_chm (iequa) % Lewis /= 1.0_rp) then
                    write(momod(modul) % lun_outpu,'(19X,A5,1X,A)')                                                                &
                                                                 & 'Le:', trim(retost(mixedEq_eqs_chm (iequa) % Lewis))
                 endif

                 !
                 ! Fix diffusivity
                 !
                 if (mixedEq_eqs_chm (iequa) % kfl_fix_diffusion /= 0_ip) then
                    write(momod(modul) % lun_outpu,'(19X,A5,1X,A)')                                                                &
                                                           & 'Diff:', trim(retost(mixedEq_eqs_chm (iequa) % diffusivity)) // " m2/s"
                 endif


                 !
                 ! POSTPROCESSING
                 !
                 if (mixedEq_eqs_chm (iequa) % kfl_do_post /= 0_ip) then
                    write(momod(modul) % lun_outpu,'(19X,A)') 'DO POSTPROCESSING'
                 endif

              endif
           enddo
        enddo
        write(momod(modul) % lun_outpu,*)'----------------------'
        !
        ! Index of typical control varibles
        !
        if (kfl_icmean_chm > 0 .or. &
            kfl_icvar_chm > 0  .or. &
            kfl_izmean_chm > 0 .or. &
            kfl_izvar_chm > 0  ) then
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'INDEX OF CONTROL VARIABLES AMONG UNKNOWNS:'
           if (kfl_icmean_chm > 0)  write(momod(modul) % lun_outpu,'(5X,A5,1X,I4)') 'CMEAN', kfl_icmean_chm
           if (kfl_icvar_chm > 0)   write(momod(modul) % lun_outpu,'(5X,A5,1X,I4)') 'CVAR', kfl_icvar_chm
           if (kfl_izmean_chm > 0)  write(momod(modul) % lun_outpu,'(5X,A5,1X,I4)') 'ZMEAN', kfl_izmean_chm
           if (kfl_izvar_chm > 0)   write(momod(modul) % lun_outpu,'(5X,A5,1X,I4)') 'ZVAR', kfl_izvar_chm
        endif


        !
        ! Source term frameworks
        !
        if ( any(kfl_fw_has_sources > 0)) then
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'LOOKUP FRAMEWORKS USED FOR SOURCE TERMS:'
           do ifw = kfl_min_srcfw_chm,kfl_max_srcfw_chm
              if (kfl_fw_has_sources(ifw) > 0) then
                  write(momod(modul) % lun_outpu,'(5X,A5,3X,I3)') 'IFRWK', ifw
                  do ii = 1,kfl_fw_has_sources(ifw)
                  write(momod(modul) % lun_outpu,'(7X,A5,1X,I3)') 'IEQUA', kfl_fw_src_equa_list(ifw,ii)
                  enddo
              endif
           enddo
        endif

        !
        ! Source term frameworks
        !
        if ( any(kfl_annfw_has_sources > 0)) then
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'ANN FRAMEWORKS USED FOR SOURCE TERMS:'
           do ifw = kfl_min_src_annfw_chm,kfl_max_src_annfw_chm
              if (kfl_annfw_has_sources(ifw) > 0) then
                  write(momod(modul) % lun_outpu,'(5X,A6,2X,I3)') 'IANNFW', ifw
                  do ii = 1,kfl_annfw_has_sources(ifw)
                  write(momod(modul) % lun_outpu,'(7X,A5,1X,I3)') 'IEQUA', kfl_annfw_src_equa_list(ifw,ii)
                  enddo
              endif
           enddo
        endif


        !
        ! Column index of properties
        !
        if (kfl_W_tab_index_chm > 0        .or. &
            kfl_k_tab_index_chm > 0        .or. &
            kfl_mu_tab_index_chm > 0       .or. &
            kfl_cpCoefLT_tab_index_chm > 0 .or. &
            kfl_cpCoefHT_tab_index_chm > 0 .or. &
            kfl_T_tab_index_chm > 0        .or. &
            kfl_DtRho_tab_index_chm > 0  ) then
           write(momod(modul) % lun_outpu,*)''
           write(momod(modul) % lun_outpu,*)'COLUMN INDECES OF PROPERTIES:'
           if (kfl_W_tab_index_chm         > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4)')    'W     ', kfl_W_tab_index_chm
           if (kfl_k_tab_index_chm         > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4)')    'K     ', kfl_k_tab_index_chm
           if (kfl_mu_tab_index_chm        > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4)')    'MU    ', kfl_mu_tab_index_chm
           if (kfl_cpCoefLT_tab_index_chm  > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4,1X,I4)') 'CPijL ', &
               kfl_cpCoefLT_tab_index_chm, kfl_cpCoefLT_end_tab_index_chm
           if (kfl_cpCoefHT_tab_index_chm  > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4,1X,I4)') 'CPijH ', &
               kfl_cpCoefHT_tab_index_chm, kfl_cpCoefHT_end_tab_index_chm
           if (kfl_T_tab_index_chm         > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4)')    'T     ', kfl_T_tab_index_chm
           if (kfl_DtRho_tab_index_chm     > 0) write(momod(modul) % lun_outpu,'(5X,A6,1X,I4)')    'Dt*rho', kfl_DtRho_tab_index_chm
        endif


        !
        ! Information about enthalpy flux model
        !
        if ( kfl_model_chm  == 3) then
           if ( kfl_htran == 1 ) then
              write(momod(modul) % lun_outpu,*) 'Detailed computation for enthalpy flux'
           else
              write(momod(modul) % lun_outpu,*) 'Simplified computation for enthalpy flux'
           end if
        end if

        flush(momod(modul) % lun_outpu)
     end if

  end if


end subroutine chm_outinf

