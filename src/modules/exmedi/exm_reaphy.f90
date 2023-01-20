!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup exmedi
!> @{
!> @file    exm_reaphy.f90
!> @author  Many
!> @date    19/11/2020
!> @brief   physical problem definition
!> @details physical problem definition
!> @}
!-----------------------------------------------------------------------
!.md<module>exmedi
!.md<input>case.exm.dat
!.md<pos>0
!.md<sec>
subroutine exm_reaphy

  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      def_exmedi
  use mod_ecoute,              only :  ecoute
  use mod_messages,            only : messages_live, livinf  
  use mod_opfpos,              only : postpr_intto8
  use mod_exm_fitzhugh_nagumo, only : exm_fitzhugh_nagumo_ReadDat
  use mod_memory,              only : memory_size
  use mod_eccoupling
  use mod_exm_cellmodel
  use mod_exm_drugs,           only : exm_drugs_read_data
  use mod_exm_activation,      only : exm_stim_allocate
  use mod_exm_ecg,             only : exm_ecg_allocate

  implicit none
  integer(ip) :: imate,iipar,iauxi,nauxi
  integer(ip) :: ivalu,istim,jstim,kstim,iall_materials,icoor,nstimtable
  real(rp)    :: stimfreqtable,offset_stimfreqtable
  character(300)           :: messa
  integer(ip) :: kfl_gemod_exm             ! General model type of ionic 
  integer(ip) :: ipara

  fisoc_exm = 0.0_rp   ! isochrones are taken when value becomes higher than threshold
  nstrb_exm     = 0_ip
  kfl_nodif_exm = 0_ip     ! compute diffusion terms 
  nvint_exm     = 2_ip*mitim
  ! 
  ! Allocate memory for kfl_cellmod
  ! 
  if( .not. kfl_exmsld_ecc )then
      call eccou_allocate_memory(2_ip)
  endif

  if( INOTSLAVE ) then

     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('exm_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('exm_reaphy')
     end do

     !
     ! Initializations (defaults)
     !
     kfl_gcoup_exm = 0_ip                     ! No geometric coupling with SOLIDZ
     kfl_genal_exm = 1_ip                     ! General algorithm type (EXPLICIT=1, decoupled 
     !    implicit, monolithic implicit, ...)

!!!!! LEAVE THIS FOR RETROCOMPATIBILITY, TO BE ELIMINATED
     kfl_gemod_exm = 0_ip                     ! Default: FHN
     !    (TT,LR,BR, ...) or 0=no-subcell or approximate models (FHN, FENTON...)
!!!!!

     kfl_cemod_exm = 1_ip                     ! Cell model (propagation model): monodomain 1 o bidomain 2

     kfl_stree_exm = 0_ip                     ! No streeter fiber field created
     kfl_voini_exm = 0_ip                     ! No initial potential imposed

     kcopeecg_exm = -1_ip 
     kfl_paced_exm = 1_ip 

     kfl_fract_diffusion_exm = -1_ip ! Fractional diffusion is deactivated by default
     kfl_ortho_diffusion_exm = .false. ! Orthotropic diffusion
     fract_diff_coef_exm = 1.0_rp        ! Default value of Fractional diffusion is 1.
     fract_diff_nintp_exm = 0_ip       ! Default value for integration points 0

     modst_exm     = 0_ip                     ! Field: Starting stimuli
     nstim_exm     = 0_ip                     ! Zero starting stimulus, to throw an error. The code is not designed to run without stimulus
     nstis_exm     = 0_ip

     nstimtable    = 1_ip                     ! local variable
     stimfreqtable = 0.0_rp                   ! in seconds

     xmccmmate_exm     = 1.0_rp
     modab_exm = 0_ip
     iall_materials= 0                                    ! Each material reads its own set of properties

     fiaxe_exm = 1.0_rp
     strbo_exm(1)=1_ip
     strbo_exm(2)=2_ip
     strbo_exm(3)=3_ip
     stran_endo_exm=60.0_rp
     stran_epi_exm=200.0_rp ! This vale is to check if there's an input

     !kfl_hfmod_exm = EXM_MYOCYTE_NONE
     kfl_hfmodmate_exm = EXM_MYOCYTE_NORMAL
     kfl_isac_exm = 0_ip
     gdiff_exm = 0.0_rp
     kfl_active_material = .TRUE.

     nrootecg_exm = 0_ip

     moneclmate_exm(1_ip,:) = 1000_ip ! Default beats
     moneclmate_exm(2_ip,:) = 750_ip  ! Default cycle length

     !
     ! Begin to read data
     !
     messa = &
          '        READING PHYSICS...'
     call livinf(0_ip,messa,one)
     messa = &
          '        MODULE EXMEDI USES THESE UNITS:'
     call livinf(0_ip,messa,one)
     messa = &
          '          [D] = cm2 / msec  ;   [C_m] = muF / cm2'
     call livinf(0_ip,messa,one)
     messa = &
          '          [I] = muA         ;   [phi] = mV'
     call livinf(0_ip,messa,one)
     messa = &
          '          [x] = cm          ;   [t]   = s'
     call livinf(0_ip,messa,one)
     messa = &
          '          Stimuli: [Current density] =  muA / mm3      '
     call livinf(0_ip,messa,one)
     messa = &
          '        WARNING: TO OBTAIN t IN SECONDS D IS INTERNALLY MULTIPLIED BY 1000.'
     call livinf(0_ip,messa,one)

     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Physical properties definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !
     do while(words(1)/='ENDPH')
        call ecoute('exm_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('exm_reaphy')
           !
           ! ADOC[1]> PROBLEM_DEFINITION
           !
           do while(words(1)/='ENDPR')
              !
              ! ADOC[2]> ALL_MATERIALS: ON | OFF
              !
              ! ADOC[d]> ALL_MATERIALS: Properties are applied to all the materials in the problem. If the flag is not present, everty 
              if(words(1)=='ALLMA') then   ! When many materials are present, all of them share properties
                 iall_materials = 1
                 if(words(2)=='OFF  ') iall_materials = 0
              else if(words(1)=='APPLI') then               ! Applied currents           

                 call runend("EXM_COMAPP: ALWAYS STARTING POTENTIAL OPTION")

                 call ecoute('exm_reaphy')

                 !
                 ! ADOC[2]> PSEUDOECG: ON
                 !
                 ! ADOC[d]> Compute the pseudo-ecg
                 !
              else if(words(1)=='PSEUD') then
                 if (words(2)=='ON') then
                    nvint_exm = 1
                    kcopeecg_exm = kcopeecg_exm + 1
                    messa = &
                         '        CALCULATING PSEUDO-ECG'
                    call livinf(0_ip,messa,one)
                    do while (words(1)/='ENDPS')
                       call ecoute('exm_reaphy')   
                       if (words(1) == 'NUMBE') then
                           nrootecg_exm = getint('NUMBE',-1_ip,'!numnber of ecg derivations') 
                           call exm_ecg_allocate() !allcate space for the table on MASTER
                       elseif (words(1) == 'FREQU') then
                          nvint_exm = getint('FREQU',-1_ip,'!Frequency of ecg saving') 
                       elseif (words(1) == 'COORD') then
                          if ( nrootecg_exm<1_ip ) then
                             call messages_live("Ecg section is present in exm.dat, however NUMBER of ecg derivations is less than 1","WARNING")
                          end if
                          do icoor = 1,nrootecg_exm
                             call ecoute('exm_reaphy')   

                             if( ( nnpar .ne. 3_ip ) .and. (nnwor > 1_ip ) ) then
                                 call runend( "PSEUDOECG: Coord line "//trim(intost(icoor))//" does not have 3 floats and at most 1 word" )
                             end if 
                             if( ( words(1) .ne. '     ' ) .or. ( words(2) .ne. '     ' ) .or. ( words(3) .ne. '     ' ) ) then
                                 call runend( "PSEUDOECG: Coord line "//trim(intost(icoor))//" has a non-number in the first 3 positions. The firstFIRST 4 words are: '//words(1)//','//words(2)//','//words(3))','//words(4))")
                             end if

                             ecg_points(icoor) % coords(1:ndime) = param(1:ndime)
                             if ( words(4) == '     ' ) then
                                 ecg_points(icoor) % label = trim(intost(icoor))
                             else
                                 ecg_points(icoor) % label = words(4)(1:5)
                             end if
                          end do
                       end if
                    end do
                 end if

                 !
                 ! ADOC[2]> GEO_COUPLING: EULERIAN | LAGRANGIAN 
                 !
                 ! ADOC[d]> GEO_COUPLING: Mesh in wich the excitable media problem is solved.
                 ! ADOC[d]> EULERIAN: The problem is solved in the original fixed mesh.
                 ! ADOC[d]> LAGRANGIAN: The problem is solved in the deformed mesh
                 ! Pseudo-ecg

              else if(words(1)=='GEOCO') then               ! Geometric coupling: 
                 if (words(2)=='LAGRA') then
                    kfl_gcoup_exm = 1_ip    !    is SOLIDZ moving my mesh?
                    messa = &
                         '        RUNNING EXMEDI ON A DEFORMED MESH'
                    call livinf(0_ip,messa,one)
                 else if (words(2)=='EULERIAN') then
                    kfl_gcoup_exm = 0_ip 
                    messa = &
                         '        RUNNING EXMEDI ON A NON-DEFORMED MESH'
                    call livinf(0_ip,messa,one)
                 end if
                 !
                 ! ADOC[2]> STARTING POTENTIAL
                 !
              else if(words(1)=='START') then               ! Starting potential
                 kfl_appty_exm = 1_ip                           ! Decay time considered, good for TT-like models
                 messa = &
                      '        STARTING IMPULSES...'
                 call livinf(0_ip,messa,one)                 
                 !
                 ! ADOC[2]> FLASH: ON | OFF
                 !
                 ! ADOC[d]> FLASH: When ON, no temporal variation of the stimuli is considered. the option LAPSE has no effect.
                 ! ADOC[d]> This option must only be used with Fitzhugh nagumo celular model
                 if(exists('FLASH')) then
                    kfl_appty_exm = 2_ip       ! No decay considered, just a flash, good for FHN
                    if(words(2)=='OFF  ') kfl_appty_exm = 1_ip 
                    messa = &
                         '          TYPE OF INITIAL IMPULSE: FLASH'
                    call livinf(0_ip,messa,one)
                 end if

                 if(exists('TABLE')) then
                    offset_stimfreqtable= 0.0_rp
                    kfl_appva_exm = 1_ip 
                    messa = &
                         '          READING STARTING POTENTIAL FROM TABLE...'
                    call livinf(0_ip,messa,one)
                    call ecoute('exm_reaphy')
                    if (words(1) == 'NSTIM') then       ! Reach from the center
                       nstim_exm = getint( 'NSTIM', 0_ip, "!Number of stimuli")
                       if (words(2) == 'BOUND') nstim_exm = -1_ip 
                       if (words(3) == 'CURDE') kfl_appva_exm = 2_ip 
                       if (words(3) == 'VOLTA') kfl_appva_exm = 1_ip
                       if (exists('CYCLE')) then
                          stimfreqtable = &
                               getrea('CYCLE',0.00_rp,'#frequency at which the initial stimuli table is repeated ') 
!!!!                          stimfreqtable = stimfreqtable / 1000.0_rp  ! read in mSecs and transformed to Secs
                       end if
                       if (exists('REPEA')) then
                          nstimtable = getint('REPEA',1_ip,'#times you want to repeat the initial stimuli table') 
                       end if
                       if (exists('OFFSE')) then
                          offset_stimfreqtable =&
                               getrea('OFFSE',0.0_rp,'#time offset for the intial stimuli table') 
                       end if
                    else
                       call runend('EXM_REAPHY: GIVE FIRST THE TOTAL NUMBER OF STARTING STIMULI')
                    end if

                    nauxi= nstim_exm
                    if (nstim_exm > 0) then
                       nstim_exm = nstim_exm * nstimtable ! the total amount of initial stimuli
                       call exm_stim_allocate()               !allocate space on master for the table
                    end if
                    
                    if (nstim_exm < 0) nauxi= 1                    
                    !
                    ! read the table
                    !
                    do istim=1,nauxi
                       ! ALWAYS IN THIS ORDER: VALUE, CENTER (DE TAMANYO NDIME), REACH, TSTAR, LAPSE
                       call ecoute('exm_reaphy')

                       !Check if there enough numeric columns in the file
                       if (nnpar .ne. ndime+4) then
                          call runend("EXMEDI: Starting potential table row "//trim(intost(istim))//" does not have "//trim(intost(ndime+4))//" numeric columns")
                       end if
                       if (nnwor > 0_ip) then
                          call runend("EXMEDI: Starting potential table row "//trim(intost(istim))//" has unexpected words (non-numbers)")
                       end if

                       apval_exm(istim) = param(1)
                       apcen_exm(1:ndime,istim) = param(2:ndime+1)
                       aprea_exm(istim) = param(ndime+2)
                       aptim(istim) = param(ndime+3) + offset_stimfreqtable
                       aplap_exm(istim) = param(ndime+4)
                       
                    end do                    
                    !
                    ! extend the table, when repeated over many beats
                    !
                    do jstim=1,nstimtable-1
                       do istim= 1,nauxi
                          kstim= jstim * nauxi 
                          apval_exm(kstim + istim)         = apval_exm(istim)
                          apcen_exm(1:ndime,kstim + istim) = apcen_exm(1:ndime,istim)  
                          aprea_exm(kstim + istim)         = aprea_exm(istim)          
                          aptim(kstim + istim)             = &
                               aptim(istim) + stimfreqtable * real(jstim)             
                          aplap_exm(kstim + istim)         = aplap_exm(istim)          
                       end do
                    end do

                 else  if(exists('FIELD')) then
                    call ecoute('exm_reaphy')                       
                    if (exists('CURDE')) then 
                       kfl_appva_exm = 2_ip 
                       modst_exm = -getint('FIELD',0_ip,'!Field Number for the current density starting stimuli') 
                       messa = &
                            '          CURRENT DENSITY INITIAL STIMULI STORED IN FIELD: '//trim(intost(-modst_exm))
                       call livinf(0_ip,messa,one)
                    else if (exists('VOLTA')) then 
                       modst_exm = -getint('FIELD',0_ip,'!Field Number for the voltage starting stimuli') 
                       kfl_appva_exm = 1_ip 
                       messa = &
                            '          VOLTAGE INITIAL STIMULI STORED IN FIELD: '//trim(intost(-modst_exm))
                       call livinf(0_ip,messa,one)
                    else 
                       call runend("EXMEDI: Type of starting potential from FIELD is not specified.")
                    end if
                 else       
                    call ecoute('exm_reaphy')
                    messa = &
                         '          READING STARTING POTENTIAL FROM PARAMETERS...'
                    call livinf(0_ip,messa,one)
                    if (words(1) == 'NSTIM') then       ! Reach from the center
                       nstim_exm = getint('NSTIM',0_ip,'#Number of initial stimuli')
                       call exm_stim_allocate() ! allocate table on master
                       
                       if (words(2) == 'BOUND') then
                          nstim_exm = -1_ip 
                          nstis_exm = getint('BOUND', 2_ip,'#Set number for the initial stimuli')
                       end if
                    else 
                       call runend('EXM_REAPHY: GIVE FIRST THE TOTAL NUMBER OF STARTING STIMULI')
                    end if
                    nauxi= nstim_exm
                    if (nstim_exm < 0) nauxi = 1_ip 
                    do while(words(1)/='ENDST')
                       if (words(1) == 'VALUE' .or. words(1) == 'VOLTA') then                 ! Value
                          kfl_appva_exm = 1_ip 
                          apval_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'CURDE') then            ! Density current            
                          kfl_appva_exm = 2_ip 
                          apval_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'CENTE') then       ! Center of application
                          do iauxi=1,nauxi
                             iipar=(iauxi-1)*ndime + 1
                             apcen_exm(1:ndime,iauxi) = param(iipar:ndime+iipar-1)
                          end do
                       else if (words(1) == 'LAPSE') then       ! Time lapse of application
                          aplap_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'TSTAR') then       ! Start time
                          aptim(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'REACH') then       ! Reach from the center
                          aprea_exm(1:nauxi) = param(1:nauxi)
                       end if
                       call ecoute('exm_reaphy')

                    end do

                 end if

              end if
              call ecoute('exm_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !conve= 0_ip
           ! Default: Ca+ concentration is not coming from a cell model, but from time spent after a wave 
           ! has passed, like FHN or Fenton.
           !kfl_voini_exm     = 0_ip           

           imate= 1_ip
           messa = &
                '        READING MATERIALS...'
           call livinf(0_ip,messa,one)
           do while(words(1)/='ENDPR')
              if( words(1) == 'NODIF' ) then
                 kfl_nodif_exm = 1_ip                  ! do not compute diffusion terms
                 messa = &
                      '        WARNING: NO DIFFUSION TERMS COMPUTED, NO MATTER THE DIFFUSION VALUE READ!!!! '
                 call livinf(0_ip,messa,one)
              else if (words(1)=='APEXT') then   
                 modab_exm = -getint('FIELD',0_ip,'!apex-to-base heterogeneity FIELD')
                 messa = &
                      '        READING APEX-TO-BASE HETEROGENEITY FROM FIELD '//trim(intost(modab_exm))
                 call livinf(0_ip,messa,one)
              else if( words(1) == 'MATER' ) then
                 imate = getint('MATER',1_ip,'#Current material')
                 messa = &
                      '        MATERIAL: '//trim(intost(imate))
                 call livinf(0_ip,messa,one)

                 if( imate > nmate_exm ) then
                    call runend('EXM_REAPHY: THIS MODULE HAS MORE MATERIALS THAN THOSE DEFINED IN DOM.DAT')
                 end if

              else if( words(1)=='INACT' ) then
                 kfl_active_material(imate) = .FALSE.
                 call messages_live("MATERIAL "//trim(intost(imate))//" MARKED AS INACTIVE","WARNING")

              else if (words(1)=='CONTI') then         ! continuum model properties 
                 messa = &
                      '        CONTINUUM MODEL '
                 call livinf(0_ip,messa,one)

                 ! ORTHOTROPIC or TRANSVERSAL ISOTROPIC diffusion
                 if( exists('ORTHO') ) kfl_ortho_diffusion_exm(imate) = .true.
                 if( exists('TRANS') ) kfl_ortho_diffusion_exm(imate) = .false.

                 do while(words(1)/='ENDCO')

                    !Fractional  diffusion
                    ! BCAM collaboration (ncusimano / lgerardo [at] bcamath.org )
                    if(words(1)=='FRACT') then
                       messa = '           USING FRACTIONAL DIFFUSION'  
                       call livinf(0_ip,messa,one)

                       do while(words(1)/='ENDFR')
                          call ecoute('exm_reaphy')
                          if(words(1)=='COEFF')then
                             fract_diff_coef_exm(imate) = getrea('COEFF',1.0_rp,'#Fractional diffusion coefficient')  
                          elseif(words(1)=='ITERA')then
                             kfl_fract_diffusion_exm(imate)=getint('ITERA',1_ip,'#Fractional diffusion iterations')
                          elseif(words(1)=='NUMBE')then
                             fract_diff_nintp_exm(imate)=getint('NUMBE',1_ip,'#Number of integration points')
                          endif
                       end do
                       ! Checking parameters for fractional diffusion
                       if (fract_diff_coef_exm(imate).gt.1.0_rp .or. fract_diff_coef_exm(imate).le.0.0_rp) then
                          call runend('FRACTIONAL COEFFICIENT SHOULD BE 0>=S>=1')
                       endif
                       if(kfl_fract_diffusion_exm(imate).eq.-1_ip)then
                          call runend('FRACTIONAL DIFFUSION ITERATIONS CANT BE ZERO')
                       endif
                       if (fract_diff_nintp_exm(imate).le.0_ip) then
                          call runend('NUMBER OF INTEGRATION POINTS SHOULD BE >=0')
                       endif
                    endif

                    if(words(1)=='INCON' .or. words(1)=='CONDU') then  ! we retain INCON for retrocompatibility
                       call runend('EXM_REAPHY: USE INDIF FOR THE INTRACELLULAR DIFFUSIVITY')
                    else if(words(1)=='INDIF') then  
                       gdiff_exm(1,1,imate) = param(1)
                       gdiff_exm(1,2,imate) = param(2)
                       gdiff_exm(1,3,imate) = param(3)                       
                    else if (words(1)=='INITI') then
                       kfl_voini_exm(imate) = 1_ip 
                       vminimate_exm(:,imate) = param(1)                         
                    else if (words(1)=='ORTHO') then
                       kfl_ortho_diffusion_exm(imate) = .true.
                    end if
                    call ecoute('exm_reaphy')
                 end do
                 messa = &
                      '           DIFFU 1: '//trim(retost(gdiff_exm(1,1,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           DIFFU 2: '//trim(retost(gdiff_exm(1,2,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           DIFFU 3: '//trim(retost(gdiff_exm(1,3,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           INITIAL VOLTAGE: '//trim(retost(vminimate_exm(1_ip,imate)))//' mV'
                 call livinf(0_ip,messa,one)

                 ! internally convert diffusion units from ms to 
                 gdiff_exm(1,1,imate) = 1000.0_rp * gdiff_exm(1,1,imate) 
                 gdiff_exm(1,2,imate) = 1000.0_rp * gdiff_exm(1,2,imate) 
                 gdiff_exm(1,3,imate) = 1000.0_rp * gdiff_exm(1,3,imate) 



              else if(words(1)=='CELLM') then  ! cell model properties
                 messa = &
                      '        CELL MODEL '
                 call livinf(0_ip,messa,one)

                 if (words(2)=='FITZH') then
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_FITZHUGH )

                    call exm_fitzhugh_nagumo_ReadDat(imate)

                 else if (words(2) =='FENTO') then                      
                    messa = words(2)
                    call eccou_set_cellmod(imate,  EXMSLD_CELL_FENTON)
                    !kfl_fento_exm = 1_ip 

                    call runend('EXM_REAPHYNEW: FENTON-KARMA MODELS TO BE REPROGRAMMED.')

                    !if (words(3) == 'BEELE') kfl_fento_exm = 1_ip 
                    !if (words(3) == 'MODBE') kfl_fento_exm = 2_ip 
                    !if (words(3) == 'LUORU') kfl_fento_exm = 3_ip 
                    !if (words(3) == 'GIROU') kfl_fento_exm = 4_ip 
                    !
                    ! Define some model parameters and dimensions for the fenton model
                    !

                    !if (kfl_fento_exm == 1_ip ) then

                    !   xmopa_exm( 3,imate)  = 0.13_rp    ! phi c
                    !   xmopa_exm( 4,imate)  = 3.33_rp    ! tau v +
                    !   xmopa_exm( 5,imate)  = 0.04_rp    ! phi v
                    !
                    !   xmopa_exm( 7,imate)  = 1250.0_rp  ! tau v1 -
                    !   xmopa_exm( 8,imate)  = 19.6_rp    ! tau v2 -
                    !   xmopa_exm( 9,imate)  = 33.0_rp    ! tauro             
                    !   xmopa_exm( 10,imate) = 30.0_rp    ! tausi
                    !   xmopa_exm( 11,imate) = 12.5_rp    ! tauso
                    !   xmopa_exm( 12,imate) = 0.85_rp    ! phics
                    !   xmopa_exm( 13,imate) = 4.0_rp     ! gephi
                    !   xmopa_exm( 16,imate) = 10.0_rp    ! consk
                    !   xmopa_exm( 17,imate) = 870.0_rp   ! tau w +
                    !   xmopa_exm( 18,imate) = 41.0_rp    ! tau w -
                    !
                    !else if (kfl_fento_exm == 2_ip ) then
                    !
                    !   xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                    !   xmopa_exm( 4,imate)  =  3.33_rp    ! tau v +
                    !   xmopa_exm( 5,imate)  =  0.055_rp   ! phi v
                    !
                    !   xmopa_exm( 7,imate)  =  1000.0_rp  ! tau v1 -
                    !   xmopa_exm( 8,imate)  =  19.2_rp    ! tau v2 -
                    !   xmopa_exm( 9,imate)  =  50.0_rp    ! tauro          
                    !   xmopa_exm( 10,imate) =  45.0_rp    ! tausi
                    !   xmopa_exm( 11,imate) =  8.3_rp     ! tauso
                    !   xmopa_exm( 12,imate) =  0.85_rp    ! phics
                    !   xmopa_exm( 13,imate) =  4.0_rp     ! gephi
                    !   xmopa_exm( 16,imate) =  10.0_rp    ! consk
                    !   xmopa_exm( 17,imate) =  667.0_rp   ! tau w +
                    !   xmopa_exm( 18,imate) =  11.0_rp    ! tau w -
                    !
                    !else if (kfl_fento_exm == 3_ip ) then
                    !
                    !   xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                    !   xmopa_exm( 4,imate)  =  10.0_rp    ! tau v +
                    !   xmopa_exm( 5,imate)  =  0.0_rp     ! phi v 
                    !
                    !   xmopa_exm( 7,imate)  =  18.2_rp    ! tau v1 - 
                    !   xmopa_exm( 8,imate)  =  18.2_rp    ! tau v2 -
                    !   xmopa_exm( 9,imate)  =  130.0_rp   ! tauro                
                    !   xmopa_exm( 10,imate) =  127.0_rp   ! tausi
                    !   xmopa_exm( 11,imate) =  12.5_rp    ! tauso 
                    !   xmopa_exm( 12,imate) =  0.85_rp    ! phics
                    !   xmopa_exm( 13,imate) =  5.8_rp     ! gephi
                    !   xmopa_exm( 16,imate) =  10.0_rp    ! consk
                    !   xmopa_exm( 17,imate) =  1020.0_rp  ! tau w +
                    !   xmopa_exm( 18,imate) =  80.0_rp    ! tau w -
                    !
                    !
                    !else if (kfl_fento_exm == 4_ip ) then
                    !
                    !
                    !   xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                    !   xmopa_exm( 4,imate)  =  10.0_rp    ! tau v +
                    !   xmopa_exm( 5,imate)  =  0.025_rp   ! phi v
                    !
                    !   xmopa_exm( 7,imate)  =  333.0_rp   ! tau v1 -
                    !   xmopa_exm( 8,imate)  =  40.0_rp    ! tau v2 -
                    !   xmopa_exm( 9,imate)  =  25.0_rp    ! tauro                
                    !   xmopa_exm( 10,imate) =  22.0_rp    ! tausi
                    !   xmopa_exm( 11,imate) =  12.5_rp    ! tauso 
                    !   xmopa_exm( 12,imate) =  0.85_rp    ! phics
                    !   xmopa_exm( 13,imate) =  8.7_rp     ! gephi
                    !   xmopa_exm( 16,imate) =  10.0_rp    ! consk
                    !   xmopa_exm( 17,imate) =  1000.0_rp  ! tau w +
                    !   xmopa_exm( 18,imate) =  65.0_rp    ! tau w -
                    !end if

                 else if (words(2) =='TENTU') then 
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_TENTUSCHER )
                    do while(words(1)/='ENDCE')                       

                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate) = getrea( "CMCON", 1.0_rp, "!Model constant Cm")
                       end if

                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do

!                 else if (words(2) =='TTHET') then 
!                    messa = words(2)
!                    call eccou_set_cellmod( imate, EXMSLD_CELL_TT2006 )
!                    do while(words(1)/='ENDCE')                       
!                       if(words(1)=='CMCON') then               ! Model constant Cm
!                          xmccmmate_exm(imate)=param(1)                          
!                       end if
!                       ! if needed, user-defined model parameters must go here
!                       call ecoute('exm_reaphy')
!                    end do
                 else if (words(2) =='OHARA') then 
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_OHARA )
                    if (words(3) == 'INAPA') then
                       messa = words(3)
                       call eccou_set_cellmod( imate, EXMSLD_CELL_OHARA_INAPA )
                    end if
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=getrea( "CMCON", 1.0_rp, "!Model constant Cm")
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) == 'TOROR') then
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_TORORD )
                    ! if (words(3) == 'BORDE') then
                    !    messa = words(2)//' '//words(3)
                    !    kfl_borde_exm(imate) = 1_ip
                    !    border_gkrsc_exm(imate) = param(3)
                    ! end if
                 else if (words(2) =='NOMOD') then                      ! no ionic current model
                    messa = words(2)
                    call eccou_set_cellmod( imate,  EXMSLD_CELL_NOMODEL )
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=getrea( "CMCON", 1.0_rp, "!Model constant Cm")
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='SCATR') then 
                    messa = words(2)
                    call eccou_set_cellmod( imate,  EXMSLD_CELL_SCATRIA )
                    do while(words(1)/='ENDCE')   
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=getrea( "CMCON", 1.0_rp, "!Model constant Cm")
                       end if
                       if(words(1)=='PACED') then               ! The model can be paced or it can beat on it's own
                          kfl_paced_exm = 1_ip                     
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='SCVEN') then 
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_SCVENTRI )
                    do while(words(1)/='ENDCE')    
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=getrea( "CMCON", 1.0_rp, "!Model constant Cm")
                       end if
                       if(words(1)=='PACED') then               ! The model can be paced or it can beat on it's own
                          kfl_paced_exm = 1_ip   
                       else 
                          kfl_paced_exm = 0_ip                         
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='COURT') then 
                    messa = words(2)
                    call eccou_set_cellmod( imate, EXMSLD_CELL_COURTE ) 
                    if (words(3) == 'ISACC') then
                       messa = words(3)
                       kfl_isac_exm(imate) = 1_ip
                    end if
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 end if

                 messa = &
                      '           TYPE: '//adjustl(trim(messa))
                 call livinf(0_ip,messa,one)

              else if(words(1)=='INICE') then                       ! CURRENTS PARAMETERS AND CONSTANTS
                 if (stimfreqtable > 0.0_rp) moneclmate_exm(2,imate) = int(stimfreqtable) * 1000_ip ! back to mSecs... this is the default value, if stimfreqtable is defined
                 do while(words(1)/='ENDIN')
                    call ecoute('exm_reaphy')                   
                    if(words(1)=='BEATS') then
                       moneclmate_exm( 1,imate)=getint("BEATS",1000_ip,"!Number of beats to convergence")             ! Number of Beats to simulate until steady state     
                    else if(words(1)=='TIMES') then
                       timestep_cellm(imate) = getrea('TIMES',-1.0_rp,'#TIMESTEP FOR CELLMODEL(S)')  !in seconds! 
                    else if(words(1)=='CYCLE') then
                       moneclmate_exm( 2,imate)=getint("CYCLE",750_ip,"!Heart beat cycle length")              ! Heart beat cycle length

                    else if(words(1)=='CELLT') then                        !CELLTYPES = ENDO MID EPI, which celltypes will be used. By default all. This will change for which celltypes the ODE will be solved in exm_oneohr.f90
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_ENDO) = 0_ip !reset
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_MID) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_EPI) = 0_ip

                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_RA) = 0_ip !reset
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_CTBBRA) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_BBLA) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_TVR) = 0_ip 
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_MVR) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_RAA) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_LAA) = 0_ip 
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_LA) = 0_ip
                       !kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_PV) = 0_ip
                       if ( exists('ENDO ') ) then
                          kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_ENDO) = 1_ip          
                       end if
                       if ( exists('MID  ') ) then
                          kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_MID) = 1_ip          
                       end if
                       if ( exists('EPI  ') ) then
                          kfl_user_specified_celltypes_exm(imate, EXM_CELLTYPE_EPI) = 1_ip          
                       end if

                       !IGNORE_STEADYSTATE = ENDO MID EPI. Write the names of the celltypes (ENDO MID EPI). For the specified celltypes Alya will NOT terminate if steady state  is not reached. Only warning will be displayed
                    else if(words(1)=='IGNOR') then                        
                       if ( exists('ENDO ') ) then
                          kfl_ignore_steadystate_celltypes_exm(imate, EXM_CELLTYPE_ENDO) = 1_ip          
                       end if
                       if ( exists('MID  ') ) then
                          kfl_ignore_steadystate_celltypes_exm(imate, EXM_CELLTYPE_MID) = 1_ip          
                       end if
                       if ( exists('EPI  ') ) then
                          kfl_ignore_steadystate_celltypes_exm(imate, EXM_CELLTYPE_EPI) = 1_ip          
                       end if
                       if ( exists('ALL  ') ) then
                          kfl_ignore_steadystate_celltypes_exm(imate, :) = 1_ip          
                       end if


                       !STEADY_STATE_VARIABLE = (CALCIUM|VOLTAGE) [TOLERANCE = real]
                    else if(words(1)=='STEAD') then  
                       if (exists('CALCI')) then
                          kfl_steadystate_variable(imate) = EXM_CELL_STEADY_CALCIUM
                       end if

                       if (exists('VOLTA')) then
                          kfl_steadystate_variable(imate) = EXM_CELL_STEADY_VOLTAGE
                       end if

                       if ( exists('TOLER') ) then
                          steady_tol_cellm(imate) = getrea('TOLER',-1.0_rp,'#TOLERANCE FOR STEADY STATE')
                       end if

                       !call runend("EXM_REAPHY: KEYWORD STEAD IN CELL DEFINITION FOR MATERIAL "//trim(intost(imate))//" IS NOT FOLLOWED BY ANYTHING RECOGNISEABLE.")
                       
                       !kfl_hfmodmate_exm 0 or 1 -- flag if hardcoded initial conditions should be used:
                       ! 0 - uses hardcoded initial conditions for cell, 1 - modified cell type, meaning initial conditions need to be calculated by alya
                       !kfl_hfmod_exm - type of cell: male/female/pig. 

                    else if(words(1)=='MYOCY') then               ! Normal Cell or Heart Failure simulation
                       if(words(2)=='NORMA') then
                          call runend("EXMEDI: MATERIAL "//trim(intost(imate))//", MYOCYTE NORMAL is disabled in the code for now. Use MYOCYTE MODIFIED")
                          !kfl_hfmodmate_exm(imate) = EXM_MYOCYTE_NORMAL
                          !ttparmate_exm(:, :,                           imate) = 1.0_rp
                          !ttparmate_exm(:, ohara_conductance_ikatp_row, imate) = 0.0_rp   !IKatp current should only be non-zero if ischemia
                          !!call ecoute('exm_reaphy')
                          !if(words(3)=='MALE') then
                          !   kfl_hfmod_exm(imate) = EXM_MYOCYTE_MALE
                          !else if(words(3)=='FEMAL') then
                          !   kfl_hfmod_exm(imate) = EXM_MYOCYTE_FEMALE 
                          !end if  
                       else if(words(2)=='MODIF') then
                          kfl_hfmodmate_exm(imate) = EXM_MYOCYTE_MODIFIED
                          !kfl_hfmod_exm(imate) = EXM_MYOCYTE_TABLE  !So that the parameter table is read                     
                          !if(words(3)=='MALE') then
                          !   kfl_hfmod_exm(imate) = EXM_MYOCYTE_MALE
                          ! else if(words(3)=='FEMAL') then
                          !    kfl_hfmod_exm(imate) = EXM_MYOCYTE_FEMALE 
                          ! end if                             
                        else
                           call runend( "EXMEDI: MATERIAL "//trim(intost(imate))//" HAS UNFAMILIAR MYOCYTE TYPE" )
                        end if
         
                     !The following case should execute only if NORMAL is specified in celltype, but not PIG,MALE or FEMALE
                     !also CONDU should appear after myocyte definition
                     else if(words(1)=='CONDU') then
                        !if ( kfl_hfmod_exm(imate) .eq. EXM_MYOCYTE_NONE ) then
                        !    call messages_live("NO MYOCYTE TYPE SPECIFIED, CONDUCTANCE TABLE IS IGNORED","WARNING")
                        !else
 
                            call messages_live( "           EXM_REAPHY: Cell model reading conductances for material "//trim(intost(imate)) )
                            call ecoute('exm_reaphy')
                            ivalu = 0                     
                            if(words(1)=='INAME') then  
                               call ecoute('exm_reaphy')
                               do while(words(1)/='ENDCO') 
                                  ivalu=ivalu+1
 
                                  if ( nnpar .ne. 3_ip ) then 
                                        call runend('EXM_REAPHY: MATERIAL '//trim(intost(imate))//', ROW '//trim(intost(ivalu))//&
                                        ': NUMBER OF VALUES READ IS '//trim(intost(nnpar))//' INSTEAD OF EXPECTED 3')
                                  end if  
                                  if ( nnwor > 0_ip ) then 
                                        call runend('EXM_REAPHY: MATERIAL '//trim(intost(imate))//', ROW '//trim(intost(ivalu))//&
                                        ': READ '//trim(intost(nnwor))//' WORDS. NO WORDS ARE EXPECTED. FIRST 3 WORDS ARE: '//words(1)//','//words(2)//','//words(3))
                                  end if  

                                  if ( any(param(1:3)<0.0_rp) ) then 
                                    call runend('EXM_REAPHY: MATERIAL '//trim(intost(imate))//', ROW '//trim(intost(ivalu))//' CONTAINS NEGATIVE CONDUCTANCE VALUE')
                                  end if  


                                  ttparmate_exm(1,ivalu,imate) = param(1)
                                  ttparmate_exm(2,ivalu,imate) = param(2)
                                  ttparmate_exm(3,ivalu,imate) = param(3)
 
                                  if ( ivalu .ne. ohara_conductance_ikatp_row ) then !row 14, Ik_atp has [0 0 0] by default for normal celltype
                                        do ipara=1_ip,3_ip
                                              if( (abs(param(ipara)) < 1.0e-15_rp) ) then
                                                 !PLEASE, we need to allow blocking some of the currents as not all the experimental models have all the currents
                                                 call messages_live('EXM_REAPHY: CURRENTS ARE BLOCKED IN CELL TYPE '//trim(intost(ipara))//&
                                                       ', MATERIAL '//trim(intost(imate))//', ROW '//trim(intost(ivalu)),'WARNING')
                                              endif
                                        enddo
                                     end if
 
                                  call ecoute('exm_reaphy')
                               end do
 
                               if (ivalu .ne. nvars_ttparmate_exm) then !basic check that we read all the variables, not more not less
                                  call messages_live( "EXM_REAPHY: Material "//trim(intost(imate))//". Read "//trim(intost(ivalu))//" conductances instead of "//trim(intost(nvars_ttparmate_exm)),'WARNING') 
                               end if 
 
                            else
                               call runend("INAME NOT FOUND, CONDUCTANCES NOT READ")
                            end if     
                        ! end if                
                     end if
 
                     call exm_drugs_read_data(imate)
                 end do

              else if(words(1)=='ISOCH') then
                  if(words(2)=='HIGHE') then
                     fisoc_exm(1)=  1.0_rp         ! isochrones are taken when value becomes higher than threshold
                     fisoc_exm(2)= param(2)        ! isochrones threshold
                  else if(words(2)=='LOWER') then  ! TODO: remove this or see how to do it properly, otherwise it will trigger saving instantly at t=0
                     call runend("ISOCH LOWER in exm.dat is not supported")
                     !fisoc_exm(1)= -1.0_rp         ! isochrones are taken when value becomes lower than threshold
                     !fisoc_exm(2)= param(2)        ! isochrones threshold
                  else
                     call runend("ISOCH inrecognised parameters found in exm.dat")
                  end if

              else if(words(1)=='IONIZ') then               ! Ionization current function (FHN)
                 continue
              end if

              call ecoute('exm_reaphy')

           end do ! PROPE

        else if (words(1)=='FIBER') then              
            call messages_live('EXMEDI: Fibers are now defined in ker.dat','WARNING')

        else if (words(1)=='ORTHO') then              
            call runend("EXMEDI: Orthotropic fiber directions are now defined in ker.dat")

        else if (words(1)=='CREAT') then
           call ecoute('exm_reaphy')
           do while(words(1).ne.'ENDCR')
              if(words(1)=='FUNCT') then
                 if (words(2) == 'CUBIC') then
                    kfl_stree_exm=3
                 elseif (words(2) == 'LINEA') then
                    kfl_stree_exm=1
                 else
                    kfl_stree_exm=1 ! Default option is linear
                 endif
              elseif (words(1) == 'VAXIS') then
                 fiaxe_exm(1:3)= param(1:3)

              elseif (words(1) == 'EPICA') then
                 strbo_exm(1) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'LEFTE') then
                 strbo_exm(2) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'RIGHT') then
                 strbo_exm(3) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'ANGLE') then
                 stran_endo_exm = param(1)
                 stran_epi_exm = param(2)
                 if(stran_epi_exm .eq. 0.0_rp) stran_epi_exm = - stran_endo_exm
              endif
              call ecoute('exm_reaphy')
           enddo
        else if (words(1)=='APEXT') then   
           !TODO: This is left for backwards compatilibity, remove when convenient
           modab_exm = -getint('FIELD',0_ip,'#apex-to-base heterogeneity')
           call messages_live("APEX-BASE HETEROGENEITY MUST BE INSIDE THE PROPERTIES SECTION. THE CURRENT SYNTAX WILL BE DISABLED IN THE FUTURE","WARNING")
        end if
     end do ! ENDPH

     !
     ! Unknowns element-wise dimensions
     !
     ndofn_exm = 1 
     if (kfl_cemod_exm == 2) ndofn_exm = 2_ip 
     ndof2_exm = ndofn_exm*ndofn_exm

     !  FOR NO CELL IONIC CURRENTS MODELS (FHN, ...)
     nevat_exm = mnode    ! <--- explicit (default value)
     if (kfl_genal_exm == 2) nevat_exm = mnode                   ! <--- decoupled implicit
     if (kfl_genal_exm == 3) nevat_exm = ndofn_exm*mnode         ! <--- monolithic implicit

     !
     ! Initialization
     !
     dtinv_exm=0.0_rp

     tcardiac_cycle= 0.75_rp ! default value
     if (stimfreqtable > 0.0_rp) tcardiac_cycle= stimfreqtable
     
     ! WARNING WARNING: Esta opcion ahora solo tiene sentido para gdiff.
     if (iall_materials == 1) then

        do imate= 1,nmate
           gdiff_exm(1,1,imate) = gdiff_exm(1,1,1) 
           gdiff_exm(1,2,imate) = gdiff_exm(1,2,1) 
           gdiff_exm(1,3,imate) = gdiff_exm(1,3,1) 
        end do

     end if


     iauxi= -1
     do imate= 1,nmate
        if (kfl_cellmod(imate) .ge. 0_ip ) iauxi= 1
     end do
     if (iauxi == -1) then

        call runend('EXM_REAPHY: CELL MODEL WAS NOT DEFINED FOR ANY MATERIAL')

     end if

     if( (kfl_appva_exm .ne. 1_ip) .and. (kfl_appty_exm==2_ip) ) then
         call runend("EXMEDI: FLASH CAN BE ONLY USED WITH VOLTAGE")
     end if


     if (nstim_exm==0_ip .and. modst_exm >= 0_ip) then
         call runend("EXMEDI: AT LEAST ONE STIMULUS HAS TO BE SPECIFIED")
     end if
     !if (kfl_gemod_exm == 1 .or. conve == 1) then
     ! TEN TUSSCHER model  
     ! kfl_cellmod(imate) = 1 --> TenTusscher model 
     !       (Paper 'A model for human ventricular tissue', 2003)
     ! kfl_cellmod(imate) = 2 --> LuoRudyII  model (Paper ..................... , 1994)
     ! kfl_cellmod(imate) = 3 --> BeelerReuter model (Paper ................., 1977)
     !gdiff_exm(1,:,imate) = gdiff_exm(1,:,imate) * 8991997122.56092078050535_rp
     !(1.0_rp / (4.0_rp * 3.1415_rp * 0.00000000000885_rp)) !conversion from SI to CGS

  end if      !! kfl_paral


end subroutine exm_reaphy


!.md# Physical Properties Definition
!.md<code>
!.md<0>PHYSICAL_PROBLEM
!.md<1>PROBLEM_DEFINITION
!.md<2>ALL_MATERIALS: ON | OFF
!.md<field>ALL_MATERIALS
!.md<com>Properties are applied to all the materials in the problem. <b>TODO: EXPLAIN</b>
!.md<2>PSEUDO_ECG= ON | OFF
!.md<field>PSEUDO_ECG
!.md<com>Integrates the potentials at given points. The values are stored in the exm.vin file. 
!.md<com>To generate ECG the data needs to be postprocessed (i.e. subtract the potentials)
!.md<com>There is an optional 4th column with the string label for the point. By default the label is row number. Alya reads only 5 letters! The labels are used as column names in .vin
!.md<3>NUMBER_ROOT=  int !number of leads
!.md<3>FREQUENCY=    int !frequency of saving the ECG
!.md<3>COORDS            !coordinates, one triple per row, number of rows equal to NUMBER_ROOT
!.md<4>real real real [label]
!.md<4>...
!.md<2>END_PSEUDO
!.md<2>
!.md<2>GEO_COUPLING: EULERIAN | LAGRANGIAN 
!.md<field>GEO_COUPLING
!.md<com>Mesh in wich the excitable media problem is solved.
!.md<com>EULERIAN: The problem is solved in the original fixed mesh.
!.md<com>LAGRANGIAN: The problem is solved in the deformed mesh

!.md<2>STARTING_POTENTIAL [FLASH: ON | OFF] 
!.md<field>STARTING_POTENTIAL
!.md<com>Define the starting potential. There are 3 ways to define the structure, above in the code you see all three, you need to chose one that best suits you.<br/>
!.md<com>FLASH: When ON, no temporal variation of the stimuli is considered (For Fitzhugh Nagumo cell model). the option LAPSE has no effect. <b>TODO: TEST</b><br/>
!.md<com>If "STARTING_POTENTIAL, FIELD" is used, the field should be defined (in dom.dat) as follows:<br/>
!.md<com>FIELD= int, DIMENSION= 4, NODES<br/>
!.md<com>FIELD= int <br/>
!.md<com>   node_number  magnitude time lapse<br/>
!.md<com>   node_number  magnitude time lapse<br/>
!.md<com>   ...<br/>
!.md<com>END_FIELD

!.md<3>NSTIM=    N, BOUND = X                    ! Numnber of stimuli to apply, option BOUND is unclear. <b>TODO: DESCRIBE, TEST</b>
!.md<3>CURDENSITY | VOLTA=   -80 -80 ...         ! Chose either voltage or currdensity, followed by the values, one per stimulus
!.md<3>LAPSE=    0.0001 0.0001 ...               ! Duration of the stimulus, one per stimulus
!.md<3>CENTER=   x y z x y z ...                 ! Coordinates of the stimulus, one per stimulus
!.md<3>REACH=    0.3 0.3 ...                     ! Readius of the stimulus (assigned to all points witih the radius)
!.md<3>T_START=  0.0 0.0 ...                     ! Stimulus start time in s.
!.md<2>END_STARTING_POTENTIAL
!.md<2>
!.md<2>STARTING_POTENTIAL, FIELD    ! read starting potential from field defined on nodes. You cannot use FLASH or other options
!.md<3>CURDENSITY|VOLTAGE:  FIELD= int
!.md<2>END_STARTING_POTENTIAL

!.md<2>
!.md<2>STARTING_POTENTIAL, [FLASH: ON | OFF] TABLE ! read starting potential from a table
!.md<3>NSTIM= int, NOBOU, CURDE ! The options to put here are: (BOUND|NOBOU), (CURDE|VOLTA),  <b>TODO: TEST BOUND, CURDE, VOLTA</b>
!.md<3>     CURDE   COORDX  COORDY  COORDZ RADIUS  TSTART DURATION
!.md<3>     CURDE   COORDX  COORDY  COORDZ RADIUS  TSTART DURATION
!.md<3>     ....
!.md<2>END_STARTING_POTENTIAL
!.md<1>END_PROBLEM_DEFINITION
!.md<1>  
!.md<1>PROPERTIES
!.md<3>NODIF             ! Optional. Do not compute diffusion terms. <b>TODO: DESCRIBE, TEST</b>
!.md<3>MATERIAL int      !if only 1 material is used, words MATERIAL...END_MATERIAL are not necessary
!.md<4>INACTIVE          ! Optional. Inactive material
!.md<4>
!.md<field>MATERIAL INACTIVE
!.md<com>For the inactive material you can still specify all the parameters like voltage, cell model, etc, that will define initial values at the nodes.
!.md<com>However the elements with this material will not participate in matrix assembly and ECG calculation
!.md<4>CONTINUUM_MODEL
!.md<5>FRACT !Fractional  diffusion BCAM collaboration (ncusimano / lgerardo [at] bcamath.org ) <b>TODO: DESCRIBE, TEST</b>
!.md<6>COEFF real #Fractional diffusion coefficient
!.md<6>ITERA int #Fractional diffusion iterations
!.md<6>NUMBE int #Number of integration points
!.md<5>ENDFR
!.md<5>IN_DIFFUSION= real real real                 ! Diffusion coefficients (cm/s if CGS is used for the mesh): along fibers and in 2 transversal directions
!.md<5>INITIAL_VOLTAGE = real                       ! Intial voltage to put at timestep 0
!.md<4>ENDCONTINNUM_MODEL
!.md<4>
!.md<4>CELL_MODEL:   NOMODEL  ! no ionic current model
!.md<field>CELL_MODEL
!.md<com>There are several types of cell models, each has its unique parameters. In the code above you see all of them, you need to choose one.<br/>
!.md<com>CELL_MODEL:   OHARA INA_PASSINI: INA_PASSINI is an optoinal keyword, it changes the sodium current as described in E Passini, J Mol Cell Cardiol 2016(http://dx.doi.org/10.1016/j.yjmcc.2015.09.003). Recommended to use, it add stability to the cell model.
!.md<5>CM_CONST=  1.0 $[muF/cm2]
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   OHARA INA_PASSINI 
!.md<5>CM_CONST=  1.0 $[muF/cm2]
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   TOROR    ! Tor cell model (collaboration with Blanca Rodriguez)
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   FENTO [BEELE | MODBE | LUORU | GIROU]  ! <b>TODO: DESCRIBE, TEST</b>
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   TENTU    !Tentuscher model, <b>TODO: DESCRIBE, TEST</b>
!.md<5>CM_CONST=  1.0 $[muF/cm2]
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   SCATR | SCVEN  ! <b>TODO: DESCRIBE, TEST</b>
!.md<5>CM_CONST=  1.0 $[muF/cm2]
!.md<5>PACED   
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:   COURTEMANCHE  ! Cell model for the atria (M.Courtemanche, R.J Ramirez and S. Nateel, , 1998)
!.md<5>CM_CONST=  1.0 $[muF/cm2]
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>CELL_MODEL:  FITZH     !Fitzhugh Nagumo cell model 
!.md<field>CELL_MODEL: FITZH
!.md<com>Fitzhugh Nagumo cell model 
!.md<com>As per  Rogers & McCulloch 1994 <br/>
!.md<com>If the stimulus is CURDE, it must be a negative value (e.g. -30) <br/>
!.md<com>Calcium model: cai = ca0 + (cam-ca0) * ( (t-CADEL)/TAU )*exp(1.0_rp-((t-CADEL)/TAU)) from Hunter-McCulloh 1998 (Hunter model)<br/>
!.md<com>Calcium model parameters (cam and Ca0) are in units of Exmedi (1000 times smaller than the values published in the paper). 
!.md<5>CM_CONST=   1.0   $[microF/cm2]
!.md<5>A          = real  
!.md<5>B          = real
!.md<5>C1         = real
!.md<5>C2         = real
!.md<5>D          = real
!.md<5>TIME_SCALE = real, Time scale. Normal FHN lasts 1s. This will change the duration of the pulse, if 0.5, duration is 0.5
!.md<5>VMIN       = real, min voltage, to convert FNH voltage range [0,1] to [VMIN,VMAX]
!.md<5>VMAX       = real, max voltage, to convert FNH voltage range [0,1] to [VMIN,VMAX]
!.md<5>TAU        = real, tau for the calcium model
!.md<5>CAM        = real, cam for the calcium model, exmedi untis
!.md<5>CA0        = real, Ca0 for the calcium model, exmedi units
!.md<5>CADEL      = real, delay in seconds after activation when calcium model kicks in
!.md<5>CATHR      = real, voltage threshold (in FNH scale [0,1]) to decide if activation happened
!.md<4>END_CELL_MODEL
!.md<4>
!.md<4>INI_CELL       ! CURRENTS PARAMETERS AND CONSTANTS
!.md<5>BEATS int      ! Number of beats to run to convergence. 1000 is a good value
!.md<5>TIMESTEP float ! timestep for the 0d cell model initialization run, in seconds. If unspecified, the model uses whatever is hardcoded
!.md<5>CYCLE int      ! Cardiac cycle duration in miliseconds, should be integer
!.md<5>CELLTYPES = ENDO MID EPI           !For which celltypes calcuate initial conditions. By default all (better to use all)
!.md<5>IGNORE_STEADYSTATE = ENDO MID EPI ALL !For the specified celltypes Alya will NOT terminate if steady state  is not reached. Only warn
!.md<5>STEADY_STATE_VARIABLE = (CALCIUM|VOLTAGE) [TOLERANCE = real]
!.md<5>MYOCYTE: MODIFIED 
!.md<field>MYOCYTE
!.md<com>Defines the myocyte type. For now this field is not used. Keep it at "MODIFIED" for now
!!!!<5>MYOCYTE: (NORMAL|MODIFIED) [ |MALE|FEMALE|PIG]
!!!!<field>MYOCYTE
!!!!<com>Defines the myocyte type. This is used essentially to load precalcuakted initial conditions. <br/>
!!!!<com>To force calculation of the initial conditions, just use "MYOCYTE MODIFIED" <br/>
!!!!<com>There are precalculated initial comnditions (that do not trigger calculation.):
!!!!<com>
!!!!<com> - CYCLE 1000, MYOCYTE NORMAL, CELL_MODEL not INA_PASSINI -- uses initial conditions from the O'Hara paper
!!!!<com> - CYCLE 857, MYOCYTE NORMAL, CELL_MODEL INA_PASSINI -- uses initial conditions from for a normal myocyte with 70bpm
!!!!<com> - CYCLE 600, MYOCYTE NORMAL MALE, CELL_MODEL INA_PASSINI -- uses initial conditions from for a normal male  myocyte with 600ms cardiac cycle
!!!!<com> - CYCLE 600, MYOCYTE NORMAL FEMALE, CELL_MODEL INA_PASSINI -- uses initial conditions from for a normal male  myocyte with 600ms cardiac cycle
!!!!<com>
!.md<5>CONDUCTANCES  
!.md<6>INAME str   ! Name for the conductance table. 
!.md<6>float, float, float $---1 endo, mid, epi conductance scale ito channel (GIto) wrt to the one hardcoded in Alya
!.md<6>float, float, float $---2 endo, mid, epi conductance scale slow potassium channel (GKs) 
!.md<6>float, float, float $---3 endo, mid, epi conductance scale potassium channel (GK1)  
!.md<6>float, float, float $---4 endo, mid, epi conductance scale rapid potassium channel (GKr) 
!.md<6>float, float, float $---5 endo, mid, epi conductance scale sodium channel (GNa)  
!.md<6>float, float, float $---6 endo, mid, epi conductance scale late sodium channel (GNaL)  
!.md<6>float, float, float $---7 endo, mid, epi conductance scale sodium calcium exchanger (GNaCa) 
!.md<6>float, float, float $---8 endo, mid, epi conductance scale background potassium background current conductance (gKb) 
!.md<6>float, float, float $---9 endo, mid, epi conductance scale L-type calcium current permeability (pCa equiv to gCaL)  
!.md<6>float, float, float $---10 endo, mid, epi conductance scale sodium potassium pump current permeability (pNaK)
!.md<6>float, float, float $---11 endo, mid, epi conductance scale Calmodulin
!.md<6>float, float, float $---12 endo, mid, epi conductance scale Serca Pump (Jup)
!.md<6>float, float, float $---13 scale for Sodium, Calcium, Potassium (wrt to Na=140 mmol/l, Ca=1.8 mmol/l, K=5.4 mmol/l hardcoded in alya)
!.md<6>float, float, float $---14 endo, mid, epi value (not scale) for Ik_atp. 
!.md<6>float, float, float $---16 endo, mid, epi HFjleak ( scale factor for jleak )
!.md<6>float, float, float $---15 endo, mid, epi scale for JrelNP aproximate Ca sensitivity
!.md<6>float, float, float $---17 endo, mid, epi scale factor for tauHL (wrt to tauHL=200)
!.md<5>END_CONDUCTANCES
!.md<5>DRUGS NUMBER int !number of drugs below
!.md<6>DOSIS NAME char1  !name of the drug, 5 characters max, drug effect: g = g_control ( 1 + (dosis/ic50)^h )^(-1)
!.md<7>float float float ! triplet: dosis, ic50, h for L-type Calcium current 
!.md<7>float float float ! triplet: dosis, ic50, h for Rapid Potassium Current 
!.md<7>float float float ! triplet: dosis, ic50, h for Sodium Channel 
!.md<7>float float float ! triplet: dosis, ic50, h for IK1 Channel 
!.md<7>float float float ! triplet: dosis, ic50, h for INaL Channel 
!.md<7>float float float ! triplet: dosis, ic50, h for IKs Channel 
!.md<7>float float float ! triplet: dosis, ic50, h for Ito Channel 
!.md<6>END_DOSIS 
!.md<6>DOSIS NAME char2  !name of the drug, 5 characters max
!.md<6>...
!.md<6>END_DOSIS 
!.md<6>...
!.md<5>END_DRUGS
!.md<4>END_INI_CELL
!.md<4>
!.md<3>END_MATERIAL
!.md<3>
!.md<3>ISOCH HIGHER real   ! Define how to save the isochrones. 
!.md<3>REPOL THRESHOLD=real DELAY=real  ! Define how to save the repolarisation time. 
!.md<field>ISOCH
!.md<com>Defines the INTRA threshold for the isochrones calculation, records the time INTRA crosses threshold in an upstroke. The iosochrones are reset and recalculated after their postprocessing<br/>
!.md<com>Normal use: ISOCH HIGHER 0.0<br/> 
!.md<field>ISOCH
!.md<com>Defines the REPOL threshold and delay for the isochrones calculation, records the time INTRA crosses threshold in an downstroke after DELAY seconds after the upstroke. The times are reset and recalculated after their postprocessing<br/>
!.md<3>
!.md<3>APEXT FIELD = int  ! apex base heterogeneity, a scalar field defined on nodes, each value in [0,1]
!.md<3>
!.md<1>END_PROPERTIES
!.md<1>
!.md<1>CREATE_STREETER ! see the test for streeter fiber generation, generates streeter fibers
!.md<3>FUNCT CUBIC | LINEA  (default linear)
!.md<3>VAXIS real real real
!.md<3>EPICA int
!.md<3>LEFTE int
!.md<3>RIGHT int
!.md<3>ANGLE value1 value2
!.md<1>END_CREATE_STREETER
!.md<1>
!.md<0>END_PHYSICAL_PROBLEM
!.md</code>



