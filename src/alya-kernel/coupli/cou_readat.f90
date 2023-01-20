!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    ker_readat.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Read coupling data
!> @details Read coupling data
!>          Three types of Dirichlet-type coupling, which lead
!>          to the same solution if meshes i and j coincide.
!>          \verbatim
!>
!>          Dirichlet explicit   Dirichlet implicit   Unknown
!>
!>          => x_i = x_j 
!>
!>          +-- do iter          +-- do iter          +-- do iter         
!>          |                    |                    |
!>          | Ax, x.y            ! A.x                ! A.x
!>          ! x.y                ! x.y                ! x.y
!>          |                    | x_i = x_j          | A.x|_i = A.x|_j
!>          |                    |                    |
!>          +-- end do           +-- end do           +-- end do
!>
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------

subroutine cou_readat()
  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_coupli
  use mod_parall, only :  par_code_zone_subd_to_color
  use mod_ecoute, only :  ecoute
  use mod_ecoute, only :  ecoute_set_read_unit
  use mod_ecoute, only :  ecoute_set_write_unit

  implicit none
  integer(ip)          :: icoup,ipara,kblok,ielty,jelty
  character(5)         :: where_type_char

  if( INOTSLAVE .and. lun_coupl_dat /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Read/write unit
     !
     !-------------------------------------------------------------------

     call ecoute_set_read_unit (lun_coupl_dat) ! Reading file
     call ecoute_set_write_unit(lun_coupl_res) ! Writing file

     kfl_timco_cou = 2
     ngaus_elem    = 0_ip
     
     !-------------------------------------------------------------------
     !
     ! Physical problem
     !
     !-------------------------------------------------------------------

     call ecoute('cou_readat')
     do while( words(1) /= 'PHYSI' )
        call ecoute('cou_readat')
     end do
     if( exists('NUMBE') ) mcoup = getint('NUMBE',1_ip,'#NUMBER OF COUPLINGS')
     call ecoute('cou_readat')

     do while( words(1) /= 'ENDPH' )
        !
        ! Physical problem
        !
        if( words(1) == 'NUMBE' ) then

           mcoup = getint('NUMBE',1_ip,'#NUMBER OF COUPLINGS')
           call cou_memory(1_ip)

        else if( words(1) == 'LOSTW' ) then
           !
           ! Lost wet points strategy
           !
           if(      words(2) == 'ZERO ' ) then
              kfl_lost_wet_point_cou = IMPOSE_ZERO
           else if( words(2) == 'DONTD' ) then
              kfl_lost_wet_point_cou = DONT_DO_ANYTHING 
           else if( words(2) == 'STOPA' ) then
              if( exists('WARNI') ) then
                 kfl_lost_wet_point_cou = STOP_ALYA_WITH_WARNINGS
              else
                 kfl_lost_wet_point_cou = STOP_ALYA
              end if
           end if
           
        else if( words(1) == 'ABSOL' ) then
           !
           ! Absolute tolerance for bin bounding boxes
           ! >= 0, take toler_absolute_cou
           ! <  0, take -toler_absolute_cou * (Vol_ave)^{1/ndime}
           ! where Vol_ave is the average mesh volume
           !
           toler_absolute_cou = getrea('ABSOL',-1.0_rp,'#ABSOLUTE TOLERANCE')
           if( toler_absolute_cou < 0.0_rp ) then
              kfl_absolute_cou = -1_ip
           else if( exists('OFF  ') ) then
              kfl_absolute_cou = 0_ip
           else
              kfl_absolute_cou = 1_ip
           end if
           toler_absolute_cou = abs(toler_absolute_cou)
           
        else if( words(1) == 'RELAT' ) then
           !
           ! Relative tolerance for bin bounding boxes
           !
           toler_relative_cou = getrea('RELAT',0.01_rp,'#RELATIVE TOLERANCE')
       
        else if( words(1) == 'PARTR' ) then

           if( mcoup == 0 ) call runend('COU_READAT: SET FIRST THE NUMBER OF COUPLINGS TO MODIFY A TODA COSTA OPTION')
           if(      words(2) == 'OFF  ' ) then
              coupling_type(1:mcoup) % kfl_par_transm = 0
           else if( words(2) == 'ON   ' ) then
              coupling_type(1:mcoup) % kfl_par_transm = 1
           end if

        else if( words(1) == 'CHECK' ) then

           if( mcoup == 0 ) call runend('COU_READAT: SET FIRST THE NUMBER OF COUPLINGS TO MODIFY A TODA COSTA OPTION')
           if(      words(2) == 'OFF  ' ) then
              coupling_type(1:mcoup) % kfl_check_exha = 0
           else if( words(2) == 'ON   ' ) then
              coupling_type(1:mcoup) % kfl_check_exha = 1
           end if

       else if( words(1) == 'ATODA' ) then

           if( mcoup == 0 ) call runend('COU_READAT: SET FIRST THE NUMBER OF COUPLINGS TO MODIFY A TODA COSTA OPTION')
           if(      words(2) == 'OFF  ' ) then
              coupling_type(1:mcoup) % kfl_toda_costa = 0
           else if( words(2) == 'ON   ' ) then
              coupling_type(1:mcoup) % kfl_toda_costa = 1
           end if

        else if( words(1) == 'COUPL' ) then

           icoup = getint('COUPL',1_ip,'#Number of the set or the field')
           if( icoup < 1 .or. icoup > mcoup )    call runend('COU_READAT: WRONG COUPLING NUMBER')
           if( .not. associated(coupling_type) ) call runend('COU_READAT: NUMBER_COUPLING TYPE IS MISSING')
           coupling_type(icoup) % number = icoup
           
           call ecoute('cou_readat')
           do while( words(1) /= 'ENDCO' )

              if(      words(1) == 'TARGE' .and. words(2) == 'MODUL' ) then
                 !
                 ! MODULE_TARGET
                 !
                 coupling_type(icoup) % module_target = idmod(words(3))
                 if( coupling_type(icoup) % module_target < 0 ) call runend('COU_READAT: WRONG TARGET MODULE')
                 coupling_type(icoup) % zone_target   = lzone(coupling_type(icoup) % module_target)

              else if( words(1) == 'SOURC' .and. words(2) == 'MODUL' ) then
                 !
                 ! MODULE_SOURCE 
                 !
                 coupling_type(icoup) % module_source = idmod(words(3))
                 if( coupling_type(icoup) % module_source < 0 ) call runend('COU_READAT: WRONG TARGET MODULE')

              else if( words(1) == 'VARIA' ) then
                 !
                 ! VARIABLE
                 !
                 coupling_type(icoup) % variable = trim(words(2))
                 
              else if( words(1) == 'MULTI' ) then
                 !
                 ! MULTIPLICITY OF SOURCE
                 !
                 if( option('MULTI') ) coupling_type(icoup) % kfl_multi_source =  1
                 if( exists('EXCLU') ) coupling_type(icoup) % kfl_multi_source = -1
                 
              else if( words(1) == 'LOSTW' ) then
                 !
                 ! LOST WET POINTS
                 !
                 if(      words(2) == 'STOP ' ) then
                    coupling_type(icoup) % kfl_lost_wet_points = 0
                 else if( words(2) == 'CONTI' .or. words(2) ==  'DISCA' ) then
                    coupling_type(icoup) % kfl_lost_wet_points = 1
                 end if
                 
              else if( words(1) == 'CONSE' ) then
                 !
                 ! CONSERVATION
                 !
                 if( words(2) == 'LOCAL' .or. words(2) == 'INTER' ) then
                    coupling_type(icoup) % conservation = INTERFACE_MASS
                 else if( words(2) == 'GLOBA' .or. words(2) == 'TOTAL' ) then
                    coupling_type(icoup) % conservation = GLOBAL_MASS
                 end if

              else if( words(1) == 'OVERL' ) then
                 !
                 ! VARIABLE
                 !
                 if( words(2) == 'DISJO' .or. words(2) == 'NO' ) then
                    coupling_type(icoup) % overlap = 0
                 else
                    coupling_type(icoup) % overlap = getint('OVERL',0_ip,'#OVERLAP FOR CHIMERA-TYPE COUPLING')
                 end if

              else if(      words(1) == 'TARGE' .and. words(2) == 'SUBDO' ) then
                 !
                 ! SUBDOMAIN_TARGET
                 !
                 coupling_type(icoup) % subdomain_target = getint('SUBDO',1_ip,'#TARGET SUBDOMAIN')

              else if( words(1) == 'SOURC' .and. words(2) == 'SUBDO' ) then
                 !
                 ! SUBDOMAIN_SOURCE
                 !
                 coupling_type(icoup) % subdomain_source = getint('SUBDO',1_ip,'#SOURCE SUBDOMAIN')

              else if( words(1) == 'TARGE' .and. words(2) == 'CODE ' ) then
                 !
                 ! CODE_TARGET
                 !
                 coupling_type(icoup) % code_target   = int(param(2),ip)

              else if( words(1) == 'SOURC' .and. words(2) == 'CODE ' ) then
                 !
                 ! CODE_SOURCE
                 !
                 coupling_type(icoup) % code_source   = int(param(2),ip)

              else if( words(1) == 'WHERE' ) then
                 !
                 ! WHERE
                 !
                 where_type_char = getcha('WHERE','     ','#Where')
                 if( exists('SOURC') ) then
                    coupling_type(icoup) % where_number_source = getint('NUMBE',1_ip,'#Number of the set or the field')
                    if( where_type_char == 'SET  ' ) then
                       coupling_type(icoup) % where_type_source = ON_SET
                    else if( where_type_char == 'FIELD' ) then
                       coupling_type(icoup) % where_type_source = ON_FIELD
                    else if( where_type_char == 'CODE ' ) then
                       coupling_type(icoup) % where_type_source = ON_CODE
                    else if( where_type_char == 'WHOLE' ) then
                       coupling_type(icoup) % where_type_source = ON_WHOLE_MESH
                    else if( where_type_char == 'CHIME' ) then
                       coupling_type(icoup) % where_type_source = ON_CHIMERA_MESH
                    else if( where_type_char == 'IMMER' ) then
                       coupling_type(icoup) % where_type_source = ON_IMMERSED_MESH
                    else if( where_type_char == 'EMBED' ) then
                       coupling_type(icoup) % where_type_source = ON_EMBEDDED_MESH
                    else 
                       call runend('COU_READAT: WRONG WHERE IN SOURCE')
                    end if
                 else
                    coupling_type(icoup) % where_number = getint('NUMBE',1_ip,'#Number of the set or the field')
                    if( where_type_char == 'SET  ' ) then
                       coupling_type(icoup) % where_type = ON_SET
                    else if( where_type_char == 'FIELD' ) then
                       coupling_type(icoup) % where_type = ON_FIELD
                    else if( where_type_char == 'CODE ' ) then
                       coupling_type(icoup) % where_type = ON_CODE
                    else if( where_type_char == 'WHOLE' ) then
                       coupling_type(icoup) % where_type = ON_WHOLE_MESH
                    else if( where_type_char == 'CHIME' ) then
                       coupling_type(icoup) % where_type = ON_CHIMERA_MESH
                    else if( where_type_char == 'MIRRO' ) then
                       coupling_type(icoup) % where_type = ON_MIRROR
                    else if( where_type_char == 'IMMER' ) then
                       coupling_type(icoup) % where_type = ON_IMMERSED_MESH
                    else if( where_type_char == 'EMBED' ) then
                       coupling_type(icoup) % where_type = ON_EMBEDDED_MESH
                    else if( where_type_char == 'FLOAT' ) then
                       coupling_type(icoup) % where_type = ON_FLOATING_POINTS
                    else 
                       call runend('COU_READAT: WRONG WHERE IN TARGET')
                    end if
                 end if
                 
              else if( words(1) == 'WHAT ' ) then
                 !
                 ! WHAT
                 !
                 if(      words(2) == 'UNKNO' ) then
                    coupling_type(icoup) % what = UNKNOWN
                 else if( words(2) == 'DIRIC' ) then
                    if( words(3) == 'IMPLI' ) then
                       coupling_type(icoup) % what = DIRICHLET_IMPLICIT
                    else
                       coupling_type(icoup) % what = DIRICHLET_EXPLICIT
                    end if
                 else if( words(2) == 'RESID' ) then
                    coupling_type(icoup) % what = RESIDUAL
                 else if( words(2) == 'SCALA' ) then
                    coupling_type(icoup) % what = SCALAR
                 else
                    call runend('COU_READAT: WRONG COUPLING WHAT')
                 end if

              else if( words(1) == 'TYPE ' ) then
                 !
                 ! TYPE
                 !
                 if(      words(2) == 'ELEME' ) then
                    coupling_type(icoup) % itype = ELEMENT_INTERPOLATION
                 else if( words(2) == 'NEARE' ) then
                    if( words(3) == 'ELEME' ) then
                       coupling_type(icoup) % itype = NEAREST_ELEMENT_NODE
                    else
                       coupling_type(icoup) % itype = NEAREST_BOUNDARY_NODE
                    end if
                 else if( words(2) == 'BOUND' ) then
                    coupling_type(icoup) % itype = BOUNDARY_INTERPOLATION
                 else if( words(2) == 'STRES' ) then
                    coupling_type(icoup) % itype = STRESS_PROJECTION
                    if( exists('GAUSS') ) then
                       if( exists('DEFAU') ) then
                          coupling_type(icoup) % ngaus = 0
                       else if( exists('AUTOM') ) then
                          coupling_type(icoup) % ngaus = -1
                          call runend('COU_READAT: AUTOMATIC NOT CODED YET')
                       else
                          coupling_type(icoup) % ngaus = getint('GAUSS',1_ip,'#Number of Gauss points')
                       end if
                    end if
                 else if( words(2) == 'PROJE' ) then
                    coupling_type(icoup) % itype = PROJECTION
                    if( exists('GAUSS') ) then
                       if( exists('DEFAU') ) then
                          coupling_type(icoup) % ngaus = 0
                       else if( exists('AUTOM') ) then
                          coupling_type(icoup) % ngaus = -1
                          call runend('COU_READAT: AUTOMATIC NOT CODED YET')
                       else
                          coupling_type(icoup) % ngaus = getint('GAUSS',1_ip,'#Number of Gauss points')
                       end if
                    end if
                 else if( words(2) == 'TRANS' ) then
                    coupling_type(icoup) % itype = TRANSPOSE_MIRROR
                 else if( words(2) == 'GLOBA' ) then
                    coupling_type(icoup) % itype = GLOBAL_NUMBERING
                 else if( words(2) == 'SAMEC' ) then
                    coupling_type(icoup) % itype = SAME_COORDINATE
                  else
                    call runend('COU_READAT: NON-EXISTING TYPE OF COUPLING')
                 end if

              else if ( words(1) == 'GAUSS' ) then
                 !
                 ! Number of Gauss points for FLOATING
                 !
                 jelty = 0_ip
                 do ielty = 1, nelty
                    if( lexis(ielty) > 0 ) then
                       jelty = jelty + 1
                       ngaus_elem(ielty) = int(param(jelty))
                    end if
                 end do

              else if( words(1) == 'RELAX' ) then
                 !
                 ! RELAX
                 !
                 coupling_type(icoup) % relax = getrea('RELAX',1.0_rp,'#Relaxation') 

              else if( words(1) == 'SCHEM' ) then
                 !
                 ! SCHEME
                 !
                 if( words(2) == 'RELAX' ) then
                    coupling_type(icoup) % scheme = RELAXATION_SCHEME
                 else if( words(2) == 'AITKE' ) then
                    coupling_type(icoup) % scheme = AITKEN_SCHEME
                 else if( words(2) == 'BROYD' ) then
                    coupling_type(icoup) % scheme = BROYDEN_SCHEME
                 else if( words(2) == 'IQNLS' ) then
                    coupling_type(icoup) % scheme = IQNLS_SCHEME
                    coupling_type(icoup) % ranku_iqnls = getint('RANKU',3_ip,'#Number of past iterations stored')
                    coupling_type(icoup) % history_iqnls = getint('HISTO',0_ip,'#Number of past time iterations (history) stored')
                    coupling_type(icoup) % efilter_iqnls = getrea('FILTER',-1.0_rp,'#Epsilon for the CIQN filter')
                    coupling_type(icoup) % scaling_iqnls = getrea('SCALI',1.0_rp,'#scaling for the INQLS')

                 else
                    call runend('COU_READAT: NOT CODED')
                 end if
              else if( words(1) == 'TEMPO' ) then
                 !
                 ! TEMPORAL PREDICTION
                 !
                 if( words(2) == 'OFF' )then
                    coupling_type(icoup) % temporal_predictor = 0_ip
                    coupling_type(icoup) % temporal_predictor_order = -1_ip ! Zeroth order prediction
                 else if( words(2) == 'ZEROT' )then
                    coupling_type(icoup) % temporal_predictor = 1_ip
                    coupling_type(icoup) % temporal_predictor_order = 0_ip ! Zeroth order prediction

                 else if( (words(2) == 'ON') .or. (words(2)=='SECOND') )then
                    coupling_type(icoup) % temporal_predictor = 1_ip
                    coupling_type(icoup) % temporal_predictor_order = 2_ip ! Second order prediction
                 else
                   call runend('COU_READAT: TEMPORAL PREDICTION ORDER NOT CODED')
                 end if

              else if( words(1) == 'SENDA' )then
                 !
                 ! COMPUTE & SEND
                 !
                 coupling_type(icoup) % frequ_send = 1_ip ! Initialization of frequency of send
                 if(      words(2) == 'INIUN' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_INIUNK
                 else if( words(2) == 'DOITE' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_DOITER
                 else if( words(2) == 'BEGIT' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_BEGITE
                 else if( words(2) == 'ENDIT' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_ENDITE
                 else if( words(2) == 'TURNO' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_TURNON
                 else if( words(2) == 'TURNF' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_TURNOF
                 else if( words(2) == 'CONCO' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_CONCOU
                 else if( words(2) == 'BEGZO' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_BEGZON
                 else if( words(2) == 'ENDZO' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_ENDZON
                 else if( words(2) == 'BEGST' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_BEGSTE
                 else if( words(2) == 'ENDST' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_ENDSTE
                 else if( words(2) == 'TIMST' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_TIMSTE
                 else if( words(2) == 'ENDTI' ) then
                    coupling_type(icoup) % task_compute_and_send = ITASK_ENDTIM
                 end if

                 if( words(3) == 'BEFOR' ) then
                    coupling_type(icoup) % when_compute_and_send = ITASK_BEFORE
                 else if( words(3) == 'AFTER' ) then
                    coupling_type(icoup) % when_compute_and_send = ITASK_AFTER
                 end if

                 if( words(4) == 'FREQU' )then
                    !
                    ! Frequency of send for coupling
                    !
                    coupling_type(icoup) % frequ_send = getint('FREQU',1_ip,'#Frequency of send')
                 end if

              else if( words(1) == 'RECEI' ) then
                 !
                 ! RECEIVE & ASSEMBLE
                 !
                 coupling_type(icoup) % frequ_recv = 1_ip ! Initialization of frequency of recv
                 if(      words(2) == 'INIUN' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_INIUNK
                 else if( words(2) == 'DOITE' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_DOITER
                 else if( words(2) == 'BEGIT' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_BEGITE
                 else if( words(2) == 'ENDIT' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_ENDITE
                 else if( words(2) == 'TURNO' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_TURNON
                 else if( words(2) == 'TURNF' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_TURNOF
                 else if( words(2) == 'BEGZO' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_BEGZON
                 else if( words(2) == 'ENDZO' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_ENDZON
                 else if( words(2) == 'BEGST' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_BEGSTE
                 else if( words(2) == 'ENDST' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_ENDSTE
                 else if( words(2) == 'TIMST' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_TIMSTE
                 else if( words(2) == 'ENDTI' ) then
                    coupling_type(icoup) % task_recv_and_assemble = ITASK_ENDTIM
                 end if

                 if( words(3) == 'BEFOR' ) then
                    coupling_type(icoup) % when_recv_and_assemble = ITASK_BEFORE
                 else if( words(3) == 'AFTER' ) then
                    coupling_type(icoup) % when_recv_and_assemble = ITASK_AFTER
                 end if
                 if( words(4) == 'FREQU' )then
                    !
                    ! Frequency of receive for coupling
                    !
                    coupling_type(icoup) % frequ_recv = getint('FREQU',1_ip,'#Frequency of recv')
                 end if

              else if( words(1) == 'FIXIT' ) then
                 !
                 ! FIXITY 
                 !
                 if(      words(2) == 'OFF  ' ) then
                     coupling_type(icoup) % kfl_fixity = 0
                 else 
                     coupling_type(icoup) % kfl_fixity = 1
                 end if

              end if
              call ecoute('cou_readat')
           end do

        end if
        call ecoute('cou_readat')

     end do

     !-------------------------------------------------------------------
     !
     ! Numerical treatment
     !
     !-------------------------------------------------------------------

     do while( words(1) /= 'NUMER' )
        call ecoute('cou_readat')
     end do
     call ecoute('cou_readat')
     do while( words(1) /= 'ENDNU' )

        if( words(1) == 'BLOCK' ) then
           !
           ! Block iterations
           !
           if( exists('NUMBE') ) then
              kblok = getint('NUMBE',1_ip,'#Number of the block')
           else
              kblok = int(param(1),ip)
           end if
           if( kblok < 1 .or. kblok > nblok ) call runend('COU_READAT: WRONG BLOCK NUMBER')
           if( kblok > max_block_cou )        call runend('COU_READAT: WRONG BLOCK NUMBER')

           call ecoute('cou_readat')
           do while( words(1) /= 'ENDBL' )

              if(      words(1) == 'COUPL' ) then            
                 if( nnpar > max_coupl_cou ) call runend('COU_READAT: WRONG NUMBER OF COUPLING FOR BLOCK DEFINITION')
                 coupling_driver_number_couplings(kblok) = nnpar
                 do ipara = 1,nnpar
                    coupling_driver_couplings(ipara,kblok) = int(param(ipara),ip)
                 end do
              else if( words(1) == 'ITERA' ) then
                 coupling_driver_max_iteration(kblok) = getint('ITERA',1_ip,'#Number of iterations')
              else if( words(1) == 'TOLER' ) then
                 coupling_driver_tolerance(kblok) = getrea('TOLER',1.0e-2_rp,'#Tolerance')
              end if

              call ecoute('cou_readat')
           end do

        else if( words(1) == 'TIMEC' ) then
           !
           ! Time coupling strategy
           !
           if(     words(2) == 'GLOBA' ) then
              kfl_timco_cou = 0
           else if(words(2) == 'LOCAL' ) then
              kfl_timco_cou = 2                                   ! dt=dt(module)
           else
              call runend('UNKNOWN COUPLING TIME STRATEGY')
           end if
           
        end if
        call ecoute('cou_readat')
        
     end do
     !
     ! Block ordering has not been defined
     !
     if( coupling_driver_number_couplings(1) == 0 ) then
        coupling_driver_number_couplings(1) = mcoup
        do icoup = 1,mcoup
           coupling_driver_couplings(icoup,1) = icoup
        end do
     end if
     
  end if
  
end subroutine cou_readat
 
