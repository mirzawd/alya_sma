!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    readat.f90
!> @date    16/11/1966
!> @author  Mariano Vazquez
!> @brief   Reading data
!> @details Reading data
!> @}
!------------------------------------------------------------------------
subroutine readat()
  use def_kintyp
  use def_master
  use def_inpout
  use def_domain
  use def_solver
  use mod_memory, only : memory_alloca
  use mod_ecoute, only : ecoute
  use mod_ecoute, only : ecoute_set_read_unit
  use mod_ecoute, only : ecoute_set_write_unit
  implicit none
  integer(ip) :: imodu,ipara,jmodu,itime,kblok

  !--><group>
  !-->        <groupName>run_data</groupName>
  !-->        <subGroup>
  !-->            <!--inputLine se pinta con un checkbox donde se puede seleccionar o no la opcion-->
  !-->            <inputLine>
  !-->                <inputLineName>ALYA</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <!--inputs que tiene la linea, se tienen que rellenar en caso de chequear el inputLine-->
  !-->                <inputElement>
  !-->                    <inputElementType>edit2</inputElementType>
  !-->                    <inputLineEditValue>cavhex</inputLineEditValue>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->        </subGroup>
  !-->    </group>
  !-->    <group>
  !-->        <groupName>problem_data</groupName>
  !-->        <subGroup>
  !-->            <inputLine>
  !-->                <inputLineName>Time_coupling</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>combo</inputElementType>
  !-->                    <item>
  !-->                        <itemName>Global</itemName>
  !-->                    </item>
  !-->                    <item>
  !-->                        <itemName>Default</itemName>
  !-->                    </item>
  !-->                </inputElement>
  !-->                <inputElement>
  !-->                    <inputElementType>combo</inputElementType>
  !-->                    <item>
  !-->                        <itemName>from_critical</itemName>
  !-->                    </item>
  !-->                    <item>
  !-->                        <itemName>Default</itemName>
  !-->                    </item>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->            <inputLine>
  !-->                <inputLineName>Time_interval</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>edit2</inputElementType>
  !-->                    <inputLineEditValue>cavhex</inputLineEditValue>
  !-->                </inputElement>
  !-->                <inputElement>
  !-->                    <inputElementType>edit2</inputElementType>
  !-->                    <inputLineEditValue>cavhex</inputLineEditValue>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->            <inputLine>
  !-->                <inputLineName>Time_step_size</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>edit2</inputElementType>
  !-->                    <inputLineEditValue>cavhex</inputLineEditValue>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->            <inputLine>
  !-->                <inputLineName>Number_of_steps</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>edit2</inputElementType>
  !-->                    <inputLineEditValue>cavhex</inputLineEditValue>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->        </subGroup>
  !-->        <subGroup>
  !-->            <subGroupName>nastin_module</subGroupName>
  !-->            <subGroupType>module</subGroupType>
  !-->        </subGroup>
  !-->        <subGroup>
  !-->            <subGroupName>parall_service</subGroupName>
  !-->            <subGroupType>service</subGroupType>
  !-->            <inputLine>
  !-->                <inputLineName>Output_file</inputLineName>
  !-->            </inputLine>
  !-->            <inputLine>
  !-->                <inputLineName>Postprocess</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>combo</inputElementType>
  !-->                    <item>
  !-->                        <itemName>Master</itemName>
  !-->                    </item>
  !-->                    <item>
  !-->                        <itemName>Default</itemName>
  !-->                    </item>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->            <inputLine>
  !-->                <inputLineName>Partition_type</inputLineName>
  !-->                <inputLineHelp>Falta introducir la ayuda</inputLineHelp>
  !-->                <inputElement>
  !-->                    <inputElementType>combo</inputElementType>
  !-->                    <item>
  !-->                        <itemName>Faces</itemName>
  !-->                    </item>
  !-->                    <item>
  !-->                        <itemName>Default</itemName>
  !-->                    </item>
  !-->                </inputElement>
  !-->            </inputLine>
  !-->       </subGroup>
  !-->    </group>
  if( INOTSLAVE ) then
     !
     ! Time coupling
     !
     timei         = 0.0_rp                                 ! Initial time
     timef         = 1.0e9_rp                               ! Final time
     !mitim         = 0                                      ! Max. number of steps( initialized elsewhere as this can be given as a command option)
     dtime         = 0.0_rp                                 ! Time step size dt
     kfl_timco     = 1                                      ! Time coupling strategy
     kfl_timei     = 0                                      ! Problem is stationary
     kfl_timef     = 0                                      ! Time step function
     kfl_dtfun     = 0                                      ! Time step defined as a piecewise function
     micou         = 1                                      ! Do at least one global iteration
     kfl_lumped    = 0                                      ! Consistent time evolution matrix (not lumped)
     !
     ! Mesh refinement
     !
     mitrf         = 0                                      ! Maximum # of refinement iterations
     mitsm         = 0                                      ! Maximum # of smoothing iterations
     !
     ! Modules
     !
     kfl_modul(1:mmodu-1) = 0                               ! Modules are off
     lmord         = 0                                      ! Order of problems solutions
     !
     ! Others
     !
     kfl_wwork     = 1                                      ! Master works

     !-------------------------------------------------------------------
     !
     ! Read/write unit
     !
     !-------------------------------------------------------------------

     call ecoute_set_read_unit (lun_pdata) ! Reading file
     call ecoute_set_write_unit(lun_outpu) ! Writing file
     !
     ! Reach the section
     !
     call ecoute('readat')
     if ( words(1) /= 'PROBL')&
          call runend('readat: WRONG PROBLEM_DATA CARD')
     call ecoute('readat')
     !
     ! Read data
     !
     do while( words(1) /= 'ENDPR')

        !-------------------------------------------------------------
        !
        ! Main problem data
        !
        !-------------------------------------------------------------

        do modul = 1,mmodu                                      ! Do not read modules
           if(words(1) == namod(modul)(1:5)) then
              do while( words(1) /= 'END'//namod(modul)(1:2))
                 call ecoute('readat')
              end do
           end if
        end do
        modul = 0

        ! NEW WAY OF DESCRIBING TIME EVOLUTION (THE OLD ONE IS STILL HERE FOR RETROCOMPATIBILITY)

        if( words(1) == 'TIMEE' ) then    ! TIME_EVOLUTION
           ! Default values
           mitim= 1   ! do one time step
           do while( words(1) /= 'ENDTI')
              if( words(1) == 'INTER' ) then         ! INTERVAL
                 timei=param(1)
                 timef=param(2)
              else if( words(1) == 'TIMES' ) then    ! TIMESTEP
                 dtime=param(1)
                 if(dtime>0.0_rp) dtinv = 1.0_rp/dtime
                 if(exists('PRESC')) then
                    kfl_timco=0                                ! prescribed dt
                 else if(exists('FROMC')) then
                    kfl_timco=1                                ! dt=min(dt,f_c*dt_c)
                 else if(exists('PERMO')) then
                    kfl_timco=1                                ! dt=min(dt,f_c*dt_c)
                 end if
              else if( words(1) == 'STEPS' ) then    ! STEPS
                 mitim = int(param(1))
              end if
              call ecoute('readat')
           end do
           
        else if( words(1) == 'PARAL' ) then                         
           !
           ! Read Parall service options
           !
           call par_reapro()

        else if( words(1) == 'TIMEI' ) then                          ! Time interval
           timei=param(1)
           timef=param(2)

        else if( words(1) == 'TIMES' ) then                     ! Time step size
           dtime=param(1)
           if(dtime>0.0_rp) dtinv = 1.0_rp/dtime
           if(exists('FUNCT'))&
                kfl_timef=getint('FUNCT',0_ip,'#Time step function')

        else if( words(1) == 'MAXIM' ) then                     ! Maximum number of iterations
           micou = int(param(1:mblok))
           
        elseif(exists('FUNCT')) then
            if (kfl_timco.ne.0) call runend('kernel/master/readat: IMPOSING AN EXTERNAL FUNCTION IS PRESCRIBING THE TIME STEP. PRESCRIBED OPTION MUST BE USED')
            if (words(2).eq.'SPACET') then
               kfl_timef=getint('SPACET',0_ip,'#Time step function')
            elseif(words(2).eq.'DISCR')then
               kfl_timef=-1
            else
               call runend('READAT: TYPE OF FUNCTION NOT DEFINED')
            endif

        else if( words(1) == 'DTIME' ) then                     ! Time step defined as piecewise constant function
           call ecoute('readat')
           kfl_dtfun = 1
           nfundt = 0
           do while(words(1) /= 'ENDDT')
              if(words(1) == 'TIMES') then
                 if (words(2) == 'DISCR') then
                    kfl_timef= -1                               ! ... because the time step is defined with a table
                    call ecoute('readat')
                    if (words(1) == 'SHAPE') then
                       call ecoute('readat')
                       nfundt = int(param(1))
                       itime = 0
                       call memory_alloca(memma,'dtfun','readat',dtfun,2_ip,nfundt)
                       call ecoute('readat')
                       do while (words(1) /= 'ENDSH')
                          itime = itime + 1
                          dtfun(1,itime) = param(1) ! Time
                          dtfun(2,itime) = param(2) ! Step size
                          call ecoute('readat')
                       end do
                       dtime = dtfun(2,1)           ! First time step
                       if(dtime>0.0_rp) dtinv = 1.0_rp/dtime
                    end if
                 end if
              end if
              if ( dtfun(1,nfundt) /= timef ) call runend('READAT: FINAL TIME DEFINED IN DTIME MUST BE EQUAL TO TOTAL TIME')
              call ecoute('readat')
           end do

        else if( words(1) == 'TIMEC' ) then                     ! Time coupling strategy
           if(words(2)=='GLOBA' ) then
              kfl_timco=0
              if(exists('PRESC')) then
                 kfl_timco=0                                ! prescribed dt
              else if(exists('FROMC')) then
                 kfl_timco=1                                ! dt=min(dt,f_c*dt_c)
              end if
           else if(words(2)=='LOCAL' ) then
              kfl_timco=2                                   ! dt=dt(module)
           else if(words(2)=='PERMO' ) then
              kfl_timco=2                                   ! dt=dt(module)
           else if( words(1) == 'OFF  ' ) then
              kfl_timco=-1
           end if
           if(exists('LUMPE')) kfl_lumped=1                 ! Lumped mass matrix evolution
           if(exists('GALER')) kfl_lumped=2                 ! Galerkin lumped mass matrix evolution

        else if( words(1) == 'NUMBE' ) then
           !
           ! Number of time steps... do not read if it was given by command line option
           !
           if( mitim == 0 ) mitim = int(param(1))

        else if( words(1) == 'ZONES' ) then                     ! Zones
           call ecoute('readat')
           do while( words(1) /= 'ENDZO')
              ipara = -1
              modules3: do jmodu=0,mmodu
                 if(words(1)==namod(jmodu)(1:5)&
                      .and.namod(jmodu)/='     ' ) then
                    ipara=jmodu
                    exit modules3
                 end if
              end do modules3
              if( ipara == - 1 ) call runend('ERROR IN ZONES FIELD: WRONG MODULE NAMES')
              lzone(ipara) = int(param(1))
              call ecoute('readat')
           end do

        else if( words(1) == 'BLOCK' ) then
           nblok=getint('BLOCK',1_ip,'#Number of blocks')
           call ecoute('readat')
           do while( words(1) /= 'ENDBL')
              kblok=int(param(1))
              do imodu=1,mmodu
                 ipara=int(param(imodu+1))
                 modules1: do jmodu=1,mmodu
                    if(words(imodu+1)==namod(jmodu)(1:5)&
                         .and.namod(jmodu)/='     ' ) then
                       ipara=jmodu
                       exit modules1
                    end if
                 end do modules1
                 if(kblok==0) call runend('WRONG BLOCK DEFINITION IN DATA FILE')
                 lmord(imodu,kblok)=ipara
              end do
              call ecoute('readat')
           end do

        end if
        call ecoute('readat')

     end do

     if( kfl_timco == 2 ) timei     = 0.0_rp
     if( timef < 0.0_rp ) kfl_timei = 1

  end if

end subroutine readat
