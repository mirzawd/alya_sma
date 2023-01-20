!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_reaous.f90
!> @author  Guillaume Houzeaux
!> @brief   Read postprocess data
!> @details Read postprocess data
!!          - Array postprocess
!!          - Witness points
!!          - Element, boundary and node sets
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_reaous
  !-----------------------------------------------------------------------
  !.md<module>nastin
  !.md<input>case.nsi.dat
  !.md<pos>2
  !.md<sec>

  use def_parame
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_output_postprocess, only : output_postprocess_read
  use mod_output_postprocess, only : output_postprocess_sets_ordering
  implicit none
  integer(ip) :: ii

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     kfl_exacs_nsi = 0            ! Exact solution
     kfl_exfix_nsi = 1            ! Dirichlet condition imposed 
     kfl_inert_nsi = 0            ! Non-inertial velocity 
     kfl_psmat_nsi = 0            ! Postscript file of the matrix
     kfl_miniapp_nsi = 0          ! No miniapp output
     cloth_nsi     = -500.0_rp    ! CLO
     metab_nsi     = -500.0_rp    ! MET
     wetme_nsi     = -500.0_rp    ! WME
     ambie_nsi     = -500.0_rp    ! TA
     radia_nsi     = -500.0_rp    ! TR
     relat_nsi     = -500.0_rp    ! RH
     avtim_nsi     = 0.0_rp       ! Start averaging time
     avste_nsi     = 0_ip         ! Averaging window in steps
     !--><group>
     !-->            <groupName>OUTPUT_&amp;_POST_PROCESS</groupName>
     !-->            <groupType>oupos</groupType>
     !-->            <groupTypeDefinition>
     !-->        <![CDATA[    POSTPROCESS XXXXX, STEPS=int1 [, FILTER=int2]                                                   $ Variables computed on all nodes or filtered
     !-->                     ELEMENT_SET
     !-->                       XXXXX                                                                                         $ Variables computed in volume integrals
     !-->                     END_ELEMENT_SET
     !-->                     BOUNDARY_SET
     !-->                       XXXXX                                                                                         $ Variables computed in boundary integrals
     !-->                     END_BOUNDARY_SET
     !-->                     NODE_SET
     !-->                       XXXXX                                                                                         $ Variables computed on nodes
     !-->                     END_NODE_SET
     !-->                     WITNESS_POINTS
     !-->                       XXXXX                                                                                         $ Variables computed witness points
     !-->                     END_WITNESS_POINTS
     !-->                     OUTPUT ERROR, SOLUTION = int1                                                                   $ Manufactured solution
     !-->                  ]]>
     !-->      </groupTypeDefinition>
     !-->            <postProcess>
     !-->                <inputLineHelp><![CDATA[POSTPROCESS:
     !-->                  Postprocess variables on nodes at each int1 time steps. the name of the file
     !-->                  where the variables is stored at time step nnn is: <tt> sample-XXXXX-00000nnn.post.alyabin</tt>.
     !-->                  A fliter can be applied to the variable. Filters are defined in kermod.]]></inputLineHelp>
     !-->                <var>VELOCITY_MODULE</var>
     !-->                <var>VORTICITY</var>
     !-->                <var>KINETIC_ENERGY</var>
     !-->            </postProcess>
     !-->            <elementSet>
     !-->                <inputLineHelp><![CDATA[ELEMENT_SET:
     !-->                  Variables computed as volume integrals on element sets. The results are stored
     !-->                  in file <tt>sample-element.nsi.set</tt>. Let V be the size of the bounday set. 
     !-->                  The following variables can be postprocessed:
     !-->                  
     !-->                        -VELOCITY_MODULE: int_V u^2 
     !-->                        -VORTICITY:       int_V (dv/dx-du/dy)^2
     !-->                        -KINETIC_ENERGY:  int_V 1/2*rho*u^2 
     !-->                  ]]></inputLineHelp>
     !-->                <var>VELOCITY_MODULE</var>
     !-->                <var>VORTICITY</var>
     !-->                <var>KINETIC_ENERGY</var>
     !-->            </elementSet>
     !-->            <boundarySet>
     !-->                <inputLineHelp><![CDATA[BOUNDARY_SET:
     !-->                 Variables computed as boundary integrals on boundary sets. The results are stored
     !-->                 in file <tt>sample-boundary.nsi.set</tt>. Let S be the size of the boundary set. 
     !-->                 The following variables can be postprocessed:
     !-->                 
     !-->                      -MEAN_PRESSURE: int_S P/|S| 
     !-->                      -MASS:          int_S rho*u.n
     !-->                      -FORCE:         int_S 2 mu eps(u).n, int_S (-pI).n
     !-->                      -TORQUE:        int_S (r-rc) x ( 2 mu eps(u).n ), int_S  (r-rc) x (-pI.n) 
     !-->                      -MEAN_YPLUS:    int_S y+ / S
     !-->                      -MEAN_VELOCITY: int_S u.n / S
     !-->                      -WEAT_FORCE:    int_Sw 2 mu eps(u).n, int_Sw (-pI).n 
     !-->                      -WEAT_SURFACE:  int_Sw, Sw=wet surface
     !-->                ]]></inputLineHelp>
     !-->                <var>MEAN_PRESSURE</var>
     !-->                <var>MASS</var>
     !-->                <var>FORCE</var>
     !-->                <var>TORQUE</var>
     !-->                <var>MEAN_YPLUS</var>
     !-->                <var>MEAN_VELOCITY</var>
     !-->                <var>WEAT_FORCE</var>
     !-->                <var>WEAT_SURFACE</var>
     !-->            </boundarySet>
     !-->            <nodeSet>
     !-->                <inputLineHelp><![CDATA[NODE_SET:
     !-->                 Variables computed as boundary integrals on boundary sets. The results are stored
     !-->                 in file <tt>sample-node.nsi.set</tt>. 
     !-->                 The following variables can be postprocessed:
     !-->                 
     !-->                       -VELOX: X-velocity
     !-->                       -VELOY: Y-velocity
     !-->                       -VELOZ: Z-velocity
     !-->                       -PRESS: Pressure
     !-->                       -YPLUS: y+
     !-->                 
     !-->                 WITNESS_POINTS
     !-->                    XXXXX                                                       $ Variables computed witness points
     !-->                 END_WITNESS_POINTS
     !-->                 WITNESS_POINTS:
     !-->                 
     !-->                       -VELOX: X-velocity
     !-->                       -VELOY: Y-velocity
     !-->                       -VELOZ: Z-velocity
     !-->                       -PRESS: Pressure
     !-->                 ]]></inputLineHelp>
     !-->                <var>VELOX</var>
     !-->                <var>VELOY</var>
     !-->                <var>VELOZ</var>
     !-->                <var>PRESS</var>
     !-->                <var>YPLUS</var>
     !-->            </nodeSet>
     !-->            <witnessPoints>
     !-->              <inputLineHelp><![CDATA[
     !-->                    -VELOX: X-velocity
     !-->                    -VELOY: Y-velocity
     !-->                    -VELOZ: Z-velocity
     !-->                    -PRESS: Pressure
     !-->                ]]></inputLineHelp>
     !-->                <var>VELOX</var>
     !-->                <var>VELOY</var>
     !-->                <var>VELOZ</var>
     !-->                <var>PRESS</var>
     !-->            </witnessPoints>
     !-->            <outputError>
     !-->               <inputLineHelp><![CDATA[OUTPUT ERROR:
     !-->                Manufactured solutions used to carry out code typing checking and mesh convergence.
     !-->                The L2, L1 and Linf norm of the velocity, pressure, velocity gradients and pressure 
     !-->                gradients can be found in file <tt>sample.nsi.log</tt>. Dirhclet boundary
     !-->                conditions are automatically imposed on the velocity on all the boundary. Pressure
     !-->                should be imposed on one node to the right value as the flow is confined.
     !-->                The following manufactured solutions are available:
     !-->                
     !-->                      - SOLUTION = 5:  r=sqrt(x^2+y^2); u = 2*y*(r-0.5), v = -2*x*(r-0.5),  p = r-1 
     !-->                      - SOLUTION = 13: u=2x +  y + 3z, v=- x - 3y - 4z, w=6x + 7y - z, p=-5x + 6y + 2z 
     !-->                      - SOLUTION = 20: u=-cos(pi*x)*sin(pi*y), v=sin(pi*x)*cos(pi*y), p=sin(pi*x)*sin(pi*x) + sin(pi*y)*sin(pi*y)
     !-->                ]]></inputLineHelp>
     !-->            </outputError>
     !-->            <turbulenceAverages>
     !-->               <inputLineHelp><![CDATA[TURBULENCE AVERAGES:
     !-->                Some requiremenets so that turbulent averages work fine.
     !-->                   1) use constant time step
     !-->                   2) in .dat RUN_TYPE    Frequency has to be idem to the postprocesing frequency of average values.
     !-->                   3) The postprocesisng frequency of all average values has to be the same (my recomendation is that you use STEPO  in ker.dat)
     !-->                   4) If you postprocess averages, for each of the modules where you postprocess averages you must postprocess at least the corresponding value
     !-->                       nastin : AVVEL
     !-->                       nastal : AVVEL
     !-->                       turbul : AVTVI
     !-->                       temper : AVTEM
     !-->                       (see  'if( mod(ittim, postp(1) % npp_stepi(arrays_number('AVVEL') ) == 0 ) then'  in nsi_output and similar lines in the rest of modules.)
     !-->            </turbulenceAverages>
     !
     !.md# Output and Post-process
     !.md<code>
     !.md<0>OUTPUT_&_POST_PROCESS
     !
     !.md<1>POSTPROCESS XXXXX, STEPS=int1 [, FILTER=int2]                                                   $ Variables computed on all nodes or filtered
     !.md<field>POSTPROCESS
     !.md<com>Postprocess variables on nodes at each int1 time steps. the name of the file
     !.md<com>where the variables is stored at time step nnn is: <tt> sample-XXXXX-00000nnn.post.alyabin</tt>.
     !.md<com>A fliter can be applied to the variable. Filters are defined in kermod.
     !
     !.md<1>ELEMENT_SET
     !.md<1>  XXXXX                                                                                         $ Variables computed in volume integrals
     !.md<1>END_ELEMENT_SET
     !.md<field>ELEMENT_SET
     !.md<com>Variables computed as volume integrals on element sets. The results are stored
     !.md<com>in file <tt>sample-element.nsi.set</tt>. Let V be the size of the bounday set. 
     !.md<com>The following variables can be postprocessed:
     !.md<com>
     !.md<com>    - VELOCITY_MODULE: int_V u^2 
     !.md<com>    - VORTICITY:       int_V (dv/dx-du/dy)^2
     !.md<com>    - KINETIC_ENERGY:  int_V 1/2*rho*u^2 
     !.md<com>
     !
     !.md<1>BOUNDARY_SET
     !.md<1>  XXXXX                                                                                         $ Variables computed in boundary integrals
     !.md<1>END_BOUNDARY_SET
     !.md<field>BOUNDARY_SET
     !.md<com>Variables computed as boundary integrals on boundary sets. The results are stored
     !.md<com>in file <tt>sample-boundary.nsi.set</tt>. Let S be the size of the boundary set. 
     !.md<com>The following variables can be postprocessed:
     !.md<com>
     !.md<com>    - MEAN_PRESSURE: int_S P/|S| 
     !.md<com>    - MASS:          int_S rho*u.n
     !.md<com>    - FORCE:         int_S 2 mu eps(u).n, int_S (-pI).n
     !.md<com>    - TORQUE:        int_S (r-rc) x ( 2 mu eps(u).n ), int_S  (r-rc) x (-pI.n) 
     !.md<com>    - MEAN_YPLUS:    int_S y+ / S
     !.md<com>    - MEAN_VELOCITY: int_S u.n / S
     !.md<com>    - WEAT_FORCE:    int_Sw 2 mu eps(u).n, int_Sw (-pI).n 
     !.md<com>    - WEAT_SURFACE:  int_Sw, Sw=wet surface
     !.md<com>
     !
     !.md<1>NODE_SET
     !.md<1>  XXXXX                                                                                         $ Variables computed on nodes
     !.md<1>END_NODE_SET
     !.md<field>NODE_SET
     !.md<com>Variables computed as boundary integrals on boundary sets. The results are stored
     !.md<com>in file <tt>sample-node.nsi.set</tt>. 
     !.md<com>The following variables can be postprocessed:
     !.md<com>
     !.md<com>    - VELOX: X-velocity
     !.md<com>    - VELOY: Y-velocity
     !.md<com>    - VELOZ: Z-velocity
     !.md<com>    - PRESS: Pressure
     !.md<com>    - YPLUS: y+
     !.md<com>
     !.md<com> 
     !.md<1>WITNESS_POINTS
     !.md<1>  XXXXX                                                                                         $ Variables computed witness points
     !.md<1>END_WITNESS_POINTS
     !.md<field>WITNESS_POINTS
     !.md<com>
     !.md<com>    - VELOX: X-velocity
     !.md<com>    - VELOY: Y-velocity
     !.md<com>    - VELOZ: Z-velocity
     !.md<com>    - PRESS: Pressure
     !.md<com>
     !.md<com>
     !
     call ecoute('nsi_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('nsi_reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('nsi_reaous')

        call output_postprocess_read()

        if(words(1)=='INERT') then
           !
           ! Post process in inertial frame of reference
           !
           kfl_inert_nsi=1

        else if(words(1)=='OUTPU') then

           if(words(2)=='ERROR') then
              !.md<1>OUTPUT ERROR, SOLUTION = int1, DIRICHLET: ON | OFF                                     $ Manufactured solution
              !.md<field>OUTPUT ERROR
              !.md<com>Manufactured solutions used to carry out code typing checking and mesh convergence.
              !.md<com>The L2, L1 and Linf norm of the velocity, pressure, velocity gradients and pressure 
              !.md<com>gradients can be found in file <tt>sample.nsi.log</tt>. Dirhclet boundary
              !.md<com>conditions are automatically imposed on the velocity on all the boundary. Pressure
              !.md<com>should be imposed on one node to the right value as the flow is confined.
              !.md<com>The following manufactured solutions are available:
              !.md<com>
              !.md<com>    - SOLUTION = 5:  r=sqrt(x^2+y^2); u = 2*y*(r-0.5), v = -2*x*(r-0.5),  p = r-1 
              !.md<com>    - SOLUTION = 13: u=2x +  y + 3z, v=- x - 3y - 4z, w=6x + 7y - z, p=-5x + 6y + 2z 
              !.md<com>    - SOLUTION = 20: u=-cos(pi*x)*sin(pi*y), v=sin(pi*x)*cos(pi*y), p=sin(pi*x)*sin(pi*x) + sin(pi*y)*sin(pi*y)
              !.md<com>
              !.md<com>If DIRICHLET option is ON, Alya imposes automatically exact Dirichlet conditions on the boundary.
              !
              kfl_exacs_nsi = getint('SOLUT',1_ip,'#Exact solution')
              expar_nsi     = param(4:13)
              if( exists('DIRIC') ) then
                 if( trim(getcha('DIRIC','ON   ','#Automatic Dirichlet for exact solution')) == 'OFF  ') then
                    kfl_exfix_nsi = 0
                 end if
              else if( exists('ONLYO') ) then
                 kfl_exfix_nsi = 2
              end if

           else if(words(2)=='MATRI') then
              !
              ! Matrix profile
              !
              kfl_psmat_nsi=getint('ITERA',1_ip,'#Iteration to postprocess matrix') 

           else if( words(2)=='MINIA' ) then
              !
              ! Miniapp database
              !
              !.md<1>OUTPUT MINIAPP                                                                         $ Output for nastin-miniapp
              !.md<field>OUTPUT MINIAPP
              !.md<com>Output files to be further used by nastin-miniapp. One file per slave is generated             
              kfl_miniapp_nsi = 1

           end if

        else if(words(1)=='PARAM') then
           !
           ! Comfort factor for postprocess
           !
           cloth_nsi = getrea('CLO  ',-500.0_rp,'#CLO') 
           metab_nsi = getrea('MET  ',-500.0_rp,'#MET')  
           wetme_nsi = getrea('WME  ',-500.0_rp,'#WME')  
           ambie_nsi = getrea('TA   ',-500.0_rp,'#TA')  
           radia_nsi = getrea('TR   ',-500.0_rp,'#TR')  
           relat_nsi = getrea('RH   ',-500.0_rp,'#RH')  

        else if(words(1)=='AVERA') then
           !
           ! Averaging starting time
           !
           if(words(2)=='STEPS') then
              avste_nsi = getint('STEPS',0_ip,'#Averaging steps')
           else
              avtim_nsi = getrea('AVERA',0.0_rp,'#Averaging starting time')  
           endif

        end if
     end do
     !
     ! Force, Torque and wet force 
     !
     if( postp(1) % npp_setsb(3) == 1 ) then
        postp(1) % npp_setsb(3:8) = 1
     end if
     if( postp(1) % npp_setsb(9) == 1 ) then       
        postp(1) % npp_setsb(9:14) = 1 
        do ii = 10,14
           postp(1) % pabse(1:3,ii) = postp(1) % pabse(1:3,9)   ! Copy Torque center 
        end do
     end if
     if( postp(1) % npp_setsb(17) == 1 ) then
        postp(1) % npp_setsb(17:23) = 1
     end if
     if( postp(1) % npp_setsb(24) == 1 ) then
        postp(1) % npp_setsb(24:26) = 1
     end if
     if( postp(1) % npp_setsb(36) == 1 ) then
        postp(1) % npp_setsb(36:38) = 1
     end if
     !
     ! Reattachment: force the force to be computed
     !
     if( postp(1) % npp_setsb(27) == 1 ) then
        postp(1) % npp_setsb(27:30) = 1
        postp(1) % npp_setsb( 3)    = 1 
     end if
     !
     ! Recompute sets to be postprocessed
     !
     call output_postprocess_sets_ordering()
     !
     !-->        </group>
     !
     !.md<0>END_OUTPUT_&_POST_PROCESS
     !.md</code>
     ! 
  end if

end subroutine nsi_reaous
    
 
