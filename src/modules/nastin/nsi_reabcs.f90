!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @ingroup    Nastin
!> @{
!> @file    nsi_reabcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Read boundary conditions 
!> @details Read boundary conditions, initial conditions and parameters
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_reabcs()
!-----------------------------------------------------------------------
 !.md<module>nastin
 !.md<input>case.nsi.dat
 !.md<pos>3
 !.md<sec>
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_nastin 
  use mod_memchk
  use mod_opebcs
  use def_kermod
  use mod_ker_space_time_function
  use mod_vortex, only: VORTEX_NAME, vortex_reabcs  
  use mod_ecoute, only : ecoute
  implicit none
  integer(ip)  :: ifunc
  character(5) :: wfname
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opebcs_initialization_structure(3_ip,tncod_nsi)
     call opebcs_initialization_variable (ndime,tncod_nsi(1)) ! Velocity
     call opebcs_initialization_variable ( 1_ip,tncod_nsi(2)) ! Pressure
     call opebcs_initialization_variable ( 1_ip,tncod_nsi(3)) ! Divergence correction
  end if
  if( kfl_icodb > 0 ) then
     call opebcs_initialization_structure(1_ip,tbcod_nsi)     ! Velocity 
     call opebcs_initialization_variable (5_ip,tbcod_nsi)
   end if
  if( kfl_geome > 0 ) then
     call opebcs_initialization_structure( 1_ip,tgcod_nsi)    ! Velocity, geometrical
     call opebcs_initialization_variable (ndime,tgcod_nsi)    ! Velocity
  end if
  if( num_spare_meshes > 0 ) then
     call opebcs_initialization_structure(1_ip,tscod_nsi)
     call opebcs_initialization_variable (ndime,tscod_nsi(1)) ! Velocity
  end if
 
  if( INOTSLAVE ) then
     !
     ! Initialization global variables
     !
     kfl_confi_nsi           = -1                                   ! Flow is not confined
     kfl_local_nsi           = 0                                    ! No local system of reference
     kfl_conbc_nsi           = 1                                    ! Constant boundary conditions
     kfl_syntu_nsi           = 0                                    ! Synthetic eddy method (SEM) (OFF:      = 0; ON:      =1)
     kfl_initi_nsi           = 0                                    ! Initial condition
     kfl_inico_nsi           = 0                                    ! Coarse grid system not solved
     kfl_inipr_nsi           = 0                                    ! Initial pressure
     kfl_nopen_nsi           = 0                                    ! No penetration condition in strong form
     kfl_cadan_nsi           = 0                                    ! No coupling with ADAN for the pressure condition
     neddy_nsi               = 0                                    ! Number of eddies in the inlet box for SEM
     itebc_nsi               = 1                                    ! Initial step boundary condition
     nodpr_global_nsi        = 0                                    ! Node where to impose pressure (global numbering)
     exfpr_nsi               = 0                                    ! Extension of o fixpr one layer of elements
     kfl_imppr_nsi           = 0                                    ! Imposses pressure in nodes w/ fixpr>0
     kfl_divcorrec_nsi       = 0                                    ! No div(u) correction
     delta_nsi               = 0.0_rp                               ! Distance to the wall
     relbc_nsi               = 1.0_rp                               ! Boundary condition relaxation
     valpr_nsi               = 0.0_rp                               ! Pressure value
     mulpr_nsi               = 100.0_rp                             ! Multiplicative factor for pressure imposition
     velin_nsi               = 0.0_rp                               ! Initial velocity
     poise_nsi               = 0.0_rp                               ! Parameters for Poiseuille distribution
     divcorrec_nsi           = 0.0_rp                               ! div(u)=0 parameter (alpha)
     if( mflow_nsi > 0 ) then
        kfl_flow_rate_codes_nsi = 0                                 ! Code on which flow rates are imposed 
        kfl_flow_rate_set_nsi   = 0                                 ! Code of the set on which pressure has to be computed
        kfl_flow_rate_normal_nsi= 0                                 ! Code on which flow rates are imposed 
        kfl_flow_rate_stfn_nsi  = 0                                 ! Flow rate with space time functions
        kfl_flow_rate_tfn_nsi   = 0                                 ! Flow rate with space time functions
        kfl_flow_rate_pfn_nsi   = 0                                 ! Flow rate with discrete functions
        kfl_flow_rate_alg_nsi   = 1                                 ! Flow rate algorithm is multiplicative
        flow_rate_values_nsi    = 0.0_rp                            ! Values of flow rates
        flow_rate_normals_nsi   = 0.0_rp                            ! Values of flow rates
        flow_rate_relax_nsi     = 1.0_rp                            ! Flow rate time relaxation
        flow_rate_press_nsi     = 0.0_rp                            ! Flow rate pressure
     end if
     nflow_rates             = 0                                    ! Number of flow rates
     kfl_press_lev_nsi       = 0                                    ! Level pressure set
     val_press_lev_nsi       = 0.0_rp                               ! Level pressure value
     kfl_immersed_nsi        = 0                                    ! IB strategy
     fact_immersed_nsi       = 100.0_rp                             ! IB factor
     inflow_cosine_nsi       = 1.0e-03_rp                           ! Inflow cosine
     kfl_velom_nsi           = 1                                    ! Velom used in bc
     !
     ! Coupling with ALEFOR: bc must be non-consatnt
     !
     if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
        kfl_conbc_nsi = 0
     end if
     !
     ! Reach the nodal-wise section
     !
     call ecoute('nsi_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('nsi_reabcs')
     end do
     !--><group>
     !-->       <groupName>BOUNDARY_CONDITIONS</groupName>
     !-->       <groupType>bouco</groupType>
     !-->       <groupTypeDefinition><![CDATA[
     !-->         PARAMETERS
     !-->           INITIAL_CONDITION: OFF | COARSE | SPACE_TIME_FUNCTON=char1 | CONSTANT, VALUES=real1,real2,real3 [COARSE]  $ Type of initial condition
     !-->           FIX_PRESSURE:      OFF | AUTOMATIC | ON[, ON_NODE= int1][, VALUE= real1]       $ Fix pressure for confined flow
     !-->           VARIATION:         CONSTANT | NON_CONSTANT                                     $ Boundary conditions constant/variable in time
     !-->         END_PARAMETERS
     !-->         CODES, NODES, PRESSURE
     !-->         ...
     !-->         int1 1                                                                       $  int1=code
     !-->         ...
     !-->        END_CODES
     !-->        CODES, NODES
     !-->          ...
     !-->          int1 int2 real1 real2 real3 AXES:LOCAL | GLOBAL, SPACE_TIME_FUNCTION:char1, TIME_FUNCTION:char2, VALUE_FUNCTION=int3 $  int1=code, int2=velocity code,  char1=space/time function, char2=time function
     !-->          -int1 int2 ...                                                                                                        $ -int1=code, int2=value function, int3=time function, int4=basis 
     !-->          ...
     !-->       END_CODES
     !-->       CODES, BOUNDARIES
     !-->         ...
     !-->         int1 int2 real1                                                              $  int1=code, int2=type of condition, real1=value
     !-->         ...
     !-->       END_CODES]]></groupTypeDefinition>
     !-->       <groupHelp><![CDATA[BOUNDARY_CONDITIONS:
     !-->        In this section, the boundary conditions are defined.
     !-->        Previously, in the <tt>sample.dom.dat</tt> file, codes have been applied to
     !-->        boundaries and/or nodes, or extrapolated from boundaries to nodes.
     !-->        Nastin must translate these conditions in terms of its own degrees of freedom: velocity 
     !-->        and pressure.]]></groupHelp>
     !-->       <parameters><inputLineHelp><![CDATA[PARAMETERS:
     !-->         In this field, some parameters are defined for the Navier-Stokes equations
     !-->         concerning the initial and boundary conditions.]]></inputLineHelp></parameters>
     !-->       <initialCondition><inputLineHelp><![CDATA[<b>INITIAL_CONDITION</b>: 
     !-->         Special initial conditions of the Navier-Stokes equations. The initial condition 
     !-->         defined here is automatically overwritten by Dirichlet conditions.
     !-->         
     !-->          <li>
     !-->           OFF: initial conditions are prescribed in the field "CODES, NODES" with velocity fixity 0.
     !-->          
     !-->          <li>
     !-->            SPACE_TIME_FUNCTION: Prescribe the char function on all nodes.
     !-->          
     !-->          <li>
     !-->            CONSTANT, VALUES = real1, real2, real3: initial velocity is (real1,real2,real3)
     !-->            on nodes where no Dirichlet boundary condition has been prescribed. This option
     !-->            overwrites the initial condition prescribed in the field "CODES, NODES" with velocity fixity 0.
     !-->            If option COARSE is present, the procedure is the following: 
     !-->            an initial velocity (real1,real2,real3) is considered; 
     !-->            then the monolithic stationary Navier-Stokes system is solved on the coarse
     !-->            mesh, using the initial velocity as advection field; finally the solution is prolongated on the original mesh.
     !-->            The coarse problem is solved agglomerating the fine mesh matrix and RHS on the 
     !-->            groups using defined in the <tt>sample.dom.dat</tt> file.
     !-->          
     !-->          <li>
     !-->            COARSE: the procedure is the following:
     !-->            the monolithic stationary Navier-Stokes system is solved on the coarse
     !-->            mesh; finally the solution is prolongated on the original mesh.
     !-->            The coarse problem is solved agglomerating the fine mesh matrix and RHS on the 
     !-->            groups using defined in the <tt>sample.dom.dat</tt> file.
     !-->            This option uses the initial condition prescribed in the file CODES, NODES with code 0.
     !-->          
     !-->        ]]></inputLineHelp></initialCondition>
     !-->       <fixPressure><inputLineHelp><![CDATA[FIX_PRESSURE:
     !-->         By default the pressure is not fixed anywhere. However, if the flow is confined
     !-->         that is if the velocity is prescribed with a Dirichlet boundary condition on 
     !-->         all the boundary nodes, the pressure should be prescribed on one node int1 to a value real1.
     !-->         This is necessary because the pressure is defined up to a constant,
     !-->         The node on which the pressure is prescribed should be located on the boundary.
     !-->         Example:  FIX_PRESSURE: ON, ON_NODE= 1251, VALUE= 0.0. The pressure is prescribed
     !-->         to zero on node 1251. Note that if the mesh multiplication is used, the node number
     !-->         refers to the original mesh. Option AUTOMATIC enables to let Alya choose the node
     !-->         where to prescribe the pressure to a zero value.]]></inputLineHelp></fixPressure>
     !-->       <variation><inputLineHelp><![CDATA[VARIATION:
     !-->             If a time function is used, NON_CONSTANT option should be selected.
     !-->              Alya does not identify non-constant boundary conditions automatically.]]></inputLineHelp></variation>
     !-->       <codesNodesPressure><inputLineHelp><![CDATA[<b>CODES, NODES, PRESSURE</b>:
     !-->                  Interpret the code on nodes to impose the pressure in the Schur complement matrix.
     !-->                  Nastin imposes automatically the pressure in the Schur complement preconditioner on outflows. However, one
     !-->                  can add some nodes to this aoutmatic list by prescribing the node code int2. ]]></inputLineHelp></codesNodesPressure>
     !-->       <codesNodes><inputLineHelp><![CDATA[<b>CODES, NODES</b>:
     !-->                   Interpret the code on nodes to impose velocity degrees of freedom.
     !-->                   
     !-->                     <li> 
     !-->                        int1 is the code node, int=1,2... or 99.
     !-->                        If the node has multiple codes (e.g. if the code was extrapolated from boundary codes), the
     !-->                        Syntax is int1 = 1 & 2 & 4. It means that nodes with the three codes 1,2 and 4 are considered.
     !-->                      
     !-->                     <li> 
     !-->                       int2 has the form f1f2f3: 11, 10, 01 in 2D and 101, 111, etc in 3D. 0 means free or initial condition, 1 means prescribed.
     !-->                       Values f1, f2 and f3 are fixity codes of the velocity degrees of freedom.
     !-->                       f1 refers to the first velocity components, f2 to the second and f3 to the third.
     !-->                      
     !-->                     <li> 
     !-->                       real1, real2 and real3 are the corresponding values for each degree of freedom (first to third velocity component).
     !-->                      
     !-->                     <li> 
     !-->                       char1 is the time function to applied to the nodes with code int1. If not time function is given, the conditions will be constant in time.
     !-->                       The function should be defined in sample.ker.dat file.
     !-->                      
     !-->                     <li> 
     !-->                       char2 is the space/time function f(x,t) to applied to the nodes with code int1. The value of the nodal condition will be
     !-->                       f(x,t)*(real1,real2,real3) is the function is a scalar or (f1(x,t)*real1,f2(x,t)*real2,f3(x,t)*real3) if the funciton
     !-->                       is vectorial.
     !-->                      
     !-->                     <li> 
     !-->                       AXES is the basis in which the velocity is prescibed. GLOBAL means that the three components of the velocity
     !-->                       are the three Cartesian components. LOCAL mean the local basis computed by Alya. The first component is the
     !-->                       normal component and the other two the tangential ones.
     !-->                      
     !-->                     <li> 
     !-->                       If the option VALUE_FUNCTION is present, then the nodal condition will be taken from the corresponding field defined
     !-->                       in the BOUNDARY_CONDITIONS field of the file sample.dom.dat   
     !-->                      
     !-->                    
     !-->                    Examples:
     !-->                    
     !-->                      <li> 
     !-->                         3 100 0.0 0.0 0.0 0 AXES:LOCAL => Slip condition: zero normal velocity is prescribed and tangent velocity is free on nodes with code 3.
     !-->                         To impose a wall law, a condition on boundaries of type 3 should be prescribed.
     !-->                       
     !-->                      <li> 
     !-->                        5 111 0.0 0.0 0.0 0 FIELD=2 => Prescribe value function number 2 as a Dirichlet boundary condition on all nodes.
     !-->                       
     !-->                      <li> 
     !-->                        4 & 6 111 0.0 0.0 0.0 => No-slip conditions on nodes with codes 4 and 6.
     !-->                       
     !-->                      <li> 
     !-->                        7 000 1.0 2.0 3.0 => Initial velocity is set to (1,2,3) on nodes with code 7.
     !-->                       
     !--> ]]></inputLineHelp></codesNodes>
     !--> <codesBoundaries><inputLineHelp><![CDATA[CODES, BOUNDARIES:
     !-->                  Impose a natural boundary condition of type int2 on boundaries with code int1.
     !-->                  The different boundary conditions available are:
     !-->                  
     !-->                   <li> 
     !-->                      int2 = 3: law of the wall. It requires the distance to the wall as given in Kermod. real1 not needed.
     !-->                    
     !-->                   <li> 
     !-->                     int2 = 2: impose the traction to value real1. If the flow is uniform or the Reynolds number very high, this is 
     !-->                     equivalent, this is equivalent to imposing the pressure.
     !-->                    
     !-->                  
     !-->                  Examples:
     !-->                  
     !-->                   <li> 
     !-->                     5 3: Impose the law of the wall to boundaries with code 5.
     !-->                    
     !-->                   <li> 
     !-->                     6 2 5.0: Impose the pressure to 5.0 on boundaries with code 6.
     !-->                    
     !-->                   ]]></inputLineHelp></codesBoundaries>
     !--></group>
     !
     !.md# Boundary Conditions
     !.mdIn this section, the boundary conditions are defined.
     !.mdPreviously, in the <tt>sample.dom.dat</tt> file, codes have been applied to
     !.mdboundaries and/or nodes, or extrapolated from boundaries to nodes.
     !.mdNastin must translate these conditions in terms of its own degrees of freedom: velocity 
     !.mdand pressure.
     !.md<code>
     !.md<0>BOUNDARY_CONDITIONS
     !
     call ecoute('nsi_reabcs')
     do while(words(1)/='ENDBO')

        if (words(1) == 'PARAM' ) then

           !-------------------------------------------------------------
           !
           ! Parameters
           !
           !-------------------------------------------------------------

           call ecoute('nsi_reabcs')
           do while(words(1)/='ENDPA') 
              !
              !.md<1>PARAMETERS
              !.md<field>PARAMETERS
              !.md<com>In this field, some parameters are defined for the Navier-Stokes equations
              !.md<com>concerning the initial and boundary conditions.
              !
              !if( words(1) == 'HYDRO' ) then
              !   !
              !   ! Initial hydrostatic pressure: KFL_HYDRO_NSI, HYDRO_NSI
              !   !
              !   if( words(2) == 'ON   ') kfl_hydro_nsi = -1
              !   if(exists('Z    ')) then
              !      kfl_hydro_nsi = 1
              !      hydro_nsi     = getrea('Z    ',0.0_rp,'#Hydrostatic z-plane')
              !   else if(exists('Y    ')) then
              !      kfl_hydro_nsi = 1
              !      hydro_nsi     = getrea('Y    ',0.0_rp,'#Hydrostatic y-plane')
              !   end if

              if(      words(1) == 'INFLO' ) then
                 !
                 ! Inflow cosine
                 !
                 !
                 !.md<2>INFLOW_COSINE= val1
                 !.md<field>INFLOW_COSINE
                 !.md<com>For adaptative inflow/outflow condition. Inflow is imposed
                 !.md<com>if if( udotn <= inflow_cosine_nsi ) t
                 !
                 inflow_cosine_nsi = getrea('INFLO',1.0e-3_rp,'#Value of inflow cosine')
                 
              else if( words(1) == 'PRESS' ) then
                 !
                 ! Initial hydrostatic pressure: KFL_HYDRO_NSI, HYDRO_NSI
                 !
                 !.md<2>PRESSURE= LAPLACIAN/HYDROSTATIC
                 !.md<field>PRESSURE
                 !.md<com>When prescribing Neumann boundary condition on pressure, use
                 !.md<com>LAPLACIAN to solve for Dp=0 with Dirichlet conditions corresponding
                 !.md<com>to the prescribed pressure on boundaries.
                 !
                 if( words(2) == 'LAPLA' ) then 
                    kfl_inipr_nsi = 1
                 else if( words(2) == 'HYDRO' ) then
                    kfl_inipr_nsi = 2
                 end if

              else if( words(1) == 'IMMER' ) then
                 !
                 ! Immersed boudnary
                 !
                 !.md<2>IMMERSED_BOUNDARY: PENALTY/WEAK, VALUE=char1
                 !.md<field>IMMERSED_BOUNDARY
                 !.md<com>When using weak Dirichlet condition on a spare mesh, decide if the integral
                 !.md<com>should penalize the equation (PENALTY) or subsitute the equation (WEAK)
                 !
                 if( words(2) == 'PENAL' ) then                   
                    kfl_immersed_nsi   = 0                         
                    fact_immersed_nsi  = getrea('VALUE',100.0_rp,'#Value of pressure level')
                 else if( words(2) == 'WEAK ' ) then      
                    kfl_immersed_nsi   = 1
                 end if
                      
              else if( words(1) == 'LEVEL' ) then
                 !
                 ! Level pressure
                 !
                 !.md<2>LEVEL_PRESSURE: SET=char1, VALUE=val2
                 !.md<field>LEVEL_PRESSURE
                 !.md<com>Level pressure. It can be useful when coupling to solidz to change the level of pressure.
                 !.md<com>Set value are not leveld, while output pressure is
                 !
                 if( exists('SET  ') ) kfl_press_lev_nsi = getint('SET  ',0_ip  ,'#Set where to level pressure')
                 if( exists('VALUE') ) val_press_lev_nsi = getrea('VALUE',0.0_rp,'#Value of pressure level')

              else if(words(1) == 'NOPEN' ) then
                 !
                 ! No penetration condition: weak or strong form
                 !
                 if( words(2) == 'WEAK ' .or. words(2) == 'NATUR' .or. words(2) == 'WEAKF' ) then
                    kfl_nopen_nsi = 1
                 else if( words(2) == 'STRON' .or. words(2) == 'DIRIC' ) then
                    kfl_nopen_nsi = 0
                 end if

              else if( words(1) == 'VELOM' ) then
                 !
                 !.md<2>VELOM: ON/OFF
                 !.md<field>VELOM
                 !.md<com>Velom is added to prescribed velocity boundary condition.
                 !
                 if( option('VELOM') ) then
                    kfl_velom_nsi = 1
                 else
                    kfl_velom_nsi = 0
                 end if
                 
              else if( words(1) == 'INITI' ) then                 
                 !
                 ! Initial solution KFL_INITI: Constant, Stokes, etc.
                 !
                 !.md<2>INITIAL_CONDITION: OFF | COARSE | SPACE_TIME_FUNCTON=char1 | CONSTANT, VALUES=real1,real2,real3 [COARSE]  $ Type of initial condition
                 !.md<field>INITIAL_CONDITION
                 !.md<com>Special initial conditions of the Navier-Stokes equations. The initial condition 
                 !.md<com>defined here is automatically overwritten by Dirichlet conditions.
                 !.md<com>
                 !.md<com>
                 !.md<com>    - OFF: initial conditions are prescribed in the field "CODES, NODES" with velocity fixity 0.
                 !.md<com>    - SPACE_TIME_FUNCTION: Prescribe the char function on all nodes.
                 !.md<com>    - CONSTANT, VALUES = real1, real2, real3: initial velocity is (real1,real2,real3)
                 !.md<com>      on nodes where no Dirichlet boundary condition has been prescribed. This option
                 !.md<com>      overwrites the initial condition prescribed in the field "CODES, NODES" with velocity fixity 0.
                 !.md<com>      If option COARSE is present, the procedure is the following: 
                 !.md<com>      an initial velocity (real1,real2,real3) is considered; 
                 !.md<com>      then the monolithic stationary Navier-Stokes system is solved on the coarse
                 !.md<com>      mesh, using the initial velocity as advection field; finally the solution is prolongated on the original mesh.
                 !.md<com>      The coarse problem is solved agglomerating the fine mesh matrix and RHS on the 
                 !.md<com>      groups using defined in the <tt>sample.dom.dat</tt> file.
                 !.md<com>    - COARSE: the procedure is the following:
                 !.md<com>      the monolithic stationary Navier-Stokes system is solved on the coarse
                 !.md<com>      mesh; finally the solution is prolongated on the original mesh.
                 !.md<com>      The coarse problem is solved agglomerating the fine mesh matrix and RHS on the 
                 !.md<com>      groups using defined in the <tt>sample.dom.dat</tt> file.
                 !.md<com>      This option uses the initial condition prescribed in the file CODES, NODES with code 0.
                 !.md<com>
                 !.md<com>
                 !
                 if( exists('NONIN') .or. exists('CONST') ) then                           ! 1: Constant
                    if (kfl_initi_nsi /= 6 .and. kfl_initi_nsi /= 7) kfl_initi_nsi = 1
                    velin_nsi(1:3) = param(3:5)    
                 else if( exists('INERT') ) then                                           ! 2: Constant in inertial f.o.r     
                    velin_nsi(1:3) = param(3:5)
                    kfl_initi_nsi = 2
                 else if( exists('STOKE') ) then                                           ! 3: Stokes
                    kfl_initi_nsi = 3
                 else if( exists('POTEN') ) then                                           ! 4: Potential flow
                    kfl_initi_nsi = 4
                 else if( exists('LAPLA') ) then                                           ! 5: Laplacian
                    kfl_initi_nsi = 5
                    solve(3)%kfl_solve = 1                                                 ! Output flag
                 else if( exists('POISE') ) then                                           ! 6: Poiseuille distribution 
                    if (exists('BOUND')) then
                       kfl_initi_nsi = 7
                       poise_nsi(1:6) = param(4:9)                                         ! Parameters: Flow axis, radius, max. velocity
                    else
                       kfl_initi_nsi = 6
                       poise_nsi(1:6) = param(3:8)               
                    endif
                 else if( exists('RANDO') ) then                                           ! 8: Random field
                    kfl_initi_nsi = 8
                 else if( exists(VORTEX_NAME) ) then                                       ! 9: Vortex 
                   call vortex_reabcs( param(:), kfl_initi_nsi )  

                 else if( exists('SPACE') ) then                                           ! >100: Space time function
                    wfname = getcha('SPACE','NULL ','#Space/time Function name')
                    kfl_initi_nsi = 100 + space_time_function_number(wfname)
                 else if( exists('VALUE') ) then                                           ! <0: value function
                    if( exists('PRESS') )then
                       kfl_inipr_nsi = &
                            -getint('VALUE',1_ip,'#Initial condition for pressure is from value function')
                    else   
                       kfl_initi_nsi = &
                            -getint('VALUE',1_ip,'#Initial condition is from value function')
                    end if
                 else if( exists('FIELD') ) then                                           ! <0: value function
                    if( exists('PRESS') )then
                       kfl_inipr_nsi = &
                            -getint('FIELD',1_ip,'#Initial condition for pressure is from value function')
                    else   
                       kfl_initi_nsi = &
                            -getint('FIELD',1_ip,'#Initial condition is from value function')
                    end if
                 end if
                 if( exists('COARS') ) kfl_inico_nsi = 1                                   ! Coarse grid system

              else if( words(1) == 'FIXPR' ) then
                 !
                 ! Pressure fixity
                 !
                 !.md<2>FIX_PRESSURE:      OFF | AVERAGE | AUTOMATIC | ON[, ON_NODE= int1][, VALUE= real1]       $ Fix pressure for confined flow
                 !.md<field>FIX_PRESSURE
                 !.md<com>By default the pressure is not fixed anywhere. However, if the flow is confined
                 !.md<com>that is if the velocity is prescribed with a Dirichlet boundary condition on 
                 !.md<com>all the boundary nodes, the pressure should be prescribed on one node int1 to a value real1.
                 !.md<com>This is necessary because the pressure is defined up to a constant,
                 !.md<com>The node on which the pressure is prescribed should be located on the boundary.
                 !.md<com>Example:  FIX_PRESSURE: ON, ON_NODE= 1251, VALUE= 0.0. The pressure is prescribed
                 !.md<com>to zero on node 1251. Note that if the mesh multiplication is used, the node number
                 !.md<com>refers to the original mesh. Option AUTOMATIC enables to let Alya choose the node
                 !.md<com>where to prescribe the pressure to a zero value.
                 !
                 if( words(2) == 'AUTOM' ) then
                    
                    kfl_confi_nsi = 0
                 else if( words(2) == 'NO   '.or. words(2) == 'OFF  ' ) then
                    
                    kfl_confi_nsi = -1
                 else if( words(2) == 'AVERA' ) then
                    
                    kfl_confi_nsi = -2
                 else if( words(2) == 'COORD' ) then
                    
                    kfl_confi_nsi = -3
                 else if( words(2) == 'YES  '.or. words(2) == 'ON   ' ) then
                    
                    kfl_confi_nsi =  1
                    if( exists('ONNOD') ) then
                       nodpr_global_nsi = getint('ONNOD',1_ip,'#Node where to impose pressure')
                    else if( exists('NODE ') ) then
                       nodpr_global_nsi = getint('NODE ',1_ip,'#Node where to impose pressure')
                    end if
                    valpr_nsi = getrea('VALUE',0.0_rp,'#Pressure value')
                    
                 else if( words(2) == 'WEAK ' ) then
                    kfl_confi_nsi =  2
                    if( exists('ONNOD') ) then
                       nodpr_global_nsi = getint('ONNOD',1_ip,'#Node where to impose pressure')
                    else if( exists('NODE ') ) then
                       nodpr_global_nsi = getint('NODE ',1_ip,'#Node where to impose pressure')
                    end if
                    valpr_nsi = getrea('VALUE',0.0_rp,'#Pressure value')
                    mulpr_nsi = getrea('MULTI',0.0_rp,'#Pressure multiplicative factor')
                 end if
                 
              else if( words(1) == 'VARIA' ) then
                 !
                 !.md<2>VARIATION:         CONSTANT | NON_CONSTANT                                     $ Boundary conditions constant/variable in time
                 !.md<field>VARIATION
                 !.md<com>If a time function is used, NON_CONSTANT option should be selected.
                 !.md<com>Alya does not identify non-constant boundary conditions automatically.
                 !
                 ! if( words(2) == 'NONCO' .and. kfl_conbc_nsi /= 0 ) then   
                 if( words(2) == 'NONCO' ) then
                    kfl_conbc_nsi = 0 
                 else
                    kfl_conbc_nsi = 1
                 end if

              else if( words(1) == 'SYNTH' ) then
                 !
                 !.md<2>SYNTHETIC:             $ Syhntetic eddy method (SEM) 
                 !.md<field>SYNTHETIC
                 !.md<com>Synthetic eddy method (SEM) based on the work by Jarrin et al. (2006).
                 !.md<com>A Synthetic-Eddy Method for Generating Inflow Conditions 
                 !.md<com>for LES, Jarrin, N., Benhamadouche, S., Laurence, D. and Prosser, R 
                 !.md<com>International Journal of Heat and Fluid Flow, Vol. 27, pp. 585-593, (2006). 
                 !
                 if( exists('NUMBE') ) &
                       neddy_nsi = getint('NUMBE',0_ip,'#Node where to impose pressure')
                 kfl_syntu_nsi = 1

              else if( words(1) == 'DIVER' ) then
                 !
                 !.md<2>DIVERGENCE:  ON/OFF, VALUES=a1
                 !.md<field>DIVERGECE
                 !.md<com>Correct initial velocity so that div(u)=0. 
                 !.md<com>The parameter alpha (a1) controls the relative adjustments of horizontal
                 !.md<com>and vertical wind components and it is related to the atmospheric thermal
                 !.md<com>stratification. Values of alpha > 1 are representative of unstable atmosphere,
                 !.md<com>values of alpha < 1 reflect a stable atmosphere, and alpha close to 1 indicate
                 !.md<com>a neutral situation.
                 !.md<com>See Homicz, G. F.: Three-Dimensional Wind Field Modeling: A Review, (2002)
                 !
                 if( option('DIVER') ) then
                    kfl_divcorrec_nsi = 1
                    divcorrec_nsi = param(3)
                 else if( exists('SPACE') ) then                                           ! >100: Space time function
                    wfname = getcha('SPACE','NULL ','#Space/time Function name')
                    kfl_initi_nsi = 100 + space_time_function_number(wfname)
                 end if

              end if

              call ecoute('nsi_reabcs')
              !
              !.md<1>END_PARAMETERS
              !
           end do

        else if( words(1) == 'FLOWR' ) then
           !
           ! Flow rates
           !.md<1>FLOW_RATE= real, CODE= flowr [,SPACE_TIME_FUNCTION = name, NORMALS= x, y, z], MULTIPLICATIVE/CONSTRAINT
           !.md<field>FLOW_RATE
           !.md<com>Impose a specified flow rate on the nodes with code flowr.
           !.md<com>Two algorithms are available.
           !.md<com>    - By constraint: solve the following minimzation problem with Lagrange multiplier
           !.md<com>      ```math
           !.md<com>      \begin{equation}
           !.md<com>      \left\{ \begin{array}{ll}
           !.md<com>      {\rm Minimize} & \int_S | \mathbf{u} - \mathbf{u}* |^2 \,ds  \\
           !.md<com>      {\rm Under}    & \int_S \mathbf{u} \cdot \mathbf{n} \,ds= {\rm flowr}
           !.md<com>      \end{array} \right.
           !.md<com>      \end{equation}
           !.md<com>      ```
           !.md<com>    - Multiplicative: scale velocity by multiplying by a constant such that
           !.md<com>      ```math
           !.md<com>      \int_S \mathbf{u} \cdot \mathbf{n} \,ds = {\rm flowr}
           !.md<com>      ```
           !
           nflow_rates = nflow_rates + 1
           if( nflow_rates > size(kfl_flow_rate_codes_nsi,KIND=ip) ) &
                call runend('NSI_REABCS: INCREASE FLOW RATE SIZE IN DEF_NASTIN')
           if( nflow_rates > mflow_nsi ) &
                call runend('NSI_REABCS: TOO MANY FLOW RATE CONDITIONS. INCREASE mflow_nsi IN DEF_NASTIN')
           kfl_flow_rate_codes_nsi(nflow_rates) = getint('CODE ',1_ip,  '#Imposed flow rate code')  
           flow_rate_values_nsi(nflow_rates)    = getrea('FLOWR',0.0_rp,'#Imposed flow rate value')                 
           kfl_flow_rate_stfn_nsi(nflow_rates)  = 0_ip
           kfl_flow_rate_tfn_nsi (nflow_rates)  = 0_ip
           kfl_flow_rate_alg_nsi (nflow_rates)  = 0_ip
           if (exists('SETBO ')) then
              kfl_flow_rate_set_nsi(nflow_rates) = getint('SETBO',1_ip,  '#Computed pressure in the set code')
           end if
           if ( exists('PRESS')) then
            flow_rate_press_nsi(nflow_rates) = getrea('PRESS',0.0_rp,  '#Imposed pressure in the boundary')
           end if
           if (exists('NORMA')) then
              call runend('NSI_REABCS: FLOW RATE NORMAL NOT IMPLEMENTED YET, TALK TO YOU RESPONSIBLE')
              kfl_flow_rate_normal_nsi(nflow_rates) = kfl_flow_rate_codes_nsi(nflow_rates) 
              flow_rate_normals_nsi(nflow_rates,1:ndime) = param(2:ndime+1)            ! Parameters: normal x, y [, z]
           end if
           if( exists('SPACE') ) then                                           ! >100: Space time function
              wfname = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_flow_rate_stfn_nsi(nflow_rates) = space_time_function_number(wfname)
           end if
           if( exists('RELAX') ) then
              flow_rate_relax_nsi(nflow_rates) = getrea('RELAX',1.0_rp,  '#Temporal relaxation parameter')
           end if
           !
           ! Method
           !
           if(      exists('MULTI') ) then
              kfl_flow_rate_alg_nsi(nflow_rates) = 1
           else if( exists('CONST') ) then
              kfl_flow_rate_alg_nsi(nflow_rates) = 0
           end if
           if( exists('METHO') ) then
              if(      getcha('METHO','     ','#Flow rate method') == 'MULTI' ) then
                 kfl_flow_rate_alg_nsi(nflow_rates) = 1
              else if( getcha('METHO','     ','#Flow rate method') == 'CONST' ) then
                 kfl_flow_rate_alg_nsi(nflow_rates) = 0                 
              end if
           end if
           !
           ! KFL_FUNNO_TMP > 0: Time function
           ! KFL_FUNNO_TMP < 0: Space/time function
           !
           ! Time function
           !
           wfname = '     '
           if( exists('TIMEF') ) then
              wfname  = getcha('TIMEF','     ','#Time Function name')
              do ifunc = 1,number_time_function
                 if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                    kfl_flow_rate_tfn_nsi(nflow_rates) = ifunc
                 end if
              end do
           end if
           !
           ! Pump function
           !
           wfname = '     '!
           if( exists('PUMP ') ) then
              wfname  = getcha('PUMP ','     ','#PUMP Function name')
              do ifunc = 1,number_pump_curve
                 if( trim(wfname) == trim(pump_curve(ifunc) % name) ) then
                    kfl_flow_rate_pfn_nsi(nflow_rates) = ifunc
                 end if
              end do
           end if

        else if( words(1) == 'COUAD') then           
           !
           ! Pressure prescriptions will respond to coupling with ADAN
           !           
           kfl_cadan_nsi=1
           kfl_aiobo_nsi= int(param(1)) !Boundary connected with ADAN

        else if( words(1) == 'CODES' .and. exists('NODES') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !-------------------------------------------------------------

           if( exists('PRESS') ) then
              !
              !.md<1>CODES, NODES, PRESSURE
              !.md<2>  ...
              !.md<2>  int1 1                                                                       $  int1=code
              !.md<2>  ...
              !.md<1>END_CODES
              !.md<field>CODES, NODES, PRESSURE
              !.md<com>Interpret the code on nodes to impose the pressure in the Schur complement matrix.
              !.md<com>Nastin imposes automatically the pressure in the Schur complement preconditioner on outflows. However, one
              !.md<com>can add some nodes to this aoutmatic list by prescribing the node code int2. 
              ! 
              tncod => tncod_nsi(2:)
              !tncod => momod(ID_NASTIN) % tncod(2:)
              !call reacod(READ_NODE_CODES)
              call boundary_conditions_read_node_codes('PRESSURE')

           else if( exists('DIVER') ) then
              !
              ! Codes for divergence correction
              !
              tncod => tncod_nsi(3:)
              !call reacod(READ_NODE_CODES)
              call boundary_conditions_read_node_codes('DIVERGENCE')

           else if( exists('GEOME') ) then
              !
              ! Velocity: geometrical node code
              !              
              tgcod => tgcod_nsi(1:)
              !tgcod => momod(ID_NASTIN) % tgcod(1:)
              call reacod(READ_GEOMETRICAL_CODES)

           else
              !
              !.md<1>CODES, NODES
              !.md<2>  ...

              !.md<2>  ...
              !.md<1>END_CODES
              !.md<field>CODES, NODES
              !.md<com>Interpret the code on nodes to impose velocity degrees of freedom.
              !.md<com>   - int1 is the code node, int=1,2... or 99.
              !.md<com>     If the node has multiple codes (e.g. if the code was extrapolated from boundary codes), the
              !.md<com>     Syntax is int1 = 1 & 2 & 4. It means that nodes with the three codes 1,2 and 4 are considered.
              !.md<com>   - int2 has the form f1f2f3: 11, 10, 01 in 2D and 101, 111, etc in 3D. 0 means free or initial condition, 1 means prescribed.
              !.md<com>     Values f1, f2 and f3 are fixity codes of the velocity degrees of freedom.
              !.md<com>     f1 refers to the first velocity components, f2 to the second and f3 to the third.
              !.md<com>   - real1, real2 and real3 are the corresponding values for each degree of freedom (first to third velocity component).
              !.md<com>   - char1 is the time function to applied to the nodes with code int1. If not time function is given, the conditions will be constant in time.
              !.md<com>     The function should be defined in sample.ker.dat file.
              !.md<com>   - char2 is the space/time function f(x,t) to applied to the nodes with code int1. The value of the nodal condition will be
              !.md<com>     f(x,t)*(real1,real2,real3) is the function is a scalar or (f1(x,t)*real1,f2(x,t)*real2,f3(x,t)*real3) if the funciton
              !.md<com>     is vectorial.
              !.md<com>   - AXES is the basis in which the velocity is prescibed. GLOBAL means that the three components of the velocity
              !.md<com>     are the three Cartesian components. LOCAL mean the local basis computed by Alya. The first component is the
              !.md<com>     normal component and the other two the tangential ones.
              !.md<com>   - If the option VALUE_FUNCTION is present, then the nodal condition will be taken from the corresponding field defined
              !.md<com>     in the BOUNDARY_CONDITIONS field of the file sample.dom.dat            
              !.md<com>     Examples:
              !.md<com>     - 3 100 0.0 0.0 0.0 0 AXES:LOCAL => Slip condition: zero normal velocity is prescribed and tangent velocity is free on nodes with code 3.
              !.md<com> To impose a wall law, a condition on boundaries of type 3 should be prescribed.
              !.md<com>     - 5 111 0.0 0.0 0.0 0 FIELD=2 => Prescribe value function number 2 as a Dirichlet boundary condition on all nodes.
              !.md<com>     - 4 & 6 111 0.0 0.0 0.0 => No-slip conditions on nodes with codes 4 and 6.
              !.md<com>     - 7 000 1.0 2.0 3.0 => Initial velocity is set to (1,2,3) on nodes with code 7.
              !.md<com> 
              !
              tncod => tncod_nsi(1:)
              call boundary_conditions_read_node_codes('VELOCITY')

           end if

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !          
           !-------------------------------------------------------------
           !
           !.md<1>CODES, BOUNDARIES
           !.md<2>  ...
           !.md<2>  int1 int2 real1                                                              $  int1=code, int2=type of condition, real1=value
           !.md<2>  ...
           !.md<1>END_CODES
           !.md<field>CODES, BOUNDARIES
           !.md<com>Impose a natural boundary condition of type int2 on boundaries with code int1.
           !.md<com>The different boundary conditions available are:
           !.md<com>    - int2 = 3: law of the wall. It requires the distance to the wall as given in Kermod. real1 not needed.
           !.md<com>    - int2 = 2: impose the traction to value real1. If the flow is uniform or the Reynolds number very high, this is 
           !.md<com>     equivalent, this is equivalent to imposing the pressure.
           !.md<com>
           !.md<com>       Examples:
           !.md<com>       - 5 3: Impose the law of the wall to boundaries with code 5.
           !.md<com>       - 6 2 5.0: Impose the pressure to 5.0 on boundaries with code 6.
           !.md<com>
           !
           tbcod => tbcod_nsi(1:)
           call boundary_conditions_read_boundary_codes('MOMENTUM')

        else if( words(1) == 'CODES' .and. exists('SPARE') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on spare meshes
           !          
           !-------------------------------------------------------------
           !
           !.md<1>CODES, SPARE
           !.md<2>  ...
           
           !.md<2>  ...
           !.md<1>END_CODES
           !.md<field>CODES, SPARE
           !.md<com>Impose velocity value on a spare mesh
           !
           call boundary_conditions_read_node_codes('VELOCITY',tscod_nsi)
              
        else if( words(1) == 'EXTEN' ) then
           exfpr_nsi = 1

        else if( words(1) == 'IMPOS' ) then
           kfl_imppr_nsi = 1
        end if

        call ecoute('nsi_reabcs')

     end do
     !
     !.md<0>END_BOUNDARY_CONDITIONS
     !.md</code>
     !

  end if

end subroutine nsi_reabcs
