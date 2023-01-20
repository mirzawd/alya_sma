!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup SolidzInput
!> @{
!> @file    sld_reabcs.f90
!> @author  Mariano Vazquez
!> @date    
!> @brief   Read boundary conditions
!> @details Read boundary conditions
!> @}
!-----------------------------------------------------------------------

subroutine sld_reabcs()
  !
  !.md<module>solidz
  !.md<input>case.sld.dat
  !.md<pos>3
  !.md<sec>
  !
  use def_kintyp, only : ip, rp
  use def_inpout, only : words, exists, param, getint, getrea, getcha
  use def_master, only : INOTSLAVE
  use def_domain, only : kfl_icodb, kfl_icodn
  use def_domain, only : ndime, tbcod, tncod
  use mod_opebcs, only : opnbcs, opbbcs
  use mod_opebcs, only : boundary_conditions_read_node_codes, boundary_conditions_read_boundary_codes
  use mod_ecoute, only : ecoute
  use mod_messages, only : messages_live
  use def_solidz, only : tbcod_sld, tncod_sld
  use def_solidz, only : kfl_funty_sld, nfunc_sld, mtloa_sld, tload_sld
  use def_solidz, only : fubcs_sld
  use def_solidz, only : kfl_local_sld, kfl_csysl_sld
  use def_solidz, only : csysl_sld, SLD_CSYS_CARTESIAN, SLD_CSYS_CYLINDRICAL
  use def_solidz, only : SLD_CSYS_SPHERICAL, SLD_CSYS_EXNOR
  use def_solidz, only : fubcs_sld, rtico_sld
  use def_solidz, only : kfl_conta_sld
  use def_solidz, only : contactbou_sld, neumann_relax, coupling_contact_tol, coupling_contact_its
  use def_solidz, only : vect_proje_sld, petol_sld, kfl_mrele_sld, kfl_contf_sld
  use def_solidz, only : kfl_bodyf_sld, kfl_conbc_sld, kfl_newbc_sld
  use def_solidz, only : ncoef_sld
  use def_solidz, only : kfl_follo_sld
  use def_solidz, only : kfl_insdv_sld,kfl_invel_sld,invel_sld
  use def_solidz, only : kfl_bvesv_sld
  use def_solidz, only : ncrak_sld, crkco_sld
  use def_solidz, only : kfl_conta_stent, conbo_sld, r_fin_stent, contact_friction_sld
#if defined COMMDOM && COMMDOM == 2
  use mod_sld_pdn_contact_plepp, only : sld_pdn_contact_reabcs
#endif
  use mod_sld_atm, only : sld_atm_reabcs
  implicit none

  integer(ip)   :: itloa,ifunc,nauxi,icrak,idime,dummi,ipoin,nposi
  real(rp)      :: rauxi(2)
  character(5)  :: wcsys

  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(0_ip,1_ip,ndime,0_ip,tncod_sld)
  end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,ndime,tbcod_sld)
  end if

  if( INOTSLAVE ) then
     !
     ! Initializations global variables
     !
     kfl_conbc_sld  = 0                    ! Constant boundary conditions
     kfl_newbc_sld  = 0                    ! New boundary conditions
     kfl_local_sld  = 0                    ! No local system of reference
     kfl_conta_sld  = 0                    ! PDN contact
     kfl_contf_sld  = 0                    ! PDN contact
     kfl_bodyf_sld  = 0                    ! Number of body forces (by default is OFF)
     kfl_follo_sld  = 0                    ! Follower loads
     kfl_invel_sld  = 0                    ! Initial condition: velocity
     kfl_insdv_sld  = 0                    ! Initial condition: state dependent variables
     kfl_bvesv_sld  = 0                    ! Flag for velocity prescription
     kfl_conta_stent = 0                   ! Flag for stenting: 1 is crimpring, 2 i expansion, 3 is charge
     r_fin_stent = 0                       ! Initialization of final raidus of stent
     nfunc_sld      = 0                    ! Number time-dependent boundary condition types
     mtloa_sld(:)   = 0                    ! Number of discrete times per function
     contactbou_sld = 0                    ! Identificator of contact boundary
     neumann_relax  = 0.9_rp               ! Neumann relaxation
     coupling_contact_tol = 1e-4_rp        ! Contact coupling tolerance
     coupling_contact_its = 20_ip          ! Contact coupling iterations
     vect_proje_sld  = (/ 0.0_rp,1.0_rp,0.0_rp /)  ! Projection direction
     petol_sld       = 0.0_rp              ! Penetration tolerance for PDN-contact
     contact_friction_sld = 0.0_rp         ! Contact friction
     kfl_mrele_sld  = 0_ip
     csysl_sld      = 0.0_rp               ! Parameters local coordinate system for boundary conditions
     invel_sld(:)   = 0.0_rp               ! Initial velocity values
     !
     ! Reach the section 
     !
     call ecoute('sld_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('sld_reabcs')
     end do
     !
     ! Read header
     !
     if( exists('CONST') ) then
        kfl_conbc_sld = 1
     else if( exists('TRANSI') .or. exists('NONCO') ) then
        kfl_conbc_sld = 0
     end if
     !
     ! Allocate the prescription time function vector of vectors
     !
     call sld_membcs(3_ip)
     !
     ! Begin to read data
     !
     call ecoute('sld_reabcs')
     !
     !.md# Boundary Conditions
     !.md<code>
     !.md<0>BOUNDARY_CONDTIONS, TRANSIENT | CONSTANT
     !
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'PARAM' ) then

           !-------------------------------------------------------------
           !
           ! Parameters
           !
           !-------------------------------------------------------------
           !
           !.md<1>PARAMETERS
           !.md<com>
           !.md<com><strong>Parameters</strong>
           !
           call ecoute('sld_reabcs')
           do while( words(1)/='ENDPA' )

              if( words(1) == 'INITI' ) then
                 !
                 !.md<2>INITIAL_CONDITIONS:    VELOCITY | STATE | TEMPERATURE         $ Initial condition
                 !.md<field>INITIAL_CONDITIONS
                 !.md<com>Set an initial condition for the problem. The initial conditions can be specified for particular
                 !.md<com>nodes or elements. The data can be provided directly or by using a field depending on the variable.
                 !.md<com>If the initial conditions are not specified, all of them are set to zero.
                 !.md<com>
                 !
                 if(      words(2) == 'VELOC' ) then
                    !
                    !.md<com>    - Set VELOCITY to prescribe initial velocities at the beginning of the simulation.
                    !.md<com>This option can only be used for DYNAMIC (transient) problems.
                    !.md<com>
                    !.md<com>        - Use FIELD=<tt>int</tt> to set nodal velocities to a deformable body.
                    !.md<com>          Where <tt>int</tt> is the code field of nodal velocities.
                    !.md<com>        - Use VX=<tt>real</tt>, VY=<tt>real</tt>, VZ=<tt>real</tt> to set initial
                    !.md<com>          velocities for rigid body or to all nodes for deformable bodies.
                    !
                    if ( words(3) == 'FIELD' ) then
                       kfl_invel_sld = getint('FIELD',1_ip,'#Field Number for initial velocity')
                    else
                       kfl_invel_sld = -1_ip
                       invel_sld(1) = getrea('VX   ',0.0_rp,'#x-component of velocity')
                       invel_sld(2) = getrea('VY   ',0.0_rp,'#y-component of velocity')
                       if ( ndime == 3_ip ) invel_sld(3) = getrea('VZ   ',0.0_rp,'#z-component of velocity')
                    end if

                 else if( words(2) == 'STATE' ) then
                    !
                    !.md<com>    - Set STATE to prescribe initial state dependent variables at the beginning of the simulation.
                    !.md<com>This option can only be used for material models with state dependent variables.
                    !.md<com>        - Use FIELD=<tt>int</tt> to set element values at gauss point level.
                    !.md<com>          Where <tt>int</tt> is the code number for the field of state dependent variables.
                    !
                    if ( words(3) == 'FIELD' ) then
                       kfl_insdv_sld = getint('FIELD',1_ip,'#Field Number for initial state dependent variable')
                    end if

                 else if( words(2) == 'TEMPE' ) then
                    !
                    !.md<com>    - Set TEMPERATURE to prescribe an initial temperature for all nodes
                    !.md<com>      at the beginning of the simulation.
                    !.md<com>This option can only be used for material models including temperature.
                    !.md<com>        - Use FIELD=<tt>int</tt> to set nodal values.
                    !.md<com>          Where <tt>int</tt> is the code number for the field of temperature variable.
                    !.md<com>        - Use VALUE=<tt>real</tt> to set nodal temperature values to all mesh.         
                    !
                    call sld_atm_reabcs()

                 end if
                 
              else if( words(1) == 'COORD' ) then
                 !
                 ! Type of coordinate system 
                 ! NOTE: this is only used for stents
                 if ( words(2) == 'BASIS' ) then
                    nposi = 2_ip
                    kfl_local_sld = 1_ip
                    wcsys = getcha('BASIS','     ','#Local coordinate system')
                    if (      wcsys == 'CARTE' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_CARTESIAN
                    else if ( wcsys == 'CYLIN' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_CYLINDRICAL
                    else if ( wcsys == 'SPHER' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_SPHERICAL
                    else if ( wcsys == 'EXNOR' ) then
                       kfl_csysl_sld = SLD_CSYS_EXNOR
                    end if
                    !
                    ! Read parameters
                    if ( wcsys /= 'EXNOR' ) csysl_sld(1:ndime*3) = param(nposi:ndime*3+nposi-1)
                 end if

              else if( words(1) == 'CONTA' ) then

                 !----------------------------------------------------------
                 !
                 ! Contact parameters
                 !
                 !----------------------------------------------------------
                 !
                 !.md<2>CONTACT_PDN                        $ PDN Contact
                 !.md<3>UNILATERAL | BILATERAL | RBO_DEFORMABLE
                 !.md<3>BOUNDARY=   int1
                 !.md<3>PROJECTION: X | Y | Z
                 !.md<2>END_CONTACT_PDN
                 !.md<field>CONTACT_PDN
                 !.md<com>Partial Dirichlet-Neumann contact algorithm (Rivero et. al 2018)
                 !.md<com>
                 !
                 kfl_local_sld = 1_ip ! local axes
                 !
                 call ecoute('sld_reabcs')
                 do while( words(1) /= 'ENDCO' )

                    if (words(1) == 'CODES') then
                       conbo_sld(1) = int(param(1))  !EXTERIOR SURFACE
                       conbo_sld(2) = int(param(2))  !INTERIOR SURFACE

                    else if ( words(1) == 'STENT' ) then
                       if ( words(2) == 'CRIMP' ) then
                         kfl_conta_stent = 1_ip
                         if ( words(3) == 'RMIN' ) then
                           r_fin_stent = param(3)
                         else if ( words(3) == 'RELAX' ) then
                           kfl_conta_stent = -1_ip
                         end if
                       else if ( words(2) == 'EXPAN' ) then
                         kfl_conta_stent = 2_ip
                         if ( words(3) == 'RMAX' ) then
                           r_fin_stent = param(3)
                         else if ( words(3) == 'RELAX' ) then
                           kfl_conta_stent = -2_ip
                         end if
                       else if ( words(2) == 'CHARG' ) then
                           kfl_conta_stent = 3_ip
                       else if (words(2) == 'MOVEV' ) then
                           kfl_conta_stent = 4_ip
                       end if

                    else if ( words(1) == 'FRICT' ) then
                       contact_friction_sld = param(1)

                    else if ( words(1) == 'UNILA' ) then
                       kfl_conta_sld = 1

                    else if ( words(1) == 'BILAT' ) then
                       kfl_conta_sld = 2
                       if ( exists('NEUMA') ) neumann_relax        = getrea('NEUMA',0.9_rp,'#Neumann relaxation factor')
                       if ( exists('TOLER') ) coupling_contact_tol = getrea('TOLER',1.0E-4_rp,'#Coupling contact tolerance')
                       if ( exists('ITERA') ) coupling_contact_its = getint('ITERA',20_ip,'#Coupling contact iterations')

                    else if ( words(1) == 'RBODE' ) then
                       kfl_conta_sld = 3

                    else if ( words(1) == 'BOUND' ) then
                       contactbou_sld = int(param(1),ip)

                    else if ( words(1) == 'PROJE' ) then
                       if (      words(2) == 'X    ' ) then
                          vect_proje_sld(1:3) = (/ 1.0_rp,0.0_rp,0.0_rp /)
                       else if ( words(2) == 'Y    ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,1.0_rp,0.0_rp /)
                       else if ( words(2) == 'Z    ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,0.0_rp,1.0_rp /)
                       else if ( words(2) == '-X   ' ) then
                          vect_proje_sld(1:3) = (/ -1.0_rp,0.0_rp,0.0_rp /)
                       else if ( words(2) == '-Y   ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,-1.0_rp,0.0_rp /)
                       else if ( words(2) == '-Z   ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,0.0_rp,-1.0_rp /)
                       else
                          vect_proje_sld(1:3) = (/ 0.0_rp,1.0_rp,0.0_rp /)
                       end if

                    else if ( words(1) == 'PENET' ) then
                       petol_sld = param(1)
                       
                    else if( words(1) == 'RELEA') then
                       if( words(2) == 'PREVI' ) then
                          kfl_mrele_sld = 1_ip
                       else
                          kfl_mrele_sld = 0_ip
                       end if
                    else if( words(1) == 'CONVE' ) then
                       kfl_contf_sld = 1_ip
                    end if
                    call ecoute('sld_reabcs')

                 end do

              else if( words(1) == 'NEWBO' ) then
                 !
                 !.md<2>NEW_BOUNDARY_CONDITIONS: ON | OFF  $ New boundary conditions
                 !.md<field>NEW_BOUNDARY_CONDITIONS
                 !.md<com>Set displacements, velocities and accelerations to zero after a restart.
                 !.md<com>This is useful when you want to start a new simulation (after restart) 
                 !.md<com>with new boundary conditions respect to the previous simulation. 
                 !.md<com>By default is set to <tt>OFF</tt>.
                 !.md<com>
                 !
                 if( words(2)=='ON' ) kfl_newbc_sld = 1_ip
                 
              else if( words(1) == 'FOLLO' ) then
                 !
                 !.md<2>FOLLOWER_LOADS: ON | OFF           $ Follower loads
                 !.md<field>FOLLOWER_LOADS
                 !.md<com>Activation of follower loads such as hydrostatic pressure from a fluid on a solid.
                 !.md<com>By default is set to <tt>OFF</tt>.
                 !.md<com>
                 !
                 if( words(2)=='ON' ) kfl_follo_sld = 1_ip
                 
              end if
              call ecoute('sld_reabcs')
           
           end do
           !
           !.md<1>END_PARAMETERS
           !
        else if( words(1) == 'CONTA' ) then
           
           !----------------------------------------------------------
           !
           ! Contact algorithm
           !
           !----------------------------------------------------------
           !
           !.md<1>CONTACT_ALGORITHM                  $ PDN Contact
           !.md<2>RBO_DEFORMABLE
           !.md<2>BOUNDARY=        int1
           !.md<2>PENETRATION_TOL= real
           !.md<2>CONVERGENCE:     ON | OFF
           !.md<1>END_CONTACT_ALGORITHM
           !.md<field>CONTACT_ALGORITHM
           !.md<com>Clean version of the partial Dirichlet-Neumann contact
           !.md<com>algorithm (Rivero et. al 2018) for SHERLOC project
           !.md<com>
           !
#if defined COMMDOM && COMMDOM == 2
           call sld_pdn_contact_reabcs()
           call messages_live("SOLIDZ: PDN CONTACT ALGORITHM USING PLE++")
#endif
           
        else if( words(1) == 'CODES' .and. exists('NODES') ) then

           !----------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !----------------------------------------------------------
           !
           !.md<com><strong>Codes</strong>
           !
           !.md<1>CODES, NODES
           !.md<2>int1 int2 rea1 rea2 rea3  $ int1=code, int2=fixity codes, real1-real3=DoFs 
           !.md<2>int1 int2 rea1 rea2 rea3, MODULE_FUNCTION=int3 | SPACE_&_TIME_FUNCTIONS=char1
           !.md<2>int1 int2 rea1 rea2 rea3, AXES: LOCAL_BASIS=int4  
           !.md<2>...
           !.md<1>END_CODES
           !.md<field>CODES, NODES
           !.md<com>Interpret the code on nodes to impose displacement for each degrees of freedom. <tt>int1</tt> is the node code, 
           !.md<com>e.g. 1, 2, 3 etc. If the node has multiple codes (e.g. if the code was extrapolated from boundary codes), the
           !.md<com>syntaxis is <tt>int1 = 1 & 2 & 4</tt>. It means that nodes with the three codes 1,2 and 4 are considered;
           !.md<com><tt>int2</tt> is the fixity code for each DoF (11, 10, 01 in 2D and 101, 111, etc in 3D). 0 means free or
           !.md<com>initial condition, 1 means prescribed.
           !.md<com>In solid mechanics are the displacement fixity for each degrees of freedom.
           !.md<com><tt>rea1</tt>, <tt>rea2</tt> and <tt>rea3</tt> are the corresponding values for each degree of freedom 
           !.md<com>(first to third displacement component).
           !.md<com><ul>
           !.md<com><ul>
           !.md<com><li>Set <tt>MODULE_FUNCTION=int3</tt> to use solidz discrete functions.
           !.md<com>The function has to be defined in caseName.sld.dat file as <tt>FUNCTIONS</tt> section.</li>
           !.md<com>
           !.md<com><li>Set <tt>SPACE_&_TIME_FUNCTIONS=char1</tt> to use a space and time function defined in 
           !.md<com><tt>caseName.ker.dat</tt> file.</li>
           !.md<com>
           !.md<com><li>Set <tt>AXES: LOCAL_BASIS=int4</tt> to use local boundary conditions. In this case a <tt>LOCAL_BASIS</tt>
           !.md<com>must be defined in the <tt>caseName.ker.dat file</tt>.</li>
           !.md<com></ul>
           !.md<com></ul>
           !.md<com>
           !
           if( exists('VELOC') ) then
              !
              ! Velocity
              !
              kfl_bvesv_sld = 1_ip
              tncod => tncod_sld(:)
              call boundary_conditions_read_node_codes('VELOCITY')
              
           else
              !
              ! Displacement
              !
              kfl_bvesv_sld = 0_ip
              tncod => tncod_sld(:)
              call boundary_conditions_read_node_codes('DISPLACEMENT')
              
           end if


        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !
           !-------------------------------------------------------------
           !
           !.md<1>CODES, BOUNDARIES
           !.md<2>int1 int2 int3 int4 int5   $ int1=code, int2=fixbo, int3-int5=DoFs
           !.md<2>...
           !.md<1>END_CODES
           !.md<field>CODES, BOUNDARIES:
           !.md<com>Impose a natural boundary condition (Neumaan) of type <tt>int2</tt> on boundaries with code <tt>int1</tt>. 
           !.md<com>The different boundary conditions available are:
           !.md<com><ul>
           !.md<com><ul>
           !.md<com><li> Set int2 = 3 to apply pressure</li>
           !.md<com><li> Set int2 = 8 to apply springs using F = -k x + mu dx/dt: boundary_code 8 k mu</li>
           !.md<com><li> Set int2 = 9 to apply cardiac cycle: boundary_code 9 cavity_code</li>
           !.md<com><li> Set int2 = 10 to apply cardiac cycle with springs F = -k x + mu dx/dt: boundary_code 10 
           !.md<com>cavity_code k mu</li>
           !.md<com></ul>
           !.md<com></ul>
           !.md<com><ul>
           !.md<com> The Neumaan boundary conditions can be applied using functions such as:
           !.md<com><ul>
           !.md<com><li>Set <tt>MODULE_FUNCTION=int3</tt> to use solidz discrete functions.
           !.md<com>The function has to be defined in caseName.sld.dat file as <tt>FUNCTIONS</tt> section.</li>
           !.md<com>
           !.md<com><li>Set <tt>SPACE_&_TIME_FUNCTIONS=char1</tt> to use a space and time function defined in 
           !.md<com><tt>caseName.ker.dat</tt> file.</li>
           !.md<com>
           !.md<com></ul>
           !.md<com></ul>
           !.md<com>
           !
           tbcod => tbcod_sld(1:)
           call boundary_conditions_read_boundary_codes('DISPLACEMENT')

        else if( words(1) == 'CODES' ) then
           
            call runend('CODES section without NODES or BOUNDARIES')

        else if ( words(1) == 'CRACK' ) then

           !----------------------------------------------------------
           !
           ! Cracks
           !
           !----------------------------------------------------------

           if( exists('NUMBE') ) then
              ncrak_sld = getint('NUMBE',1_ip,'#Number of cracks')
           else
              call runend('SLD_RABCS: NUMBER OF CRACKS MISSING')
           end if
           call ecoute('sld_reabcs')
           call sld_membcs(31_ip)
           if( ndime == 2 ) then
              dummi = 2
           else
              dummi = 4
           end if
           do while( words(1) /= 'ENDCR' )
              if( words(1) == 'CRACK') then
                 icrak = getint('CRACK',1_ip,'#CRACK NUMBER')
                 if( icrak > ncrak_sld .or. icrak < 1 ) call runend('SLD_REABCS: WRONG NUMBER OF CRACKS')
                 call ecoute('sld_reabcs')
                 ipoin = 0
                 do while( words(1) /= 'ENDCR' )
                    ipoin = ipoin + 1
                    if( ipoin > dummi ) call runend('SLD_REABCS: WRONG CRACK DEFINITION')
                    do idime = 1,ndime
                       crkco_sld(idime,ipoin,icrak) = param(idime)
                    end do
                    call ecoute('sld_reabcs')
                 end do
              end if
              call ecoute('sld_reabcs')
           end do

        else if ( words(1) == 'BODYF' ) then

           !----------------------------------------------------------
           !
           ! Body forces
           !
           !----------------------------------------------------------

           kfl_bodyf_sld = int(param(1),ip)

        else if ( words(1) == 'FUNCT' ) then

           !----------------------------------------------------------
           !
           ! Functions
           !
           !----------------------------------------------------------
           !
           !.md<com><strong>Solidz functions</strong>
           !
           !.md<1>FUNCTIONS                            $ Functions (transient only)
           !.md<2>TOTAL_NUMBER: int1               $ Total number of functions
           !.md<2>CONDITION
           !.md<3>FUNCTION_NUMBER: int2        $ Function number
           !.md<3>TIME_SHAPE: DISCRETE         $ Function type
           !.md<3>SHAPE_DEFINITION
           !.md<3>int3                         $ Number of lines
           !.md<3>rea1 rea2 rea3 rea4          $ rea1=time, rea2-4=DoFs
           !.md<3>...
           !.md<3>END_SHAPE_DEFINITION
           !.md<2>END_CONDITION
           !.md<1>END_FUNCTIONS
           !.md<field>FUNCTIONS
           !.md<com>Solidz time functions applied to nodes/boundaries.
           !.md<com><ul>
           !.md<com><ul>
           !.md<com><li>
           !.md<com>Set <tt>DISCRETE,LINEAR</tt>
           !.md<com>to apply a linear step (ramp) between discrete times.
           !.md<com></li>
           !.md<com><li>
           !.md<com>Set <tt>DISCRETE,SMOOTH</tt> to apply a smooth step between a time period.
           !.md<com>This function is such that first and second derivatives are zero between two discrete times.
           !.md<com>This definition is intended to ramp up or down smoothly from one amplitude value to another.
           !.md<com></li>
           !.md<com></ul>
           !.md<com></ul>
           !.md<com>
           !
           call ecoute('sld_reabcs')

           ! Initializations
           nauxi = 0
           ifunc = 0

           ! Total number of functions
           if (words(1) == 'TOTAL') nfunc_sld =  int(param(1),ip)
           if (nfunc_sld == 0) then
              call runend('SLD_REABCS: PROVIDE A TOTAL_NUMBER OF FUNCTIONS')
           else if (nfunc_sld > 9_ip) then
              call runend('SLD_REABCS: THE TOTAL_NUMBER OF FUNCTIONS MUST BE LOWER THAN 10')
           end if

           do while(words(1) /= 'ENDFU')
              if(kfl_conbc_sld == 0) then                   ! non-constant (transient) BC

                 ! Condition
                 do while(words(1) /= 'ENDCO')
                    if (words(1) == 'FUNCT') then
                       ! Function number
                       ifunc = getint('FUNCT',1_ip,'#FUNCTION NUMBER')
                       kfl_funty_sld(8,ifunc) = 1_ip        ! default: only apply function to fixity codes 1
                       if (ifunc < 0 .or. ifunc > 10_ip) then
                          call runend('SLD_REABCS: WRONG FUNCION NUMBER, MUST BE GT.0 AND LT.10')
                       else
                          nauxi = nauxi + 1
                          if (nauxi > nfunc_sld) call runend('SLD_REABCS: MORE FUNCTIONS THAN THE TOTAL_NUMBER ARE GIVEN')
                       end if

                       ! Function type
                    else if ( words(1)=='TIMES' ) then      ! time shape

                       if ( words(2)=='PARAB' .or. words(2)=='LINEA' .or. words(2)=='POLYN' ) then
                          kfl_funty_sld(1,ifunc)=1
                          if (words(2) == 'LINEA') kfl_funty_sld(1,ifunc)=10
                          kfl_funty_sld(2,ifunc)=20
                          kfl_funty_sld(3:7,ifunc)=1        ! default: apply to all equations

                       else if ( words(2) == 'PERIO' ) then
                          kfl_funty_sld(1,ifunc)=2
                          kfl_funty_sld(2,ifunc)=20
                          kfl_funty_sld(3:7,ifunc)=1        ! default: apply to all equations

                       else if ( words(2) == 'DISCR' ) then
                          kfl_funty_sld(1,ifunc)=3_ip
                          if (words(3) == 'SMOOT') then
                             kfl_funty_sld(2,ifunc)= 1_ip
                          else
                             kfl_funty_sld(2,ifunc)= 0_ip
                          end if

                          ! Shape defintion
                          call ecoute('sld_reabcs')
                          rauxi(:) = 0.0_rp
                          if ( words(1) == 'SHAPE' ) then   ! defining the shape by discrete points
                             ! Reference time and value
                             rauxi(1)= getrea('REFER', 0.0_rp, 'Reference value, sometimes useful')
                             rauxi(2)= getrea('TSTAR', 0.0_rp, 'Startint time, sometimes useful')

                             ! Save total number of tabular lines
                             call ecoute('sld_reabcs')
                             mtloa_sld(ifunc) = int(param(1),ip)

                             ! Allocate the prescription time function vector for ifunc
                             call sld_membcs(10_ip + ifunc)

                             ! Save discrete times and values
                             call ecoute('sld_reabcs')
                             itloa = 0
                             do while(words(1)/='ENDSH')
                                itloa= itloa + 1
                                tload_sld(ifunc)%a(ndime+1,itloa)= param(1) - rauxi(2)         ! time
                                tload_sld(ifunc)%a(1:ndime,itloa)= param(2:ndime+1) - rauxi(1) ! prescribed value
                                call ecoute('sld_reabcs')
                             end do
                          end if

                          kfl_funty_sld(3:7,ifunc)=1       ! default: apply to all equations

                       else
                          kfl_funty_sld(1,ifunc)=0

                       end if

                       !
                       ! OLD WAY: funcrion without an input file
                       !
                       if(kfl_funty_sld(1,ifunc)>0) then
                          if(kfl_funty_sld(1,ifunc) /= 3) then
                             fubcs_sld(1:10,ifunc)=param(2:11)
                          end if
                       end if

                    else if ( words(1)=='REFVA' ) then    ! reference values

                       fubcs_sld(15,ifunc)=param(1)
                       fubcs_sld(16,ifunc)=param(2)
                       fubcs_sld(17,ifunc)=param(3)
                       fubcs_sld(18,ifunc)=param(4)
                       fubcs_sld(19,ifunc)=param(5)

                    else if ( words(1)=='TIMEL' ) then    ! time lapse

                       fubcs_sld(11,ifunc)   =param(1)  !   start
                       fubcs_sld(12,ifunc)   =param(2)  !   end
                       fubcs_sld(13,ifunc)   =1.0e10_rp !   when repeat (default value, very high)
                       rtico_sld( 1,ifunc)    =param(1) !   initial start time for all eqs.
                       rtico_sld( 2,ifunc)    =param(2) !   initial final time for all eqs.
                       if (kfl_funty_sld(1,ifunc)==10) then   ! LINEAR: divide the load interval in equispaced steps
                          kfl_funty_sld(1,ifunc)=1
                          fubcs_sld(3:10,ifunc)= 0.0_rp
                          fubcs_sld(2,ifunc)   = 1.0_rp/(fubcs_sld(12,ifunc)-fubcs_sld(11,ifunc))
                       end if

                    else if ( words(1)=='TIMER' ) then    ! time repeat

                       fubcs_sld(13,ifunc)=param(1)     !   when repeat

                    else if ( words(1)=='FIXIT' ) then    ! special behavior controlled by the fixity code

                       if (exists('CONST')) then
                          kfl_funty_sld(8,ifunc) = 1    ! constrained in the direction with fixity code other than 1
                          call runend('SLD_REABCS: CONSTRAINED IS A DEPRECATED OPTION')
                       else if (exists('UNCON')) then   ! free in the direction with fixity code other than 1
                          kfl_funty_sld(8,ifunc) = 2
                          call runend('SLD_REABCS: UNCONSTRAINED IS A DEPRECATED OPTION')
                       else if (exists('ALLDI')) then   ! apply function to all dimensions
                          kfl_funty_sld(8,ifunc) = 3
                          call runend('SLD_REABCS: ALLDIMENSIONS IS A DEPRECATED OPTION')
                       end if

                    end if

                    call ecoute('sld_reabcs')
                 end do
                 ! End condition

              end if

              call ecoute('sld_reabcs')

           end do
           ! End functions

        end if

        call ecoute('sld_reabcs')

     end do
     !
     !.md<0>END_BOUNDARY_CONDITIONS
     !.md</code>
     ! 
     ! give a default size to non-defined tloads
     do ifunc= 1,10
        if (mtloa_sld(ifunc) == 0) then
           mtloa_sld(ifunc) = 1_ip         ! done to allocate a default memory space
           call sld_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
        end if
     end do

  end if

end subroutine sld_reabcs
