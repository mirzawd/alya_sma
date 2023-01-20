!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup SolidzInput
!> @{
!> @file    sld_reaous.f90
!> @author  Solidz Team
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Read post-process Solidz data
!>
!> @details Read post-process data
!>           - Array post-process (Field outputs)
!>           - Element, boundary and node sets
!>           - Witness points
!>           - Other outputs
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_reaous()
  !
  !.md<module>solidz
  !.md<input>case.sld.dat
  !.md<pos>2
  !.md<sec>
  !
  use def_kintyp,             only : ip, rp
  use def_inpout,             only : words, getint, getrea, getcha, exists
  use def_master,             only : INOTSLAVE
  use mod_memory,             only : memory_alloca
  use mod_ecoute,             only : ecoute
  use mod_output_postprocess, only : output_postprocess_read
  use def_solidz,             only : kfl_exacs_sld, kfl_foten_sld, kfl_rotei_sld
  use def_solidz,             only : kfl_psmat_sld, psmat_sld
  use mod_sld_post_reaction,  only : sld_post_reaction_reaous
  use mod_sld_strain,         only : strain_basis

  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations global variables
     !
     kfl_exacs_sld = 0                  ! Exact solution (FEM errors)
     kfl_foten_sld = 0                  ! Tensors of forces caust, green, lepsi (0 if none)
     kfl_rotei_sld = 0                  ! Rotate and correct sigma eigenvalues (ONLY FOR 3D PRISMATIC SHELLS)
     kfl_psmat_sld = 0                  ! Postscript file of the matrix
     psmat_sld(1)  = 1_ip
     psmat_sld(2)  = 1_ip
     !
     ! Reach the section
     !
     call ecoute('sld_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('sld_reaous')
     end do
     !
     ! Begin to read data
     !
     !.md# Output and Post-process
     !.md<code>
     !.md<0>OUTPUT_&_POST_PROCESS
     !
     !.md<com><strong>General post-process</strong>
     !
     !.md<1>START_POSTPROCESS_AT: STEP=int | TIME=int   $ Start from step number or time
     !.md<field>START_POSTPROCESS_AT:
     !.md<com>This option is used to start the post-process at a certain step or time in the simulation.
     !.md<com>Set TIME=<tt>int</tt> to specify the initial time for post-process
     !.md<com>or STEP=<tt>int</tt> to specify the first time step for post-process.
     !
     !.md<1>POSTPROCESS                                 $ Variable for post-process
     !.md<field>POSTPROCESS
     !.md<com>This option is to post-process a variable for visualitzation.
     !.md<com>Set STEPS=<tt>int</tt> equal to the post-process frequency.
     !.md<com>Results are written in <tt>caseName-XXXXX-YYYYYYYY.post.alyabin</tt> file. The following
     !.md<com>variables can be postprocess for visualization:
     !
     !.md<com><ul>
     !.md<com>Pre-process mesh checks:
     !.md<com><ul>
     !.md<com><li>BOCOD: Boundary codes for boundary conditions</li>
     !.md<com><li>NOCOD: Node codes for boundary conditions</li>
     !.md<com><li>NELEM: Element global number</li>
     !.md<com><li>NPOIN: Node global number</li>
     !.md<com><li>ELCHA: Element characteristic code</li>
     !.md<com><li>ELNOR: Element normal for <tt>HEX08</tt>,<tt>PEN06</tt> and <tt>QUA04</tt> elements</li>
     !.md<com><li>FIXNO: Fixity on nodes</li>
     !.md<com><li>FIXBO: Fixity on boundaries</li>
     !.md<com><li>FIXRS: Fixity on nodes using local basis</li>
     !.md<com><li>ROTM1: First vector rotation matrix</li>
     !.md<com><li>ROTM2: Second vector rotation matrix</li>
     !.md<com><li>ROTM3: Third vector rotation matrix</li>
     !.md<com><li>BOSET: Boundary sets</li>
     !.md<com><li>ELSET: Element sets</li>
     !.md<com><li>NOSET: Node sets</li>
     !.md<com><li>CELEN: Characteristic element lenght</li>
     !.md<com><li>PELCH: Element characteristic</li>
     !.md<com><li>PELTY: Element type</li>
     !.md<com><li>PERIO: Periodic nodes</li>
     !.md<com><li>PARTI: Subdomain partition at element level</li>
     !
     !.md<com>General solid mechanics variables:
     !.md<com><ul>
     !.md<com><li>DISPL: Displacements in global coordinate system</li>
     !.md<com><li>VELOC: Velocities in global coordinate system</li>
     !.md<com><li>ACCEL: Accelerations in global coordinate system</li>
     !.md<com><li>FRXID: Reaction forces</li>
     !.md<com><li>FEXTE: External forces</li>
     !.md<com><li>BVESS: Prescribed displacements on nodes</li>
     !.md<com><li>SIGMA: Stresses in global coordinate system</li>
     !.md<com><li>LNENR: Logarithmic strains in global coordinate system</li>
     !.md<com><li>SDVAR: State dependent variables at node level</li>
     !.md<com><li>TEMPE: Temperature</li>
     !.md<com><li>PRESS: Pressure</li>
     !.md<com><li>ENDEN: Energy density calculated from the material model</li>
     !.md<com><li>DONNA: Donnan osmosis swelling for sm400</li>
     !.md<com><li>WATER: Water content for poroelasticity, for sm400</li>
     !.md<com></ul>
     !.md<com>Variables for composites:
     !.md<com><ul>
     !.md<com><li>ORIEN: Orientation angle for composite materials</li>
     !.md<com><li>STACK: Element stacking direction for <tt>HEX08</tt>,<tt>PEN06</tt> and <tt>QUA04</tt> elements</li>
     !.md<com><li>AXIS1: Longitudinal fiber direction</li>
     !.md<com><li>AXIS2: Transversal fiber direction</li>
     !.md<com><li>AXIS3: Out-of-plane direction</li>
     !.md<com><li>DCOHE: Damage variable for interface cohesive elements</li>
     !.md<com><li>DAMAG: Group of damage variables for intra-laminar damage models</li>
     !.md<com><li>SRPRO: Strength reduction properties for sm154 model</li>
     !.md<com></ul>
     !.md<com>Variables for PDN-contact algorithm:
     !.md<com><ul>
     !.md<com><li>FCONT: Contact force for PDN-contact</li>
     !.md<com></ul>
     !.md<com></ul>
     !
     !.md<com><strong>Element sets</strong>
     !
     !.md<1>ELEMENT_SET                                    
     !.md<2>...                                     $ Variables computed in volume integrals
     !.md<1>END_ELEMENT_SET
     !.md<field>ELEMENT_SET
     !.md<com>Variables computed as volume integrals on element sets. The results are written in
     !.md<com><tt>caseName-element.sld.set</tt> file.
     !.md<com>The following variables can be post-processed:
     !.md<com><ul>
     !.md<com><ul>
     !.md<com><li>EPSRT: Averaged green strain using polar-cylindrical coordiantes (eps_rr + eps_tt)</li>
     !.md<com><li>LEPRT: Averaged logarithminc strain using polar-cylindrical coordiantes (eps_rr + eps_tt)</li>
     !.md<com><li>VELOC: Averaged velocity</li>
     !.md<com><li>TMASS: Total mass</li>
     !.md<com><li>ALLKE: All kinetic energy</li>
     !.md<com></ul>
     !.md<com></ul>
     !
     !.md<com><strong>Boundary sets</strong>
     !
     !.md<1>BOUNDARY_SET                               
     !.md<2>...                                     $ Variables computed in boundary integrals
     !.md<1>END_BOUNDARY_SET
     !.md<field>BOUNDARY_SET
     !.md<com>Variables computed as boundary integrals on boundary sets. The results are written in
     !.md<com><tt>caseName-boundary.sld.set</tt> file.
     !.md<com>The following variables can be post-processed:
     !.md<com><ul>
     !.md<com><ul>
     !.md<com><li>DIBOX, DIBOY, DIBOZ, DIBOU (Magnitude): Averaged displacement</li>
     !.md<com><li>FRBOX, FRBOY, FRBOZ, FRBOU (Magnitude): Sum of reaction forces</li>
     !.md<com><li>FCONX, FCONY, FCONZ, FCONT (Magnitude): Contact forces</li>
     !.md<com><li>FCONO, FCOT1, FCOT1: Contact forces in local system</li>
     !.md<com><li>PRESS: Pressure</li>
     !.md<com></ul>
     !.md<com></ul>
     !
     !.md<com><strong>Node sets</strong>
     !
     !.md<1>NODE_SET                                     
     !.md<2>...                                     $ Variables computed on nodes
     !.md<1>END_NODE_SET
     !.md<field>NODE_SET
     !.md<com>Variables on node sets. The results are written
     !.md<com>in <tt>caseName-node.sld.set</tt> file.
     !.md<com>The following variables can be post-processed:
     !.md<com><ul>
     !.md<com><ul>
     !.md<com><li>DISPX, DISPY, DISPZ: Displacements</li>
     !.md<com><li>VELOX, VELOY, VELOZ: Velocities</li>
     !.md<com><li>ACCEX, ACCEY, ACCEZ: Accelerations</li>
     !.md<com><li>FRXIX, FRXIY, FRXIZ: Reaction Forces</li>
     !.md<com><li>SIGXX, SIGYY, SIGZZ, SIGYZ, SIGXZ, SIGXY: Global stresses</li>
     !.md<com><li>EPSXX, EPSYY, EPSZZ, EPSYZ, EPSXZ, EPSXY: Global strains</li>
     !.md<com><li>LEPXX, LEPYY, LEPZZ, LEPYZ, LEPXZ, LEPXY: Global logarithmic strains</li>
     !.md<com><li>COORX, COORY, COORZ: Coordinates in reference configuration</li>
     !.md<com><li>SEQVM: Von Mises Stress</li>
     !.md<com><li>BVESX, BVESY, BVESZ: Dirichlet displacement</li>
     !.md<com><li>SBVNX, SBVNY, SBVNZ: Solver Neumann values</li>
     !.md<com><li>TEMPE: Temperature</li>
     !.md<com><li>ENDEN: Energy density calculated from the material model</li>
     !.md<com></ul>
     !.md<com></ul>
     !
     !.md<com><strong>Witness points</strong>
     !
     !.md<1>WITNESS_POINTS
     !.md<2>...                                     $ Variables computed witness points
     !.md<1>END_WITNESS_POINTS
     !.md<field>WITNESS_POINTS
     !.md<com>Variables on witness points. The results are stored
     !.md<com><tt>caseName.sld.wit</tt>.
     !.md<com>The following variables can be postprocessed:
     !.md<com><ul>
     !.md<com><ul>
     !.md<com><li>DISPX, DISPY, DISPZ: Displacements</li>
     !.md<com><li>VELOX, VELOY, VELOZ: Velocities</li>
     !.md<com><li>ACCEX, ACCEY, ACCEZ: Accelerations</li>
     !.md<com><li>FRXIX, FRXIY, FRXIZ: Reaction Forces</li>
     !.md<com><li>SIGXX, SIGYY, SIGZZ, SIGYZ, SIGXZ, SIGXY: Global stresses</li>
     !.md<com><li>EPSXX, EPSYY, EPSZZ, EPSYZ, EPSXZ, EPSXY: Global strains</li>
     !.md<com><li>LEPXX, LEPYY, LEPZZ, LEPYZ, LEPXZ, LEPXY: Global Logarithmic strains</li>
     !.md<com><li>COORX, COORY, COORZ: Coordinates</li>
     !.md<com><li>SEQVM: Von Mises Stress</li>
     !.md<com><li>TEMPE: Temperature</li>
     !.md<com><li>ENDEN: Energy density calculated from the material model</li>
     !.md<com></ul>
     !.md<com></ul>
     !.md<com>
     !

     call ecoute('sld_reaous')

     do while(words(1)/='ENDOU')
        call output_postprocess_read()
        !
        !.md<com><strong>Other outputs</strong>
        !
        if( words(1)=='OUTPU' ) then

           !----------------------------------------------------------
           !
           ! Other Outputs 
           !
           !----------------------------------------------------------

           if( words(2)=='ERROR' ) then
              !
              !.md<1>OUTPUT, ERROR, SOLUTION= int                $ Manufactured solution
              !.md<field>OUTPUT, ERROR
              !.md<com>Manufactured solutions used to carry out code typing checking and mesh convergence.
              !.md<com>The L2, L1 and Linf norm of displacement
              !.md<com>gradients can be found in file <tt>caseName.sld.log</tt>. Dirichlet boundary
              !.md<com>conditions are automatically imposed on all the boundaries.
              !.md<com>The manufactured solutions available are the following:
              !.md<com>
              !.md<com>   - SOLUTION = 2. Exact solution using Solidz Implicit integration scheme.
              !.md<com>
              !
              kfl_exacs_sld = getint('SOLUT',1_ip,'#Exact solution')

           else if( words(2) == 'MATRI' ) then
              !
              !.md<1>OUTPUT, MATRIX, STEP = int1, ITERA = int2   $ Matrix output
              !.md<field>OUTPUT, MATRIX
              !.md<com>Postprocess global system matrix. Set STEP=<tt>int1</tt> to choose the time step and
              !.md<com>ITERA=<tt>int2</tt> to choose the N-R iteration. The available formats for matrix output are the following:
              !.md<com>
              !.md<com>   - Postcript file (.ps) 
              !.md<com>   - Data file (.dat)
              !.md<com>
              !
              kfl_psmat_sld = 1_ip
              psmat_sld(1) = getint('STEP ',1_ip,'#Step to postprocess matrix')
              psmat_sld(2) = getint('ITERA',1_ip,'#Iteration to postprocess matrix')
              
           end if
           
        else if( words(1) == 'REACT' ) then
           !
           !.md<1>REACTION_SET, NUMBER=int, FORCE_DROP_PERCENTAGE=real, ON_REACTION_SET=int $ Reaction and displacement output
           !.md<2>REACTION_SET=int, PLANE: AXIS=<tt>int</tt>, GREATER_THAN=<tt>real</tt> LESS_THAN=<tt>real</tt>
           !.md<2>REACTION_SET=int, BOX:   XMIN=<tt>real</tt>, YMIN=<tt>real</tt>, ZMIN=<tt>real</tt>, 
           !.md<2>XMAX=<tt>real</tt>, YMAX=<tt>real</tt>, ZMAX=<tt>real</tt>
           !.md<2>...
           !.md<1>END_REACTION_SET
           !.md<field>REACTION_SET
           !.md<com>Postprocess reaction and displacement on Dirichlet nodes defined by the user. Different methods 
           !.md<com>can be used in order to select a set of nodes (internal and external). The results file name is
           !.md<com><tt>caseName-reaction.sld.res</tt>.
           !.md<com>
           !.md<com>   - Set PLANE: AXIS=<tt>int</tt>, GREATER_THAN=<tt>real</tt> or LESS_THAN=<tt>real</tt> 
           !.md<com>     to select a set of nodes using a global axis and distance greater or less than a fixed point.
           !.md<com>   - Set BOX: XMIN=<tt>real</tt>, YMIN=<tt>real</tt>, ZMIN=<tt>real</tt>, 
           !.md<com>     XMAX=<tt>real</tt>, YMAX=<tt>real</tt>, ZMAX=<tt>real</tt>
           !.md<com>     to select a set of nodes within box.
           !.md<com>
           !.md<com> Optional arguments to end the run due to a reaction force drop:
           !.md<com>
           !.md<com>   - Set FORCE_DROP_PERCENTAGE=<tt>int</tt> to define the percentage the reaction force drop. The
           !.md<com>     user have to define the code set using ON_REACTIONSET=<tt>int</tt>. By default it takes set no. 1.
           !.md<com>
           !
           call sld_post_reaction_reaous()

        else if( words(1) == 'ROTAT' ) then
           !
           ! ADOC[1]> ROTATION: ON | OFF                          $ Rotation
           !
           kfl_rotei_sld = 1
           
        else if (words(1) == 'STRAI') then
           strain_basis%on = .true.
           strain_basis%field_lng = getint('LONGI', 0_ip, '!Longitudinal strain vector field')
           strain_basis%field_rad = getint('RADIA', 0_ip, '!Radial strain vector field')
           strain_basis%field_cir = getint('CIRCU', 0_ip, '!Cirumferential strain vector field')
        end if

        call ecoute('sld_reaous')
     end do
     !
     !.md<0>END_OUTPUT_&_POST_PROCESS
     !.md</code>
     ! 
  end if

end subroutine sld_reaous
