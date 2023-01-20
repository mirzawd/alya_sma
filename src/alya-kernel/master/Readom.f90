!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    readom.f90
!> @author  houzeaux
!> @date    2018-11-16
!> @brief   Read domain data
!> @details Read all the variables and arrays of dom.dat file
!-----------------------------------------------------------------------

subroutine Readom()

  use def_kintyp,                   only : ip
  use mod_ker_timeline,             only : ker_timeline
  use mod_messages,                 only : messages_live
  use mod_mpio_par_configure,       only : mpio_init
  use mod_read_domain_arrays,       only : read_domain_arrays_reageo
  use mod_read_domain_arrays,       only : read_domain_arrays_reaset
  use mod_read_domain_arrays,       only : read_domain_arrays_reabcs
  use mod_read_domain_arrays,       only : read_domain_arrays_reafie
  use mod_read_domain_arrays,       only : read_domain_arrays_finalize
  use mod_run_config,               only : run_config
  implicit none  
  
  call messages_live('READ MESH ARRAYS','START SECTION')
  if( run_config%timeline ) call ker_timeline(0_ip,'INI_READ_MESH')
  !
  ! Open files
  !
  call opfdom(1_ip)
  !
  ! Read geometry data
  !
  call readim()                                 ! Read dimensions:      NDIME, NPOIN, NELEM, LEXIS...
  call reastr()                                 ! Read strategy:        LQUAD, NGAUS, SCALE, TRANS...
  call par_broadcast_dimensions()               ! Broadcast what have just been read
  call cderda()                                 ! Derivated param:      NTENS, HNATU, MNODE, MGAUS.
  !
  ! Determine the reading mode
  !
  call mpio_init()

  !-------------------------Arrays----------------------------------------------------------------------
  !
  ! Read geometry arrays                        ! LTYPE, LNODS, LMATE, LESUB, LELCH, LNODB,
  !                                             ! LTYPB, LBOEL, LELBO, LBOCH, COORD, LNOCH
  call read_domain_arrays_reageo()
  !
  ! Read sets                                   ! LESET, LBSET, LNSET
  !
  call read_domain_arrays_reaset()
  !
  ! Read bcs                                    ! KFL_CODNO, KFL_CODBO
  !
  call read_domain_arrays_reabcs()
  !
  ! Read fields                                 ! XFIEL.TAG1.TAG2
  !
  call read_domain_arrays_reafie()  
  !
  !-------------------------Arrays----------------------------------------------------------------------
  !
  ! Deallocate MPIO/IO arrays
  !
  call read_domain_arrays_finalize()
  !
  ! Compute LELBO if the elements connected to the boundaries are unknown (kfl_bouel=0)
  !
  call lbouel()
  !
  ! Automatic boundaries
  !
  call dombou()

  if( run_config%timeline ) call ker_timeline(0_ip,'END_READ_MESH')
  call messages_live('READ MESH ARRAYS','END SECTION')
  !
  ! Close domain unit
  !
  call opfdom(3_ip)
  !
  ! Some checks according to optimization options
  !
  call mescek(8_ip)
  
end subroutine Readom

!------------------------------------------------------------------------
!
! READIM
! ------
!
!   Variable                      | Description                                                                | Name              | Compulsory
! | ----------------------------- | -------------------------------------------------------------------------- | ----------------- | ----------    
! | <tt>ndime</tt>                | Space dimension (if not defined with an ifdef through <tt>-DNDIMEPAR</tt>) | SPACE_DIMENSIONS  | yes   
! | <tt>npoin</tt>                | Number of nodes                                                            | NODAL_POINTS      | yes   
! | <tt>nelem</tt>                | Number of elements                                                         | ELEMENTS          | yes                      
! | <tt>nboun</tt>                | Number of boundaries                                                       | BOUNDARIES        | -                      
! | <tt>nperi</tt>                | Number of periodic nodes                                                   | PERIODIC_NODES    | -              
! | <tt>mcodb</tt>                | Maximum number of boundary codes                                           | MCODB             | -              
! | <tt>lexis(:)</tt>             | List of volume elements present in the mesh                                | TYPES_OF_ELEMENTS | yes                    
! | <tt>nfiel</tt>                | Number of fields                                                           | FIELDS            | -              
! | <tt>kfl_field(:,:)</tt>       | The dimensions and characteristics of each field                           | FIELD             | -    
! | <tt>time_field(:) % a(:)</tt> | Time parameters of the fields                                              | TIMES             | -     
! 
! REASTR
! ------
! 
!   Variable                      | Description                                                                 | Name                      | Compulsory           
! | ----------------------------- | --------------------------------------------------------------------------- | ------------------------- | ---------- 
! | <tt>lquad</tt>                | Integration rule for each type of element                                   | INTEGRATION_RULE          | -   
! | <tt>ngaus</tt>                | Number of integration points                                                | DOMAIN_INTEGRATION_POINTS | -   
! | <tt>ngaus</tt>                | Number of integration points (same as previous but in more intuitve format) | GAUSS_POINTS              | -   
! | <tt>ngrou_dom</tt>            | Number of groups for deflated solvers                                       | GROUPS                    | -   
! | <tt>kfl_ngrou</tt>            | Strategy for computing groups                                               | -                         | -   
! | <tt>xscal</tt>                | Scaling of the geometry                                                     | SCALE                     | -   
! | <tt>trans</tt>                | Translation of the geometry                                                 | TRANSLATION               | -   
! | <tt>kfl_geome</tt>            | Geometrical boundary conditions                                             | GEOMETRICAL               | -   
! | <tt>kfl_bouel</tt>            | If boundary elements are given in geometry fields (LELBO)                   | BOUNDARY_ELEMENT          | -   
! | <tt>kfl_chege</tt>            | Flag to check the validity of the mesh                                      | CHECK_GEOMETRY            | -   
! | <tt>kfl_extra</tt>            | If boundary codes should be extrapolated to nodes                           | EXTRAPOLATE               | -   
! | <tt>materials_icode</tt>      | To generate materials from boundaries: codes                                | MATERIALS / CODE          | -   
! | <tt>materials_imate</tt>      | To generate materials from boundaries: material number                      | MATERIALS / MATERIAL      | -   
! | <tt>materials_nlaye</tt>      | To generate materials from boundaries: number of layers                     | MATERIALS / LAYER         | -            
! 
! REGEO
! ------
! 
!   Variable                      | Description                                                                 | Name                      | Compulsory           
! | ----------------------------- | --------------------------------------------------------------------------- | ------------------------- | ---------- 
!
!
!
!
!
! NOTES
! 
!    Flags and scalars:
!       KFL_CHEGE .................. Check geometry
!       KFL_NAXIS .................. Axi-symmetric geometry (=1)
!    Element types:
!       NELTY ...................... Total number of element types (def_elmtyp)
!       MGAUS ...................... Maximum # of element Gauss points in the mesh
!       MLAPL ...................... 1 if at least one element needs Laplacian
!       CENAL(NELTY) ............... Alternative element name (for Ensight format)
!       CENAM(NELTY) ............... Element name: 3 letters (type) + 2 figures (# nodes)
!       CETOP(NELTY) ............... Topology name
!       HNATU(NELTY) ............... Element natural length: 2 Quad/Hexa, 1 Others
!       LDIME(NELTY) ............... Dimension
!       LLAPL(NELTY) ............... If Laplacian exists
!       LNUTY(NELTY) ............... Number of elements for corresponding type
!       LORDE(NELTY) ............... Order of the interpolation (1,2,3, etc)
!       LQUAD(NELTY) ............... Integration quadrature
!       LRULE(NELTY) ............... Integration rule
!       LTOPO(NELTY) ............... Topology (line, quad, tru, penta, pyr)
!       NGAUS(NELTY) ............... # of integration points
!       NNODE(NELTY) ............... Number of nodes
!    Element arrays:
!       NELEM ...................... Number of elements
!       MNODE ...................... Maximum number of nodes per element
!       NBOPO ...................... Number of boundary nodes
!       NTENS ...................... # of independent Hessian matrix values
!       NINER ...................... # of independent inertia tensor values
!    ** COORD(NDIME,NPOIN) ......... Node coordinates
!       EXNOR(NDIME,NDIME,NBOPO) ... Exterior normal
!       LEXIS(NELEM) ............... If element type exists
!    ** LNODS(MNODE,NELEM) ......... Connectivity
!       LPOTY(NPOIN) ............... Node type: 0=interior, IBOPO=boundary
!    ** LTYPE(NELEM) ............... Type (BAR02, QUA04, etc as defined in def_elmtyp)
!       VMASS(NPOIN) ............... Diagonal mass matrix
!       VMASC(NPOIN) ............... Diagonal mass matrix with close rule
!    Boundary arrays:
!       MGAUB ...................... Maximum # of boundary Gaus points in the mesh
!       MNODB ...................... Maximum number of nodes per boundary
!       NDIMB ...................... Boundary dimension (NDIME-1)
!       LBOEL(MNODB,NBOUN) ......... Boundary to element ndoes
!    ** LELBO(NBOUN) ............... Boundary connected element
!    ** LNODB(MNODB,NBOUN) ......... Connectivity
!       LTYPB(NBOUN) ............... Boundary Type
!       NBONO ...................... Number of boundary nodes (Explicitly declared in LNODB)
!       LBONO(NBONO) ............... List of boundary nodes (Explicitly declared in LNODB)
!    Graphs:
!       R_DOM(NPOIN+1) ............. Node/node       graph: pointer to C_DOM
!       C_DOM(NZDOM) ............... Node/node       graph: list of connectivities
!       LELPO(NPOIN+1) ............. Node/element    graph: pointer to PELPO
!       PELPO(LELPO(NPOIN+1)) ...... Node/element    graph: list of connectivities
!       LELEL(NPOIN+1) ............. Element/element graph: pointer to PELEL
!       PELEL(LELEL(NPOIN+1)) ...... Element/element graph: list of connectivities
!    Witness points:
!       MWITN ...................... Max # Number of witness points
!       NWITN ...................... # of witness points
!       LEWIT(NWITN) ............... Host elements of witness points
!       COWIT(NWITN) ............... Coordinates of witness points
!       SHWIT(MNODE,MWITN) ......... Shape function in host element
!    Materials:
!       NMATE ...................... # materials
!    ** LMATE(NELEM) ............... List of materials
!       LMATN(NPOIN) ............... List of nodol materials
!    Sets:
!       NESET ...................... # element sets
!       NBSET ...................... # boundary sets
!       NNSET ...................... # node sets
!    ** LESET(NELEM) ............... List of element sets
!       LESEC(NESET) ............... Element sets numbers
!    ** LBSET(NBOUN) ............... List of boundary sets
!       LBSEC(NBSET) ............... Boundary sets numbers
!    ** LNSET(NPOIN) ............... List of node sets
!       LNSEC(NNSET) ............... Node sets numbers
!    Subdomains:
!       NSUBD ...................... Number of subdomains
!    ** LESUB(NELEM) ............... List of element subdomains (value=1:NSUBD)
!    Fields:
!    ** XFIEL(NFIEL) % A(:,NPOIN/NBOUN/nelem) ... Values of field
!    Permutation with original mesh:
!       LNLEV(NPOIN) ............... Node permutation
!       LELEV(NELEM) ............... Element permutation
!       LBLEV(NBOUN) ............... Boundary permutation
!    Groups:
!       LGROU_DOM(NPOIN) ........... List of groups (computed)
!
!    README:
!    -------
!
!     *: scalar read from geometry file
!    **: array read from geometry file
!
!------------------------------------------------------------------------
