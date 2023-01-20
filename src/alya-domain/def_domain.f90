!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_domain

  !-----------------------------------------------------------------------
  !****f* defmod/def_domain
  ! NAME
  !   def_domain
  ! DESCRIPTION
  !   This module is the header of the domain
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp_basic,               only : ip,rp,r1p,r3p,i1p,i2p
  use def_kintyp_boundary_conditions, only : bc_bound
  use def_kintyp_boundary_conditions, only : bc_nodes
  use def_kintyp_domain,              only : ompss_domain
  use def_kintyp_domain,              only : mfiel
  use def_kintyp_domain,              only : elm
  use def_kintyp_domain,              only : elmgp
  use def_kintyp_domain,              only : elm_cloud
  use def_kintyp_domain,              only : typ_element_bin
  use def_kintyp_domain,              only : typ_lobas
  use def_kintyp_domain,              only : nelty
  use def_kintyp_mesh,                only : mesh_type
  use def_elmtyp,                     only : element_max
  use mod_htable,                     only : hash_t    
  use def_kintyp_dims
  
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip),   parameter :: &
       mfree          =  200, &      ! Max # free surfaces
       mnode_max      =   64, &      ! Max # nodes per element
       interval_funno = 1000         ! kfl_funno now has diferent meaning if it is >1000 or not.
  
  !------------------------------------------------------------------------
  ! Units
  !------------------------------------------------------------------------
  
  integer(ip)              :: &
       lun_pdata_dom ,        &      ! Domain data file unit
       lun_outpu_dom,         &      ! Output domain file unit
       lun_elsta_dom,         &      ! Elsest statistics
       lun_elmsh_dom,         &      ! Elsest mesh
       lun_elres_dom                 ! Elsest results

  !------------------------------------------------------------------------
  ! Dimensions: read in readim
  !------------------------------------------------------------------------

!!$  integer(ip)              :: &
!!$       kfl_autbo,             &      ! Automatic boundaries
!!$       npoin,                 &      ! # of nodal points
!!$       nelem,                 &      ! # of elements
!!$       necnt,                 &      ! # of contact elements
!!$       nncnt,                 &      ! # of contact nodes
!!$       nboun,                 &      ! # of boundary elements
!!$       nperi,                 &      ! # periodic nodes
!!$       lexis(nelty),          &      ! List of existing elements
!!$       nexis,                 &      ! Number of different element types in the geometry file
!!$       utype,                 &      ! Value of the unique type, -1 if several types
!!$       npoib,                 &      ! # immersed nodes
!!$       nboib,                 &      ! # immersed boundaries
!!$       lexib(nelty),          &      ! List of existing IB elements
!!$       nhang,                 &      ! # hanging nodes
!!$       nimbo,                 &      ! # IB
!!$       nrbod,                 &      ! # RB
!!$       nzone,                 &      ! Number of zones
!!$       nsubd,                 &      ! Number of subdomains
!!$       nmate,                 &      ! # of materials
!!$       nfiel,                 &      ! Number of fields
!!$       kfl_field(7,mfiel),    &      ! Fields dimensions
!!$       mcodb                         ! Max # codes
!!$
!!$#ifdef NDIMEPAR
!!$  ! much more comfortable this way; you can have a forder unix2d with -DTWODIM in the config.in
!!$#ifdef TWODIM
!!$  integer(ip), parameter   :: ndime = 2  ! # of space dimensions
!!$#else
!!$  integer(ip), parameter   :: ndime = 3  ! # of space dimensions
!!$#endif
!!$#else
!!$  integer(ip)              :: ndime      ! # of space dimensions
!!$#endif

    
  !------------------------------------------------------------------------
  ! Strategy: read in reastr
  !------------------------------------------------------------------------

  integer(ip), target      :: &
       ngaus(nelty),          &      ! # of Gauss points per element
       ngaib(nelty),          &      ! IB: # of Gauss points per element
       ngaus_cloud_gp(nelty)         ! # of Gauss points per element in the cloud of points
  integer(ip)              :: &
       lquad(nelty),          &      ! List of quadrature
       kfl_ngrou,             &      ! Strategy to construct groups
       ngrou_dom,             &      ! Groups (for deflated CG)
       ngrou_dom_target,      &      ! Groups wanted (for deflated CG) (could not correspond to the one computed)
       ngrou_boxes_coarse,    &      ! Number of coarse boxes for SFC
       ngrou_boxes_fine,      &      ! Number of fines boxes for SFC
       kfl_savda,             &      ! Save element data base
       lquib(nelty),          &      ! IB: List of quadrature
       kfl_geome,             &      ! Geometrical normals should be computed
       kfl_convx,             &      ! What to do with convex nodes
       kfl_frees,             &      ! Freestream criterion
       kfl_extra,             &      ! Extrapolate from boundary to nodes
       npbcs(8),              &      ! # parameters geometrical bcs
       lsbcs(100,8),          &      ! # list of geometrical bcs
       kfl_chege,             &      ! Check geometry
       kfl_naxis,             &      ! Axi-symmetry
       kfl_spher,             &      ! Spherical
       kfl_bouel,             &      ! Boundary-element connectivity
       kfl_divid,             &      ! Divide element into TET04
       curvatureDataField,    &      ! The field that occupies the curved data for mesh division
       curvatureField,        &      ! The field that occupies the curved geometry for mesh division
       kfl_xscal_fields              ! Scaling of fields
  
  real(rp),    pointer     :: &
       xscal_fields(:)               ! Scaling for the fields

  integer(ip), pointer     :: &
       materials_nlaye(:),    &      ! Automatic generation of materials
       materials_icode(:),    &      ! Automatic generation of materials
       materials_imate(:)            ! Automatic generation of materials

  real(rp)                 :: &
       xscal(3),              &      ! Geometric scale factors
       trans(3)                      ! Geometric translation factors
  real(rp)                 :: &
       awind,                 &      ! Wind angle (for freestream condition)
       tolan,                 &      ! Tolerance used to define inflow from freestream
       geoan                         ! Geometrical angle

  !------------------------------------------------------------------------
  ! Geometry: read in reageo
  !------------------------------------------------------------------------

  integer(ip),             pointer :: &
       lnods(:,:),                    &      ! Interior element connectivity
       ltype(:),                      &      ! List of element types
       lesub(:),                      &      ! List of element subdomains
       lgaus(:),                      &      ! List of element Gauss points
       lnodb(:,:),                    &      ! Boundary element connectivity
       ltypb(:),                      &      ! List of boundary types
       lboch(:),                      &      ! List of boundary characteristics
       lelbo(:),                      &      ! List of boudnary elements (old LBOEL)
       lnnod(:),                      &      ! Element number of nodes
       lelch(:),                      &      ! Element characteristic
       lmate(:),                      &      ! Materials (elements)
       lnoch(:),                      &      ! List of node characteristics
       lmast(:),                      &      ! List of masters for periodicity
       lgrou_dom(:),                  &      ! List of groups (deflated CG)
       lperi(:,:)                            ! List of Master/Slave
  real(rp),              pointer   :: &
       coord(:,:),                    &      ! Coordinates
       coord_cloud_gp(:,:)                   ! Coordinates for cloud of points
  real(rp),    pointer             :: &
       cooin(:,:)                            ! Initial Coordinates
  type(i1p),   pointer             :: &
       lhang(:)                              ! List of hanging nodes
  type(r3p)                        :: &
       xfiel(mfiel)                          ! Fields
  type(r1p)                        :: &
       time_field(mfiel)                     ! Time for Fields for more than 1 step

  !------------------------------------------------------------------------
  ! Sets: read in reaset
  !------------------------------------------------------------------------

  integer(ip)              :: &
       neset,                 &      ! # of element sets
       nbset,                 &      ! # of boundary sets
       nnset,                 &      ! # of node sets
       neset_origi,           &      ! # of element sets
       nbset_origi,           &      ! # of boundary sets
       nnset_origi                   ! # of node sets
  integer(ip), pointer     :: &
       leset(:),&                    ! List of element sets
       lbset(:),&                    ! List of boundary sets
       lnset(:)                      ! List of node sets

  !------------------------------------------------------------------------
  ! Sets: read in reabcs
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_icodn,             &      ! If codes on nodes
       kfl_icode,             &      ! If codes on edges
       kfl_icodb,             &      ! If codes on boundaries
       mcono                         ! Max # codes per nodes
  integer(ip),   pointer   :: &
       kfl_codno(:,:),        &      ! Node codes
       kfl_codbo(:),          &      ! Boundary codes
       kfl_coded(:,:)                ! Edge codes

  !------------------------------------------------------------------------
  ! Global variables
  !------------------------------------------------------------------------
  !
  ! Dimensions
  !
#ifdef PNODE_VALUE
  integer(ip), parameter   :: mnode = PNODE_VALUE
#endif
  integer(ip)              :: &
       nedge,                 &      ! Number of edges
       medge,                 &      ! Maximum number of edges
       mecnt,                 &      ! Total number of contact elements
#ifndef PNODE_VALUE
       mnode,                 &      ! Max # of nodes per element
#endif
       mnoga,                 &      ! Max mnode and mgaus
       mnoib,                 &      ! Max # of nodes per IB (boundary)
       mnodi,                 &      ! Max # of nodes per IB (volume)
       mnodb,                 &      ! Max # of nodes per boundary
       mgaus,                 &      ! Max # of Gauss points per element
       mgaib,                 &      ! Max # of Gauss points per IB (surface)
       mgaui,                 &      ! Max # of Gauss points per IB (volume)
       mgaub,                 &      ! Max # of Gauss points per boundary
       mlapl,                 &      ! Max # of llapl(ielty), ielty=1,nelty
       nfacs,                 &      ! Number of faces
       ntens,                 &      ! # of components of symmetric tensors
       nrule,                 &      ! # of integration rule
       nbopo,                 &      ! # of boundary points
       ndimb,                 &      ! # of space dimensions-1
       nzdom,                 &      ! # of nonzero elements in the mesh graph
       nzdom_own,             &      ! # of nonzero elements in the own mesh graph
       nzbou,                 &      ! # of nonzero elements in the boundary mesh graph
       nzsky,                 &      ! # of comp. in the skyline matrix of the graph
       nzsol,                 &      ! = nzdom + 2*nslav: Components due to slaves
       nzmat,                 &      ! = max(nzsol,nzsky,nzexp): Components of A
       nzmbt,                 &      ! Components of B
       nzrhs,                 &      ! = max(nzsol,nzsky,nzexp): Components of RHS
       nzpre,                 &      ! Components of Preconditioner
       neige,                 &      ! Size of eigen value vector
       neiva,                 &      ! Number of eigenvalues
       nzerr,                 &      ! Components of the Error Estimator
       elmin,                 &      ! Element with minimum volume
       elmax,                 &      ! Element with maximum volume
       nzsym,                 &      ! # of nonzero elements in the symmetric mesh graph
       nxexa,                 &      ! # x-element for example
       nyexa,                 &      ! # y-element for example
       nzexa,                 &      ! # z-element for example
       nxexp,                 &      ! # x-partition for example
       nyexp,                 &      ! # y-partition for example
       nzexp,                 &      ! # z-partition for example
       kfl_parex,             &      ! Automatic parallelization of domain
       nbono,                 &      ! Number of boundary nodes
       gp_total_cloud_gp,     &      ! Total number of Gauss points in the cloud of points
       mgaus_cloud_gp                ! Maximum number of Gauss points in the cloud of points
  integer(ip)              :: &
       nzmax,                 &      ! = max(nzsol,nzsky,nzexp): Components of A
       nzrhx,                 &      ! = max(nzsol,nzsky,nzexp): Components of RHS
       nzprx                         ! Components of Preconditioner

  !$acc declare create(mnode)  
  
  integer(8)               :: &
       memor_dom(2)                  ! Memory counter
  integer(ip),     pointer :: &
       lpoib(:),              &      ! IB List of point IB
       nmatn(:),              &      ! Number of nodes per material
       nodemat(:),            &      ! Nodal material (computed from elements material)
       lbono(:),              &      ! List of nodes attached to boundaries
       lnnob(:),              &      ! Boundary number of nodes
       lesec(:),              &      ! Element set connectivity
       lbsec(:),              &      ! Boundary set connectivity
       lnsec(:)                      ! Node set connectivity
  type(i1p),   pointer     :: &
       lmatn(:)                      ! Materials (nodes)
  !
  ! Edge
  !
  integer(ip),     pointer :: &
       edge_to_node(:,:),     &      ! NEDGE
       ledgs(:,:),            &      ! NELEM: element edge connectivity  (see LNODS)
       ledgb(:,:),            &      ! NBOUN: boundary edge connectivity (see LNODB)
       lnned(:),              &      ! NELEM: Number of edge by element  (see LNNOD)
       lnneb(:),              &      ! NBOUN: Number of edge by element  (see LNNOB)
       c_edg(:),              &      ! Edge graph
       r_edg(:)                      ! Edge graph
  !
  ! Domain properties
  !
  integer(ip)              :: &
       bandw_dom,             &      ! Bandwidth
       naved_dom,             &      ! Average number of edges
       nmied_dom,             &      ! Min number of edges
       nmaed_dom                     ! Max number of edges
  real(rp)                 :: &
       vodom,                 &      ! Measure of the domain
       vomin,                 &      ! Minimum element volume
       vomax,                 &      ! Maximum element volume
       voave,                 &      ! Averaged element volume
       profi_dom                     ! Profile
  real(rp),        target  :: &
       xmima(2,3)                    ! Bounding box
  real(rp),       pointer  :: &
       vomat(:)                      ! Material volume
  !
  ! Reals
  !
  real(rp)                 :: &
       permx(9)                      ! Rotation matrix relating periodic faces
  !
  ! Graph
  !
  integer(ip)              :: &
       mepoi,                 &      ! Max. # of elements by node
       mpopo,                 &      ! Max. # of point-point connectivity
       kfl_crbou,             &      ! If boundary graph has been computed
       kfl_pelel,             &      ! If element graph has been computed
       kfl_domar,             &      ! If geometrical arrays should be recomputed (then goto domarr)
       kfl_domar_world,       &      ! If geometrical arrays should be recomputed (in Alya world)
       npoin_ii,              &      ! Interior nodes
       npoin_bb                      ! Boundary nodes
  integer(ip), pointer :: &
       nepoi(:),              &      ! # of neighbor elements
       pelpo(:),              &      ! Pointer node/element connectivity
       lelpo(:),              &      ! List node/element connectivities
       pelpo_2(:),            &      ! Pointer node/element extended connectivity
       lelpo_2(:),            &      ! List node/element extended connectivities
       pelel(:),              &      ! Pointer element/element connectivity
       lelel(:),              &      ! List element/element connectivities
       pelel_2(:),            &      ! Pointer element/element extended connectivity
       lelel_2(:),            &      ! List element/element extended connectivities
       lezdo(:,:,:),          &      ! Lnods to graph array
       lbzdo(:,:,:)                  ! Lnodb to graph array

  integer(ip), pointer     :: &
       r_sol(:),              &      ! Row array for the matrix CSR storage
       c_sol(:)                      ! Column array for the matrix CSR storage
  integer(ip), pointer     :: &
       r_dom(:),              &      ! Row array for the domain CSR storage
       c_dom(:),              &      ! Column array for the domain CSR storage
       r_bou(:),              &      ! Row array for the boundary domain CSR storage
       c_bou(:),              &      ! Column array for the boundary domain CSR storage
       r_dom_aii(:),          &      ! Row array for the domain CSR storage
       permr_aii(:),          &      ! Permutation for Aii
       invpr_aii(:),          &      ! Inverse permutation for Aii
       c_dom_aii(:),          &      ! Column array for the domain CSR storage
       r_dom_aib(:),          &      ! Row array for the domain CSR storage
       c_dom_aib(:),          &      ! Column array for the domain CSR storage
       r_dom_abi(:),          &      ! Row array for the domain CSR storage
       c_dom_abi(:),          &      ! Column array for the domain CSR storage
       r_dom_abb(:),          &      ! Row array for the domain CSR storage
       c_dom_abb(:),          &      ! Column array for the domain CSR storage
       permr_abb(:),          &      ! Permutation for Abb
       invpr_abb(:),          &      ! Inverse permutation for Abb
       r_dom_prec(:),         &      ! Graph for Schur preconditioner
       c_dom_prec(:),         &      ! Graph for Schur preconditioner
       r_dom_own(:),          &      ! Full row graph
       c_dom_own(:),          &      ! Full row graph
       permr_prec(:),         &      ! Permutation for Schur preconditioner
       invpr_prec(:),         &      ! Inverse permutation for Schur preconditioner
       r_sym(:),              &      ! Row array for the matrix CSR symmetric storage
       c_sym(:)                      ! Column array for the matrix CSR symmetric storage

  !------OMPS DOMAIN--------
  type(ompss_domain),   pointer :: ompss_domains(:)
  type(ompss_domain),   pointer :: ompss_boundaries(:)
  !------END OMPSS DOMAIN------

  !
  ! Integer Arrays
  !
  integer(ip), pointer     :: &
       lpoty(:),              &      ! List of point types
       lfacs(:,:),            &      ! List faces
       lboel(:,:),            &      ! List of boundary elements
       lgaib(:),              &      ! IB: number of Gauss points on IB
       lrenn(:),              &      ! Node renumbering
       lfcnt(:,:),            &      ! List of contact element faces
       lncnt(:,:),            &      ! List of contact element nodes
       lessl(:,:)                    ! Groups in parallel
  integer(ip)              :: &
       lnuty(element_max),    &      ! Number of elements for each ielty
       lnuib(nelty)                  ! Number of IB types
  !
  ! Node and boundary codes
  !
  integer(ip), pointer     :: &
       kfl_fixno(:,:),        &      ! General node fixity
       kfl_fixbo(:),          &      ! General boundary fixity
       kfl_funno(:),          &      ! Function number for nodes
       kfl_funbo(:),          &      ! Function number for boundaries
       kfl_funtb(:),          &      ! Function type on boundaries
       kfl_funtn(:),          &      ! Function type on nodes
       kfl_fixrs(:),          &      ! Axes
       kfl_geobo(:),          &      ! Geometrical boundary b.c
       kfl_geono(:),          &      ! Geometrical node b.c.
       wallo(:),              &      ! Order of the Index of the nearest wall node
       lpoin(:)                      ! Type of point (for geometrical arrays)
  integer(ip)              :: &
       iffun,                 &      ! If function should be read
       ifloc,                 &      ! If axis should be read
       ifbop,                 &      ! If bc are imposed on boundary nodes
       ifbes                         ! If value should be assigned
  real(rp),    pointer     :: &
       bvess(:,:),            &      ! General node fixity
       bvnat(:,:)                    ! General boundary fixity
  type(bc_nodes), pointer  :: &
       tncod(:)                      ! Node code type
  type(bc_nodes), pointer  :: &
       tgcod(:)                      ! Geometrical Node code type
  type(bc_bound), pointer  :: &
       tbcod(:)                      ! Boundary code type
  !
  ! Real Arrays
  !
  real(rp),    contiguous, pointer     :: &
       vmass(:),              &      ! Lumped mass matrix
       vmasc(:),              &      ! Mass matrix with close rule
       dmass(:),              &      ! Lumped mass with diagonal scaling
       cmass(:)                      ! Consistent mass matrix
  real(rp),    pointer     :: &
       cmass_weighted(:),     &      ! Consistent weighted mass matrix
       mass_right_fact(:),    &      ! Matrix for approx mass
       walld(:),              &      ! Distance to the wall
       wallcoor(:,:),         &      ! Coordinates of the to the nearest wall points
       walln(:,:),            &      ! Normal to the wall
       rough(:),              &      ! Roughness
       canhe(:),              &      ! Canopy height
       heiov(:),              &      ! Height over terrain
       canla(:),              &      ! Canopy Leaf Area Density
       ywalb(:),              &      ! Distance to the wall at each boundary for variable wall distance
       ywalp(:),              &      ! Projection of boundary wall distance onto boundary nodes
       ywale(:),              &      ! Minimum of ywalb for each element
       yscab(:),              &      ! Variable scaling distance for machine learning
       yscap(:)                      ! Projection of yscap on to nodes
  !
  ! Coordinate systems and local basis
  !
  real(rp),        pointer :: &
       exnor(:,:,:),          &      ! Exterior normal
       skcos(:,:,:)                  ! Cosine matrices of skew systems
   integer(ip)             :: &
       num_lobas                     ! Number of coordinate systems
  type(typ_lobas), pointer :: &      ! Coordinate system
       lobas(:) 
  !
  ! Element shape functions and derivatives
  !
  type(elm),   pointer     :: &
       elmar(:)                      ! Element data base
  type(elmgp), pointer     :: &
       elmda(:)                      ! Element Gauss point data base
  real(rp), pointer        :: &
       elmda_gpvol(:,:),      &
       elmda_gpcar(:,:,:,:)
  type(elm_cloud), pointer :: &
       elmar_cloud_gp(:)             ! Element data base for cloud of points

  real(rp)                 :: &
       hnatu(nelty)                  ! Natural element length
  integer(ip)              :: &
       lenex(mnode_max+1,nelty),&    ! List of next element node
       ldime(nelty),          &      ! List of element dimensions
       ltopo(nelty),          &      ! List of element topology
       llapl(nelty),          &      ! List of element Laplacian
       lrule(nelty),          &      ! List of element integration rules
       linte(nelty),          &      ! List of element interpolation methods
       lruib(nelty),          &      ! List of IB integration rules
       lorde(nelty),          &      ! List of element order
       nnode(-nelty:nelty),   &      ! List of element # of nodes
       iesta_dom,             &      ! Where element starts
       iesto_dom,             &      ! Where element stops
       ibsta_dom,             &      ! Where boundary starts
       ibsto_dom,             &      ! Where boundary stops
       kfl_horde,             &      ! IF high order element exist
       kfl_elcoh,             &      ! Cohesive elements
       kfl_elint,             &      ! Interface elements
       mface,                 &      ! Maximum number of faces
       lrule_cloud_gp(nelty)         ! List of element integration rules for the cloud of points
  !
  ! File name
  !
  character(150)           :: &
       fil_outpu_dom                 ! Output domain mesh file
  !
  ! Old mesh data
  !
  integer(ip)              :: &
       npoin_old,             &      ! # of nodal points
       nelem_old,             &      ! # of elements
       nboun_old                     ! # of boundary elements
  !
  ! Special arrays for Level Set reinitialization
  !
  integer(ip)              :: &
       nelwh                         ! # of total elements in the whole mesh
  integer(ip), pointer     :: &
       pefpo(:),              &      ! Pointer node/element (where interface must be sought) connectivity
       lefpo(:),              &      ! List node/element (where interface must be sought) connectivities
       lnuew(:)                      ! Numeration of an element in the whole mesh

  integer(ip), pointer     :: &
       leldo(:,:)                    ! Fringe elements: subdomains/local numbering
  !
  ! Mesh multiplication
  !
  integer(ip),  pointer    :: &
       facel(:,:,:)                  ! List of faces

  !------------------------------------------------------------------------
  !
  ! Mesh structures
  !
  !------------------------------------------------------------------------

  integer(ip),   pointer  :: lnlev(:)        ! List of node level
  integer(ip),   pointer  :: lelev(:)        ! List of element level
  integer(ip),   pointer  :: lblev(:)        ! List of boundary level

  type(mesh_type), pointer :: meshe(:)
  integer(ip),     pointer :: lpmsh(:)
  integer(ip),     pointer :: lemsh(:)
  integer(ip),     pointer :: lbmsh(:)
  !
  ! Cut elements structure
  !
  type subel_type
     integer(ip)               :: inout      ! Type of elemen: -1 if is inside, 1 1 if is outside.
     real(rp),pointer          :: elcod(:,:) ! Coordinates of a subelement
  end type subel_type
  type subbo_type
     real(rp),pointer          :: bocod(:,:) ! Coordinates of the boundary
  end type subbo_type

  type cutel_type
     integer(ip)               ::  nelem     ! List of cut elements
     integer(ip)               ::  iimbo     ! Id of the cut particle
     type(subel_type), pointer ::  l(:)      ! List of subelements inside each cut element
     integer(ip)               ::  nboun     ! List of boundaries
     type(subbo_type), pointer ::  lb(:)     ! List of boundary elements that cut an element
     integer(ip),      pointer ::  linou(:)  ! Determina if a gauss point is inside or outside the particle
  end type cutel_type

  type(cutel_type),    pointer :: cutel(:)
  !
  ! Element bin
  !
  integer(ip)                    :: element_bin_boxes(3)
  type(typ_element_bin), pointer :: element_bin(:)
  !
  ! Additional vectors recalculated every time step for transient fields
  !
  real(rp)                  :: x_tran_fiel(mfiel)      ! interpolation for transient field
  integer(ip)               :: k_tran_fiel(mfiel)      ! begining of interval for transient field. Always = 1 for fields loaded on demand
  integer(ip)               :: k_tran_fiel_real(mfiel) ! real value of the k_tran_fiel. For checking if the timestep changed and files need to be loaded
  integer(ip), parameter    :: nsteps_fiel_ondemand=2  ! how many timesteps allocate for the fields loaded on demand (kfl_field(6,ifiel) == 1)
  integer(ip)               :: kexist_tran_fiel        ! do transient fields exist?
  !
  ! Finite volume arrays
  !
  real(rp),         pointer ::   &
       fv_center_coord(:,:) ,    &                     ! Coordinates of centroides
       fv_cell_volume(:)    ,    &                     ! Volume of cells
       fv_face_area(:)      ,    &                     ! Face area
       fv_face_normal(:,:)  ,    &                     ! Face normals
       fv_face_orientation(:),   &                     ! Face orientation
       fv_center_distance(:),    &                     ! Distance bewteen centroids
       fv_center_vector(:,:),    &                     ! Vector between centroids
       fv_center_face(:,:,:)                           ! Vector from centroide to face centroide
  integer(ip),      pointer ::   &
       fv_face_boundary(:),      &                     ! Correspondance face booundary
       fv_face_graph(:,:),       &                     ! Correspondance face to graph position
       fv_graph_diag(:)                                ! Diagonal position of element in graph
  type(hash_t)              :: htable_lninv_loc        ! Hash table for global->local

  
end module def_domain
