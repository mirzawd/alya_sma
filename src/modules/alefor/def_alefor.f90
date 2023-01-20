!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_alefor
  !------------------------------------------------------------------------
  !****f* Alefor/def_alefor
  ! NAME 
  !    def_alefor
  ! DESCRIPTION
  !    Heading for the Alefor routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip),   parameter              :: &
       mfunc_ale = 50                        ! Maximum number of functions
  integer(ip),   parameter              :: &
       nprev_ale = 3                         ! Maximum number of functions

  integer(ip), parameter :: &
       lun_outpu_ale = 714,                &
       lun_resta_ale = 715
  ! -----------------------------------------------------------------------
  ! 
  ! Variable to choose if the ALE module is to be solved
  !
  ! -----------------------------------------------------------------------
  integer(ip) :: &
       kfl_solve_ale                         ! Variable to choose if the ALE module is to be solved
  !------------------------------------------------------------------------
  !
  ! Physical problem : rigid body
  !
  !------------------------------------------------------------------------
  !--BEGIN REA GROUP
  integer(ip)                           :: &
       kfl_rigid_ale,                      & ! Solve rigid body problem
       kfl_sensi_ale,                      & ! Field number for the sensitivities
       nstli_ale(3),                       & ! Step to start linear motion for each dimension
       nstro_ale(3),                       & ! Step to start rotation motion for each dimension
       kfl_mvext_ale,                      & ! Move external boundaries of the domain in x & y idem CoG 
       kfl_ralei_ale,                      & ! Raleight damping
       nstra_ale,                          & ! Start eliminating Raleight damping
       nenra_ale,                          & ! End eliminating Raleight damping
       kfl_catfo_ale,                      & ! Add external forces for the catamaran
       kfl_genco_ale,                      & ! Solve the rigid body equations in generalized coordinates
       kfl_sprin_ale,                      & ! Use the spring model for the rigid body with generalized coordinates
       kfl_pertu_ale,                      & ! Use perturbation for generalized forces
       kfl_topmo_ale,                      & ! Use the non-spinning top model for the rigid body with generalized coordinates
       kfl_grafo_ale                         ! Add gravity force

  real(rp)                              :: &
       ralei_ale,                          & ! Raleight damping parameter
       sprin_ale(3_ip),                    & ! Spring properties, 1: k (stiffness), 2: damping. 3:scale for non-dimen. force
       defor_param_ale                       ! Parameter for deformation
  
  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in ale_reanut
  !
  !------------------------------------------------------------------------
  
  integer(ip)                           :: &
       kfl_smoot_ale,                      & ! No mesh smoothing
       kfl_timef_ale,                      & ! Timefunction (only has sense when ALEFOR runs alone)
       nsmoo_ale,                          & ! Number of smoothing steps (loading steps)
       kfl_defor_ale,                      & ! No mesh deformation
       ndefo_ale,                          & ! Number of deformation steps (loading steps)
       kfl_smobo_ale,                      & ! Boundary smoothing
       kfl_fixsm_ale,                      & ! Fix boundary nodes by default
       nsmob_ale,                          & ! Number of boundary smoothing iterations
       moddi_ale,                          & ! Initial Displacement field number
       modvi_ale,                          & ! Initial Velocity field number
       kfl_crist_ale,                      & ! Cristobal 1 , new 0  ! RIGID BODY
       kfl_foexo_ale,                      & ! Force ( & Torque) extrapolation order  ! RIGID BODY
       kfl_disor_ale,                      & ! Integration order for the RB displacements  ! RIGID BODY
       kfl_nforc_ale                         ! Use force N or average N,N-1 for kfl_crist_ale

  real(rp)                              :: &                        
       ansmo_ale,                          & ! Sharp edge detection
       resmo_ale                             ! Relaxation factor
       
  integer(ip)                           :: &
       npoin_ad,                           & ! Support geometry for MM: number of nodes
       nboun_ad,                           & ! Support geometry for MM: number of boundaries
       mnodb_ad                              ! Support geometry for MM: max number of nodes per boundary
   integer(ip),    pointer              :: &
       lnodb_ad(:,:),                      & ! Support geometry for MM: connectivity
       ltypb_ad(:)  ,                      & ! Support geometry for MM: connectivity
       kfl_fixrs_ale(:)                      ! Reference system for the BV
  real(rp),        pointer              :: &
       coord_ad(:,:)                         ! Support geometry for MM: coordinates
  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in ale_reabcs
  !
  !------------------------------------------------------------------------
 
  type(bc_nodes), pointer               :: &     
       tncod_ale(:)                          ! Node code type
  type(bc_nodes), pointer               :: &     
       tgcod_ale(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &     
       tbcod_ale(:)                          ! Boundary code type
  
  !------------------------------------------------------------------------
  !
  ! Boundary conditions
  ! kfl_fixno_ale is in master because it is used as coupling variable
  !
  !------------------------------------------------------------------------

!--END REA GROUP
  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------
  !
  ! Boundary conditions
  !
  integer(ip), pointer                  :: &
       kfl_funno_ale(:),                   & ! Nodal function 
       kfl_funbo_ale(:),                   & ! Boundary function 
       kfl_funtn_ale(:),                   & ! Nodal # function 
       kfl_funtb_ale(:),                   & ! Nodal type function 
       kfl_fixbo_ale(:)                      ! Boundary fixity
  !
  ! General
  !
  real(rp),    pointer                  :: &
       coord_ale(:,:,:),                   &  ! Update coordinates
       bvess_ref(:,:),                     &  ! Coeficients of directions for displacements set *.ale.dat to allow for motion in multiple directions
       coord_ori(:,:)                         ! Original reference

  real(rp)                              :: &                        
       xrota_ale(3),                       & ! Rigid body rotation multiplier (obtained in ale_begste) 
       xline_ale(3)                          ! Rigid body linear motion multiplier (obtained in ale_begste) 

end module def_alefor
