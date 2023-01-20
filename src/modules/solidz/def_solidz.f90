!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    def_solidz.f90
!> @author  Mariano Vázquez
!> @author  Matias Rivero
!> @author  Gerard Guillamet
!> @author  Eva Casoni
!> @brief   Solidz global variables
!> @details
!> @}
!------------------------------------------------------------------------

module def_solidz

  use def_kintyp

  !----------------------------------------------------------------------
  !
  ! Parameters
  !
  !----------------------------------------------------------------------

  integer(ip), parameter              :: &
       SLD_EXPLICIT_SCHEME         =  1, &
       SLD_IMPLICIT_SCHEME         =  2, &
       SLD_DYNAMIC_PROBLEM         =  1, &
       SLD_STATIC_PROBLEM          =  0, &
       SLD_TANGENT_ANALYTICAL      =  0, &
       SLD_TANGENT_CHECK           =  1, &
       SLD_TANGENT_NUMERICAL       =  2, &
       SLD_SECANT                  =  3, &
       SLD_INFINITESIMAL           =  0, &
       SLD_GREEN                   =  1, &
       SLD_CSYS_CARTESIAN          =  0, &
       SLD_CSYS_CYLINDRICAL        =  1, &
       SLD_CSYS_SPHERICAL          =  2, &
       SLD_CSYS_EXNOR              =  3, &
       SLD_PDN_UNILATERAL          =  1, &
       SLD_PDN_BILATERAL           =  2, &
       SLD_PDN_RBO_DEFORMABLE      =  3, &
       SLD_CSHEL_EAS               =  1, &
       SLD_CSHEL_ANS_SHEAR         =  1, &
       SLD_CSHEL_ANS_TRAPEZOIDAL   =  1
  
  integer(ip), parameter               :: &
       lun_psmat_sld               = 114, &
       lun_carcy_res_sld           = 115, &
       lun_carcy_cvg_sld           = 116, &
       ! same value as the carcy file because alya opens one or the other
       lun_sysnet_heart_res_sld    = 115, & 
       ! same value as the carcy file because alya opens one or the other
       lun_sysnet_system_res_sld   = 116, & 
       lun_react_sld               = 117

  integer(ip),   parameter             :: &
       ncoef_sld                   =  50    ! # coefficient for material laws

  !--BEGIN REA GROUP
  !----------------------------------------------------------------------
  !
  ! Physical problem: read in sld_reaphy
  !
  !----------------------------------------------------------------------
  
  integer(ip)                         :: &
       kfl_timei_sld,                    & ! Existance of dd/dt (static or dynamic problem)
       kfl_rigid_sld,                    & ! Solve rigid body problem
       kfl_fiber_sld,                    & ! Anisotropic media with fibers
       kfl_rmate_sld,                    & ! Material axesrotation
       kfl_csysm_sld,                    & ! Existence of a material coordinate system
       kfl_restr_sld,                    & ! Residual Stress in voigt notation
       kfl_indis_sld(2),                 & ! Initial displacements
       kfl_prdef_sld,                    & ! Fibers field nodal- or element- wise
       kfl_sdvar_sld,                    & ! State dependent variables flag
       kfl_damag_sld,                    & ! Damage model
       kfl_cohes_sld,                    & ! Cohesive model
       kfl_cohft_sld,                    & ! Friction in cohesive model (at the momento in always 0)
       kfl_vofor_sld,                    & ! Applied external volume force
       kfl_moduf_sld(5),                 & ! Modulator fields for properties
       kfl_tange_sld,                    & ! Calculation of the tangent moduli (only for implicit computations)
       kfl_strai_sld,                    & ! Type of deformations (Large (finite) or Small strains)
       kfl_plast_sld,                    & ! Plasticity model
       kfl_donna_sld,                    & ! Donnan Osmosis for porous hyperelastic swelling
       kfl_excou_sld,                    & ! Exmedi coupling model (WITH EXMEDI)
       kfl_nopassive_sld,                & ! Only active model (only works for 134implicit!!)
       kfl_cshel_sld(3),                 & ! Parameters continuum shell
       kfl_compressible,                 & ! Comrpessible modification for HGO model
       ncalcio_sld,                      & ! Calcium values given from MUSICO
       oripa_sld(3)                        ! Orientation parameters

  real(rp)                            :: &
       gravi_sld(3),                     & ! Gravity vector
       grnor_sld,                        & ! Gravity norm
       invel_sld(3),                     & ! Initial velocity vector
       csysm_sld(9),                     & ! Material coordinate system points
       thick_sld,                        & ! Out-of-plane thickness (2-d problems)
       calmusico_sld(2,100),             & ! Ca read from MUSICO values
       gfibe_sld(3)                        ! Global fiber orientation

  integer(ip),  pointer               :: &
       modfi_sld(:),                     & ! Type of model used to define the fibers
       modor_sld(:,:),                   & ! Fibers with sheet-like structure (orthotropy)
       kfl_dampi_sld(:),                 & ! Rayleigh damping
       kfl_coupt_sld(:),                 & ! Coupling tensor (transversally zero, trans isotropic o trans aniso)
       lawco_sld(:),                     & ! Constitutive law
       lawst_sld(:),                     & ! Law for the stress model
       lawch_sld(:),                     & ! Law for cohesive model
       lawho_sld(:),                     & ! Law for Holzapfel model
       lawpl_sld(:),                     & ! Law for plasticity model
       lawta_sld(:),                     & ! Calculation of the tangent stiffness matrix for each material law.
       nocoh_sld(:),                     & ! List of nodes where the tractions exceed the cohesive limit
       lawmo_sld(:)                        ! Auxiliar when some versions of one constitutive model are defined
  real(rp),     pointer               :: &
       pres0_sld(:),                     & ! Pre-stress for cardiac
       parsp_sld(:,:),                   & ! Parameter of the spheroid-geometry model
       densi_sld(:,:),                   & ! Reference Density (rho)
       dampi_sld(:,:),                   & ! Rayleigh damping coefficients
       parco_sld(:,:),                   & ! Parameters built-in constitutive law
       parch_sld(:,:),                   & ! Parameters cohesive law
       parcf_sld(:,:),                   & ! Parameters contact and friction
       velas_sld(:,:)                      ! Reference velocity of sound

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in sld_reanut
  !
  !------------------------------------------------------------------------
  
  logical(lg)                           :: &
       kfl_vecto_sld                         ! Vectorised assembly
  
  integer(ip)                           :: &
       kfl_stabi_sld,                      & ! Stabilization
       kfl_ellen_sld,                      & ! =0,1 for min/max element length
       kfl_xfeme_sld,                      & ! Enrichment strategy (XFEM)
       kfl_xfcra_sld,                      & ! Crack/discontinuity definition (XFEM)
       kfl_resid_sld,                      & ! Types of residuals calculation
       kfl_timet_sld,                      & ! Time treatment
       kfl_ninex_sld,                      & ! Inexact Newton (approximated tangent stiffness matrix)
       kfl_tisch_sld,                      & ! Time integration scheme
       kfl_serei_sld,                      & ! Selective reduced integration
       kfl_limit_sld,                      & ! Limiter on
       kfl_plane_sld,                      & ! Plane-stress 2D assumption
       kfl_pseud_sld,                      & ! Pseudo time step
       kfl_penal_sld,                      & ! Penalization
       kfl_savdt_sld,                      & ! Save dt at the beginning of the simulation
       kfl_dttyp_sld,                      & ! Strategy for dt calculation
       kfl_celen_sld,                      & ! Type of calculation characteristic lenght
       kfl_isoch_sld,                      & ! Isochrones for solid
       kfl_safet_table_sld,                & !
       minex_sld,                          & ! Inexact Newton Counter
       nisaf_sld,                          &
       miinn_sld,                          & ! Max # of iterations
       last_iters_sld,                     & ! Solver iterations for the current sub-iteration
       miinn_pseud_sld                       ! Max # pseudo-time of iterations

  real(rp)                              :: &
       dafac_sld,                          & ! Damping factor for stabilize
       tifac_sld(5),                       & ! Time integration factors (for CN, newmark, etc.)
       cotol_sld,                          & ! Convergence tolerance
       safet_sld,                          & ! Safety factor for time step
       safet_table_sld(2,10),              & ! Safety factor for time step
       safet_pseud_sld,                    & ! Safety factor for pseudo time step
       factor_penal_sld,                   & ! Penalization factor
       safex_sld,                          & ! Safety factor for time step
       epsex_sld,                          & ! Perturbation to compute the secant stiffness matrix
       safma_sld,                          & ! Safety factor for time step
       sstol_sld,                          & ! Steady state tolerance
       masss_sld                             ! Mass scaling factor

  !----------------------------------------------------------------------
  !
  ! Boundary conditions: read in sld_reabcs
  !
  !----------------------------------------------------------------------

  integer(ip)                         :: &
       kfl_conbc_sld,                    & ! Constant b.c.
       kfl_newbc_sld,                    & ! New boundary conditons after restart
       kfl_local_sld,                    & ! Local system of reference
       kfl_csysl_sld,                    & ! Existance of a local coordinate system for boundary conditions
       kfl_bodyf_sld,                    & ! Body force function id
       kfl_follo_sld,                    & ! Follower loads
       kfl_bvesv_sld,                    & ! Flag for velocity prescription
       kfl_invel_sld,                    & ! Initial condition: velocity
       kfl_insdv_sld,                    & ! Initial condition: state dependent variables
       kfl_funty_sld(20,20),             & ! Function type and number of paremeters
       kfl_windk_sld,                    & ! Windkessel functions
       kfl_conta_stent,                  & ! Stent options: crimping, expansion, charge
       mtloa_sld(10),                    & ! Number of data points for the transient boundary functions
       nfunc_sld,                        & ! Number of transient boundary conditions function
       ncrak_sld                           ! Number of cracks
  real(rp)                            :: &
       csysl_sld(9),                     & ! Local coordinate system parameters
       rtico_sld(2 ,20),                 & ! Time counters for the transient boundary conditions
       r_fin_stent,                      & ! Final radius imposed to stent in crimpring or expansion
       fubcs_sld(20,20)
  type(bc_nodes), pointer             :: &
       tncod_sld(:)                        ! Node code type
  type(bc_nodes), pointer             :: &
       tgcod_sld(:)                        ! Geometrical node code type
  type(bc_bound), pointer             :: &
       tbcod_sld(:)                        ! Boundary code type
  type(r2p),      pointer             :: &
       tload_sld(:)                        ! time dependent boundary functions (to be read in a file)
  real(rp),       pointer             :: &
       crkco_sld(:,:,:)                    ! Initial crack coordinates

  !----------------------------------------------------------------------
  !
  ! Output and Postprocess: read in sld_reaous
  !
  !----------------------------------------------------------------------

  integer(ip)                         :: &
       kfl_exacs_sld,                    & ! Exact solution
       kfl_psmat_sld,                    & ! Postprocess global system matrix
       kfl_rotei_sld,                    & ! Rotate and correct sigma in prismatic shells
       kfl_foten_sld                       ! Tensors of forces are used (caust_sld, green_sld and lepsi_sld)
  integer(ip)                         :: &
       psmat_sld(2)                        ! Parameters postprocess matrix
  real(rp)                            :: &
       rorig_sld(3),                     & ! Postprocess rotations origin
       thiso_sld(2)                        ! Threashold for isochrones
  real(rp),     pointer               :: &
       isoch_sld(:,:)                      ! Isochrones for solidz

  !--END REA GROUP
  !----------------------------------------------------------------------
  !
  ! Others
  !
  !----------------------------------------------------------------------
  !
  ! Boundary conditions
  !
  integer(ip),  pointer               :: &
       kfl_fixbo_sld(:),                 & ! Element boundary fixity
       kfl_fixno_sld(:,:),               & ! Nodal fixity
       kfl_immer_sld(:,:),               & ! Nodal fixity for immersed boundary
       kfl_fixrs_sld(:),                 & ! Reference system when local axes is used
       kfl_funbo_sld(:),                 & ! Functions for boundary bc
       kfl_funno_sld(:),                 & ! Functions for node bc
       kfl_funtb_sld(:),                 & ! Functions type in BCs
       kfl_funtn_sld(:)                    ! Functions type in BCs

  real(rp)               :: &
       pressEndoIntegralNew_sld(3),         & ! Next step integral of the endocardiac pressure    
       pushForwardIntegralNew_sld,          & ! Next step Integral of the norm of the push forward of the exterior normal
       pressEndoIntegralOld_sld(3),         & ! Current step Integral of the endocardiac pressure    
       pushForwardIntegralOld_sld             ! Current step Integral of the norm of the push forward of the exterior normal
       
  !
  ! Stent
  !
  integer(ip), pointer                :: &
       kfl_contn_stent(:)                  ! Flag to set the nodes of the contact surface

  real(rp),     pointer               :: &
       bvess_sld(:,:,:),                 & ! Essential bc values (on dd/dt)
       bvnat_sld(:,:,:),                 & ! Natural bc values   (on traction)
       jacrot_du_dq_sld(:,:,:),          & ! Jacobian+Rotation matrix  U -> Q, defined on nbopo
       jacrot_dq_du_sld(:,:,:)             ! Transpose of the Jacobian+Rotation matrix  Q -> U, defined on nbopo

  integer(ip),  pointer               :: &
       lmate_sld(:),                     & ! Materials (element-wise)
       lcrkf_sld(:),                     & ! List of cracked faces
       resnode_sld(:),                   & ! List of nodes for measuring the residual force
       iswav_sld(:,:)                      ! Displacement and stress wave

  real(rp),     pointer                :: &
       calcium_sld(:,:),                  & ! Calcium concentration for solidz module
       fiber_sld(:,:),                    & ! Fiber corientation field
       fiemo_sld(:,:),                    & ! Modulator fields
       restr_sld(:,:),                    & ! Residual Internal stresses in voigt notation
       vofor_sld(:,:),                    & ! Volume external force in the reference frame
       disep_sld(:,:)                       ! Perturbed displacements
  !
  ! XFEM
  !
  real(rp),     pointer                 :: &
       dxfem_sld(:,:,:),                   & ! XFEM Displacement enrichment global function (alpha in denny's report)
       vxfem_sld(:,:,:),                   & ! XFEM Velocity enrichment global function (alpha in denny's report)
       axfem_sld(:,:,:),                   & ! XFEM Acceleration enrichment global function (alpha in denny's report)
       crapx_sld(:,:),                     & ! XFEM Crack/discontinuity position (per element)
       cranx_sld(:,:),                     & ! XFEM Crack/discontinuity normal (per element)
       sgmax_sld(:),                       & ! XFEM Maximum tensile stress (per element)
       skcos_sld(:,:,:),                   & ! Cosine matrices of skew systems
       cockf_sld(:,:),                     & ! List of cracked faces coordinates
       crtip_sld(:,:),                     & ! List of actual crack tip faces
       frxid_sld(:),                       & ! Global/assembled reaction force
       fexte_sld(:),                       & ! Global/assembled external force vector
       finte_sld(:),                       & ! Global/assembled internal force vector
       macce_sld(:),                       & ! Global/assembled Mass*accel force vector
       cohnx_sld(:,:)                        ! Normal to the midplane of cohesive element (per element)

  integer(ip)                           :: &
       nmate_sld,                          & ! # of materials
       kfl_goite_sld,                      & ! Keep iterating
       kfl_stead_sld,                      & ! Steady-state has been reached
       kfl_gdepo,                          & ! If GDEPO is needed
       nvgij_sld(6,2),                     & ! Voigt - indicial conversion table
       nvgij_inv_sld(3,3),                 & ! Voigt - inverse of the indicial conversion table
       ncomp_sld,                          & ! Number components temporal
       nprev_sld,                          & ! # previous time step or iteration
       ndofn_sld,                          & ! # of d.o.f. of the problem
       nvoig_sld,                          & ! voigt dimension
       nsvar_sld,                          & ! Number of state variables
       nzmat_sld,                          & ! Matrix size
       nzrhs_sld,                          & ! RHS size
       nderi_sld,                          & ! Number of ODES for RK4
       lasti_sld,                          & ! use to calculate the isovolumetric phases of the cardiac cycle
       ifase_sld(4),                       & ! flag for the cardiac cycle phases (detect the first time the phase is entered)
       kfase_sld,                          & ! Phase ID for the cardiac cycle
       iwave_sld,                          & ! Polarization wave
       numnodfor_sld                         ! Total number of nodes where the residual force is read

  real(rp), pointer                     :: &
       dunkn_sld(:),                       & ! corrector for Newton--Raphson
       ddisp_sld(:,:,:),                   & ! Displacement increment
       accel_sld(:,:,:),                   & ! Acceleration
       veloc_sld(:,:,:),                   & ! Velocity
       fextt_sld(:,:,:),                   & ! External forces (Temporal)
       fintt_sld(:,:,:),                   & ! Internal forces (Temporal)
       allie_sld(:),                       & ! All internal energy
       allwk_sld(:),                       & ! All external work
       allke_sld(:),                       & ! All kinetic energy
       etota_sld(:),                       & ! Total energy
       ddism_sld(:,:),                     & ! FSI movemente variation (ALEFOR COUPLING)
       gpsl0_sld(:,:,:),                   & ! Initial Logarithmic strain (BRIDGE)
       rstr0_sld(:,:,:),                   & ! Rotation matrix (F=RU) (BRIDGE)
       dfric_sld(:,:,:,:),                 & ! Friction force for contact/friction
       dslip_sld(:,:,:,:),                 & ! Crack tangential jump/surface slip
       dceff_sld(:,:,:,:),                 & ! Effective scalar opening for cohesive law
       dcmax_sld(:,:,:),                   & ! Maximum scalar opening for cohesive law
       treff_sld(:),                       & ! Effective traction for cohesive law
       unknotmp_sld(:),                    &
       veloctmp_sld(:,:)
  !
  ! Rigid body variables
  !
  integer(ip)                           :: &
       kfl_rbfix_sld(3)
  real(rp)                              :: &
       rbfor_sld(3),                       & ! Rigid body force
       rbmas_sld                             ! Rigid body mass

  real(rp), pointer                     :: &
       crkpo_sld(:,:),                     & ! Crack positions
       crkno_sld(:,:),                     & ! Crack normals
       srpro_sld(:,:),                     &
       celen_sld(:)                          ! Characteristic element length

  integer(ip), pointer                  :: &
       lnenr_sld(:),                       & ! List of enriched nodes
       lecoh_sld(:),                       & ! List of cohesive elements
       leenr_sld(:),                       & ! List of enriched elements
       ledam_sld(:)                          ! List of damaged elements

  real(rp),     pointer               :: &
       fibts_sld(:,:),                   & ! Orthotropic fiber direction s
       fibtn_sld(:,:),                   & ! Orthotropic fiber direction n
       nopio_sld(:,:),                   & ! First Piola--Kirchhoff Stress tensor on nodes
       caust_sld(:,:),                   & ! Cauchy Stress tensor
       roloc_sld(:,:,:),                 & ! Local rotation tensor
       sigei_sld(:),                     & ! Largest eigenvector of the Cauchy Stress tensor (sigma)
       lepsi_sld(:,:),                   & ! Log strain in a local system define in the const. law (ex:fiber-sheet-normal)
       fibde_sld(:,:),                   & ! Fiber postproc vector in the current configuration (relative to ini. lenght)
       seqvm_sld(:),                     & ! Combined stress: Von Mises
       caunn_sld(:),                     & ! n . Cauchy . n: different than zero in the surface
       green_sld(:,:),                   & ! Green strain tensor
       eprin_sld(:),                     & ! Principal strains at ipoin
       dttau_sld(:),                     & ! Local pseudo-time step (tau)
       vmass_sld(:),                     & ! Weighted mass matrix
       vmasx_sld(:),                     & ! Enriched mass matrix
       grlst_sld(:,:),                   & ! green lagrangian strain tensor per element
       epsee_sld(:,:),                   & ! elastic strains in voig notain at ipoin
       ebfil_sld(:),                     & ! Integral of strain projection onto biofiber longitudinal
       sbfil_sld(:)                        ! Integral of stress projection onto biofiber longitudinal


  type(r3p),   pointer                :: &
       cause_sld(:)                        ! Element-wise Cauchy Stress tensor

  real(rp)                              :: &
       dtinv_sld,                          & ! 1/dt
       dtcri_sld,                          & ! Critical time step
       dtmin_sld,                          & ! Min. critical time step set by the user
       sgult_sld,                          & ! Ultimate maximum level of tensile stress
       pabdf_sld(10),                      & ! BDF paremeters
       resid_sld,                          & ! Residual outer iterations (block-coupling)
       volvr_sld(2),                       & ! Reference ventricle volumes
       ejrat_sld(2),                       & ! Ejection rate
       dltap_sld(12,2),                    & ! dP for calculation of pression in the heart
       timst_sld(12),                      & ! time for calculation of pression in the heart
       sidev_sld(2,2),                     & ! ...for cardiac isovolumetric phase
       forcm_sld,                          & ! Max force for EC coupling
       velov_sld,                          & ! rate of change of volume - for cardiac isovolumetric phase
       bodyf_sld(3),                       & ! Body force vector
       ptota_sld(4),                       & ! Pressure from the windkessel bc (cardiac)
       tactv_sld,                          & ! Volumetric active force along the fibers
       fint2_sld,                          & ! L2-norm assembled/global internal force vector
       fext2_sld,                          & ! L2-norm assembled/global external force vector
       fine2_sld,                          & ! L2-norm assembled/global M*a force vector
       fnatu_sld,                          & ! L2-norm of the natural force imposed in the solver
       eener_sld,                          & ! Error in energy
       volum_sld(2),                       & ! Reference (1) and deformed configuration (2) volumes
       cpu_ass_sol_sld(4)                    ! CPU time assembly and solver at each iteration

  real(rp),  pointer :: &
       du_sld(:,:)

  ! Global vectors defined per element and per gauss point
  type(r3p),   pointer                  :: &
       epsel_sld(:),                       &
       lepse_sld(:),                       &
       gppio_sld(:),                       & ! 1st PK stress tensor on elements
       gpgdi_sld(:),                       & ! Deformation gradient tensor on elements
       fibeg_sld(:)                          ! Fiber orientation field

  type(r1p),   pointer                  :: &
       dedef_sld(:),                       & ! Determinant of the deformation gradient tensor grdef
       enden_sld(:),                       & ! Helmholz free energy, currently only for sm400, could be expanded
       water_sld(:)                          ! Water content of element for sm400

  !
  ! Exact solution
  !
  real(rp)                              :: &
       err01_sld(2),                       &
       err02_sld(2),                       &
       err0i_sld(2),                       &
       err11_sld(2),                       &
       err12_sld(2),                       &
       err1i_sld(2)

  !
  ! Aitken relaxation factors, computed by the master and passed to the slaves
  !
  real(rp)                              :: &
       relpa_sld                             ! Solver
  !
  ! PDN Contact stuff
  !
  real(rp)                              :: &
       neumann_relax,                      & ! Relaxation for Neumann part (input)
       contact_friction_sld,               & ! Contact friction coefficient, no frictions is default
       coupling_contact_tol,               & ! Coupling tolerance
       relax_bilateral,                    & ! Relaxation (fixed or computed)
       vect_proje_sld(3),                  & ! Projection vector
       petol_sld                             ! Penetration tolerance for PDN contact

  integer(ip)                           :: &
       kfl_conta_sld,                      & ! PDN Contact flag
       kfl_contf_sld,                      & ! PDN Contact convergence file
       kfl_minco_sld = 0_ip,               & ! Was here. Has nothing to do with PDN contact
       contactbou_sld,                     & ! Identificator of contact boundary. Input keyword: CONTACT_BOUNDARY
       conbo_sld(2),                       & ! Identificator of outer or inner surface of contact in stent
       new_to_release,                     & ! PDN contact. Flag. If there are new nodes to be released
       kfl_mrele_sld,                      & ! Method release nodes
       coupling_contact_its                  ! Number of coupling iterations for contact

  integer(ip), pointer                  :: &
       release_nodes(:)                      ! PDN contact. Marked nodes to be released

  real(rp),  pointer                    :: &
       fcont_sld(:,:),                     & ! Contact force
       rhsid1st_it(:),                     & ! PDN contact. Save rhsid for each first iteration of the NR.
       saved_rhsid(:,:),                   & ! ndime x nre:q
       reaction_ant(:,:,:)                   ! Bi-lateral contact. Save reaction from previous iteration.
  !
  ! Composite stuff ( *AQU* )
  !
  real(rp), pointer                     :: &
       stiff0_sld(:, :, :),                & ! Stiffness matrix without damage
       axis1_sld(:,:),                     & ! Local CSYS: material orientation (axis-1)
       axis2_sld(:,:),                     & ! Local CSYS: material orientation (axis-2)
       axis3_sld(:,:),                     & ! Local CSYS: material orientation (axis-3)
       exm_lambda_sld(:,:,:),              & ! Stretch field (for Land)
       orien_sld(:)                          ! Orientation angle for material CSYS

  type(r3p), pointer                    :: &
       caunp_sld(:),                       & ! Recovered stresses at nodes, element-wise Cauchy Stress tensor
       caulo_sld(:)                          ! Recovered stresses at nodes, element-wise Cauchy Stress tensor in local csys

  type(r3p), pointer                    :: &
       rmate_sld(:),                       & ! Rotation matrix for material orientations
       svegm_sld(:)                          ! Internal state variables per element at gauss point (*AQU*)
 
  
end module def_solidz
