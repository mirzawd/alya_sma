!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_magnet
  !-----------------------------------------------------------------------
  !****f* magnet/def_magnet
  ! NAME
  !    def_magnet
  ! DESCRIPTION
  !    
  ! USES
  ! USED BY
  !
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame, only: pi
  use def_parall, only: nproc_par, iproc_par
  use def_master, only: zeror
  use mod_mag_quadra
  use mod_mag_bdfode
  use mod_mag_interp
  
  implicit none

  integer :: &
    ioun1_mag = 28, &
    ioun2_mag = 29, &
    ioun3_mag = 30, &
    ioun4_mag = 31, &
    ioun5_mag = 32, &
    ioun6_mag = 33, &
    ioun7_mag = 34, &
    ioun8_mag = 35, &
    ioun9_mag = 36

  !--------------------------------------------------------------
  ! Structures Definition
  !--------------------------------------------------------------
  type stcntr
    real(rp), dimension(:), pointer ::    &
      cnstr,    &
      unkno
    real(rp) ::    &
      valuc,    & ! External Constrain Value
      valup,    &
      valud,    &
      valui       ! Internal Constrain Value
    real(rp) ::    &
      multc = 0.0_rp,    & ! Multiplier Value
      multp = 0.0_rp,    &
      multd = 0.0_rp,    &
      rmult = 0.0_rp       ! Multiplier Residual
    real(rp) ::    &
      normc,    &
      normu
    integer(ip) ::    &
      lsite
  end type stcntr
  !--------------------------------------------------------------

  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in mag_reabcs
  !
  !------------------------------------------------------------------------
 
  type(bc_nodes), pointer :: &
    tecod_mag(:)    !  code type  
  type(bc_bound), pointer :: &
    tbcod_mag(:)    ! Boundary code type

  !----------------------------------------------
  ! Geometric data
  !----------------------------------------------
  real(rp), pointer :: &
    edglen_mag(:), &
    elevol_mag(:), &
    elesig_mag(:,:), &
    edgsig_mag(:), &
    locboucen_mag(:,:), &
    globoucen_mag(:,:), &
    locboutan_mag(:,:), &
    globoutan_mag(:,:)
  
  integer(ip), pointer :: &
    eleedg_mag(:,:),    &
    diredg_mag(:),    &
    edgnod_mag(:,:),    &
    bouedg_mag(:)

  integer(ip) :: &
    nedgdir_mag

  logical(lg), pointer :: &
    edgdir_mag(:), &
    edgbou_mag(:)

  logical(lg) :: &
    struct_mag = .false. ! Mesh made of quadrangles/hexaedral is structured/unstructured
  !----------------------------------------------


  !----------------------------------------------
  ! More Geometric Data
  !----------------------------------------------
  integer(ip), allocatable :: &
    node_to_nedge_mag(:),    &
    r_node_to_edge_mag(:),    &
    c_node_to_edge_mag(:),    &
    i_node_to_edge_mag(:),    &
    edge_to_nedge_mag(:),    &
    r_edge_to_edge_mag(:),    &
    i_edge_to_edge_mag(:),    &
    c_edge_to_edge_mag(:)

  integer(ip) :: &
    ne2e_mag
  !----------------------------------------------


  !----------------------------------------------
  ! Field variable & Auxiliary Vectors
  !----------------------------------------------
  real(rp), pointer :: &
    He_mag(:)

  real(rp), pointer :: &
    Ha_mag(:),    &
    Hp_mag(:),    &
    dH_mag(:),    &
    bhsid(:),    &
    fhsid(:),    &
    hhsid(:)
  !----------------------------------------------


  !----------------------------------------------
  ! Axisymmetric Analysis
  !----------------------------------------------
  logical(lg) :: &
    kfl_axsym_mag = .false.
  !----------------------------------------------


  !----------------------------------------------
  ! Post-processor
  !----------------------------------------------
  real(rp), pointer :: &
    Hc_mag(:,:),    &
    Hn_mag(:,:),    &
    Jcz_mag(:),    &
    Jnz_mag(:),    &
    Bc_mag(:,:),    &
    Bn_mag(:,:),    &
    Jc_mag(:,:),    &
    Jn_mag(:,:),    &
    Fn_mag(:,:),    &
    Fc_mag(:,:),    &
    Gc_mag(:),    &
    Gn_mag(:)

  type(r3p), pointer :: &
    Hgp_mag(:),    &
    Jgp_mag(:),    &
    Bgp_mag(:),    &
    Fgp_mag(:),    &
    Ggp_mag(:)

  logical :: &
    postev_mag = .false.,    &
    postce_mag = .false.
  !----------------------------------------------


  !----------------------------------------------
  ! Numerical Parameters
  ! Input through mag_reanut
  !----------------------------------------------
  !
  ! Time Integration
  !
  real(rp) :: &
    dtmax_mag, &
    dtmin_mag, &
    theta_mag
  !
  ! Newton Iteration
  !
  real(rp) :: &
    nltol_mag,  &
    abstol_mag,  &
    reltol_mag

  integer(ip) :: &
    nlide_mag, &
    nlite_mag 
  !
  ! Quadrature
  !
  integer(ip) :: &
    gslin_mag, &
    gstri_mag, &
    gsqua_mag, &
    gstet_mag, &
    gshex_mag
  !----------------------------------------------


  !----------------------------------------------
  ! Quadrature
  !----------------------------------------------
  type (stcqua) ::    &
    quadLin,    &
    quadTri,    &
    quadQua,    &
    quadTet,    &
    quadHex
  !----------------------------------------------


  !----------------------------------------------
  ! BDF
  !----------------------------------------------
  type (stcbdf) :: &
    bdfode_mag

  real(rp) :: &
    a0_mag

  real(rp), pointer :: &
    Hpav_mag(:)

  logical(lg) :: &
    kfl_reset_mag = .false.
  !----------------------------------------------


  !----------------------------------------------
  ! Contraints
  !----------------------------------------------
  logical(lg) :: &
    kfl_lagr_mag = .false. 

  type (stcntr), allocatable :: &
    constr_mag(:)

  integer(ip) :: &
    constr_total = 0_ip

  integer(ip), allocatable :: &
    constrlist_mag(:)

  real(rp), allocatable :: &
    cntmat_mag(:,:),    &
    cntrhs_mag(:),    &
    cntunk_mag(:)
  !----------------------------------------------


  !----------------------------------------------
  ! Interpolation Files
  !----------------------------------------------
  type (stcntp), allocatable :: &
    intp1_mag(:)

  integer(ip) :: &
    nintp1_mag = 0_ip
  !----------------------------------------------


  !----------------------------------------------
  ! Self field
  !----------------------------------------------
  real(rp), pointer :: &
    Hxtl_mag(:), &
    Hsfl_mag(:), &	! total self field in local boundary
    Hsfg_mag(:), &    	! partial self field in global boundary
    Hbs_mag(:,:)

  logical(lg) :: &
    kfl_self_mag = .false.,    &
    kfl_prev_mag = .false.

  logical(lg), allocatable :: selfList_mag(:)

  integer(4), pointer :: &
    proc_nbedg(:)
  !----------------------------------------------


  !----------------------------------------------
  ! Arrays Allocation
  !----------------------------------------------
  integer(ip) :: &
    maxdof_mag, &
    maxgau_mag, &
    maxmat_mag, &
    mxndof_mag
  !----------------------------------------------

  integer(ip) :: &
    lsiter0_mag

  real(rp) :: &
    lsresi1_mag, &
    lsresi2_mag

  !---------------------------------------------------------
  ! Material properties
  ! Input through mag_reaphy
  !---------------------------------------------------------  
  integer(ip), allocatable, dimension(:) :: &
    resistOpt_mag,    &
    scalinOpt_mag

  real(rp), allocatable, dimension(:) :: &
    nc0_mag,    &
    Ec0_mag,    &
    Jc0_mag,    &
    mur_mag,    &
    rho_mag,    &
    B0k_mag

  ! Number of tensor components
  ! Choose 3 to have access to XX, YY, ZZ
  ! Choose 6 to have access to XX, YY, ZZ, XY, XZ, YZ 
  integer(ip), parameter :: &
    ncomp_mag = 6_ip

  ! Dimensions ncomp_mag x maxmat_mag
  integer(ip), allocatable, dimension(:,:) :: &
    resmat_mag,    &    ! Resistivity function
    Jcrmat_mag          ! Scaling Law 

  ! Dimensions ncomp_mag x maxmat_mag
  real(rp), allocatable, dimension(:,:) :: &
    ncrmat_mag,    &
    Ecrmat_mag,    &
    Jc0mat_mag,    &
    murmat_mag,    &
    rhomat_mag,    &
    Bc0mat_mag,    &
    Tc0mat_mag

  logical(lg), allocatable, dimension(:) :: &
    kfl_resiso_mag,    &
    kfl_ncriso_mag,    &
    kfl_Ecriso_mag,    &
    kfl_Jcriso_mag,    &
    kfl_muriso_mag,    &
    kfl_rhoiso_mag,    &
    kfl_Bc0iso_mag,    &
    kfl_Tc0iso_mag,    &
    kfl_Jc0iso_mag

  ! Origin for magnetic moment calculation
  real(rp), allocatable, dimension(:,:) :: &
    momori_mag(:,:)
 
  real(rp), parameter :: &
    mu0_mag = pi * 4.0e-7_rp
  !---------------------------------------------------------


  !---------------------------------------------------------
  ! Norms
  !---------------------------------------------------------
  real(rp) :: &
    absRes_mag,    &
    absRes0_mag,    &
    relRes_mag,     &
    solVar_mag,     &
    Hnorm_mag,    &
    gloRes_mag,   &
    dHnorm_mag

  real(rp), allocatable :: &
    refres_mag(:),    &
    residu_mag(:),    &
    refval_mag(:),    &
    delval_mag(:)

  real(rp), allocatable, dimension(:) :: &
    magnen_mag,    &
    joulen_mag
  !---------------------------------------------------------

  logical(lg) :: &
    kfl_nrj_mag = .false.,    &
    kfl_mtz_mag = .false.,    &
    kfl_crn_mag = .false.,    &
    kfl_vlm_mag = .false.

  real(rp), allocatable, dimension(:,:) :: &
    magtiz_mag,    &
    cursum_mag,    &
    magsum_mag

  real (rp), allocatable, dimension(:) :: &
    volume_mag
  
  !---------------------------------------------------------
  ! Non-linear loop
  !---------------------------------------------------------
  logical(lg) :: &
    kfl_goiter_mag

  integer(ip) :: &
    nonlCount_mag
  !---------------------------------------------------------


  !---------------------------------------------------------
  ! Time loop
  !---------------------------------------------------------
  real(rp) :: &
    dt_mag, &
    dtp_mag

  real(rp) :: &
    myClockTime_mag

  integer(ip) :: &
    timeStep_mag
  !---------------------------------------------------------


  !--------------------------------------------------------------
  integer(ip) :: &
    kfl_edge_mag

  !--------------------------------------------------------------
  !
  ! Boundary conditions data structures: read in 'mag_reabcs' file
  !
  integer(ip), pointer ::  &
    kfl_fixno_mag(:,:),    & ! Nodal fixity
    kfl_fixbo_mag(:)         ! Element boundary fixity
  
  real(rp), pointer ::   & 
    bvess_mag(:,:,:),    & ! Essential BC values
    bvnat_mag(:,:,:)       ! Neumann BC values
  !--------------------------------------------------------------
  !
  ! Commonly used data structures
  !
  integer(ip), pointer :: &
    edgind_mag(:),    &
    nodind_mag(:)

  real(rp), allocatable :: &
    elstif(:,:),    &
    elstdf(:,:),    &
    elmass(:,:),    &
    elsrc(:),    &
    elheat(:),    &
    Aelem(:,:),    &
    Anelem(:,:),    &
    belem(:),    &
    rnelem(:),    &
    felem(:),    &
    helem(:),    &
!    nodi(:,:),    &
!    nodq(:,:),    &
!    bq(:,:),    &
    sigval_mag(:),    &
    lenval_mag(:)

  real(rp) :: &
    INF_MAG = huge(1.0_rp),    &
    ZER_MAG = tiny(1.0_rp),    &
    SCL_MAG = 10.0_rp

end module def_magnet

