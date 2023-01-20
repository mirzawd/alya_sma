!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup ADR_Toolbox
!> ToolBox for assembling advection-diffusion-reaction equation
!> @{
!> @name    ToolBox for advection-diffusion-reaction equation
!> @file    mod_ADR.f90
!> @date    21/09/2015
!> @author  Guillaume Houzeaux
!> @brief   ADR equation
!> @details ADR equation
!
!-----------------------------------------------------------------------

module mod_ADR

  use def_kintyp,              only : ip,rp,lg,r3p
  use def_domain,              only : mesh_type
  use def_domain,              only : elm,lezdo
  use def_domain,              only : ngaus,lorde,lezdo
  use def_domain,              only : llapl
  use def_master,              only : INOTMASTER
  use def_master,              only : ittim
  use def_master,              only : INOTSLAVE
  use def_kermod,              only : kfl_element_to_csr
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_communications,      only : PAR_MAX
  use mod_communications,      only : PAR_SUM
  use mod_elmgeo,              only : element_type
  use mod_elmgeo,              only : elmgeo_jacobian_matrix
  use mod_elmgeo,              only : elmgeo_cartesian_derivatives
  use mod_elmgeo,              only : elmgeo_element_characteristic_length
  use mod_elmgeo,              only : elmgeo_element_length
  use mod_matrix,              only : matrix_assemble_element_matrix_to_CSR,matrix_assemble_element_RHS
  use def_inpout,              only : words
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use mod_exchange,            only : exchange_add
  implicit none 

  private
  
  real(rp),    parameter :: zeror         =  epsilon(1.0_rp)
  real(rp),    parameter :: pi            =  3.141592653589793238462643383279502884197_rp
  !
  ! Maximum dimensions
  !
  integer(ip), parameter :: mreac_adr     =  4 ! Maximum number of reaction terms
  !
  ! On/Off 
  !
  integer(ip), parameter :: ON            =  1
  integer(ip), parameter :: OFF           =  0
  !
  ! Time scheme
  !
  integer(ip), parameter :: TRAPEZOIDAL    =  1
  integer(ip), parameter :: BDF            =  2
  integer(ip), parameter :: ADAMS_BASHFORD =  3
  integer(ip), parameter :: RUNGE_KUTTA    =  4
  !
  ! Linearization
  !
  integer(ip), parameter :: RHS           = -1
  integer(ip), parameter :: PICARD        =  0
  integer(ip), parameter :: NEWTON        =  1
  !
  ! Time strategy
  !
  integer(ip), parameter :: PRESCRIBED    =  0
  integer(ip), parameter :: FROM_CRITICAL =  1
  integer(ip), parameter :: LOCAL         =  2
  !
  ! Stabilization strategy
  !
  integer(ip), parameter :: GALERKIN      = -2
  integer(ip), parameter :: SU            = -1
  integer(ip), parameter :: SUPG          = -3
  integer(ip), parameter :: BUBBLE        = -4
  integer(ip), parameter :: ASGS          =  0
  integer(ip), parameter :: FULL_OSS      =  1 ! Full residual is projected
  integer(ip), parameter :: A_OSS         =  2 ! Only advection is projected
  integer(ip), parameter :: AR_OSS        =  3 ! Advection & Reaction are projected
  !
  ! Stabilization parameter
  !
  integer(ip), parameter :: TAU_OFF       =  0
  integer(ip), parameter :: TAU_CODINA    =  1
  integer(ip), parameter :: TAU_EXACT_1D  =  2
  integer(ip), parameter :: TAU_SHAKIB    =  3
  integer(ip), parameter :: TAU_TEST      =  4
  integer(ip), parameter :: TAU_INCLU_DT  =  5
  integer(ip), parameter :: TAU_DT        =  6
  !
  ! What to do
  !
  integer(ip), parameter :: ELEMENT_ASSEMBLY             = 1
  integer(ip), parameter :: PROJECTIONS_AND_SGS_ASSEMBLY = 4
  integer(ip), parameter :: BUBBLE_ASSEMBLY              = 5
  !
  ! Mesh 
  !
  integer(ip)            :: ndime
  integer(ip)            :: ntens
  integer(ip)            :: npoin
  integer(ip)            :: nelem
  integer(ip)            :: mnode
  integer(ip)            :: mgaus
  integer(ip), pointer   :: lnods(:,:)
  integer(ip), pointer   :: ltype(:)
  integer(ip), pointer   :: lnnod(:)
  real(rp),    pointer   :: coord(:,:)
  type(elm),   pointer   :: elmar_loc(:)

  !----------------------------------------------------------------------
  !
  ! ADR type
  !
  !----------------------------------------------------------------------

  type ADR_typ
     !
     ! Read variables
     !
     integer(ip)        :: kfl_time_integration    ! Time flag              
     integer(ip)        :: kfl_time_step_strategy  ! Time strategy
     integer(ip)        :: kfl_stabilization       ! Stabilization strategy
     integer(ip)        :: kfl_shock               ! Shock capturing flag            
     integer(ip)        :: kfl_time_lumped         ! Lumped time evolution
     integer(ip)        :: kfl_linearization       ! Linearization of reaction
     integer(ip)        :: kfl_tau_strategy        ! Tau strategy         
     integer(ip)        :: kfl_laplacian           ! If Laplacian should be computed
     integer(ip)        :: kfl_time_sgs            ! SGS tracking in time
     integer(ip)        :: kfl_nonlinear_sgs       ! SGS non-linear tracking
     integer(ip)        :: kfl_time_bubble         ! If bubble is tracked in time
     integer(ip)        :: kfl_time_scheme         ! Time scheme
     integer(ip)        :: kfl_time_order          ! Time scheme order
     integer(ip)        :: kfl_manufactured        ! Manufactured solution
     integer(ip)        :: kfl_length              ! Characteristic length
     integer(ip)        :: kfl_first_order_sgs     ! If SGS should be computed at first order
     integer(ip)        :: kfl_discretization      ! Discretization method (FE=0,FV=1)
     integer(ip)        :: kfl_skewsymm            ! Skew-symmetric convective term
     integer(ip)        :: number_euler_steps      ! Number of Euler time steps

     real(rp)           :: bemol                   ! Integration by parts convective term
     real(rp)           :: tau_parameters(3)       ! Stability constants             
     real(rp)           :: shock                   ! Shock capturing coefficients    
     real(rp)           :: relax_sgs               ! Relaxation of SGS
     !
     ! Memory counter      
     !
     integer(4)         :: lun_output4             ! Output file unit 
     integer(8)         :: memor(2)                ! Memory counter
     !
     ! Arrays
     ! 
     real(rp),  pointer :: proje1(:)               ! First projection 
     real(rp),  pointer :: proje2(:)               ! Second projection  
     real(rp),  pointer :: proje1_tmp(:)           ! temporary first projection 
     real(rp),  pointer :: proje2_tmp(:)           ! Temporary second projection  
     type(r3p), pointer :: sgs(:)                  ! Subgrid scale
     real(rp),  pointer :: bubble(:,:)             ! Bubble
     !
     ! Derived variables
     !
     integer(ip)        :: number_components       ! Number of components            
     integer(ip)        :: ntime                   ! Number of time values
     integer(ip)        :: kfl_time_order_save     ! Time scheme order saved value
     integer(ip)        :: nunkn                   ! Number of unknowns
     real(rp)           :: time_parameters(10)     ! BDF coefficients
     real(rp)           :: dtinv                   ! 1/dt of current time step n+1
     real(rp)           :: dtinv_old(10)           ! 1/dt of old time steps n,n-1,n-2, etc.
  end type ADR_typ

  !----------------------------------------------------------------------
  !
  ! ADR type
  !
  !----------------------------------------------------------------------

  !type ADR_elm_mgaus_type
  !   real(rp)    :: gprea(mreac_adr)               ! r ................. reaction
  !   real(rp)    :: gpvel(mdime_adr)               ! a ................. advection
  !   real(rp)    :: gpcod(mdime_adr)               ! x ................. position
  !   real(rp)    :: gpdif                          ! k+kt .............. diffusion
  !   real(rp)    :: gpgrd(mdime_adr)               ! grad(k+kt) ........ diffusion gradient
  !   real(rp)    :: gpdiv                          ! div(a) ............ divergence of convection
  !   real(rp)    :: gprhs                          ! f ................. all terms
  !   real(rp)    :: gpden                          ! rho ............... density
  !   real(rp)    :: gpunk(mcomp_adr)               ! u, u^n-1, etc ..... previous time step unknowns
  !   real(rp)    :: grunk(mdime_adr)               ! grad(u) ........... unknown gradients
  !   real(rp)    :: gpsgs(mcomp_adr)               ! u', u'^n-1, etc ... previous time step SGS
  !end type ADR_elm_mgaus_type

  interface ADR_begin_time_step
     module procedure ADR_begin_time_step_s,&
          &           ADR_begin_time_step_1
  end interface ADR_begin_time_step
  interface ADR_begin_inner_iteration
     module procedure ADR_begin_inner_iteration_s,&
          &           ADR_begin_inner_iteration_1
  end interface ADR_begin_inner_iteration
  interface ADR_end_inner_iteration
     module procedure ADR_end_inner_iteration_s,&
          &           ADR_end_inner_iteration_1
  end interface ADR_end_inner_iteration
  interface ADR_end_time_step
     module procedure ADR_end_time_step_s,&
          &           ADR_end_time_step_1
  end interface ADR_end_time_step
  interface ADR_after_restart
     module procedure ADR_after_restart_s,&
          &           ADR_after_restart_1
  end interface ADR_after_restart 

  public :: ADR_typ
  public :: ADR_element_assembly
  public :: ADR_projections_and_sgs_assembly
  public :: ADR_bubble_assembly
  public :: ADR_update_projections
  public :: ADR_load_mesh
  public :: ADR_initialize_projections
  public :: ADR_manufactured_error
  public :: ADR_add_sgs_or_bubble
  public :: ADR_check_and_compute_data
  public :: ADR_manufactured_nodal_error
  public :: ADR_allocate_projections_bubble_sgs
  public :: ADR_end_time_step
  public :: ADR_begin_inner_iteration
  public :: ADR_end_inner_iteration
  public :: ADR_begin_time_step
  public :: ADR_critical_time_step
  public :: ADR_manufactured_nodal_solution
  public :: ADR_read_data
  public :: ADR_parallel_data
  public :: ADR_arrays
  public :: ADR_time_strategy
  public :: ADR_initialize_type
  public :: ADR_after_restart
  public :: ADR_assemble_laplacian
  public :: ADR_assemble_convective
  public :: ADR_assemble_extension
  
  public :: BDF
  public :: TRAPEZOIDAL
  public :: ADAMS_BASHFORD
  public :: RUNGE_KUTTA
  public :: FROM_CRITICAL
  public :: mreac_adr
  public :: BUBBLE_ASSEMBLY
  public :: PROJECTIONS_AND_SGS_ASSEMBLY
  public :: ELEMENT_ASSEMBLY
  public :: BUBBLE
  public :: AR_OSS
  public :: A_OSS
  public :: FULL_OSS

contains

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Read data
  !> @details Read data
  !
  !-----------------------------------------------------------------------

  subroutine ADR_parallel_data(ADR)
    type(ADR_typ), intent(inout) :: ADR !< ADR type
    integer(ip)                  :: ii
 
    call exchange_add(ADR % kfl_time_integration   )
    call exchange_add(ADR % kfl_time_step_strategy )
    call exchange_add(ADR % kfl_stabilization      )
    call exchange_add(ADR % kfl_shock              )
    call exchange_add(ADR % kfl_time_lumped        )
    call exchange_add(ADR % kfl_linearization      )
    call exchange_add(ADR % kfl_tau_strategy       )
    call exchange_add(ADR % kfl_laplacian          )
    call exchange_add(ADR % kfl_time_sgs           )
    call exchange_add(ADR % kfl_nonlinear_sgs      )
    call exchange_add(ADR % kfl_time_bubble        )
    call exchange_add(ADR % kfl_time_scheme        )
    call exchange_add(ADR % kfl_time_order         )
    call exchange_add(ADR % kfl_manufactured       )
    call exchange_add(ADR % kfl_length             )
    call exchange_add(ADR % kfl_first_order_sgs    )
    call exchange_add(ADR % kfl_discretization     )
    call exchange_add(ADR % kfl_skewsymm           )
    call exchange_add(ADR % number_euler_steps     )

    call exchange_add(ADR % bemol                  ) 
    do ii = 1,size(ADR % tau_parameters,KIND=ip)
       call exchange_add(ADR % tau_parameters(ii)  )
    end do
    call exchange_add(ADR % shock                  ) 
    call exchange_add(ADR % relax_sgs              ) 

    call exchange_add(ADR % ntime) 
    do ii = 1,size(ADR % time_parameters,KIND=ip)
       call exchange_add(ADR % time_parameters     )    
    end do
    call exchange_add(ADR % number_components      ) 
    call exchange_add(ADR % dtinv                  ) 
    do ii = 1,size(ADR % dtinv_old,KIND=ip)
       call exchange_add(ADR % dtinv_old(ii)       ) 
    end do

  end subroutine ADR_parallel_data

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Read data
  !> @details Read data
  !
  !-----------------------------------------------------------------------

  subroutine ADR_read_data(ADR)
    type(ADR_typ), intent(inout) :: ADR !< ADR type

    if( words(1) == 'STABI' ) then
       !
       ! Stabilization strategy
       !
       if(      words(2) == 'GALER' ) then
          ADR % kfl_stabilization = GALERKIN
       else if( words(2) == 'SU   ' ) then
          ADR % kfl_stabilization = SU
       else if( words(2) == 'SUPG ' ) then
          ADR % kfl_stabilization = SUPG
       else if(words(2) == 'BUBBL' ) then
          ADR % kfl_stabilization = BUBBLE    
       else if( words(2) == 'ASGS ' ) then
          ADR % kfl_stabilization = ASGS
       else if(words(2) == 'FULOS' ) then
          ADR % kfl_stabilization = FULL_OSS    
       else if(words(2) == 'AOSS ' ) then
          ADR % kfl_stabilization = A_OSS    
       else if(words(2) == 'AROSS' ) then
          ADR % kfl_stabilization = AR_OSS    
       else
          call runend('ADR_READ_DATA: UNKNOWN STABIIZATION STRATEGY')
       end if

    else if( words(1) == 'TAUST' ) then
       !
       ! Tau strategy
       !
       if(      words(2) == 'OFF  ' ) then
          ADR % kfl_tau_strategy = TAU_OFF
       else if( words(2) == 'CODIN' ) then
          ADR % kfl_tau_strategy = TAU_CODINA
       else if( words(2) == 'EXACT' ) then
          ADR % kfl_tau_strategy = TAU_EXACT_1D
       else if( words(2) == 'SHAKI' ) then
          ADR % kfl_tau_strategy = TAU_SHAKIB
       else if( words(2) == 'TEST ' ) then
          ADR % kfl_tau_strategy = TAU_TEST
       else if( words(2) == 'INCLU' ) then
          ADR % kfl_tau_strategy = TAU_INCLU_DT
       else if( words(2) == 'TIMES' ) then
          ADR % kfl_tau_strategy = TAU_DT
       else
          call runend('ADR_READ_DATA: UNKNOWN TAU STRATEGY')
       end if

    else if( words(1) == 'TIMEI' .or. words(1) == 'TEMPO' ) then
       !
       ! Time integration
       !
       if(      words(2) == 'OFF  ' ) then
          ADR % kfl_time_integration = OFF
       else if( words(2) == 'ON   ' ) then
          ADR % kfl_time_integration = ON
       else
          call runend('ADR_READ_DATA: UNKNOWN TIME INTEGRATION OPTION')
       end if
      
    end if

  end subroutine ADR_read_data

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Load mesh
  !> @details Load mesh arrays and parameters
  !
  !-----------------------------------------------------------------------

  subroutine ADR_assemble_laplacian(meshe_in,elmar_in,amatr,rhs_value,rhsid,kfl_fixno,bvess)
    
    type(mesh_type), intent(in)            :: meshe_in
    type(elm),       intent(in)            :: elmar_in(*)
    real(rp),        intent(out)           :: amatr(*)
    real(rp),        intent(in),  optional :: rhs_value(*)
    real(rp),        intent(out), optional :: rhsid(*)
    integer(ip),     intent(in),  optional :: kfl_fixno(:,:)
    real(rp),        intent(in),  optional :: bvess(:,:)
    integer(ip)                            :: ndime,pnode,pgaus,inode
    integer(ip)                            :: mnode,pelty,mgaus,ielem
    integer(ip)                            :: igaus,ipoin,izdom
    real(rp),        allocatable           :: elcod(:,:)
    real(rp),        allocatable           :: gpvol(:)
    real(rp),        allocatable           :: gpsha(:,:)
    real(rp),        allocatable           :: gpcar(:,:,:)
    
    if( INOTMASTER ) THEN
       !
       ! Initialization
       !
       if( present(rhsid) ) then
          do ipoin = 1,meshe_in % npoin
             rhsid(ipoin) = 0.0_rp
          end do
       end if
       do izdom = 1,meshe_in % r_dom(meshe_in % npoin+1)-1
          amatr(izdom) = 0.0_rp
       end do

       ndime = meshe_in % ndime
       mnode = meshe_in % mnode
       mgaus = maxval(ngaus)

       allocate( elcod(ndime,mnode) )
       allocate( gpvol(mgaus) )
       allocate( gpsha(mnode,mgaus) )
       allocate( gpcar(ndime,mnode,mgaus) )

       do ielem = 1,meshe_in % nelem
          pnode = meshe_in % lnnod(ielem)
          pelty = meshe_in % ltype(ielem)
          pgaus = ngaus(pelty)
          do inode = 1,pnode
             ipoin = meshe_in % lnods(inode,ielem)
             elcod(1:ndime,inode) = meshe_in % coord(1:ndime,ipoin)
          end do
          do igaus = 1,pgaus
             call elmgeo_cartesian_derivatives(&
                  ndime,pnode,elcod,elmar_in(pelty) % deriv(:,:,igaus),gpcar(1,1,igaus),gpvol(igaus))
             gpvol(igaus) = gpvol(igaus) * elmar_in(pelty) % weigp(igaus)
          end do
          call ADR_elemental_laplacian(&
               ielem,ndime,pnode,pgaus,mgaus,meshe_in % lnods(:,ielem),&
               meshe_in % r_dom,meshe_in % c_dom,gpcar,gpvol,amatr,kfl_fixno,bvess)
          !
          ! If a constant RHS is present
          !
          if( present(rhs_value) .and. present(rhsid) ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   ipoin = meshe_in % lnods(inode,ielem)
                   rhsid(ipoin) = rhsid(ipoin) + rhs_value(ipoin) * &
                        gpvol(igaus) * elmar_in(pelty) % shape(inode,igaus)
                end do
             end do
          end if
       end do

       deallocate(elcod)
       deallocate(gpvol)
       deallocate(gpsha)
       deallocate(gpcar)

    end if

  end subroutine ADR_assemble_laplacian

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Load mesh
  !> @details Load mesh arrays and parameters
  !
  !-----------------------------------------------------------------------

  subroutine ADR_assemble_extension(ndofn,meshe_in,elmar_in,veloc,amatr,rhsid)
    integer(ip),     intent(in)            :: ndofn
    type(mesh_type), intent(in)            :: meshe_in
    type(elm),       intent(in)            :: elmar_in(*)
    real(rp),        intent(in)            :: veloc(*)
    real(rp),        intent(out)           :: amatr(ndofn,ndofn,*)
    real(rp),        intent(out), optional :: rhsid(*)
    integer(ip)                            :: ndime,pnode,pgaus,inode
    integer(ip)                            :: mnode,pelty,mgaus,ielem
    integer(ip)                            :: igaus,ipoin,izdom,ntens

    integer(ip)                            :: porde,idime,plapl

    real(rp)                               :: gpsha_bub(2)
    real(rp)                               :: gpder_bub(2)
    real(rp)                               :: cutim
    real(rp)                               :: chale(2)
    real(rp)                               :: hleng(3)

    real(rp),        allocatable           :: elcod(:,:)
    real(rp),        allocatable           :: elvel(:,:)
    real(rp),        allocatable           :: elunk(:)
    real(rp),        allocatable           :: elmat(:,:)
    real(rp),        allocatable           :: elrhs(:)

    real(rp),        allocatable           :: gpvol(:)
    real(rp),        allocatable           :: gpsha(:,:)
    real(rp),        allocatable           :: gpcar(:,:,:)
    real(rp),        allocatable           :: gphes(:,:,:)

    real(rp),        allocatable           :: gpder(:,:,:) 

    real(rp),        allocatable           :: chave(:,:)
    real(rp),        allocatable           :: tragl(:,:)

    real(rp),        allocatable           :: gpden(:) 
    real(rp),        allocatable           :: gpvel(:,:)
    real(rp),        allocatable           :: gpdif(:) 
    real(rp),        allocatable           :: gpgrd(:,:) 
    real(rp),        allocatable           :: gprea(:,:) 
    real(rp),        allocatable           :: gprhs(:) 
    real(rp),        allocatable           :: gpunk(:) 

    type(ADR_typ)                          :: ADR 

    if( INOTMASTER ) THEN
       !
       ! Initialization
       !
       call ADR_initialize_type(ADR)

       ADR % kfl_time_integration   = OFF           ! Time flag              
       ADR % kfl_time_step_strategy = PRESCRIBED    ! Time step strategy

       ADR % kfl_stabilization      = SUPG          ! Stabilization strategy
       ADR % kfl_shock              = OFF           ! Shock capturing flag            
       ADR % kfl_time_lumped        = OFF           ! Lumped time evolution
       ADR % kfl_linearization      = NEWTON        ! Newton-Raphson for reaction
       ADR % kfl_tau_strategy       = TAU_CODINA    ! Tau strategy         
       ADR % kfl_laplacian          = OFF           ! If Laplacian should be computed
       ADR % kfl_time_sgs           = OFF           ! SGS tracking in time
       ADR % kfl_nonlinear_sgs      = OFF           ! SGS tracking in non-linear term
       ADR % kfl_time_bubble        = OFF           ! If bubble is tracked in time
       ADR % kfl_time_scheme        = TRAPEZOIDAL   ! If bubble is tracked in time
       ADR % kfl_time_order         = 1             ! First order
       ADR % kfl_manufactured       = 0             ! Manufactured solution
       ADR % kfl_length             = 0             ! Characteristic length
       ADR % kfl_first_order_sgs    = 0             ! SGS time integration like grid scale
       ADR % kfl_discretization     = 0             ! Discretization method
       ADR % number_euler_steps     = 0             ! Number of Euler time steps
       
       ADR % bemol                  = 0.0_rp        ! Integration by parts convective term
       ADR % tau_parameters         = 1.0_rp        ! Stability constants             
       ADR % shock                  = 0.0_rp        ! Shock capturing coefficients       
       ADR % relax_sgs              = 1.0_rp        ! No under relaxation   
       
       ADR % lun_output4            = 6_4           ! Output file unit
       ADR % memor                  = 0_8           ! Memory counter
       
       ADR % ntime                  = 2             ! Number of time values
       ADR % time_parameters        = 0.0_rp        ! BDF coefficients        
       ADR % number_components      = 2             ! Number of components  
       ADR % dtinv                  = 0.0_rp        ! 1/dt of n+1
       ADR % dtinv_old              = 0.0_rp        ! 1/dt of n,n-1,n-2, etc.

       if( present(rhsid) ) then
          do ipoin = 1,meshe_in % npoin * ndofn
             rhsid(ipoin) = 0.0_rp
          end do
       end if
       do izdom = 1,meshe_in % r_dom(meshe_in % npoin+1)-1
          amatr(1:ndofn,1:ndofn,izdom) = 0.0_rp
       end do

       ndime = meshe_in % ndime
       mnode = meshe_in % mnode
       mgaus = maxval(ngaus)
       ntens = 3 * ndime - 3

       allocate( elcod(ndime,mnode) )
       allocate( elvel(ndime,mnode) )
       allocate( elunk(mnode) )
       allocate( elmat(mnode,mnode) )
       allocate( elrhs(mnode) )

       allocate( gpvol(mgaus) )
       allocate( gpsha(mnode,mgaus) )
       allocate( gpcar(ndime,mnode,mgaus) )
       allocate( gphes(ntens,mnode,mgaus) )

       allocate( gpder(ndime,mnode,mgaus) )
       
       allocate( tragl(ndime,ndime) )
       allocate( chave(ndime,2) )

       allocate( gpden(mgaus) )
       allocate( gpvel(ndime,mgaus) )
       allocate( gpdif(mgaus) )
       allocate( gpgrd(ndime,mgaus) )
       allocate( gprea(mgaus,mreac_adr) )
       allocate( gprhs(mgaus) )
       allocate( gpunk(mgaus) )

       gpsha_bub = 0.0_rp
       gpder_bub = 0.0_rp 
       cutim     = 0.0_rp

       gpden     = 1.0_rp
       gpvel     = 0.0_rp
       gpgrd     = 0.0_rp
       gprea     = 0.0_rp
       gprhs     = 0.0_rp
       gpunk     = 0.0_rp
       gpdif     = 0.0_rp
       
       do ielem = 1,meshe_in % nelem
          pnode = meshe_in % lnnod(ielem)
          pelty = meshe_in % ltype(ielem)
          pgaus = ngaus(pelty)
          porde = lorde(pelty)
          plapl = llapl(pelty)
          do inode = 1,pnode
             ipoin = meshe_in % lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = meshe_in % coord(idime,ipoin)
                elvel(idime,inode) = veloc((ipoin-1)*ndime+idime)
             end do
             elunk(inode) = 0.0_rp
          end do

          call element_shape_function_derivatives_jacobian(&
               pnode,pgaus,plapl,elmar_in(pelty) % weigp,elmar_in(pelty) % shape,&
               elmar_in(pelty) % deriv,elmar_in(pelty) % heslo,&
               elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
          
          do igaus = 1,pgaus
             do idime = 1,ndime
                gpvel(idime,igaus) = dot_product(elvel(idime,1:pnode),elmar_in(pelty) % shape(1:pnode,igaus))
             end do
          end do

          call elmgeo_element_characteristic_length(&
               ndime,pnode,elmar_in(pelty) % dercg(:,:),elcod,hleng,element_type(pelty) % natural_length,tragl)
          call elmgeo_element_length(&
               ndime,pnode,porde,tragl,hleng,elcod,elvel,chave,chale,&
               element_type(pelty) % natural_length,1_ip,ADR % kfl_length)
          call ADR_element_assembly(&
               ielem,pnode,pgaus,elcod,elmar_in(pelty) % shape,gpcar,elmar_in(pelty) % deriv,gphes,&
               gpvol,chale,gpsha_bub,gpder_bub,ADR,cutim,gpden,gpvel,gpdif,gpgrd,&
               gprea,gprhs,gpunk,elunk,elmat,elrhs)
          !
          ! Assembly element matrix and RHS
          !
          call matrix_assemble_element_matrix_to_CSR(&
               kfl_element_to_csr,1_ip,pnode,pnode,&
               ielem,lnods(:,ielem),elmat,meshe_in % r_dom,meshe_in % c_dom,amatr,lezdo)
       end do

       deallocate( elcod )
       deallocate( elvel )
       deallocate( elunk )
       deallocate( elmat )
       deallocate( elrhs )

       deallocate( gpvol )
       deallocate( gpsha )
       deallocate( gpcar )
       deallocate( gphes )

       deallocate( gpder )
       deallocate( tragl )
       deallocate( chave )

       deallocate( gpden )
       deallocate( gpvel )
       deallocate( gpdif )
       deallocate( gpgrd )
       deallocate( gprea )
       deallocate( gprhs )
       deallocate( gpunk )

    end if

  end subroutine ADR_assemble_extension

    subroutine ADR_assemble_convective(ndofn,meshe_in,elmar_in,veloc,amatr,force, rhsid)
    !*******************************************************************************
    ! Matrix and rhs assembly for a pure convective equation using SUPG stabilization
    ! veloc grad u = rhs
    !********************************************************************************
    integer(ip),     intent(in)            :: ndofn
    type(mesh_type), intent(in)            :: meshe_in
    type(elm),       intent(in)            :: elmar_in(*)
    real(rp),        intent(in)            :: veloc(*)
    real(rp),        intent(in)            :: force(*)
    real(rp),        intent(out)           :: amatr(ndofn,ndofn,*)
    real(rp),        intent(out)           :: rhsid(*)
    integer(ip)                            :: ndime,pnode,pgaus,inode
    integer(ip)                            :: mnode,pelty,mgaus,ielem
    integer(ip)                            :: igaus,ipoin,ntens
!    integer(ip)                            :: izdom
    integer(ip)                            :: porde,idime, plapl

    real(rp)                               :: gpsha_bub(2)
    real(rp)                               :: gpder_bub(2)
    real(rp)                               :: cutim
    real(rp)                               :: chale(2)
    real(rp)                               :: hleng(3)

    real(rp),        allocatable           :: elcod(:,:)
    real(rp),        allocatable           :: elrhs(:), elfor(:)
    real(rp),        allocatable           :: elvel(:,:)
    real(rp),        allocatable           :: elunk(:)
    real(rp),        allocatable           :: elmat(:,:)

    real(rp),        allocatable           :: gpvol(:)
    real(rp),        allocatable           :: gpsha(:,:)
    real(rp),        allocatable           :: gpcar(:,:,:)
    real(rp),        allocatable           :: gphes(:,:,:)
    real(rp),        allocatable           :: gpder(:,:,:)
    
    real(rp),        allocatable           :: chave(:,:)
    real(rp),        allocatable           :: tragl(:,:)

    real(rp),        allocatable           :: gpden(:) 
    real(rp),        allocatable           :: gpvel(:,:)
    real(rp),        allocatable           :: gpdif(:) 
    real(rp),        allocatable           :: gpgrd(:,:) 
    real(rp),        allocatable           :: gprea(:,:) 
    real(rp),        allocatable           :: gprhs(:) 
    real(rp),        allocatable           :: gpunk(:) 

    type(ADR_typ)                          :: ADR 

    if( INOTMASTER ) THEN
       !
       ! Initialization
       !
       call ADR_initialize_type(ADR)

       ADR % kfl_time_integration   = OFF           ! Time flag              
       ADR % kfl_time_step_strategy = PRESCRIBED    ! Time step strategy
       ADR % kfl_stabilization      = SUPG          ! Stabilization strategy
       ADR % kfl_shock              = OFF           ! Shock capturing flag            
       ADR % kfl_time_lumped        = OFF           ! Lumped time evolution
       ADR % kfl_linearization      = PICARD        ! Newton-Raphson for reaction
       ADR % kfl_tau_strategy       = TAU_CODINA    ! Tau strategy         
       ADR % kfl_laplacian          = OFF           ! If Laplacian should be computed
       ADR % kfl_time_sgs           = OFF           ! SGS tracking in time
       ADR % kfl_nonlinear_sgs      = OFF           ! SGS tracking in non-linear term
       ADR % kfl_time_bubble        = OFF           ! If bubble is tracked in time
       ADR % kfl_time_scheme        = TRAPEZOIDAL   ! If bubble is tracked in time
       ADR % kfl_time_order         = 1             ! First order
       ADR % kfl_manufactured       = 0             ! Manufactured solution
       ADR % kfl_length             = 4             ! Characteristic length (0=minimum)
       ADR % kfl_first_order_sgs    = 0             ! SGS time integration like grid scale
       ADR % kfl_discretization     = 0             ! Discretization method (FE=0)
       ADR % number_euler_steps     = 0             ! Number of Euler time steps
       
       ADR % bemol                  = 0.0_rp        ! Integration by parts convective term
       ADR % tau_parameters         = 1.0_rp        ! Stability constants             
       ADR % shock                  = 0.0_rp        ! Shock capturing coefficients       
       ADR % relax_sgs              = 1.0_rp        ! No under relaxation   
       
       ADR % lun_output4            = 6_4           ! Output file unit
       ADR % memor                  = 0_8           ! Memory counter
       
       ADR % ntime                  = 2             ! Number of time values
       ADR % time_parameters        = 0.0_rp        ! BDF coefficients        
       ADR % number_components      = 2             ! Number of components  
       ADR % dtinv                  = 0.0_rp        ! 1/dt of n+1
       ADR % dtinv_old              = 0.0_rp        ! 1/dt of n,n-1,n-2, etc.

!!$       do ipoin = 1,meshe_in % npoin * ndofn
!!$          rhsid(ipoin) = 0.0_rp
!!$       end do
!!$       do izdom = 1,meshe_in % r_dom(meshe_in % npoin+1)-1
!!$          amatr(1:ndofn,1:ndofn,izdom) = 0.0_rp
!!$       end do
       ndime = meshe_in % ndime
       mnode = meshe_in % mnode
       mgaus = maxval(ngaus)
       ntens = 3 * ndime - 3

       allocate( elcod(ndime,mnode) )
       allocate( elvel(ndime,mnode) )
       allocate( elunk(mnode) )
       allocate( elmat(mnode,mnode) )
       allocate( elrhs(mnode) )
       allocate( elfor(mnode) )

       allocate( gpvol(mgaus) )
       allocate( gpsha(mnode,mgaus) )
       allocate( gpcar(ndime,mnode,mgaus) )
       allocate( gphes(ntens,mnode,mgaus) )
       allocate( gpder(ndime,mnode,mgaus) )
       
       allocate( tragl(ndime,ndime) )
       allocate( chave(ndime,2) )

       allocate( gpden(mgaus) )
       allocate( gpvel(ndime,mgaus) )
       allocate( gpdif(mgaus) )
       allocate( gpgrd(ndime,mgaus) )
       allocate( gprea(mgaus,mreac_adr) )
       allocate( gprhs(mgaus) )
       allocate( gpunk(mgaus) )

       gpsha_bub = 0.0_rp
       gpder_bub = 0.0_rp 
       cutim     = 0.0_rp

       gpden     = 1.0_rp
       gpvel     = 0.0_rp
       gpdif     = 0.0_rp
       gpgrd     = 0.0_rp
       gprea     = 0.0_rp
       gprhs     = 0.0_rp
       gpunk     = 0.0_rp

       do ielem = 1,meshe_in % nelem
          pnode = meshe_in % lnnod(ielem)
          pelty = meshe_in % ltype(ielem)
          pgaus = ngaus(pelty)
          porde = lorde(pelty)
          plapl = llapl(pelty)
          do inode = 1,pnode
             ipoin = meshe_in % lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = meshe_in % coord(idime,ipoin)
                elvel(idime,inode) = veloc((ipoin-1)*ndime+idime)
             end do            
             elfor(inode) = force(ipoin)
             elunk(inode) = 0.0_rp
          end do
          call element_shape_function_derivatives_jacobian(&
               pnode,pgaus,plapl,elmar_in(pelty) % weigp,elmar_in(pelty) % shape,&
               elmar_in(pelty) % deriv,elmar_in(pelty) % heslo,&
               elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
          
          do igaus = 1,pgaus
             do idime = 1,ndime
                gpvel(idime,igaus) = dot_product(elvel(idime,1:pnode),elmar_in(pelty) % shape(1:pnode,igaus))
             end do
             gprhs(igaus) = dot_product(elfor(1:pnode),elmar_in(pelty) % shape(1:pnode,igaus))
          end do
          call elmgeo_element_characteristic_length(&
               ndime,pnode,elmar_in(pelty) % dercg(:,:),elcod,hleng,element_type(pelty) % natural_length,tragl)
          call elmgeo_element_length(&
               ndime,pnode,porde,tragl,hleng,elcod,elvel,chave,chale,&
               element_type(pelty) % natural_length,1_ip,ADR % kfl_length)
 
          call ADR_element_assembly(&
               ielem,pnode,pgaus,elcod,elmar_in(pelty) % shape,gpcar,elmar_in(pelty) % deriv,gphes,&
               gpvol,chale,gpsha_bub,gpder_bub,ADR,cutim,gpden,gpvel,gpdif,gpgrd,&
               gprea,gprhs,gpunk,elunk,elmat,elrhs)
          !
          ! Assembly element matrix and RHS
          !
          call matrix_assemble_element_matrix_to_CSR(&
               kfl_element_to_csr,1_ip,pnode,pnode,&
               ielem,lnods(:,ielem),elmat,meshe_in % r_dom,meshe_in % c_dom,amatr,lezdo)
          call matrix_assemble_element_RHS(&
               1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
       end do

       deallocate( elcod )
       deallocate( elvel )
       deallocate( elunk )
       deallocate( elmat )
       deallocate( elrhs )

       deallocate( gpvol )
       deallocate( gpsha )
       deallocate( gpcar )
       deallocate( gphes )
       deallocate( gpder )
       
       deallocate( tragl )
       deallocate( chave )

       deallocate( gpden )
       deallocate( gpvel )
       deallocate( gpdif )
       deallocate( gpgrd )
       deallocate( gprea )
       deallocate( gprhs )
       deallocate( gpunk )

    end if

  end subroutine ADR_assemble_convective
  !-------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Elemental Laplacian
  !> @details Elemental Laplacian
  !
  !-------------------------------------------------------------------

  subroutine ADR_elemental_laplacian(&
       ielem,ndime,pnode,pgaus,mgaus,lnods,ia,ja,gpcar,&
       gpvol,amatr,kfl_fixno,bvess)
    integer(ip), intent(in)           :: ielem
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: pgaus
    integer(ip), intent(in)           :: mgaus
    integer(ip), intent(in)           :: lnods(pnode)
    integer(ip), intent(in)           :: ia(*)
    integer(ip), intent(in)           :: ja(*)
    real(rp),    intent(in)           :: gpcar(ndime,mnode,pgaus)
    real(rp),    intent(in)           :: gpvol(pgaus)
    real(rp),    intent(inout)        :: amatr(*)
    integer(ip), intent(in), optional :: kfl_fixno(:,:)
    real(rp),    intent(in), optional :: bvess(:,:)
    real(rp)                          :: elmat(pnode,pnode),fact1
    integer(ip)                       :: inode,igaus,jpoin,ipoin,jnode,izsol,jcolu

    elmat = 0.0_rp

    do igaus = 1,pgaus
       do jnode = 1,pnode
          do inode = 1,pnode
             elmat(inode,jnode) = elmat(inode,jnode) &
                  + gpvol(igaus) * dot_product(gpcar(1:ndime,inode,igaus),gpcar(1:ndime,jnode,igaus))
          end do
       end do
    end do
    !
    ! Prescribe Laplacian 
    !
    if( present(kfl_fixno) ) then
       if( present(bvess) ) then
          call runend('NOT CODED')
          !do inode = 1,pnode
          !   ipoin = lnods(inode)
          !   if( kfl_fixno(1,ipoin) == 1 ) then
          !      fact1 = elmat(inode,inode)
          !      do jnode = 1,pnode
          !         elmat(inode,jnode) = 0.0_rp
          !         elrhs(jnode)       = elrhs(jnode)-elmat(jnode,inode) * bvess(1,ipoin)
          !         elmat(jnode,inode) = 0.0_rp
          !      end do
          !      elmat(inode,inode) = fact1
          !      elrhs(inode)       = fact1*bvess
          !   end if
          !end do
       else
          do inode = 1,pnode
             ipoin = lnods(inode)
             if( kfl_fixno(1,ipoin) == 1 ) then
                fact1 = elmat(inode,inode)
                do jnode = 1,pnode
                   elmat(inode,jnode) = 0.0_rp
                   elmat(jnode,inode) = 0.0_rp
                end do
                elmat(inode,inode) = fact1
             end if
          end do          
       end if
    end if
    !
    ! Assembly
    !
    if( kfl_element_to_csr == 1 ) then
       do inode = 1,pnode
          do jnode = 1,pnode
             izsol = lezdo(inode,jnode,ielem)         
             amatr(izsol) = amatr(izsol) + elmat(inode,jnode)
          end do
       end do
    else
       do inode = 1,pnode
          ipoin = lnods(inode)
          do jnode = 1,pnode
             jpoin = lnods(jnode)
             izsol = ia(ipoin)
             jcolu = ja(izsol)
             do while( jcolu /= jpoin .and. izsol < ia(ipoin+1)-1)
                izsol = izsol + 1
                jcolu = ja(izsol)
             end do
             if( jcolu == jpoin ) then             
                amatr(izsol) = amatr(izsol) + elmat(inode,jnode)             
             end if
          end do
       end do
    end if

  end subroutine ADR_elemental_laplacian 
 
  !-------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Check data compatibility  
  !> @details Check data compatibility and compute some derived
  !>          parameters
  !
  !-------------------------------------------------------------------

  subroutine ADR_check_and_compute_data(ADR)
    type(ADR_typ), intent(inout) :: ADR !< ADR type
    !
    ! Stabilization compatibility
    !
    if( ADR % kfl_stabilization == BUBBLE .and. ( ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) ) then
       call runend('ADR_check_and_compute_data: BUBBLE AND SGS ARE INCOMPATIBLE')
    end if
    if( ADR % kfl_stabilization == AR_OSS .and. ( ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) ) then
       call runend('ADR_check_and_compute_data: AR-OSS AND SGS ARE INCOMPATIBLE')
    end if
    if( ADR % kfl_time_step_strategy == LOCAL .and. ADR % kfl_time_order /= 1 ) then
       call runend('ADR_check_and_compute_data: LOCAL TIME STEP REQUIRES FIRST ORDER TIME INTEGRATION')
    end if
    !
    ! Time integration
    !
    if( ADR % kfl_time_integration == 0 ) then
       ADR % kfl_time_order    = 1
       ADR % number_components = 2
       ADR % ntime             = 1 
    else
       if(ADR % kfl_time_scheme == TRAPEZOIDAL ) then
          ADR % number_components = 3                              ! Trapezoidal rule
          ADR % ntime             = 2
       else if(ADR % kfl_time_scheme == BDF) then
          ADR % number_components = 2 + ADR % kfl_time_order       ! BDF scheme
          ADR % ntime             = 1 + ADR % kfl_time_order
       else if(ADR % kfl_time_scheme == ADAMS_BASHFORD) then
          ADR % number_components = 3  ! AB scheme
          ADR % ntime             = 2
       else if(ADR % kfl_time_scheme == RUNGE_KUTTA) then
          ADR % number_components = 3  ! AB scheme
          ADR % ntime             = 2
       else
          call runend('ADR_CHECK_AND_COMPUTE_DATA: UNKNOWN TIME SCHEME')
       end if
    end if
    !
    ! Save original time integration order as it can be modified
    ! along the run
    !
    ADR % kfl_time_order_save = ADR % kfl_time_order
    !
    ! Number of unknowns
    !
    if(      ADR % kfl_discretization == 0 ) then
       ADR % nunkn = npoin
    else if( ADR % kfl_discretization == 1 ) then
       ADR % nunkn = nelem       
    end if

  end subroutine ADR_check_and_compute_data

  subroutine ADR_arrays(wtask,ADR,wbubble,wproje1,wproje2,wsgs,TAG1)

    use mod_arrays, only : arrays
    use mod_arrays, only : arrays_number

    character(len=*), intent(in)             :: wtask
    type(ADR_typ),    intent(inout)          :: ADR
    character(len=*), intent(in)             :: wbubble
    character(len=*), intent(in)             :: wproje1
    character(len=*), intent(in)             :: wproje2
    character(len=*), intent(in)             :: wsgs
    integer(ip),      intent(in),   optional :: TAG1
    integer(ip)                              :: ielem,pgaus

    if( ADR % kfl_stabilization == BUBBLE ) then
       !
       ! Bubble
       !
       call arrays(arrays_number(trim(wbubble)),wtask,ADR % bubble,nelem,ADR % number_components)

    else if( ADR % kfl_stabilization > 0 ) then 
       !
       ! Projections
       !
       call arrays(arrays_number(trim(wproje1)),wtask,ADR % proje1,ADR % nunkn)

       if( ADR % kfl_stabilization == AR_OSS ) then
          call arrays(arrays_number(trim(wproje2)),wtask,ADR % proje2,ADR % nunkn)             
       end if

    end if
    !
    ! SGS
    !
    if( ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) then
       call arrays(arrays_number(trim(wsgs)),wtask,ADR % sgs,nelem)
       if( trim(wtask) == 'ALLOCATE') then
          if( ADR % kfl_first_order_sgs == 1 ) then
             do ielem = 1,nelem
                pgaus = ngaus(abs(ltype(ielem)))
                call memory_alloca(ADR % memor,'ADR % sgs','ADR_allocate_projections_and_bubble',ADR % sgs(ielem)%a,1_ip,pgaus,3_ip)
             end do
          else
             do ielem = 1,nelem
                pgaus = ngaus(abs(ltype(ielem)))
                call memory_alloca(ADR % memor,'ADR % sgs','ADR_allocate_projections_and_bubble',ADR % sgs(ielem)%a,1_ip,pgaus,ADR % number_components)
             end do
          end if
       end if
    end if

  end subroutine ADR_arrays
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Allocate memory for bubble, projections and SGS
  !> @details Allocate memory for bubble, projections and SGS 
  !
  !-----------------------------------------------------------------------

  subroutine ADR_allocate_projections_bubble_sgs(ADR)
    type(ADR_typ), intent(inout) :: ADR
    integer(ip)                  :: ielem,pgaus

    if( INOTMASTER ) then

       if( ADR % kfl_stabilization == BUBBLE ) then
          !
          ! Bubble
          !
          if( .not. associated(ADR % bubble) ) call memory_alloca(ADR % memor,'ADR % bubble','ADR_allocate_projections_and_bubble',ADR % bubble,nelem,ADR % number_components)

       else if( ADR % kfl_stabilization > 0 ) then 
          !
          ! Projections
          !
          if( .not. associated(ADR % proje1) ) call memory_alloca(ADR % memor,'ADR % proje1','ADR_allocate_projections_and_bubble',ADR % proje1,ADR % nunkn)

          if( ADR % kfl_stabilization == AR_OSS ) then
             if( .not. associated(ADR % proje2) ) call memory_alloca(ADR % memor,'ADR % proje2','ADR_allocate_projections_and_bubble',ADR % proje2,ADR % nunkn)
          end if

       end if
       !
       ! SGS
       !
       if( ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) then
          if( .not. associated(ADR % sgs) ) call memory_alloca(ADR % memor,'ADR % sgs','ADR_allocate_projections_and_bubble',ADR % sgs,nelem)
             if( ADR % kfl_first_order_sgs == 1 ) then
                do ielem = 1,nelem
                   pgaus = ngaus(abs(ltype(ielem)))
                   call memory_alloca(ADR % memor,'ADR % sgs','ADR_allocate_projections_and_bubble',ADR % sgs(ielem)%a,1_ip,pgaus,3_ip)
                end do
             else
                do ielem = 1,nelem
                   pgaus = ngaus(abs(ltype(ielem)))
                   call memory_alloca(ADR % memor,'ADR % sgs','ADR_allocate_projections_and_bubble',ADR % sgs(ielem)%a,1_ip,pgaus,ADR % number_components)
                end do
             end if
       end if

    end if

  end subroutine ADR_allocate_projections_bubble_sgs

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Initialize projections
  !> @details Initialize projections according to stabilization method
  !
  !-----------------------------------------------------------------------

  subroutine ADR_initialize_projections(ADR)
    type(ADR_typ), intent(out) :: ADR
    integer(ip)                :: ii

    if( INOTMASTER .and. ADR % kfl_stabilization > 0 ) then
       call memory_alloca(ADR % memor,'ADR % proje1_tmp','ADR_initialize_projections',ADR % proje1_tmp,ADR % nunkn)
       do ii = 1,ADR % nunkn
          ADR % proje1_tmp(ii) = 0.0_rp
       end do
       if( ADR % kfl_stabilization == AR_OSS ) then
          call memory_alloca(ADR % memor,'ADR % proje2_tmp','ADR_initialize_projections',ADR % proje2_tmp,ADR % nunkn)
          do ii = 1,ADR % nunkn
             ADR % proje2_tmp(ii) = 0.0_rp
          end do
       end if
    end if

  end subroutine ADR_initialize_projections

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Load mesh
  !> @details Load mesh arrays and parameters
  !
  !-----------------------------------------------------------------------

  subroutine ADR_load_mesh(meshe_in,elmar_in)
    type(mesh_type), intent(in)          :: meshe_in(1)
    type(elm),       intent(in), pointer :: elmar_in(:)

    ndime     =  meshe_in(1) % ndime 
    ntens     =  meshe_in(1) % ntens
    npoin     =  meshe_in(1) % npoin
    nelem     =  meshe_in(1) % nelem
    mnode     =  meshe_in(1) % mnode
    mgaus     =  meshe_in(1) % mgaus
    lnods     => meshe_in(1) % lnods
    ltype     => meshe_in(1) % ltype
    lnnod     => meshe_in(1) % lnnod
    coord     => meshe_in(1) % coord
    elmar_loc => elmar_in

  end subroutine ADR_load_mesh

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Initialize ADR type
  !> @details Initialize ADR type
  !
  !-----------------------------------------------------------------------

  subroutine ADR_initialize_type(ADR)
    type(ADR_typ), intent(inout) :: ADR

    ADR % kfl_time_integration   = OFF           ! Time flag              
    ADR % kfl_time_step_strategy = PRESCRIBED    ! Time step strategy
    ADR % kfl_stabilization      = GALERKIN      ! Stabilization strategy
    ADR % kfl_shock              = OFF           ! Shock capturing flag            
    ADR % kfl_time_lumped        = OFF           ! Lumped time evolution
    ADR % kfl_linearization      = NEWTON        ! Newton-Raphson for reaction
    ADR % kfl_tau_strategy       = TAU_CODINA    ! Tau strategy         
    ADR % kfl_laplacian          = OFF           ! If Laplacian should be computed
    ADR % kfl_time_sgs           = OFF           ! SGS tracking in time
    ADR % kfl_nonlinear_sgs      = OFF           ! SGS tracking in non-linear term
    ADR % kfl_time_bubble        = OFF           ! If bubble is tracked in time
    ADR % kfl_time_scheme        = TRAPEZOIDAL   ! If bubble is tracked in time
    ADR % kfl_time_order         = 1             ! First order
    ADR % kfl_manufactured       = 0             ! Manufactured solution
    ADR % kfl_length             = 0             ! Characteristic length
    ADR % kfl_first_order_sgs    = 0             ! SGS time integration like grid scale
    ADR % kfl_discretization     = 0             ! Discretization method
    ADR % kfl_skewsymm           = 0             ! Skew symmetric discretization of convective term
    ADR % number_euler_steps     = 0             ! Number of Euler time steps

    ADR % bemol                  = 0.0_rp        ! Integration by parts convective term
    ADR % tau_parameters         = 1.0_rp        ! Stability constants             
    ADR % shock                  = 0.0_rp        ! Shock capturing coefficients       
    ADR % relax_sgs              = 1.0_rp        ! No under relaxation   

    ADR % lun_output4            = 6_4           ! Output file unit
    ADR % memor                  = 0_8           ! Memory counter

    ADR % ntime                  = 2             ! Number of time values
    ADR % time_parameters        = 0.0_rp        ! BDF coefficients        
    ADR % number_components      = 2             ! Number of components  
    ADR % dtinv                  = 0.0_rp        ! 1/dt of n+1
    ADR % dtinv_old              = 0.0_rp        ! 1/dt of n,n-1,n-2, etc.

    nullify(ADR % proje1)                        ! First projection 
    nullify(ADR % proje2)                        ! Second projection  
    nullify(ADR % proje1_tmp)                    ! Temporary first projection 
    nullify(ADR % proje2_tmp)                    ! Temporary second projection  
    nullify(ADR % sgs)                           ! Subgrid scale
    nullify(ADR % bubble)                        ! Bubble


  end subroutine ADR_initialize_type

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and Matias Avila
  !> @brief   Elemental assembly of the ADR equation
  !> @details Elemental assembly of the ADR equation:
  !>                                                        4
  !>          rho*du/dt + rho*a.grad(u) - div[k*grad(u)] + sum ri*u^i = f 
  !>                                                       i=1
  !>
  !>          where Gauss points values are:
  !>
  !>          1/dt .................... DTINV
  !>          rho ..................... GPDEN(:)
  !>          a ....................... GPVEL(:)
  !>          ri ...................... GPREA(:,i)
  !>          k ....................... GPDIF
  !>          grad(k) ................. GPGRD(:) 
  !>          f ....................... GPRHS(:)
  !>          u^n+1,u^n,u^n-1,etc ..... GPUNK(:,:)
  !>          Same at element nodes ... ELUNK(:,:)
  !>
  !>
  !>          o---------o pnode nodes
  !>          | x     x | pgaus Gauss points
  !>          |    e    | e is the bubble
  !>          | x     x |
  !>          o---------o
  !>
  !>
  !>          * Linearization:
  !>          ----------------
  !>
  !>          According to the linearization strategy, the reaction
  !>          terms are recomputed in REACT(:) and GPRHS(:) <= RHSIT(:).
  !>          Previous time terms are thrown to the RHS, RHSIT, as well.
  !>
  !>          GPDEN         GPVEL          GPDIF         REACT RHSIT
  !>          rho*u/dt + rho*a.grad(u) - div[k*grad(u)] + r*u = f' 
  !>          
  !>          * Stabilization:
  !>          ----------------
  !>          
  !>          The reaction term in the adjoint L*(v) and in the 
  !>          computation of tau use SREAC(:)
  !>
  !>          * Bubble array:
  !>          ---------------
  !>
  !>          Bubble array ADR % bubble(1:nelem,:) is defined just as the
  !>          primary unknown:
  !>          ADR % bubble(1:nelem,1) = current value at time n+1
  !>          ADR % bubble(1:nelem,2) = last block value
  !>          ADR % bubble(1:nelem,3) = value at time n
  !>          ADR % bubble(1:nelem,4) = value at time n-1 
  !>          ...
  !>
  !
  !-----------------------------------------------------------------------

  subroutine ADR_element_assembly(&
       ielem,pnode,pgaus,&
       elcod,gpsha,gpcar,gpder,gphes,gpvol,chale,gpsha_bub,gpder_bub,ADR,&
       cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpunk,elunk,elmat,elrhs,&
       messa, eladv)
    !
    ! Element dimensions
    !
    integer(ip), intent(in)            :: ielem                         !< Current element
    integer(ip), intent(in)            :: pnode                         !< # nodes
    integer(ip), intent(in)            :: pgaus                         !< # Gauss points
    !
    ! Element characteristics at Gauss point
    !
    real(rp),    intent(in)            :: elcod(ndime,pnode)            !< Element node coordinates
    real(rp),    intent(in)            :: gpsha(pnode,pgaus)            !< Shape function Nk
    real(rp),    intent(in)            :: gpcar(ndime,mnode,pgaus)      !< Shape function Cartesian derivatives dNk/dxi
    real(rp),    intent(in)            :: gpder(ndime,pnode,pgaus)      !< Shape function derivatives DNk/dsi
    real(rp),    intent(in)            :: gphes(ntens,mnode,pgaus)      !< Hessian dNk/dxidxj
    real(rp),    intent(in)            :: gpvol(pgaus)                  !< Element Jacobian
    real(rp),    intent(in)            :: chale(2)                      !< Element characteristic length
    real(rp),    intent(in)            :: gpsha_bub(pgaus)              !< Bubble shape function Ne
    real(rp),    intent(in)            :: gpder_bub(ndime,pgaus)        !< Bubble shape function derivatives dNe/dxi
    ! 
    ! Numerical strategy
    !
    type(ADR_typ), intent(inout)          :: ADR                           !< ADR type 
    !
    ! Equation coefficients
    !
    real(rp),    intent(in)            :: cutim                         !< Current time
    real(rp),    intent(in)            :: gpden(pgaus)                  !< Density
    real(rp),    intent(in)            :: gpvel(ndime,pgaus)            !< Advection vector
    real(rp),    intent(in)            :: gpdif(pgaus)                  !< Diffusion 
    real(rp),    intent(in)            :: gpgrd(ndime,pgaus)            !< Diffusion gradient
    real(rp),    intent(in)            :: gprea(pgaus,4)                !< Reaction
    real(rp),    intent(in)            :: gprhs(pgaus)                  !< RHS
    real(rp),    intent(in)            :: gpunk(pgaus,*)                !< Unknown at Gauss point
    real(rp),    intent(in)            :: elunk(pnode,2)                !< Element unknown
    !
    ! Output
    !
    real(rp),    intent(out)           :: elmat(pnode,pnode)            !< Element matrix
    real(rp),    intent(out)           :: elrhs(pnode)                  !< Element RHS
    character(*),intent(in),  optional :: messa                         !< Message
    real(rp),    intent(in),  optional :: eladv(ndime, pnode)           !< Element velocities

    integer(ip)                        :: inode,jnode,idime
    integer(ip)                        :: itime,igaus,ipoin 
    integer(ip)                        :: jdime
    integer(ip)                        :: pnode_SU_A_OSS
    integer(ip)                        :: pnode_AR_OSS
    real(rp)                           :: fact1,fact2,fact3,gpdiv
    real(rp)                           :: fact4,fact5,fact6
    real(rp)                           :: gppe1(mnode),rhsit(mgaus)
    real(rp)                           :: xmuit,tau,gpadv(pnode),gpad1(pnode)
    real(rp)                           :: grvgr,resi2(pnode)
    real(rp)                           :: resi1(pnode),xmui3
    real(rp)                           :: SD(mgaus),gplap,CD(mgaus)
    real(rp)                           :: sreac(pgaus),gppe2(pnode)
    real(rp)                           :: galer(pnode),bemo1,fact7
    real(rp)                           :: gptau(pgaus)
    real(rp)                           :: gptau_time(pgaus)
    real(rp)                           :: exunk(1),exgra(1,1)
    real(rp)                           :: alpha,alpha1,sgs(mgaus)
    !
    ! Numerical strategy
    !
    integer(ip)                        :: kfl_stabilization             ! Stabilization strategy
    integer(ip)                        :: kfl_shock                     ! Shock capturing strategy
    integer(ip)                        :: kfl_time_lumped               ! Time term lumped
    integer(ip)                        :: kfl_tau_strategy              ! Tau strategy
    integer(ip)                        :: kfl_laplacian                 ! If Laplacian should be computed
    integer(ip)                        :: kfl_time_bubble               ! If bubble is tracked in time
    integer(ip)                        :: kfl_time_scheme               ! Time scheme
    integer(ip)                        :: kfl_time_order                ! Time order
    integer(ip)                        :: kfl_manufactured              ! Manufactured solution
    real(rp)                           :: tau_parameters(3)             ! Tau multiplication factors
    real(rp)                           :: shock                         ! Shock capturing constant
    real(rp)                           :: bemol,d                       ! Convective term, integration by parts
!    real(rp)                           :: tragl(3,3)
!    real(rp)                           :: hleng(3)
    real(rp)                           :: dtinv                         ! Inverse of time step
    !
    ! Time integration
    !
    integer(ip)                        :: ntime        
    !
    ! Bubble stuffs
    !
    real(rp)                           :: elmat_ue(pnode,1)
    real(rp)                           :: elmat_eu(1,pnode)
    real(rp)                           :: elmat_ee(1,1)
    real(rp)                           :: elrhs_ee(1)   
    real(rp)                           :: galer_bub,gpadv_bub
    real(rp)                           :: xtime,gpdet,react(mgaus)
    real(rp)                           :: xjaci(ndime,ndime)
    real(rp)                           :: gpcar_bub(ndime,pgaus)
    real(rp)                           :: elbub(10),c1sgs
    logical(lg)                        :: update_bubble
    !
    ! Projections stuffs
    !
    real(rp)                           :: gppr1(mgaus)
    real(rp)                           :: gppr2(mgaus)

    !--------------------------------------------------------------------
    !
    ! Initial computations 
    !
    !--------------------------------------------------------------------
    ! 
    ! Numerical treatment
    !
    kfl_stabilization = ADR % kfl_stabilization 
    kfl_shock         = ADR % kfl_shock
    kfl_time_lumped   = ADR % kfl_time_lumped
    kfl_tau_strategy  = ADR % kfl_tau_strategy
    kfl_laplacian     = ADR % kfl_laplacian
    kfl_time_bubble   = ADR % kfl_time_bubble
    kfl_time_scheme   = ADR % kfl_time_scheme
    kfl_time_order    = ADR % kfl_time_order
    tau_parameters    = ADR % tau_parameters
    kfl_manufactured  = ADR % kfl_manufactured
    ntime             = ADR % ntime
    shock             = ADR % shock
    bemol             = ADR % bemol
    bemo1             = 1.0_rp - bemol
    xtime             = 0.0_rp 
    sgs(1:pgaus)      = 0.0_rp
    dtinv             = ADR % dtinv
    update_bubble     = .false.
    !
    ! Initialization 
    !
    elrhs(1:pnode)         = 0.0_rp
    elmat(1:pnode,1:pnode) = 0.0_rp
    gppr1(1:pgaus)         = 0.0_rp
    gppr2(1:pgaus)         = 0.0_rp
    !
    !----------------------------------------------------------------------
    !
    ! Explicit integration: time derivatives is computed using mass matrix
    ! later on, when updating the unknown using AB or RK schemes
    !
    !----------------------------------------------------------------------
    if (kfl_time_scheme == ADAMS_BASHFORD .or. kfl_time_scheme == RUNGE_KUTTA) then
       dtinv = 0.0_rp
    end if
    !----------------------------------------------------------------------
    !
    ! Assemble projections
    !
    if( kfl_stabilization == AR_OSS ) then 
       do igaus = 1,pgaus
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             gppr1(igaus) = gppr1(igaus) + ADR % proje1(ipoin) * gpsha(inode,igaus)
             gppr2(igaus) = gppr2(igaus) + ADR % proje2(ipoin) * gpsha(inode,igaus)
          end do
       end do
    else if( kfl_stabilization > 0 ) then 
       do igaus = 1,pgaus
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             gppr1(igaus) = gppr1(igaus) + ADR % proje1(ipoin) * gpsha(inode,igaus)
          end do
       end do
    end if
    !
    ! Time difference scheme to RHS for all schemes
    ! ELBUB(1) = bubble at current time n+1
    ! ELBUB(2) = bubble at time n
    ! ELBUB(3) = bubble at time n-1 
    !
    rhsit(1:pgaus) = gprhs(1:pgaus)

    if( ADR % kfl_time_integration /= 0 ) then

       if( kfl_stabilization == BUBBLE .and. kfl_time_bubble /= 0 ) then

          elbub(1) = ADR % bubble(ielem,1)
          elbub(2) = ADR % bubble(ielem,3)
          if( kfl_time_scheme == BDF ) then
             do itime = 3,size(ADR % bubble,KIND=ip)-1
                elbub(itime) = ADR % bubble(ielem,itime+1)
             end do
          end if
          xtime = 1.0_rp
          do igaus = 1,pgaus
             fact1 = gpden(igaus) * dtinv
             do itime = 2,ntime
                rhsit(igaus) = rhsit(igaus) - fact1 * ADR % time_parameters(itime) * ( gpunk(igaus,itime) + gpsha_bub(igaus) * elbub(itime) )
             end do
          end do

       else

          do igaus = 1,pgaus
             fact1 = gpden(igaus) * dtinv
             do itime = 2,ntime
                rhsit(igaus) = rhsit(igaus) - fact1 * ADR % time_parameters(itime) * gpunk(igaus,itime)
             end do
          end do

       end if

    end if
    !
    ! Linearization of reaction term
    ! LHS + r1*u + r2*u^2 = RHS
    !        
    if( ADR % kfl_linearization == RHS ) then
       !
       ! LHS + r1*u = RHS - r2*u_i^2
       !
       react(1:pgaus) = gprea(1:pgaus,1) 
       sreac(1:pgaus) = gprea(1:pgaus,1)
       rhsit(1:pgaus) = rhsit(1:pgaus) - gprea(1:pgaus,2) * gpunk(1:pgaus,1) ** 2

    else if( ADR % kfl_linearization == PICARD ) then
       !
       ! LHS + (r1+r2*u_i)*u = RHS 
       !
       react(1:pgaus) = gprea(1:pgaus,1) + gprea(1:pgaus,2) * gpunk(1:pgaus,1) 
       sreac(1:pgaus) = react(1:pgaus)
       rhsit(1:pgaus) = rhsit(1:pgaus)

    else if( ADR % kfl_linearization == NEWTON ) then
       !
       ! LHS + u * (r1+2*r2*u_i) = RHS + r2*u_i^2
       !
       react(1:pgaus) = gprea(1:pgaus,1) + 2.0_rp * gprea(1:pgaus,2) * gpunk(1:pgaus,1)
       sreac(1:pgaus) = gprea(1:pgaus,1) + gprea(1:pgaus,2) * gpunk(1:pgaus,1)
       rhsit(1:pgaus) = rhsit(1:pgaus)   + gprea(1:pgaus,2) * gpunk(1:pgaus,1) ** 2

    end if

    if( kfl_stabilization == SUPG ) sreac = 0.0_rp
    !
    ! Element characteristic length 
    !
    !call elmgeo_element_characteristic_length(ndime,pnode,deriv,elcod,hleng,hnatu,tragl)
    !call ADR_element_lengths(tragl,hleng,elcod,chale,pnode,1_ip,hnatu,ADR % kfl_length)
    !
    ! Stabilization parameter without reactive (nonlinear) term :
    ! tau = 1.0_rp /( 4.0_rp*gpdif(igaus)/chale(2)/chale(2) + 2.0_rp*rhnve/chale(1) ) 
    !
    if( ( kfl_stabilization == GALERKIN .or. kfl_stabilization == BUBBLE ) .and. kfl_shock == 0 ) then
       gptau      = 0.0_rp
       gptau_time = 0.0_rp
       c1sgs      = 0.0_rp
    else
       call ADR_tau(&
            pgaus,kfl_tau_strategy,tau_parameters,gpden,gpvel,gpdif,sreac,chale(1),chale(2),gptau, dtinv)
       if( ADR % kfl_time_sgs /= 0 ) then
          if( ADR % kfl_first_order_sgs == 1 ) then
             c1sgs = 1.0_rp
             gptau_time(1:pgaus) = 1.0_rp / ( gpden(1:pgaus) * dtinv + 1.0_rp / gptau(1:pgaus) )
             sgs(1:pgaus) = ADR % sgs(ielem) % a(1,1:pgaus,3)
          else
             c1sgs = ADR % time_parameters(1)
             gptau_time(1:pgaus) = 1.0_rp / ( c1sgs * gpden(1:pgaus) * dtinv + 1.0_rp / gptau(1:pgaus) )
             do itime = 2,ntime
                sgs(1:pgaus) = sgs(1:pgaus) - ADR % time_parameters(itime) * ADR % sgs(ielem) % a(1,1:pgaus,itime+1)
             end do
          end if
       else
          gptau_time = gptau
          c1sgs = 0.0_rp
       end if
    end if
    !
    ! Shock capturing
    !
    if( kfl_shock /= 0 ) then
       call ADR_shock_capturing(&
            ndime,mnode,ntens,pgaus,pnode,kfl_laplacian,kfl_shock,ADR % time_parameters,shock,&
            dtinv,rhsit,gpden,gpvel,gpdif,gpgrd,gprea,gpunk,gptau,chale,gpcar,gphes,elunk,SD,CD)
    else
       CD = 0.0_rp
       SD = 0.0_rp
    end if
    !
    ! Manufactured solution
    !
    if( kfl_manufactured /= 0 ) then
       call ADR_manufactured_solution_and_rhs(&
            2_ip,kfl_manufactured,pnode,pgaus,cutim,elcod,exunk,&
            exgra,gpsha,gpden,gpdif,react,gpgrd,gpvel,rhsit)
    end if

    if( kfl_stabilization == GALERKIN .or. kfl_stabilization == SU .or. kfl_stabilization == A_OSS .or. kfl_stabilization == AR_OSS ) then

       !-------------------------------------------------------------------
       !
       ! Galerkin, SU, A_OSS, AR_OSS
       !
       ! k (grad(u),grad(v)) + (1-b)*(rho*a.grad(u),v) -b*(u,rho*a.grad(v))
       ! + (s*u,v) + SU(u,v) = (f,v)
       !
       ! SU(u,v) = (tau*rho*a.grad(v),rho*a.grad(u)-P_adv) 
       ! where P_adv = P(rho*a.grad(u)) for A_OSS
       !
       !-------------------------------------------------------------------

       if( kfl_stabilization == GALERKIN ) then
          pnode_SU_A_OSS = 0
       else
          pnode_SU_A_OSS = pnode
       end if
       if( kfl_stabilization == AR_OSS ) then
          pnode_AR_OSS = pnode
       else
          pnode_AR_OSS = 0
       end if

       do igaus = 1,pgaus
          !
          ! Coefficients
          !
          fact2 = gpvol(igaus) * gpdif(igaus)              ! k
          fact3 = gpvol(igaus) * gpden(igaus)              ! rho
          fact4 = gpvol(igaus) * rhsit(igaus)              ! f
          fact5 = gpvol(igaus) * react(igaus)              ! r
          fact6 = fact3 * dtinv * ADR % time_parameters(1) ! rho/dt
          fact7 = bemol * fact3                            ! rho * bemol
          !
          ! Advections
          ! a.grad(Ni)
          !
          gpadv(1:pnode) = 0.0_rp
          do idime = 1,ndime
             gpadv(1:pnode) = gpadv(1:pnode) + gpvel(idime,igaus) * gpcar(idime,1:pnode,igaus)           
          end do
          !
          ! Galerkin operator
          ! Nj/dt + r*Nj + rho*a.grad(Nj)
          !        
          galer(1:pnode) = ( fact6 + fact5 ) * gpsha(1:pnode,igaus) + fact3 * bemo1 * gpadv(1:pnode)
          !
          ! Galerkin operator + diffusion + convection
          ! Galerkin + k grad(Ni).grad(Nj) - b*( u, rho*a.grad(v) ) 
          !        
          do inode = 1,pnode
             do jnode = 1,pnode
                elmat(inode,jnode) = elmat(inode,jnode) + gpsha(inode,igaus) * galer(jnode)
                elmat(inode,jnode) = elmat(inode,jnode) - fact7 * gpadv(inode) * gpsha(jnode,igaus)
                do idime = 1,ndime
                   elmat(inode,jnode) = elmat(inode,jnode) + fact2 * gpcar(idime,inode,igaus) * gpcar(idime,jnode,igaus)
                end do
             end do
             elrhs(inode) = elrhs(inode) + gpsha(inode,igaus) * fact4
          end do
          !
          ! Skew symmetric advection
          !
          if (ADR % kfl_skewsymm == 1) then  ! Adds + 0.5*rho div(u)*scalar
             gpdiv = 0.0_rp
             ! calculates velocity divergence 
             do inode=1, pnode
                do idime=1, ndime
                   gpdiv = gpdiv + gpcar(idime, inode, igaus)*eladv(idime, inode)
                end do
             end do
             !  Adds the term + 0.5*rho div(u)*scalar
             do inode=1, pnode
                fact1 = gpsha(inode, igaus)*gpvol(igaus)*0.5_rp*gpden(igaus)*gpdiv
                do jnode =1, pnode
                   elmat(inode, jnode) = elmat(inode, jnode) + gpsha(jnode, igaus)*fact1
                end do
             end do
             ! Another way of skew symmetric terms is to use bemol =0.5,  needing boundary terms            
          end if
          !
          ! SU and A_OSS terms: ( tau*rho*a.grad(Ni) ) * ( rho*a.grad(Nj) + P(-rho*a.grad(u) )
          !
          tau   = gptau_time(igaus)                         ! tau'
          fact4 = tau * fact3                               ! tau'*rho*|dv|
          fact5 = gpden(igaus) * fact4                      ! tau'*rho^2*|dv|
          d     = gpden(igaus) * dtinv * c1sgs * sgs(igaus) ! rho*u'n/dt

          do inode = 1,pnode_SU_A_OSS
             fact6 = fact5 * gpadv(inode)       ! tau'*rho^2*|dv|*a.grad(Ni)
             do jnode = 1,pnode
                elmat(inode,jnode) = elmat(inode,jnode) + fact6 * gpadv(jnode)
             end do
             elrhs(inode) = elrhs(inode) + fact4 * gpadv(inode) * ( d - gppr1(igaus) )
          end do
          !
          ! Time tracking of subgrid scale for A_OSS
          ! rho*du'/dt = -rho*a.grad(u)-P - tau'/tau*(-rho*a.grad(u)-P+d)
          !            = (tau'/tau-1)*rho*a.grad(u) + (tau'/tau-1)*P + - tau'/tau*d
          ! => The (rho*du'/dt,v) gives
          !
          ! LHS + ( (tau'/tau-1)*rho*a.grad(u) , v ) = RHS + ( -(tau'/tau-1)*P , v ) + ( tau'/tau*d , v )
          !
          if( ADR % kfl_time_sgs /= 0 ) then
             xmuit  = gptau_time(igaus) / gptau(igaus)
             alpha  = xmuit * gpvol(igaus)               ! tau'/tau
             alpha1 = ( xmuit - 1.0_rp ) * gpvol(igaus)  ! tau'/tau - 1
             do inode = 1,pnode
                do jnode = 1,pnode
                   elmat(inode,jnode) = elmat(inode,jnode) &
                        + alpha1 * gpsha(inode,igaus) * gpden(igaus) * gpadv(jnode) 
                end do
                elrhs(inode) = elrhs(inode) + gpsha(inode,igaus) &
                     * ( alpha * d - alpha1 * gppr1(igaus) )
             end do
          end if
          !
          ! AR_OSS terms: ( tau*r*v ) * ( -f+r*u + P(f-r*u) )
          !
          do inode = 1,pnode_AR_OSS
             ! use gppr2 which is P(f-r*u)
             call runend('NOT CODED YET: MATIAS, PONTE LAS PILAS')
          end do

       end do

    else if( kfl_stabilization == BUBBLE ) then 

       !----------------------------------------------------------------------
       !
       ! Bubble
       !
       ! Auu: Ni * [ Nj/dt + r*Nj + rho*a.grad(Nj) ] + k*grad(Ni).grad(Nj)
       ! Aue: Ni * [ Ne/dt + r*Ne + rho*a.grad(Ne) ] + k*grad(Ni).grad(Ne)
       ! Aeu: Ne * [ Ni/dt + r*Ni + rho*a.grad(Ni) ] + k*grad(Ni).grad(Ne)
       ! Aee: Ne * [ Ne/dt + r*Ne + rho*a.grad(Ne) ] + k*grad(Ne).grad(Ne)
       ! bu:  Ni * ( f + time )
       ! be:  Ne * ( f + time )
       !
       ! +-       -+ +-  -+   +-  -+
       ! | Auu Aue | | u  |   | b  |
       ! |         | |    | = |    | 
       ! | Aeu Aee | | ue |   | be |
       ! +-       -+ +-  -+   +-  -+
       !
       ! Upon condensation we obtain:
       !
       ! [ Auu - Aue.Aee^-1.Aeu ] u = b - Aue.Aee^-1 be
       ! 
       !----------------------------------------------------------------------
       !
       ! If bubble should be actualized
       !
       if( present(messa) ) then
          if( trim(messa) == 'UPDATE BUBBLE' ) update_bubble = .true.
       end if
       !
       ! Cartesian derivatives of bubble shape function
       !
       do igaus = 1,pgaus
          call elmgeo_jacobian_matrix(&
               ndime,pnode,elcod,gpder(1,1,igaus),&
               gpdet,xjaci)
          do idime = 1,ndime
             gpcar_bub(idime,igaus) = 0.0_rp
             do jdime = 1,ndime
                gpcar_bub(idime,igaus) = gpcar_bub(idime,igaus) &
                     + xjaci(jdime,idime) * gpder_bub(jdime,igaus)
             end do
          end do
       end do

       elmat_ue = 0.0_rp
       elmat_eu = 0.0_rp
       elmat_ee = 0.0_rp
       elrhs_ee = 0.0_rp

       do igaus = 1,pgaus
          !
          ! Assembly
          !
          fact2 = gpvol(igaus) * gpdif(igaus)              ! |dv| * k       
          fact3 = gpvol(igaus) * gpden(igaus)              ! |dv| * rho     
          fact4 = gpvol(igaus) * rhsit(igaus)              ! |dv| * ( f + time terms )       
          fact5 = gpvol(igaus) * react(igaus)              ! |dv| * r       
          fact6 = fact3 * dtinv * ADR % time_parameters(1) ! |dv| * rho / dt
          fact7 = bemol * fact3                            ! |dv| * rho * bemol 
          !
          ! Advections
          ! a.grad(Ni)
          ! a.grad(Ne)
          !
          gpadv(1:pnode) = 0.0_rp
          gpadv_bub      = 0.0_rp
          do idime = 1,ndime
             gpadv(1:pnode) = gpadv(1:pnode) + gpvel(idime,igaus) * gpcar(idime,1:pnode,igaus)           
             gpadv_bub      = gpadv_bub      + gpvel(idime,igaus) * gpcar_bub(idime,igaus)           
          end do
          !
          ! Galerkin operator
          ! Nj/dt + r*Nj + rho*a.grad(Nj)
          ! Ne/dt + r*Ne + rho*a.grad(Ne)
          !        
          galer(1:pnode) = (         fact6 + fact5 ) * gpsha(1:pnode,igaus) + fact3 * bemo1 * gpadv(1:pnode)
          galer_bub      = ( xtime * fact6 + fact5 ) * gpsha_bub(igaus)     + fact3 * bemo1 * gpadv_bub
          !
          ! Diffusion operator convection and Galerkin
          ! Auu: k * grad(Ni).grad(Nj) - b * rho * a.grad(Ni).Nj + Ni * Galerkin(Nj) 
          ! Aue: k * grad(Ne).grad(Nj) - b * rho * a.grad(Ni).Ne + Ni * Galerkin(Ne)  
          ! Aeu: k * grad(Ne).grad(Ne) - b * rho * a.grad(Ne).Ni + Ne * Galerkin(Nj) 
          !        
          do inode = 1,pnode
             do jnode = 1,pnode
                xmuit =  dot_product(gpcar(1:ndime,inode,igaus),gpcar(1:ndime,jnode,igaus))
                elmat(inode,jnode) = elmat(inode,jnode)          &
                     + fact2 * xmuit                             & ! Diffusion
                     + gpsha(inode,igaus) * galer(jnode)         & ! Galerkin
                     - fact7 * gpsha(jnode,igaus) * gpadv(inode)   ! Convection integrated by parts
             end do
             xmuit             = fact2 * dot_product(gpcar(1:ndime,inode,igaus),gpcar_bub(1:ndime,igaus))
             elrhs(inode)      = elrhs(inode)      + gpsha(inode,igaus) * fact4
             elmat_ue(inode,1) = elmat_ue(inode,1) + xmuit + gpsha(inode,igaus) * galer_bub    - fact7 * gpsha_bub(igaus)   * gpadv(inode)
             elmat_eu(1,inode) = elmat_eu(1,inode) + xmuit + gpsha_bub(igaus)   * galer(inode) - fact7 * gpsha(inode,igaus) * gpadv_bub
          end do
          elrhs_ee(1)   = elrhs_ee(1) + gpsha_bub(igaus) * fact4 
          xmuit         = dot_product(gpcar_bub(1:ndime,igaus),gpcar_bub(1:ndime,igaus))
          elmat_ee(1,1) = elmat_ee(1,1) + fact2 * xmuit + gpsha_bub(igaus) * galer_bub - fact7 * gpsha_bub(igaus) * gpadv_bub
       end do
       !
       ! Condensation
       ! [ Auu - Aue.Aee^-1.Aeu ] u = b - Aue.Aee^-1 be
       !
       do inode = 1,pnode
          do jnode = 1,pnode
             elmat(inode,jnode) = elmat(inode,jnode) - elmat_ue(inode,1) * elmat_eu(1,jnode) / elmat_ee(1,1)
          end do
          elrhs(inode) = elrhs(inode) - elmat_ue(inode,1) * elrhs_ee(1) / elmat_ee(1,1)
       end do
       !
       ! Update bubble if required
       !
       if( update_bubble ) then
          ADR % bubble(ielem,1) = elrhs_ee(1) 
          do inode = 1,pnode
             ADR % bubble(ielem,1) = ADR % bubble(ielem,1) - elmat_eu(1,inode) * elunk(inode,1)
          end do
          ADR % bubble(ielem,1) = ADR % bubble(ielem,1) / elmat_ee(1,1)
       end if

    else  

       !-------------------------------------------------------------------
       !
       ! Other stabilization methods
       !
       ! ASGS and FULL OSS
       !
       ! k*(grad(u),grad(v)) + (1-b)*(rho*a.grad(u),v) - b*(rho*a.grad(v),u)
       ! + (ru,v) + (rho*du/dt,v) + (tau'*L*(v),R(u)-P) = (f,v) 
       !
       ! with L*(v) = s*v - rho*a.grad(v)
       !      R(u)  = R1(u) + R2(u)
       !      R1(u) = f - rho*du/dt - rho*a.grad(u) - r*u
       !      R2(u) = div(k*grad(u))
       !      P(u)  = Pi( R1(u) + R2(u) ) 
       !
       ! Then:
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u)
       ! + (rho*a.grad(u),v) + (r*u,v) + (rho*du/dt,v) - (f,v) 
       ! + (tau'*L*(v),R1(u)) + (tau'*L*(v),R2(u)-P) = 0
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u) 
       ! + (tau'*L*(v)-v,R1(u)) + (tau'*L*(v),R2(u)-P) = 0
       ! 
       ! So that
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u) 
       ! +  (v-tau'*L*(v),-R1(u)) + (-tau'*L*(v),-R2(u)) = (tau'*L*(v),P)
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u) 
       ! +  (v-tau'*L*(v), rho*du/dt + rho*a.grad(u) + r*u) 
       ! +  ( -tau'*L*(v),-div(k*grad(u))) =  (v-tau'*L*(v),f+time) 
       ! +  ( -tau'*L*(v),-P)
       !
       ! Let 
       !
       ! p1(v) = v-tau'*L*(v) = (1-tau'*s)*v + tau'*rho*a.grad(v)
       ! p2(v) =  -tau'*L*(v) =    -tau'*s*v + tau'*rho*a.grad(v)
       ! Beware the contribution to p1 from the viscous term is not present
       ! Therefore this is not ASGS!!!!
       !
       ! r1(u) = rho*du/dt + rho*a.grad(u) + r*u
       ! r2(u) = -div(k*grad(u))
       !
       ! Final formulation:
       ! ------------------
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u) 
       ! + (p1(v),r1(u)) + (p2(v),r2(u)) = (p1(v),f+time) + (p2(v),-P)
       !
       ! Adding time tracking of SGS: (without time tracking tau' = tau)
       ! ----------------------------
       !
       ! LHS + (rho*du'/dt,v) = RHS
       ! We have
       ! 
       !     rho*du'/dt + u'/tau = R(u)-P(u)       
       !
       ! =>  rho*du'/dt = R(u)-P(u) -  u'/tau           (1)
       !
       ! Discretizing in time:
       !
       !           u'-u'n    1
       !     rho * ------ + --- u' = R(u)-P(u)
       !             dt     tau
       !
       ! =>  u' = tau' * [ R(u)-P(u)+d ] with           (2)
       !
       !                 1                        u'n
       !     tau' = ---------    and    d = rho * ---  
       !            rho    1                       dt
       !            --- + ---
       !            dt    tau
       !         
       ! using (1) and (2), the additional term are:
       !
       !
       ! (rho*du'/dt,v) = (  R(u)-P(u) - tau'/tau * ( R(u)-P(u)+d ) , v )
       !                = ( (1-a) * (R(u)-P(u)) - a*d , v )
       ! and 
       !
       ! LHS = RHS - (d,tau'*L*(v) ) = RHS + ( d, p2(v) )
       ! 
       ! with a = tau'/tau 
       !
       ! LHS + ( (1-a)*R(u) , v ) = RHS + ( a*d , v ) + ( (1-a)*P(u) , v )
       ! R(u) = f-r1(u)-r(2)
       ! 
       ! LHS + ( (a-1)*(r1(u)+r2(u)) , v ) = RHS + ( a*d , v ) + ( (1-a)*P(u) , v )
       !                                     + ( (a-1) * f , v )
       !
       ! Therefore:
       !
       ! k*(grad(u),grad(v))  -b*(rho*a.grad(u),v) - b*(rho*a.grad(v),u) 
       ! + ( p1(v) , r1(u) ) + ( p2(v) , r2(u) )
       ! + ( (a-1)*(r1(u)+r2(u)) , v )
       ! = ( p1(v) , f+time ) + ( p2(v) , d-P )
       !   + ( a*d , v ) + ( -(a-1)*P , v )+ ( (a-1) * f , v )
       ! = ( p1(v) , f+time ) + ( p2(v) , d-P )
       !   + ( a*d + (a-1)*(f-P) , v )
       !
       !-------------------------------------------------------------------

       do igaus = 1,pgaus
          !
          ! Calculus of residual resid and perturbation function Pi=gppre
          !
          ! RESI1 = r1(u) =  rho/dt*Nj + rho*a.grad(Nj) + r1*Nj
          ! RESI2 = r2(u) =  -grad(k).grad(Nj) - k*Lap(Nj)
          ! GPPE1 = p1(v) =  [  Ni*(1-tau*s) + tau*rho*a.grad(Ni) ] 
          ! GPPE2 = p2(v) =  [    -Ni*tau*s  + tau*rho*a.grad(Ni) ] = p1(v) - v * |dv|
          !
          tau   = gptau_time(igaus)
          fact1 = gpden(igaus) * dtinv
          d     = fact1 * sgs(igaus)  ! rho*u'n/dt

          do inode = 1,pnode
             resi1(inode) = fact1 * ADR % time_parameters(1) * gpsha(inode,igaus)
             gpad1(inode) = 0.0_rp
             grvgr        = 0.0_rp
             gplap        = 0.0_rp
             do idime = 1,ndime
                gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)           
                grvgr        = grvgr        + gpgrd(idime,igaus) * gpcar(idime,inode,igaus)          
                gplap        = gplap        + gphes(idime,inode,igaus)
             end do
             gpadv(inode) = gpad1(inode) * gpden(igaus)
             resi1(inode) = resi1(inode) + gpadv(inode) + react(igaus) * gpsha(inode,igaus)
             resi2(inode) = - grvgr - gplap * gpdif(igaus)
             gppe1(inode) = (  gpsha(inode,igaus)*(1.0_rp-tau*sreac(igaus)) + tau * gpadv(inode) ) * gpvol(igaus)   
             !gppe2(inode) = ( -gpsha(inode,igaus)**sreac(igaus) + gpadv(inode) ) * tau * gpvol(igaus)   
             gppe2(inode) = gppe1(inode) - gpsha(inode,igaus) * gpvol(igaus)
          end do
          ! 
          ! Diffusion term
          !
          fact2 = gpvol(igaus) * ( gpdif(igaus) + CD(igaus) )  ! Diffusion
          fact3 = SD(igaus) * gpvol(igaus)                     ! Streamline negative diffusion
          do inode = 1,pnode
             do jnode = 1,inode-1
                xmuit = 0.0_rp
                do idime = 1,ndime
                   xmuit = xmuit + gpcar(idime,jnode,igaus) * gpcar(idime,inode,igaus)
                end do
                xmuit              = xmuit * fact2
                xmui3              = gpad1(inode) * gpad1(jnode) * fact3
                elmat(inode,jnode) = elmat(inode,jnode) + xmuit + xmui3
                elmat(jnode,inode) = elmat(jnode,inode) + xmuit + xmui3
             end do
             xmuit = 0.0_rp
             do idime = 1,ndime
                xmuit = xmuit + gpcar(idime,inode,igaus) * gpcar(idime,inode,igaus)
             end do
             elmat(inode,inode) = elmat(inode,inode) + xmuit * fact2 + gpad1(inode) * gpad1(inode) * fact3
          end do
          !
          ! bemol
          !          
          do inode = 1,pnode
             fact1 = gpsha(inode,igaus) * bemol * gpvol(igaus)*gpden(igaus)
             do jnode = 1,pnode
                fact2 = gpadv(jnode) * fact1
                elmat(inode,jnode) = elmat(inode,jnode) - fact2
                elmat(jnode,inode) = elmat(jnode,inode) - fact2
             end do
          end do
          !
          ! Assembly of the matrix and rhs
          !       
          ! GPPR1 <=> f - sum_i ri*u^i
          ! GPPR1 <=> - [ rho * a - grad(k) ] . grad(u) 
          ! GPPR1 <=> f - rho*(u - u^n)/dt - sum_i ri*u^i - [rho*u - grad(k)].grad(u) - k*lapl(u)
          !
          do inode = 1,pnode
             do jnode = 1,pnode
                elmat(inode,jnode) = elmat(inode,jnode) &
                     &             + resi1(jnode) * gppe1(inode) &        
                     &             + resi2(jnode) * gppe2(inode)
             end do
             elrhs(inode)         = elrhs(inode) &
                  &                + rhsit(igaus) * gppe1(inode) &
                  &                + ( d - gppr1(igaus) ) * gppe2(inode)
          end do
          !
          ! Skew symmetric advection
          !
          if (ADR % kfl_skewsymm ==1 ) then ! adds + 0.5*rho div(u)*scalar
             gpdiv = 0.0_rp
             ! Calculates velocity divergence
             do inode=1, pnode
                do idime=1, ndime
                   gpdiv = gpdiv + gpcar(idime, inode, igaus)*eladv(idime, inode)
                end do
             end do
             ! Adds the term + 0.5*rho div(u)*scalar
             do inode=1, pnode
                fact1 = gpsha(inode, igaus)*gpvol(igaus)*0.5_rp*gpden(igaus)*gpdiv
                do jnode =1, pnode
                   elmat(inode, jnode) = elmat(inode, jnode) + gpsha(jnode, igaus)*fact1
                end do
             end do
             ! Another way for skewsymm would be to use bemol =0.5, needing for boundary terms
          end if
          !
          ! Time tracking of subgrid scale
          !
          if( ADR % kfl_time_sgs /= 0 ) then
             xmuit  = gptau_time(igaus) / gptau(igaus)
             alpha  = xmuit * gpvol(igaus)               ! tau'/tau
             alpha1 = ( xmuit - 1.0_rp ) * gpvol(igaus)  ! tau'/tau - 1
             do inode = 1,pnode
                do jnode = 1,pnode
                   elmat(inode,jnode) = elmat(inode,jnode) &
                        + alpha1 * gpsha(inode,igaus) * ( resi1(jnode) + resi2(jnode) )
                end do
                elrhs(inode) = elrhs(inode) + gpsha(inode,igaus) &
                     * ( alpha1 * ( rhsit(igaus) - gppr1(igaus) ) + alpha * d )
             end do
          end if         
          ! 
          ! Lumped mass evolution matrix
          !
          if( kfl_time_lumped == 1 ) then
             do inode = 1,pnode
                fact1 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv
                elmat(inode,inode) = elmat(inode,inode) + fact1
                elrhs(inode) = elrhs(inode) - fact1 * gpunk(igaus,3)
                elrhs(inode) = elrhs(inode) + fact1 * elunk(inode,2)
                do jnode =1, pnode
                   elmat(inode, jnode) = elmat(inode, jnode) - fact1 * gpsha(jnode,igaus)
                end do
             end do
          end if
       end do

    end if

  end subroutine ADR_element_assembly

  !-------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   residuals at Gauss points
  !> @details Compute residuals for orthogonal projection
  !>          Reaction   = f - r*u
  !>          Convection = - [ rho*a - grad(k) ].grad(u) 
  !>          Residual   = f - rho*(u - u^n)/dt - r*u 
  !>                       - [ rho*a - grad(k) ].grad(u) + k*Lapl(u)
  !
  !-------------------------------------------------------------------

  subroutine ADR_compute_residual_projections(&
       ADR,ndime1,ntens1,mnode1,pnode,pgaus,elcod,gpsha,gpcar,gphes,&
       cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpunk,elunk,&
       reaction,convection,residual)

    type(ADR_typ), intent(in)    :: ADR                           !< ADR type

    integer(ip),   intent(in)    :: ndime1                        !< Dimension
    integer(ip),   intent(in)    :: ntens1                        !< Size of Hessian matrix
    integer(ip),   intent(in)    :: mnode1                        !< Maximum element nodes 
    integer(ip),   intent(in)    :: pnode                         !< # element nodes
    integer(ip),   intent(in)    :: pgaus                         !< # Gauss points
    real(rp),      intent(in)    :: elcod(ndime,pnode)            !< Element coordinates
    real(rp),      intent(in)    :: gpsha(pnode,pgaus)            !< Shape function
    real(rp),      intent(in)    :: gpcar(ndime1,mnode1,pgaus)    !< Shape function Cartesian derivatives
    real(rp),      intent(in)    :: gphes(ntens1,mnode1,pgaus)    !< Shape function Hessian
    !
    ! Equation coefficients
    !
    real(rp),    intent(in)      :: cutim                         !< Current time
    real(rp),    intent(in)      :: gpden(pgaus)                  !< Density
    real(rp),    intent(in)      :: gpvel(ndime1,pgaus)           !< Advection vector
    real(rp),    intent(in)      :: gpdif(pgaus)                  !< Diffusion 
    real(rp),    intent(in)      :: gpgrd(ndime1,pgaus)           !< Diffusion gradient
    real(rp),    intent(in)      :: gprea(pgaus,4)                !< Reaction
    real(rp),    intent(in)      :: gprhs(pgaus)                  !< RHS
    real(rp),    intent(in)      :: gpunk(pgaus,*)                !< Unknown at Gauss point
    real(rp),    intent(in)      :: elunk(pnode,*)                !< Unknown at element nodes
    !
    ! Projection residuals
    !
    real(rp),    intent(out)     :: reaction(pgaus)
    real(rp),    intent(out)     :: convection(pgaus)
    real(rp),    intent(out)     :: residual(pgaus)

    integer(ip)                  :: igaus,inode,idime
    integer(ip)                  :: itime,ii
    real(rp)                     :: rhsit,grunk(3),fact1
    real(rp)                     :: gplap
    !
    ! Local
    !
    real(rp)                     :: exunk(1),exgra(1,1)
    !
    ! Time parameters: du/dt = 1/dt * sum_i p(i) * u^{n+2-i}
    !
    do igaus = 1,pgaus
       !
       ! Convection + diffusion gradient: - [ rho*a - grad(k) ].grad(u)
       !
       convection(igaus) = 0.0_rp
       do idime = 1,ndime1
          grunk(idime) = 0.0_rp
          do inode = 1,pnode
             grunk(idime) = grunk(idime) + elunk(inode,1) * gpcar(idime,inode,igaus)
          end do
          convection(igaus) = convection(igaus) &
               - ( gpden(igaus) * gpvel(idime,igaus) - gpgrd(idime,igaus) ) * grunk(idime)                      
       end do
       !
       ! Time residual
       !
       fact1 = gpden(igaus) * ADR % dtinv
       rhsit = gprhs( igaus) 
       do itime = 1,ADR % ntime
          rhsit = rhsit - fact1 * ADR % time_parameters(itime) * gpunk(igaus,itime)
       end do
       !
       ! Total residual: f - time  - u/dt - sum_i ri*u^i - [ rho*a - grad(k) ].grad(u) 
       !
       reaction(igaus) = gprea(igaus,1) * gpunk(igaus,1)
       do ii = 2,mreac_adr
          reaction(igaus) = reaction(igaus) + gprea(igaus,ii) * gpunk(igaus,1) ** ii
       end do
       residual(igaus) = rhsit - reaction(igaus) + convection(igaus)
       !
       ! Total residual <= Total residual +  k*Lapl(u)
       !
       fact1 = 0.0_rp
       do inode = 1,pnode
          gplap = 0.0_rp
          do idime = 1,ndime1 
             gplap = gplap + gphes(idime,inode,igaus)
          end do
          fact1 = fact1 + gplap * elunk(inode,1)
       end do
       residual(igaus) = residual(igaus) + gpdif(igaus) * fact1
       !
       ! Reactive term: f - sum_i ri*u^i
       !
       reaction(igaus) = gprhs(igaus) - reaction(igaus)

    end do
    !
    ! Manufactured solution
    !
    if( ADR % kfl_manufactured /= 0 ) then
       call ADR_manufactured_solution_and_rhs(&
            2_ip,ADR % kfl_manufactured,pnode,pgaus,cutim,elcod,exunk,&
            exgra,gpsha,gpden,gpdif,gprea,gpgrd,gpvel,residual)
    end if

  end subroutine ADR_compute_residual_projections

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and Matias Avila
  !> @brief   Add SGS or Bubble to unknown
  !> @details u <= u + u'
  !
  !-----------------------------------------------------------------------

  subroutine ADR_add_sgs_or_bubble(ielem,pgaus,gpsha_bub,ADR,gpunk)

    integer(ip),   intent(in)    :: ielem            !< Current element
    integer(ip),   intent(in)    :: pgaus            !< Number of Gauss points
    real(rp),      intent(in)    :: gpsha_bub(pgaus) !< Bubble shape function
    type(ADR_typ), intent(in)    :: ADR              !< ADR type 
    real(rp),      intent(inout) :: gpunk(pgaus,*)   !< Unknown

    if( ADR % kfl_stabilization == BUBBLE ) then
       !
       !      n
       ! u = sum Ni*ui + Ne*ue 
       !     i=1
       !
       gpunk(1:pgaus,1) = gpunk(1:pgaus,1) + ADR % bubble(ielem,1) * gpsha_bub(1:pgaus)

    else if( ADR % kfl_nonlinear_sgs /= 0 ) then
       !
       !      n
       ! u = sum Ni*ui + u'
       !     i=1
       !
       gpunk(1:pgaus,1) = gpunk(1:pgaus,1) + ADR % sgs(ielem) % a(1,1:pgaus,1)

    end if

  end subroutine ADR_add_sgs_or_bubble

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and Matias Avila
  !> @brief   Critical time step
  !> @details Compute critical time step for ADR equation
  !
  !-----------------------------------------------------------------------

  subroutine ADR_critical_time_step(ADR,gpden,gpvel,gpdif,gprea,dtcri,h1,h2)
    type(ADR_typ), intent(in)           :: ADR           !< ADR type
    real(rp),      intent(in)           :: gpden(*)      !< Density
    real(rp),      intent(in)           :: gpvel(ndime)  !< Velocity
    real(rp),      intent(in)           :: gpdif(*)      !< Diffusion
    real(rp),      intent(in)           :: gprea(*)      !< Reaction
    real(rp),      intent(out)          :: dtcri(*)      !< Time step
    real(rp),      intent(in), optional :: h1            !< Length 1
    real(rp),      intent(in), optional :: h2            !< Length 2

    integer(ip)    :: tau_stra

    !time step calculated using stabilization tau_strategy, it should not depend on time step size
    if (ADR%kfl_tau_strategy.le.4 ) then
       tau_stra = ADR%kfl_tau_strategy
    else  ! use codina strategy when tau strategy depends on time step.
       tau_stra = TAU_CODINA
    end if
    
    if( present(h1) .and. present(h2) ) then 
       call ADR_tau(1_ip,tau_stra, ADR % tau_parameters,gpden,gpvel,gpdif,gprea,h1,h2,dtcri)
    else
       call runend('ADR_critical_time_step: NOT CODED')
    end if
    !
    ! Renormalize
    !
    dtcri(1) = gpden(1) * dtcri(1)

  end subroutine ADR_critical_time_step

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and Matias Avila
  !> @brief   Stabilization parameter
  !> @details Stabilization parameter for the advection-diffusion-reaction 
  !>          equation:
  !>          \verbatim
  !>          -k*Lapl(u) + a.grad(u) + s*u = Q 
  !>          \endverbatim
  !>           Accordining to the possible following strategies:
  !>           IMETH = 1 ... Codina
  !>                 = 2 ... Average of exact 1D equation with constant residual 
  !>                 = 3 ... Shakib
  !>                 = 4 ... Directional (not implemented)
  !>                 = 5 ... Codina with time step ( 1/tau = 1/tau_codina + rho/dt) 
  !>                 = 6 ... Time step  (tau = dt/rho)
  !
  !-----------------------------------------------------------------------

  subroutine ADR_tau(pgaus,kfl_tau_strategy,tau_parameters,gpden,gpvel,gpdif,react,h1,h2,gptau, dtinv)

    integer(ip), intent(in)  :: pgaus
    integer(ip), intent(in)  :: kfl_tau_strategy
    real(rp),    intent(in)  :: tau_parameters(3)
    real(rp),    intent(in)  :: gpden(pgaus)
    real(rp),    intent(in)  :: gpvel(ndime,pgaus)
    real(rp),    intent(in)  :: gpdif(pgaus)
    real(rp),    intent(in)  :: react(pgaus)
    real(rp),    intent(in)  :: h1
    real(rp),    intent(in)  :: h2
    real(rp),    intent(out) :: gptau(pgaus)
    real(rp),    intent(in),optional  :: dtinv
    integer(ip)              :: igaus
    real(rp)                 :: Dah,Dainv,Peh,Kh,Ah,alpha,PehDah
    real(rp)                 :: freq1,freq2,freq3,k,a,s,tau,gpnve

    if( h1 == 0.0_rp .and. h2 == 0.0_rp ) then
       gptau = 0.0_rp
       return
    end if

    do igaus = 1,pgaus

       gpnve = dot_product(gpvel(1:ndime,igaus),gpvel(1:ndime,igaus))
       gpnve = gpden(igaus) * sqrt(gpnve+zeror) 

       k     = tau_parameters(1) * gpdif(igaus) ! Diffusion 
       a     = tau_parameters(2) * gpnve        ! Advection
       s     = tau_parameters(3) * react(igaus) ! Reaction

       select case ( kfl_tau_strategy )

       case ( 0_ip )
          !
          ! No stabilization
          !
          tau = 0.0_rp

       case ( TAU_CODINA )
          !
          ! Codina
          !
          freq1 = 4.0_rp*k/(h2*h2)
          freq2 = 2.0_rp*a/h1
          freq3 = abs(s)
          tau   = freq1+freq2+freq3
          if(tau/=0.0_rp) tau=1.0_rp/tau

       case ( TAU_EXACT_1D )
          !
          ! Average of exact 1D equation with constant residual 
          !
          if(a/=0.0_rp.and.k/=0.0_rp.and.s==0.0_rp) then          ! AD
             Peh = a*h1/(2.0_rp*k)
             if(Peh<1.0e-03_rp) then
                alpha = Peh/3.0_rp
                tau   = h1*h1/(12.0_rp*k)
             else if(Peh>1.0e3_rp) then
                alpha = 1.0_rp
                tau   = h1/(2.0_rp*a)
             else
                alpha = 1.0_rp/tanh(Peh)-1.0_rp/(Peh)
                tau   = h1/(2.0_rp*a)*alpha
             end if

          else if(a==0.0_rp.and.k/=0.0_rp.and.s==0.0_rp) then     ! D
             tau  = h2*h2/(12.0_rp*k)

          else if(a/=0.0_rp.and.k==0.0_rp.and.s/=0.0_rp) then     ! AR

             Dah  = s*h1/a
             if(Dah>1.0e3_rp) then
                tau = 1.0_rp/s
             else if(Dah<-1.0e3_rp) then
                tau = (2.0_rp*exp(-Dah)/Dah)/Dah
             else if(abs(Dah)<1e-3) then
                tau = h1/(2.0_rp*a)
             else
                Dainv = 1.0_rp/Dah
                tau   = h1/(2.0_rp*a)*(2.0_rp*Dainv*&
                     &  (1.0_rp+Dainv*(exp(-Dah)-1.0_rp)))    
             end if

          else if(a/=0.0_rp.and.k==0.0_rp.and.s==0.0_rp) then     ! A
             tau   = h1/(2.0_rp*a)

          else if(a==0.0_rp.and.k==0.0_rp.and.s/=0.0_rp) then     ! R
             tau   = 1.0_rp/s

          else if(a==0.0_rp.and.k/=0.0_rp.and.s/=0.0_rp) then     ! DR
             Ah    = sqrt(s/k)*h2
             tau   = 2.0_rp*h2*(1.0_rp-cosh(Ah*h2))&
                  &  /(Ah*s*sinh(Ah*h2))
             Kh    = sqrt(2.0_rp*s*h2*h2/(2.0_rp*k))
             tau   = 1.0_rp/s*(1.0_rp+2.0_rp*(1.0_rp-cosh(Kh))&
                  &  /(Kh*sinh(Kh)))

          else if(a/=0.0_rp.and.k/=0.0_rp.and.s/=0.0_rp) then     ! ADR

             Peh    = a*h1/(2.0_rp*k)
             Dah    = s*h1/a
             PehDah = s*h1*h1/(2.0_rp*k)

             if(Dah<1.0e-6_rp) then
                alpha = 1.0_rp/tanh(Peh)-1.0_rp/(Peh)
                tau   = h1/(2.0_rp*a)*alpha
             else if(Dah<1.0e-6_rp*Peh) then
                freq2 = ( exp(-Dah)+exp(-2.0_rp*Peh-Dah) ) / (1.0_rp-exp(-2.0_rp*(Peh+Dah))) 
                freq3 = (Peh+Dah)/PehDah
                tau   = 1.0_rp/s*( 1.0_rp + (Peh+Dah)/PehDah*(freq2 - 1.0_rp/tanh(Peh+Dah)) )
             else
                Ah    = sqrt(Peh*Peh+2.0_rp*PehDah)
                freq2 = ( exp(Peh-Ah)+exp(-Peh-Ah) ) / (1.0_rp-exp(-2.0_rp*Ah)) 
                freq3 = Ah/PehDah
                tau   = 1.0_rp/s*( 1.0_rp + Ah/PehDah*(freq2 - 1.0_rp/tanh(Ah)) )
             end if

             if(s>0.0_rp) tau=max(0.0_rp,tau)

          else 
             tau   = 0.0_rp

          end if

       case( TAU_SHAKIB )
          !
          ! Shakib
          !
          tau = 9.0_rp*(4.0_rp*k/(h2*h2))**2+(2.0_rp*a/h1)**2+s*s
          if(tau>0.0_rp) tau=1.0_rp/sqrt(tau)

       case ( TAU_TEST )
          !
          ! Codina
          !
          freq1 = 4.0_rp*k/(h2*h2)
          freq2 = 2.0_rp*a/h1
          freq3 = abs(s)
          tau   = freq1+freq2+freq3
          if(tau/=0.0_rp) tau = 1.0_rp / ( 1.0_rp / h1 + abs(s) )
       case( TAU_INCLU_DT)
          !
          ! Include Dt in Tau - the rest idem Codina 
          ! 
          if( present(dtinv) ) then
             freq1 = 4.0_rp*k/(h2*h2)
             freq2 = 2.0_rp*a/h1
             freq3 = abs(s)
             tau   = freq1+freq2+freq3+dtinv*gpden(igaus)
             if(tau/=0.0_rp) tau=1.0_rp/tau
          else
             call runend('ADR_tau:.not.present(dtinv)') 
          end if
          
       case( TAU_DT)
          !
          ! Include Dt in Tau - the rest idem Codina 
          ! 
          if( present(dtinv) ) then
             tau=1.0_rp/(gpden(igaus)*dtinv)
          else
             call runend('modADR:.not.present(dtinv)') 
          end if

       case default

          call runend('ADR_tau: tau not coded')

       end select

       gptau(igaus) = tau

    end do

  end subroutine ADR_tau

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Time coefficients
  !> @details Initialize projections according to stabilization method.
  !>          \verbatim
  !>          
  !>          du         n+1         n         n-1
  !>          -- = p(1)*u    + p(2)*u  + p(3)*u     + ...
  !>          dt
  !>
  !>          This routine sets the coefficients for:
  !>
  !>          Trapezoidal rule:
  !>          
  !>          k |    n+1     n 
  !>          --+--------------
  !>          1 |      1    -1
  !>
  !>          k |  n+1/2     n 
  !>          --+--------------
  !>          2 |      2    -2
  !>          
  !>
  !>          Backward Differentation Formula (BDF):
  !>     
  !>          k |    n+1     n   n-1    n-2   n-3   n-4  n-5
  !>          --+-------------------------------------------
  !>          1 |      1    -1
  !>          2 |    3/2    -2   1/2
  !>          3 |   11/6    -3   3/2   -1/3
  !>          4 |  25/12    -4     3   -4/3   1/4
  !>          5 | 137/60    -5     5  -10/3   5/4  -1/5
  !>          6 | 147/60    -6  15/2  -20/3  15/4  -6/5  1/6
  !>          --+-------------------------------------------
  !>
  !>          \endverbatim
  !
  !-----------------------------------------------------------------------

  subroutine ADR_time_coefficients(kfl_time_scheme,kfl_time_order,time_parameters,dtinv,dtinv_old)

    integer(ip), intent(in)  :: kfl_time_scheme
    integer(ip), intent(in)  :: kfl_time_order
    real(rp),    intent(out) :: time_parameters(*)
    real(rp),    optional    :: dtinv
    real(rp),    optional    :: dtinv_old(*)
    real(rp)                 :: ratio

    if( present(dtinv) ) then
       if( dtinv < zeror ) then
          time_parameters(1) = 0.0_rp
          return
       end if
    end if

    if( kfl_time_scheme == TRAPEZOIDAL ) then

       !-----------------------------------------------------------------
       !
       ! Trapezoidal rule
       !
       !-----------------------------------------------------------------

       if( kfl_time_order == 1 ) then
          !
          ! 1/dt*( u^n - u^{n-1} ) 
          !
          time_parameters(1) =  1.0_rp
          time_parameters(2) = -1.0_rp

       else if( kfl_time_order == 2 ) then
          !
          ! ( u^n - u^{n-1} ) 
          ! -----------------
          !     1/2 * dt
          !
          time_parameters(1) =  2.0_rp
          time_parameters(2) = -2.0_rp

       end if

    else

       !-----------------------------------------------------------------
       !
       ! BDF scheme
       !
       !-----------------------------------------------------------------

       if(kfl_time_scheme == BDF) then

          if( kfl_time_order == 1 ) then
             !
             ! 1/dt*( u^n - u^{n-1} )
             !
             time_parameters(1) =  1.0_rp
             time_parameters(2) = -1.0_rp

          else if( kfl_time_order == 2 ) then

             ratio = dtinv_old(1) / dtinv
             if ( abs(ratio - 1.0_rp) < zeror ) then 
                !
                ! 1/(2*dt)*( 3*u^n - 4*u^{n-1} + u^{n-2} )
                !
                time_parameters(1) =  3.0_rp/2.0_rp
                time_parameters(2) = -4.0_rp/2.0_rp
                time_parameters(3) =  1.0_rp/2.0_rp

             else
                time_parameters(1) =  (1.0_rp + 2.0_rp*ratio) / (1.0_rp+ratio) 
                time_parameters(2) = -(1.0_rp+ratio) 
                time_parameters(3) =  (ratio*ratio) / (1.0_rp+ratio)
             end if

          else if( kfl_time_order == 3 ) then
             !
             ! 1/(6*dt)*( 11* u^n - 18*u^{n-1} + 9*u^{n-2} - 2*u^{n-3} )
             !
             call runend('ADR_time_coefficients: BDF3 FOR VARIABLE TIME STEP NOT CODED')

             time_parameters(1) =  11.0_rp/6.0_rp
             time_parameters(2) = -18.0_rp/6.0_rp
             time_parameters(3) =   9.0_rp/6.0_rp
             time_parameters(4) =  -2.0_rp/6.0_rp

          else if( kfl_time_order == 4 ) then
             !
             ! 1/(12*dt)*( 25*u^n - 48*u^{n-1} + 36*u^{n-2} - 16*u^{n-3} + 3*u^{n-4} )
             !
             time_parameters(1) =  25.0_rp/12.0_rp
             time_parameters(2) = -48.0_rp/12.0_rp
             time_parameters(3) =  36.0_rp/12.0_rp
             time_parameters(4) = -16.0_rp/12.0_rp
             time_parameters(5) =   3.0_rp/12.0_rp

          else if( kfl_time_order == 5 ) then
             !
             ! 1/(60*dt)*( 137*u^n - 300*u^{n-1} + 300*u^{n-2} - 200*u^{n-3} + 75*u^{n-4}
             !             -12*u^{n-5} )
             !
             time_parameters(1) =  137.0_rp/60.0_rp
             time_parameters(2) = -300.0_rp/60.0_rp
             time_parameters(3) =  300.0_rp/60.0_rp
             time_parameters(4) = -200.0_rp/60.0_rp
             time_parameters(5) =   75.0_rp/60.0_rp
             time_parameters(6) =  -12.0_rp/60.0_rp

          else if( kfl_time_order == 6 ) then
             !
             ! 1/(12*dt)*( 147*u^n - 360*u^{n-1} + 450*u^{n-2} - 400*u^{n-3} + 225*u^{n-4} 
             !             -72*u^{n-5} + 10*u^{n-6} ) 
             !
             time_parameters(1) =  147.0_rp/60.0_rp
             time_parameters(2) = -360.0_rp/60.0_rp
             time_parameters(3) =  450.0_rp/60.0_rp
             time_parameters(4) = -400.0_rp/60.0_rp
             time_parameters(5) =  225.0_rp/60.0_rp
             time_parameters(6) =  -72.0_rp/60.0_rp
             time_parameters(7) =   10.0_rp/60.0_rp

          end if
       end if

    end if

  end subroutine ADR_time_coefficients

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @date    5/10/2015
  !> @brief   Update bubble
  !> @details Update bubble
  !
  !-----------------------------------------------------------------------

  subroutine ADR_bubble_assembly(&
       ielem,pnode,pgaus,&
       elcod,gpsha,gpcar,gpder,gphes,gpvol,chale,gpsha_bub,gpder_bub,ADR,&
       cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpunk,&
       elunk,elmat,elrhs)
    !
    ! Element dimensions
    !
    integer(ip), intent(in)    :: ielem                         !< Element number
    integer(ip), intent(in)    :: pnode                         !< # nodes
    integer(ip), intent(in)    :: pgaus                         !< # Gauss points
    !
    ! Element characteristics at Gauss point
    !
    real(rp),    intent(in)    :: elcod(ndime,pnode)            !< Element node coordinates
    real(rp),    intent(in)    :: gpsha(pnode,pgaus)            !< Shape function Nk
    real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)      !< Shape function Cartesian derivatives dNk/dxi
    real(rp),    intent(in)    :: gpder(ndime,pnode,pgaus)      !< Shape function derivatives DNk/dsi
    real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)      !< Hessian dNk/dxidxj
    real(rp),    intent(in)    :: gpvol(pgaus)                  !< Element Jacobian
    real(rp),    intent(in)    :: chale(2)                      !< Element characteristic length
    real(rp),    intent(in)    :: gpsha_bub(pgaus)              !< Bubble shape function Ne
    real(rp),    intent(in)    :: gpder_bub(ndime,pgaus)        !< Bubble shape function derivatives dNe/dxi
    ! 
    ! Numerical strategy
    !
    type(ADR_typ), intent(inout)  :: ADR                           !< ADR type 
    !
    ! Equation coefficients
    !
    real(rp),    intent(in)    :: cutim                         !< Current time
    real(rp),    intent(in)    :: gpden(pgaus)                  !< Density
    real(rp),    intent(in)    :: gpvel(ndime,pgaus)            !< Advection vector
    real(rp),    intent(in)    :: gpdif(pgaus)                  !< Diffusion 
    real(rp),    intent(in)    :: gpgrd(ndime,pgaus)            !< Diffusion gradient
    real(rp),    intent(in)    :: gprea(pgaus,4)                !< Reaction
    real(rp),    intent(in)    :: gprhs(pgaus)                  !< RHS
    real(rp),    intent(in)    :: gpunk(pgaus,*)                !< Unknown at Gauss point
    real(rp),    intent(in)    :: elunk(pnode,2)                !< Element unknown
    !
    ! Output
    !
    real(rp),    intent(out)   :: elmat(pnode,pnode)            !< Element matrix
    real(rp),    intent(out)   :: elrhs(pnode)                  !< Element RHS

    call ADR_element_assembly(&
         ielem,pnode,pgaus,elcod,gpsha,gpcar,gpder,gphes,gpvol,chale,&
         gpsha_bub,gpder_bub,ADR,&
         cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,&
         gpunk,elunk,elmat,elrhs,'UPDATE BUBBLE')

  end subroutine ADR_bubble_assembly

  !-------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux 
  !> @brief   Update projections
  !> @details Update projections
  !>
  !-------------------------------------------------------------------

  subroutine ADR_update_projections(ADR,vmass)

    type(ADR_typ), intent(inout) :: ADR       !< ADR type
    real(rp),      intent(in) :: vmass(:)  !< mass matrix
    integer(ip)               :: ii

    if( ADR % kfl_stabilization == AR_OSS ) then

       call rhsmod(1_ip,ADR % proje1_tmp)
       call rhsmod(1_ip,ADR % proje2_tmp)
       do ii = 1,ADR % nunkn
          ADR % proje1(ii) = ADR % proje1_tmp(ii) / vmass(ii)
          ADR % proje2(ii) = ADR % proje2_tmp(ii) / vmass(ii)
       end do
       call memory_deallo(ADR % memor,'ADR % proje1_tmp','ADR_initialize_projections',ADR % proje2_tmp)

    else if( ADR % kfl_stabilization > 0 ) then

       call rhsmod(1_ip,ADR % proje1_tmp)
       do ii = 1,ADR % nunkn
          ADR % proje1(ii) = ADR % proje1_tmp(ii) / vmass(ii)
       end do
       call memory_deallo(ADR % memor,'ADR % proje1_tmp','ADR_initialize_projections',ADR % proje1_tmp)

    end if

  end subroutine ADR_update_projections

  !-------------------------------------------------------------------
  !>
  !> Compute and assemble residuals for orthogonal projection
  !>
  !>                +-
  !>          ri =  | r * Ni dv
  !>               -+
  !>               V
  !>
  !>          where r is:
  !>
  !>          1. Total residual for FULL OSS:
  !>             f - rho*(u - u^n)/dt - r*u - [ rho*a - grad(k) ].grad(u) + k*Lapl(u)
  !>          2. Convection residual for A_OSS (split OSS):
  !>             - [ rho*a - grad(k) ].grad(u)
  !>
  !-------------------------------------------------------------------

  subroutine ADR_projections_and_sgs_assembly(&
       ielem,pnode,pgaus,elcod,gpsha,gpcar,gphes,gpvol,chale,ADR,&
       cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpunk,elunk)

    integer(ip),   intent(in)            :: ielem                         !< # element number
    integer(ip),   intent(in)            :: pnode                         !< # element nodes
    integer(ip),   intent(in)            :: pgaus                         !< # Gauss points
    real(rp),      intent(in)            :: elcod(ndime,pnode)            !< Element coordinates
    real(rp),      intent(in)            :: gpsha(pnode,pgaus)            !< Shape function
    real(rp),      intent(in)            :: gpcar(ndime,mnode,pgaus)      !< Shape function Cartesian derivatives
    real(rp),      intent(in)            :: gphes(ntens,mnode,pgaus)      !< Shape function Hessian
    real(rp),      intent(in)            :: gpvol(pgaus)                  !< Element volume
    real(rp),      intent(in)            :: chale(2)                      !< Element characteristic length
    !
    ! Numerical strategy
    !
    type(ADR_typ)                        :: ADR                          !< ADR type
    !
    ! Equation coefficients
    ! 
    real(rp),      intent(in)            :: cutim                         !< Current time
    real(rp),      intent(in)            :: gpden(pgaus)                  !< Density
    real(rp),      intent(in)            :: gpvel(ndime,pgaus)            !< Advection vector
    real(rp),      intent(in)            :: gpdif(pgaus)                  !< Diffusion 
    real(rp),      intent(in)            :: gpgrd(ndime,pgaus)            !< Diffusion gradient
    real(rp),      intent(in)            :: gprea(pgaus,4)                !< Reaction
    real(rp),      intent(in)            :: gprhs(pgaus)                  !< RHS
    real(rp),      intent(in)            :: gpunk(pgaus,*)                !< Unknown at Gauss point
    real(rp),      intent(in)            :: elunk(pnode,*)                !< Unknown at element nodes
    !
    ! Local variables
    !
    real(rp)                             :: reaction(mgaus)             ! Reaction residual
    real(rp)                             :: convection(mgaus)           ! Convection residual
    real(rp)                             :: residual(mgaus)             ! Full residual

    integer(ip)                          :: igaus,inode,ipoin,itime
    real(rp)                             :: elre1(pnode)
    real(rp)                             :: elre2(pnode)
    real(rp)                             :: fact1,fact2
    real(rp)                             :: gppr1(mgaus)
    real(rp)                             :: gptau(mgaus)
    real(rp)                             :: gptau_time(mgaus)
    real(rp)                             :: sreac(mgaus)
    real(rp)                             :: sgs(mgaus)

    !--------------------------------------------------------------------
    !
    ! Residuals
    !
    !--------------------------------------------------------------------

    if( ADR % kfl_stabilization > 0 .or. ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) then

       call ADR_compute_residual_projections(& 
            ADR,ndime,ntens,mnode,pnode,pgaus,elcod,gpsha,gpcar,gphes,&
            cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpunk,elunk,&
            reaction,convection,residual) 

    end if

    !--------------------------------------------------------------------
    !
    ! Projections
    !
    !--------------------------------------------------------------------

    if( ADR % kfl_stabilization == AR_OSS ) then 
       !
       ! Assemble convection and reaction residuals
       !
       elre1(1:pnode) = 0.0_rp
       elre2(1:pnode) = 0.0_rp
       do igaus = 1,pgaus
          fact1 = convection(igaus) * gpvol(igaus)
          fact2 = reaction(igaus)   * gpvol(igaus)
          elre1(1:pnode) = elre1(1:pnode) + fact1 * gpsha(1:pnode,igaus)
          elre2(1:pnode) = elre2(1:pnode) + fact2 * gpsha(1:pnode,igaus)
       end do
       call assrhs(1_ip,pnode,lnods(1,ielem),elre1,ADR % proje1_tmp)             
       call assrhs(1_ip,pnode,lnods(1,ielem),elre2,ADR % proje2_tmp)      

    else if( ADR % kfl_stabilization == A_OSS ) then 
       !
       ! Assemble convection residual
       !
       elre1(1:pnode) = 0.0_rp
       do igaus = 1,pgaus
          fact1          = convection(igaus) * gpvol(igaus)
          elre1(1:pnode) = elre1(1:pnode) + fact1 * gpsha(1:pnode,igaus)
       end do
       call assrhs(1_ip,pnode,lnods(1,ielem),elre1,ADR % proje1_tmp)

    else if( ADR % kfl_stabilization == FULL_OSS ) then 
       !
       ! Assemble full residual
       !
       elre1(1:pnode) = 0.0_rp
       do igaus = 1,pgaus
          fact1          = residual(igaus) * gpvol(igaus)
          elre1(1:pnode) = elre1(1:pnode) + fact1 * gpsha(1:pnode,igaus)
       end do
       call assrhs(1_ip,pnode,lnods(1,ielem),elre1,ADR % proje1_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Compute SGS
    !
    ! Without time tracking:
    ! ----------------------
    !
    ! u' = tau* ( R(u) - P(u) )
    !
    !
    ! With time tracking:
    ! -------------------
    !
    !       u'-u'^n    1
    ! rho * ------- + --- u' = R(u) - P(u)
    !         dt      tau 
    !
    !             1 
    ! => u' = --------- [ R(u) - P(u) + rho/dt * u'^n ]
    !         rho    1 
    !         --- + ---
    !          dt   tau
    !
    !--------------------------------------------------------------------

    if( ADR % kfl_time_sgs /= 0 .or. ADR % kfl_nonlinear_sgs /= 0 ) then  

       gppr1(1:pgaus) = 0.0_rp
       if( ADR % kfl_stabilization == FULL_OSS .or. ADR % kfl_stabilization == A_OSS ) then 
          do igaus = 1,pgaus
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                gppr1(igaus) = gppr1(igaus) + ADR % proje1(ipoin) * gpsha(inode,igaus)
             end do
          end do
       end if

       if( ADR % kfl_linearization == RHS ) then
          !
          ! LHS + r1*u = RHS - r2*u_i^2
          !
          sreac(1:pgaus) = gprea(1:pgaus,1)

       else if( ADR % kfl_linearization == PICARD ) then
          !
          ! LHS + (r1+r2*u_i)*u = RHS 
          !
          sreac(1:pgaus) = gprea(1:pgaus,1) + gprea(1:pgaus,2) * gpunk(1:pgaus,1) 

       else if( ADR % kfl_linearization == NEWTON ) then
          !
          ! LHS + u * (r1+2*r2*u_i) = RHS + r2*u_i^2
          !
          sreac(1:pgaus) = gprea(1:pgaus,1) + gprea(1:pgaus,2) * gpunk(1:pgaus,1)

       end if

       call ADR_tau(&
            pgaus,ADR % kfl_tau_strategy,ADR % tau_parameters,gpden,gpvel,gpdif,sreac,&
            chale(1),chale(2),gptau, ADR%dtinv)   

       if( ADR % kfl_time_sgs /= 0 ) then          

          if( ADR % kfl_first_order_sgs == 1 ) then
             gptau_time(1:pgaus) = 1.0_rp / ( gpden(1:pgaus) * ADR % dtinv + 1.0_rp / gptau(1:pgaus) )
             sgs(1:pgaus) = ADR % sgs(ielem) % a(1,1:pgaus,3)
          else
             gptau_time(1:pgaus) = 1.0_rp / ( ADR % time_parameters(1) * gpden(1:pgaus) * ADR % dtinv + 1.0_rp / gptau(1:pgaus) )
             sgs(1:pgaus) = 0.0_rp
             do itime = 2,ADR % ntime
                sgs(1:pgaus) = sgs(1:pgaus) - ADR % time_parameters(itime) * ADR % sgs(ielem) % a(1,1:pgaus,itime+1)
             end do
          end if

          if( ADR % kfl_stabilization == A_OSS ) then 
             ADR % sgs(ielem) % a(1,1:pgaus,1) = &
                  gptau_time(1:pgaus) * ( gpden(1:pgaus) * ADR % dtinv * sgs(1:pgaus) + convection(1:pgaus) - gppr1(1:pgaus) ) 
          else             
             ADR % sgs(ielem) % a(1,1:pgaus,1) = &
                  gptau_time(1:pgaus) * ( gpden(1:pgaus) * ADR % dtinv * sgs(1:pgaus) + residual(1:pgaus)   - gppr1(1:pgaus) ) 
          end if

       else
          if( ADR % kfl_stabilization == A_OSS ) then 
             ADR % sgs(ielem) % a(1,1:pgaus,1) = gptau(1:pgaus) * ( convection(1:pgaus) - gppr1(1:pgaus) ) 
          else
             ADR % sgs(ielem) % a(1,1:pgaus,1) = gptau(1:pgaus) * ( residual(1:pgaus)   - gppr1(1:pgaus) ) 
          end if
       end if

    end if

  end subroutine ADR_projections_and_sgs_assembly

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux 
  !> @brief   Shcok capturing
  !> @details Shock capturing for the ADR equation
  !>
  !>          \verbatim
  !>               du
  !>           rho -- + rho*a.grad(u) - div[k*grad(u)] + r*u = f
  !>               dt
  !>          
  !>           k    = Diffusion                  [M/(L*T)]    
  !>           R    = Residual of the equation   [M*U/(L^3*T)]
  !>           C    = Shock capturing constant 
  !>           tau  = Stabilization parameter: 
  !>                  Its units are h/(rho*a)=   [L^3*T/M]
  !>                  so that rho*tau is in [T]
  !>           kiso = Isotropic SC diffusion     [M/(L*T)]
  !>           k'   = Anisotropic SC diffusion   [M/(L*T)]
  !>          
  !>                  1    |R|    h  
  !>           Pe   = - --------- - 
  !>                  2 |grad(u)| k 
  !>          
  !>                      +          2k   +           R
  !>           ac   = max | 0 , C - ----- | , a*= ----------- grad(u), therefore
  !>                      +         |a*|h +       |grad(u)|^2
  !>          
  !>                      +           1   +
  !>           ac   = max | 0 , C - ----- |
  !>                      +          Pe   +
  !>                  1          |R|
  !>           kiso = - ac*h  --------- , k'=(rho*tau)*rho*a^2
  !>                  2       |grad(u)|
  !>          
  !>           KFL_SHOCK == 1 
  !>           --------------
  !>           Isotropic:     kiso*grad(u).grad(v)  
  !>                                           
  !>           KFL_SHOCK == 2
  !>           --------------
  !>                                                   a x a
  !>           Anisotropic:   (<kiso-k'>-kiso)*grad(u).-----.grad(v)
  !>                                                    a^2
  !>          \endverbatim
  !>
  !------------------------------------------------------------------------

  subroutine ADR_shock_capturing(&
       ndime1,mnode1,ntens1,pgaus,pnode,kfl_laplacian,kfl_shock,time_parameters,shock,&
       dtinv,rhsit,gpden,gpvel,gpdif,gpgrd,react,gpunk,gptau,chale,gpcar,gphes,elunk,SD,CD)

    integer(ip), intent(in)  :: ndime1
    integer(ip), intent(in)  :: mnode1
    integer(ip), intent(in)  :: ntens1
    integer(ip), intent(in)  :: pgaus
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: kfl_laplacian
    integer(ip), intent(in)  :: kfl_shock
    real(rp),    intent(in)  :: time_parameters(*)
    real(rp),    intent(in)  :: shock
    real(rp),    intent(in)  :: dtinv
    real(rp),    intent(in)  :: rhsit(pgaus)
    real(rp),    intent(in)  :: gpden(pgaus)
    real(rp),    intent(in)  :: gpvel(ndime1,pgaus)
    real(rp),    intent(in)  :: gpdif(pgaus)
    real(rp),    intent(in)  :: gpgrd(ndime1,pgaus)
    real(rp),    intent(in)  :: react(pgaus)
    real(rp),    intent(in)  :: gpunk(pgaus,*)
    real(rp),    intent(in)  :: gptau(pgaus)
    real(rp),    intent(in)  :: chale(2)
    real(rp),    intent(in)  :: gpcar(ndime1,mnode1,pgaus)
    real(rp),    intent(in)  :: gphes(ntens1,mnode1,pgaus)
    real(rp),    intent(in)  :: elunk(pnode,*)
    real(rp),    intent(out) :: SD(pgaus)
    real(rp),    intent(out) :: CD(pgaus)
    integer(ip)              :: idime,inode,igaus
    real(rp)                 :: resid,factt,umbra
    real(rp)                 :: xfact,gplap,grnor
    real(rp)                 :: grunk(3),gpnve,ugrau
    real(rp)                 :: rhnv2,Dsupg,uscoe

!!$    real(rp) :: gpve2,vepar,gpres,cdv2,F1,F2,F3,facta
!!$    if( kfl_shock /= 0 ) then
!!$       factt = 1.0_rp                ! Other elements
!!$       factt = 0.5_rp*factt
!!$       umbra = 1.0e-6_rp
!!$       facta = real(kfl_shock-1_ip)  ! Isotropic/anisotropic shock capturing
!!$       facta = 1.0_rp
!!$
!!$        do igaus=1,pgaus                              ! Residual
!!$           F1       = -gpden(igaus) * dtinv - react(igaus)
!!$           F2       = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)
!!$           F3       = -gpden(igaus) * gpvel(2,igaus) + gpgrd(2,igaus)   
!!$
!!$           resid = rhsit(igaus) - ( xfact*gpden(igaus) + react(igaus) ) * gpunk(igaus,1) 
!!$          ugrau = 0.0_rp
!!$          grunk = 0.0_rp
!!$          grnor = 0.0_rp
!!$          do idime = 1,ndime1
!!$             do inode = 1,pnode
!!$                grunk(idime) = grunk(idime) + elunk(inode,1) * gpcar(idime,inode,igaus)
!!$             end do
!!$             ugrau = ugrau + ( gpden(igaus) * gpvel(idime,igaus) - gpgrd(idime,igaus) ) * grunk(idime)
!!$             grnor = grnor + grunk(idime) * grunk(idime)
!!$          end do
!!$          grnor = sqrt(grnor+zeror)
!!$          resid = resid - ugrau
!!$
!!$          if( kfl_laplacian > 0 ) then 
!!$             do inode = 1, pnode
!!$                gplap = 0.0_rp
!!$                do idime = 1,ndime1           
!!$                   gplap = gplap + gphes(idime,inode,igaus)
!!$                end do
!!$                resid = resid + elunk(inode,1) * gplap * gpdif(igaus)
!!$             end do
!!$          end if
!!$
!!$
!!$           gpve2 =   gpvel(1,igaus)*gpvel(1,igaus) &              ! a^2
!!$                &  + gpvel(2,igaus)*gpvel(2,igaus)
!!$           grnor = sqrt( grunk(1)*grunk(1) + grunk(2)*grunk(2) + zeror )  ! | grad(u) | 
!!$
!!$           if( grnor > umbra ) then
!!$              vepar     = resid / grnor
!!$              CD(igaus) = factt * shock * chale(2) * vepar - gpdif(igaus)
!!$              if( CD(igaus) > 0.0_rp .and. gpve2 > umbra ) then
!!$                 cdv2      = CD(igaus) / gpve2
!!$                 SD(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gptau(igaus)*gpden(igaus)) - cdv2
!!$                 SD(igaus) = facta * SD(igaus)
!!$              end if
!!$           end if
!!$
!!$        end do
!!$     end if
!!$
!!$     return

    if( kfl_shock /= 0 ) then

       factt = 0.75_rp
       umbra = 1.0e-6_rp
       xfact = dtinv * time_parameters(1)

       do igaus = 1,pgaus         
          !
          ! RESID = R(u)
          ! GRUNK = grad(u)
          ! GRNOR = |grad(u)|
          !
          resid = rhsit(igaus) - ( xfact*gpden(igaus) + react(igaus) ) * gpunk(igaus,1) 
          ugrau = 0.0_rp
          grunk = 0.0_rp
          grnor = 0.0_rp
          do idime = 1,ndime1
             do inode = 1,pnode
                grunk(idime) = grunk(idime) + elunk(inode,1) * gpcar(idime,inode,igaus)
             end do
             ugrau = ugrau + ( gpden(igaus) * gpvel(idime,igaus) - gpgrd(idime,igaus) ) * grunk(idime)
             grnor = grnor + grunk(idime) * grunk(idime)
          end do
          grnor = sqrt(grnor+zeror)
          resid = resid - ugrau

          if( kfl_laplacian > 0 ) then 
             do inode = 1, pnode
                gplap = 0.0_rp
                do idime = 1,ndime1           
                   gplap = gplap + gphes(idime,inode,igaus)
                end do
                resid = resid + elunk(inode,1) * gplap * gpdif(igaus)
             end do
          end if

          if( grnor > umbra ) then
             !
             ! |R(u)|/|grad(u)|
             !
             uscoe = abs(resid/grnor) 
             !
             ! Cross diffusion
             ! 
             CD(igaus) = max(factt*0.5_rp*shock*chale(2)*uscoe - gpdif(igaus),0.0_rp)  
             rhnv2     = dot_product(gpvel(1:ndime1,igaus),gpvel(1:ndime1,igaus))
             gpnve     = sqrt(rhnv2)
             !            
             ! Streamline diffusion tau*rho^2*u^2
             !
             if( gpnve > umbra ) then
                rhnv2     = gpden(igaus) * gpden(igaus) * rhnv2                
                Dsupg     = gptau(igaus)*rhnv2                                    ! Supg  
                SD(igaus) = max(CD(igaus)-Dsupg,0.0_rp) - CD(igaus)               ! Streamline diffusion
                SD(igaus) = -CD(igaus)                                            ! Only crosswind
                SD(igaus) = SD(igaus)/(gpnve*gpnve)   
                !SD(igaus) =  0.0_rp
             else
                SD(igaus) = 0.0_rp
             end if

          else 
             CD(igaus) = 0.0_rp
             SD(igaus) = 0.0_rp
          end if
       end do

    end if

  end subroutine ADR_shock_capturing

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    7/10/2015
  !> @brief   Compute manufactured nodal solution
  !> @details Compute manufactured nodal solution
  !
  !-----------------------------------------------------------------------

  subroutine ADR_manufactured_nodal_solution(ADR,cutim,unkno)

    type(ADR_typ), intent(in)  :: ADR       !< ADR type
    real(rp),      intent(in)  :: cutim
    real(rp),      intent(out) :: unkno(*)
    integer(ip)                :: ii
    real(rp)                   :: exunk(2)

    do ii = 1,ADR % nunkn
       call ADR_manufactured_solution_and_rhs(&
            3_ip,ADR % kfl_manufactured,1_ip,1_ip,cutim,coord(1:ndime,ii),exunk)
       unkno(ii) = exunk(1)
    end do

  end subroutine ADR_manufactured_nodal_solution

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    7/10/2015
  !> @brief   Compute nodal error wrt manufactured solution
  !> @details Compute:
  !>          ITASK = 1 ... manufactured solution and gradient
  !>                = 2 ... force term to be added to RHS 
  !
  !-----------------------------------------------------------------------

  subroutine ADR_manufactured_nodal_error(ADR,cutim,unkno,error)

    type(ADR_typ), intent(in)  :: ADR       !< ADR type
    real(rp),      intent(in)  :: cutim
    real(rp),      intent(in)  :: unkno(*)
    real(rp),      intent(out) :: error(*)
    integer(ip)                :: ii
    real(rp)                   :: exunk(2)

    do ii = 1,ADR % nunkn
       call ADR_manufactured_solution_and_rhs(&
            3_ip,ADR % kfl_manufactured,1_ip,1_ip,cutim,coord(1:ndime,ii),exunk)
       error(ii) = unkno(ii) - exunk(1)
    end do

  end subroutine ADR_manufactured_nodal_error

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    7/10/2015
  !> @brief   Compute force term or manufactured solution
  !> @details Compute:
  !>          ITASK = 1 ... manufactured solution and gradient
  !>                = 2 ... force term to be added to RHS 
  !>                = 3 ... nodal manufactured solution
  !
  !-----------------------------------------------------------------------

  subroutine ADR_manufactured_solution_and_rhs(&
       itask,kfl_manufactured,pnode,pgaus,cutim,elcod,exunk,&
       exgra,gpsha,gpden,gpdif,react,gpgrd,gpvel,rhsit)

    integer(ip), intent(in)              :: itask
    integer(ip), intent(in)              :: kfl_manufactured
    integer(ip), intent(in)              :: pnode
    integer(ip), intent(in)              :: pgaus
    real(rp),    intent(in)              :: cutim
    real(rp),    intent(in)              :: elcod(ndime,*)
    real(rp),    intent(out)             :: exunk(pgaus)
    real(rp),    intent(out),   optional :: exgra(ndime,pgaus)
    real(rp),    intent(in),    optional :: gpsha(pnode,pgaus)
    real(rp),    intent(in),    optional :: gpden(pgaus)
    real(rp),    intent(in),    optional :: gpdif(pgaus)
    real(rp),    intent(in),    optional :: react(pgaus)
    real(rp),    intent(in),    optional :: gpgrd(ndime,pgaus)
    real(rp),    intent(in),    optional :: gpvel(ndime,pgaus)
    real(rp),    intent(inout), optional :: rhsit(pgaus)
    integer(ip)                          :: inode,igaus
    real(rp)                             :: x,y,z,t,u,v,w,r,k,a,dtdx,dtdy,dtdz
    real(rp)                             :: d2tdy2,d2tdx2,d2tdz2,dtdt,Q,dkdx,dkdy,dkdz
    real(rp)                             :: rmayo,rmeno,f,sigma,x0,gpcod(3),rho
    real(rp)                             :: r2,r3,r4

    if( kfl_manufactured /= 0 ) then

       do igaus = 1,pgaus
          !
          ! Coordinates
          !
          if( itask /= 3 ) then
             gpcod = 0.0_rp
             do inode = 1,pnode
                gpcod(1:ndime) = gpcod(1:ndime) + gpsha(inode,igaus) * elcod(1:ndime,inode)
             end do
          else
             gpcod(1:ndime) = elcod(1:ndime,1)
          end if
          !
          ! Initializations
          ! 
          x = gpcod(1)
          y = 0.0_rp
          z = 0.0_rp
          if( ndime >= 2 ) y = gpcod(2)
          if( ndime == 3 ) z = gpcod(3)

          if( itask == 2 .and. present(rhsit) ) then
             dtdt    = 0.0_rp
             dtdx    = 0.0_rp
             dtdy    = 0.0_rp
             dtdz    = 0.0_rp
             d2tdx2  = 0.0_rp
             d2tdy2  = 0.0_rp
             d2tdz2  = 0.0_rp
             dkdx    = gpgrd(1,igaus)
             dkdy    = gpgrd(2,igaus)
             dkdz    = 0.0_rp
             rho     = gpden(igaus)
             u       = gpvel(1,igaus)
             v       = gpvel(2,igaus)
             w       = 0.0_rp
             k       = gpdif(igaus)
             r       = react(igaus)
             r2      = 0.0_rp
             r3      = 0.0_rp
             r4      = 0.0_rp
             if( ndime == 3 ) dkdy = gpgrd(3,igaus)
             if( ndime == 3 ) w    = gpvel(3,igaus)
          end if

          if( kfl_manufactured == 1 ) then
             !
             ! T=sin(pi*x)*sin(pi*y)*exp(x*y) in [0,1]x[0,1]
             !
             if( ndime == 2 ) then
                t      = sin(pi*x)*sin(pi*y)*exp(x*y)
                dtdx   = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
                dtdy   = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))
                d2tdx2 = sin(pi*y)*exp(x*y)*(2.0_rp*y*pi*cos(pi*x)&
                     &   -pi*pi*sin(pi*x)+y*y*sin(pi*x))
                d2tdy2 = sin(pi*x)*exp(x*y)*(2.0_rp*x*pi*cos(pi*y)&
                     &   -pi*pi*sin(pi*y)+x*x*sin(pi*y))
             else
                t      = sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
                dtdx   = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*x)+y*z*sin(pi*x))
                dtdy   = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))
                dtdz   = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))

                d2tdx2 = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(2.0_rp*y*z*pi*cos(pi*x)&
                     &   -pi*pi*sin(pi*x)+y*y*z*z*sin(pi*x))
                d2tdy2 = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(2.0_rp*x*z*pi*cos(pi*y)&
                     &   -pi*pi*sin(pi*y)+x*x*z*z*sin(pi*y))                 
                d2tdz2 = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(2.0_rp*x*y*pi*cos(pi*z)&
                     &   -pi*pi*sin(pi*z)+x*x*y*y*sin(pi*z)) 
             end if

          else if( kfl_manufactured == 2 ) then

             a = 1.0_rp
             k = 2.0_rp
             r = 3.0_rp
             Q = 4.0_rp
             !
             ! Solution of Advection-Diffusion-Reaction equation in [0,1]
             !
             if(r<=zeror .and. k>zeror .and. a>zeror) then
                !
                ! AD   :  T = 1/a*(exp(a/k)-1)*[(1-exp(a*x/k))+x/a]*Q
                !
                t=((1.0_rp-cosh(0.5_rp*a*x/k)-sinh(0.5_rp*x/k))/               &
                     &   (a*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))+x/a)*Q
                dtdx=-(0.5_rp*a/k)*(cosh(0.5_rp*a*x/k)+sinh(0.5_rp*a*x/k))  &
                     &    /(a*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))+1.0_rp/a
                d2tdx2=-0.25_rp*a*(cosh(0.5_rp*a*x/k)+sinh(0.5_rp*a*x/k))&
                     &    /(k*k*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))

             else if(r<=zeror .and. k>zeror .and. a<=zeror) then
                !
                ! D    : T = [ x/(2*k)*(1-x) ]*Q
                !
                t      = (x/(2.0_rp*k)*(1.0_rp-x))*Q
                dtdx   = 1.0_rp/(2.0_rp*k)-x/k 
                d2tdx2 = -1.0_rp/k

             else if(r<=zeror .and. k<=zeror .and. a>zeror) then
                !
                ! A   : T = [x/a]*Q
                !
                t      = (x/a)*Q
                dtdx   = 1.0_rp/a
                d2tdx2 = 0.0_rp

             else if(k<=zeror .and. r>zeror .and. a>zeror ) then
                !
                ! AR  : T = [1/r*(1-exp(-r*x/a))]*Q
                !
                t=((1.0_rp-cosh(r*x/a)+sinh(r*x/a))/r)*Q
                dtdx=(cosh(r*x/a)-sinh(r*x/a))/a
                d2tdx2=(sinh(r*x/a)-cosh(r*x/a))*r/a

             else if(k<=zeror .and. r>zeror .and. a<=zeror) then
                !
                ! R   : T = Q/s
                !
                t    = (1.0_rp/r)*Q
                dtdx = 0.0_rp
                d2tdx2=0.0_rp

             else if(a<=zeror .and. k>zeror .and. r>zeror) then
                !
                ! DR  : T = [1/(r*(exp(r/k)^(1/2)-exp(r/k)^(1/2))*(exp((r/k)^(1/2)*x)-exp(-(r/k)^(1/2)*x))-1/r]*Q
                !
                t=((sinh(sqrt(r/k)*x)/(r*sinh(sqrt(r/k)))-1.0_rp/r))*Q
                dtdx=(sqrt(r/k)*cosh(sqrt(r/k)*x)/(r*sinh(sqrt(r/k))))
                d2tdx2=sinh(sqrt(r/k)*x)/(k*sinh(sqrt(r/k)))

             else
                !
                ! ADR : T= 1/(s*a)*[(exp(r+)-1)/(exp(r-)-exp(r+)) * exp((r-)*x) -(exp(r-)-1)/(exp(r-)-exp(r+))*exp((r+)*x)+1]*Q
                !
                rmayo = a/(k*2.0_rp)+sqrt( (a/(k*2.0_rp))**2.0_rp + r/k )
                rmeno = a/(k*2.0_rp)-sqrt( (a/(k*2.0_rp))**2.0_rp + r/k )
                t=((((cosh(rmayo)+sinh(rmayo)-1.0_rp)*(cosh(rmeno*x)+sinh(rmeno*x))-       &
                     &  (cosh(rmeno)+sinh(rmeno)-1.0_rp)*(cosh(rmayo*x)+sinh(rmayo*x)))/       &
                     &  (cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo))+1.0_rp)*(1.0_rp/(r*a)))*Q
                dtdx=(((cosh(rmayo)+sinh(rmayo)-1.0_rp)*rmeno*(cosh(rmeno*x)+sinh(rmeno*x))-  &
                     &   (cosh(rmeno)+sinh(rmeno)-1.0_rp)*rmayo*(cosh(rmayo*x)+sinh(rmayo*x)))/    &
                     &  ( cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo)))*(1.0_rp/(r*a))
                d2tdx2=(1.0_rp/(r*a))*(((cosh(rmayo)+sinh(rmayo)-1.0_rp)*((cosh(rmeno*x)+sinh(rmeno*x))*rmeno*rmeno) &
                     &   -(cosh(rmeno)+sinh(rmeno)-1.0_rp)*((cosh(rmayo*x)+sinh(rmayo*x))*rmayo*rmayo))/      &
                     &   (cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo)))

             end if

          else if( kfl_manufactured == 3 ) then
             !
             ! T=2*x+3*y in [0,1]x[0,1]
             !
             t      = 2.0_rp*x+3.0_rp*y+4.0_rp*z
             dtdx   = 2.0_rp
             dtdy   = 3.0_rp
             if( ndime == 3 ) dtdz = 4.0_rp

          else if( kfl_manufactured == 4 ) then
             !
             ! T=2.0*exp(-(20.0*(x-0.5))**2)*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)   
             !
             t      =  2.0_rp*exp(-(20.0_rp*(x-0.5_rp))**2.0_rp)*sin(x)*sin(pi*x)*sin(y)*sin(pi*y) 

             dtdx   =  2.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)&
                  &   +2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*cos(x)*sin(pi*x)*sin(y)*sin(pi*y)    &
                  &   +2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)

             dtdy   =  2.0_rp*exp(-(20.0_rp*(x-0.5_rp))**2.0_rp)*sin(x)*sin(pi*x)&
                  &    *(cos(y)*sin(pi*y)+pi*sin(y)*cos(pi*y))

             d2tdx2 = -1602.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)                            &
                  &   +2.0_rp*((-800.0_rp*x+400.0_rp))**2 *exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)  &
                  &   +4.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5_rp))**2))*cos(x)*sin(pi*x)*sin(y)*sin(pi*y)        &
                  &   +4.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)     &
                  &   +4.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*cos(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)                            &
                  &   -2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*pi*pi*sin(y)*sin(pi*y)

             d2tdy2 = -2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)          &
                  &   +4.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))* sin(x)*sin(pi*x)*cos(y)*cos(pi*y)            &
                  &   *pi-2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*pi**2*sin(y)* sin(pi*y)

          else if( kfl_manufactured == 5 ) then
             !
             ! T=2x+3
             !
             t      = 2.0_rp*x*x+1.0_rp
             dtdx   = 4.0_rp*x

          else if( kfl_manufactured == 6 ) then
             !
             ! T=t*(2x+3y+4)
             !
             t      = cutim*(2.0_rp*x+3.0_rp*y+4.0_rp)
             dtdx   = cutim*2.0_rp
             dtdy   = cutim*3.0_rp
             dtdt   = 2.0_rp*x+3.0_rp*y+4.0_rp

          else if( kfl_manufactured == 5 ) then
             !
             ! T=2x^2+1
             !
             t      = 2.0_rp*x*x+1.0_rp
             dtdx   = 4.0_rp*x
             d2tdx2 = 4.0_rp

          else if( kfl_manufactured == 7 ) then
             !
             ! T=2x^2 + 3y^2 + 4z^2
             !
             t      = 2.0_rp*x*x+3.0_rp*y*y
             dtdx   = 4.0_rp*x
             d2tdx2 = 4.0_rp
             dtdy   = 6.0_rp*y
             d2tdy2 = 6.0_rp
             if( ndime == 3 ) then
                t      = t+4.0_rp*z*z
                dtdz   = 8.0_rp*z
                d2tdz2 = 8.0_rp
             end if

          else if( kfl_manufactured == 8 ) then
             !
             ! T=2x^3+3y^3
             !
             t      = 2.0_rp*x*x*x+3.0_rp*y*y*y
             dtdx   = 6.0_rp*x*x
             d2tdx2 = 12.0_rp*x
             dtdy   = 9.0_rp*y*y
             d2tdy2 = 18.0_rp*y

          else if( kfl_manufactured == 9 ) then
             !
             ! T=x^3+y^3
             !
             t      = x*x*x+y*y*y
             dtdx   = 3.0_rp*x*x
             d2tdx2 = 6.0_rp*x
             dtdy   = 3.0_rp*y*y
             d2tdy2 = 6.0_rp*y

          else if( kfl_manufactured == 10 ) then
             !
             ! T=x+2*y+3*z
             !
             if( ndime == 2 ) then
                t      = x+2.0_rp*y
                dtdx   = 1.0_rp
                dtdy   = 2.0_rp
             else
                t      = x+2.0_rp*y+3.0_rp*z
                dtdx   = 1.0_rp
                dtdy   = 2.0_rp
                dtdz   = 3.0_rp
             end if

          else if( kfl_manufactured == 11 ) then
             !
             ! T=sin*sin
             !
             f      =  4.0_rp
             t      =  sin(2.0_rp*pi*x*f)*sin(2.0_rp*pi*y*f)
             dtdx   =  2.0_rp*pi*f*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
             dtdy   =  2.0_rp*pi*f*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
             d2tdx2 = -4.0_rp*pi*pi*f*f*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
             d2tdy2 = -4.0_rp*pi*pi*f*f*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)

          else if( kfl_manufactured == 12 ) then
             !
             ! T=no forcing term
             !
             sigma  =  0.01_rp
             x0     =  0.05_rp
             f      =  (x-x0-1.0_rp*cutim)**2 / (  sigma * sigma + 4.0_rp*k*cutim ) 
             f      =  min(f,200.0_rp)
             t      =  1.0_rp / sqrt( 1.0_rp + 4.0_rp * k * cutim / (sigma*sigma) ) * exp(-f)

          else if( kfl_manufactured == 13 ) then
             !
             ! T=1
             !
             t      = 1.0_rp

          else if( kfl_manufactured == 14 ) then
             !
             ! T=x+2*y+3*z
             !
             t      = y
             dtdy   = 1.0_rp

          else if( kfl_manufactured == 15 ) then
             !
             ! T=x+2*y+3*z
             !
             f      =  4.0_rp
             t      =  cos(2.0_rp*pi*y*f)
             dtdy   = -sin(2.0_rp*pi*y*f)*2.0_rp*pi*f
             d2tdy2 = -cos(2.0_rp*pi*y*f)*2.0_rp*pi*f*2.0_rp*pi*f

          else if( kfl_manufactured == 16 ) then
             !
             ! T=2*x^2+3*y^2+4*z^2+5*x*y+6*x*z+8*y*z in [0,1]x[0,1]
             !
             t      = 2.0_rp*x*x + 3.0_rp*y*y + 4.0_rp*z*z + 5.0_rp*x*y + 6.0_rp*x*z + 7.0_rp*y*z
             dtdx   = 4.0_rp*x + 5.0_rp*y + 6.0_rp*z
             dtdy   = 6.0_rp*y + 5.0_rp*x + 7.0_rp*z
             d2tdx2 = 4.0_rp             
             d2tdy2 = 6.0_rp           
             if( ndime == 3 ) then
                dtdz   = 8.0_rp*z + 6.0_rp*x + 7.0_rp*y
                d2tdz2 = 8.0_rp             
             end if

             !t      = 2.0_rp*x*x + 3.0_rp*y*y 
             !dtdx   = 4.0_rp*x
             !dtdy   = 6.0_rp*y 
             !d2tdx2 = 4.0_rp             
             !d2tdy2 = 6.0_rp           
             !if( ndime == 3 ) then
             !   dtdz   = 8.0_rp*z + 6.0_rp*x + 7.0_rp*y
             !   d2tdz2 = 8.0_rp             
             !end if

          end if

          if( itask == 1 ) then
             !
             ! Exact unknown and gradients
             !
             exunk(igaus)   = t
             exgra(1,igaus) = dtdx
             exgra(2,igaus) = dtdy
             if( ndime == 3 ) exgra(3,igaus) = dtdz

          else if( itask == 3 ) then
             !
             ! Exact unknown
             !
             exunk(1) = t

          else if( itask == 2 .and. present(rhsit) ) then
             !
             ! Force term RHSIT= rho*[du/dt+a.grad(u)]-div[k*grad(u)]+r*u
             !
             rhsit(igaus) = rhsit(igaus)                                   &
                  &        + rho*( dtdt + u * dtdx + v * dtdy + w * dtdz ) &
                  &        - dkdx * dtdx - dkdy * dtdy - dkdz * dtdz       &
                  &        - k * ( d2tdx2 + d2tdy2 + d2tdz2 )              &
                  &        + r*t + r2*t*t + r3*t*t*t + r4*t*t*t*t 
          end if

       end do

    end if

  end subroutine ADR_manufactured_solution_and_rhs

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Compute the error
  !> @details Compute the error with respect to a manufactured solution
  !
  !-----------------------------------------------------------------------

  subroutine ADR_manufactured_error(ADR,ittim,cutim,unkno)

    type(ADR_typ),          intent(in) :: ADR
    integer(ip),            intent(in) :: ittim
    real(rp),               intent(in) :: cutim
    real(rp),      pointer, intent(in) :: unkno(:,:)
    integer(4)                         :: lun_output4
    integer(ip)                        :: ielem,inode,ipoin,igaus
    integer(ip)                        :: pelty,pgaus,idime,pnode
    real(rp)                           :: gpcar(ndime,mnode),xjaci(9),xjacm(9) 
    real(rp)                           :: elcod(ndime,mnode),elunk(mnode)       
    real(rp)                           :: diunk,digrt,abunk,abgrt
    real(rp)                           :: gpunk,gpgra(3)
    real(rp)                           :: exunk(mgaus),exgra(ndime,mgaus)
    real(rp)                           :: gpvol,gpdet,gpcod(3)
    real(rp)                           :: err01(2),err02(2),err0i(2)
    real(rp)                           :: err11(2),err12(2),err1i(2),xvolu,xinte

    err01 = 0.0_rp 
    err02 = 0.0_rp
    err0i = 0.0_rp 
    err11 = 0.0_rp 
    err12 = 0.0_rp
    err1i = 0.0_rp
    xvolu = 0.0_rp
    xinte = 0.0_rp

    elements: do ielem = 1,nelem

       pelty = ltype(ielem)         
       if( pelty > 0 ) then

          pnode = lnnod(ielem)
          pgaus = ngaus(pelty)
          !
          ! Gather operations
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             elunk(inode)         = unkno(ipoin,1)
          end do
          !
          ! Manufactured solution at Gauss points
          !
          call ADR_manufactured_solution_and_rhs(&
               1_ip,ADR % kfl_manufactured,pnode,pgaus,cutim,elcod,exunk,&
               exgra,elmar_loc(pelty) % shape)

          gauss_points: do igaus = 1,pgaus
             !
             ! Cartesian derivatives and Jacobian 
             !
             call elmder(&
                  pnode,ndime,elmar_loc(pelty) % deriv(1,1,igaus),&
                  elcod,gpcar,gpdet,xjacm,xjaci)
             gpvol = elmar_loc(pelty) % weigp(igaus) * gpdet  
             !
             ! Gauss point values
             !
             gpunk = 0.0_rp
             gpgra = 0.0_rp
             gpcod = 0.0_rp
             do inode = 1,pnode
                gpunk          = gpunk          + elunk(inode)         * elmar_loc(pelty) % shape(inode,igaus)
                gpcod(1:ndime) = gpcod(1:ndime) + elcod(1:ndime,inode) * elmar_loc(pelty) % shape(inode,igaus)
                gpgra(1:ndime) = gpgra(1:ndime) + elunk(inode)         * gpcar(1:ndime,inode)
             end do

             xvolu    = xvolu + gpvol
             xinte    = xinte + gpvol * gpcod(1)**4.0_rp

             diunk    = abs(gpunk-exunk(igaus))
             abunk    = abs(exunk(igaus))
             err01(1) = err01(1) + diunk * gpvol
             err02(1) = err02(1) + diunk * diunk * gpvol
             err0i(1) = max(err0i(1),diunk)
             err01(2) = err01(2) + abunk * gpvol
             err02(2) = err02(2) + exunk(igaus) * exunk(igaus) * gpvol
             err0i(2) = max(err0i(2),abunk)

             do idime = 1,ndime
                digrt    = abs(gpgra(idime)-exgra(idime,igaus))
                abgrt    = abs(exgra(idime,igaus))
                err11(1) = err11(1) + digrt*gpvol
                err12(1) = err12(1) + digrt*digrt*gpvol
                err1i(1) = max(err1i(1),digrt)
                err11(2) = err11(2) + abgrt*gpvol
                err12(2) = err12(2) + exgra(idime,igaus)*exgra(idime,igaus)*gpvol
                err1i(2) = max(err1i(2),abgrt)
             end do

          end do gauss_points
       end if
    end do elements

    call PAR_SUM(2_ip,err01,'IN MY CODE')
    call PAR_SUM(2_ip,err02,'IN MY CODE')
    call PAR_SUM(2_ip,err11,'IN MY CODE')
    call PAR_SUM(2_ip,err12,'IN MY CODE')
    call PAR_MAX(2_ip,err0i,'IN MY CODE')
    call PAR_MAX(2_ip,err1i,'IN MY CODE')

    err02(1) = sqrt(err02(1))
    err12(1) = sqrt(err12(1)) 
    err02(2) = sqrt(err02(2))
    err12(2) = sqrt(err12(2))

    if( err01(2) > zeror ) err01(1) = err01(1) / err01(2)  ! L1(u)
    if( err02(2) > zeror ) err02(1) = err02(1) / err02(2)  ! L2(u)
    if( err0i(2) > zeror ) err0i(1) = err0i(1) / err0i(2)  ! Li(u)
    if( err11(2) > zeror ) err11(1) = err11(1) / err11(2)  ! L1(grad(u))
    if( err12(2) > zeror ) err12(1) = err12(1) / err12(2)  ! L2(grad(u))
    if( err1i(2) > zeror ) err1i(1) = err1i(1) / err1i(2)  ! Li(grad(u))

    if( ADR % lun_output4 <= 0 ) then
       lun_output4 = 6_4
    else
       lun_output4 = ADR % lun_output4
    end if

    if( INOTSLAVE ) then       
       write( lun_output4,100 ) &
            ittim,cutim,&
            err01(1),err02(1),err0i(1),&
            err11(1),err12(1),err1i(1)
    end if

100 format(///,&
         & 5x,'FINITE ELEMENT ERRORS',/,&
         & 5x,'=====================',//,&
         & 5x,'TIME STEP= ',i5,', CURRENT TIME= ',es16.8e3,//,&
         & '     NORM       VALUE  ',/,5x,23('-'),/,&
         & '     W(0,1) ',es16.8e3,/,&
         & '     W(0,2) ',es16.8e3,/,&
         & '     W(0,i) ',es16.8e3,/,&
         & '     W(1,1) ',es16.8e3,/,&
         & '     W(1,2) ',es16.8e3,/,&
         & '     W(1,i) ',es16.8e3,/,5x,23('-'))

  end subroutine ADR_manufactured_error

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Element length
  !> @details This routine computes the characteristic element lengths CHALE 
  !>          according to a given strategy. CHALE is divided by two for
  !>          quadratic elements:
  !>          KFL_LENGTH = 0 ... CHALE(1) = Minimum element length
  !>                         ... CHALE(2) = Minimum element length
  !>          KFL_LENGTH = 1 ... CHALE(1) = Maximum element length
  !>                         ... CHALE(2) = Maximum element length
  !>          KFL_LENGTH = 2 ... CHALE(1) = Average element length
  !>                         ... CHALE(2) = Average element length
  !>          KFL_LENGTH = 3 ... IF KFL_ADVEC = 1:
  !>                             CHALE(1) = Flow direction
  !>                             CHALE(2) = Flow direction
  !>                             ELSE IF KFL_ADVEC =0:
  !>                             CHALE(1) = Minimum element length
  !>                             CHALE(2) = Minimum element length
  !>          KFL_LENGTH = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
  !>                         ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
  !>          KFL_LENGTH = 5 ... CHALE(1) = Length in flow direction
  !>                         ... CHALE(2) = Minimum element kength
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_element_lengths(&
       tragl,hleng,elcod,chale,pnode,porde,hnatu,kfl_length)
    implicit none
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: porde
    integer(ip), intent(in)  :: kfl_length
    real(rp),    intent(in)  :: hnatu
    real(rp),    intent(out) :: chale(2)
    real(rp),    intent(in)  :: tragl(ndime,ndime)
    real(rp),    intent(in)  :: hleng(ndime)
    real(rp),    intent(in)  :: elcod(ndime,pnode)
    integer(ip)              :: idime
!    integer(ip)              :: inode
!    real(rp)                 :: elno1,elno2

    if( kfl_length == 0 ) then 
       !
       ! Minimum element length
       !
       chale(1) = hleng(ndime) 
       chale(2) = chale(1)

    else if( kfl_length == 1 ) then   
       !
       ! Maximum element length
       !     
       chale(1) = hleng(1) 
       chale(2) = chale(1)

    else if( kfl_length == 2 ) then 
       !
       ! Average length
       !
       chale(1) = 0.0_rp
       do idime = 1,ndime
          chale(1) = chale(1) + hleng(idime)
       end do
       chale(1) = chale(1) / real(ndime,rp) 
       chale(2) = chale(1)

    else if( kfl_length == 3 ) then 
       !
       ! Length in flow direction
       !
       call runend('ADR_element_lengths: NOT CODED')
       !if(kfl_advec/=0) then 
       !   !
       !   ! Characteristic element velocity (average)
       !   !
       !   chave=0.0_rp
       !   do idime=1,ndime
       !      do inode=1,pnode
       !         chave(idime,1)=chave(idime,1)+elvel(idime,inode)
       !      end do
       !      chave(idime,1)=chave(idime,1)/real(pnode,rp)
       !   end do
       !   !
       !   ! Characteristic element length u^l = J^(-t) u^g
       !   !
       !   call mbvab1(chave(1,2),tragl,chave(1,1),ndime,ndime,elno2,elno1)
       !   if(elno2>1.0e-16.and.elno1>1.0e-16) then
       !      chale(1)=hnatu*elno1/elno2
       !   else
       !      chale(1)=hleng(ndime)
       !   end if
       !   chale(2)=chale(1)
       !   chale(2)=hleng(ndime)
       !   if (ndime ==3 ) then
       !      chale(2)=(hleng(ndime)*hleng(2)*hleng(1))**(1.0_rp/3.0_rp)
       !   else if (ndime==2) then
       !      chale(2)=sqrt(hleng(2)*hleng(1))
       !   end if
       !else
       !   chale(1)=hleng(ndime)       
       !   chale(2)=chale(1)
       !end if

    else if( kfl_length == 4 ) then 
       !
       ! sqrt(hmin*hmax)
       !
       chale(1) = sqrt(hleng(1)*hleng(ndime))     
       chale(2) = chale(1)

    else if( kfl_length == 5 ) then 
       !
       ! Along velocity direction
       !
       call runend('ADR_element_lengths: NOT CODED')
       !call velchl(pnode,elcod,elvel,chale,hleng)

    end if
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(1) = chale(1) / real(porde,rp)
    chale(2) = chale(2) / real(porde,rp)

  end subroutine ADR_element_lengths

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   End a time step
  !> @details End a time step (:,3) <= (:,1)
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_end_time_step_1(nsize,ADR,unkno)

    integer(ip),   intent(in)              :: nsize
    type(ADR_typ), intent(inout)           :: ADR(nsize)
    real(rp),      intent(inout), optional :: unkno(:,:,:)
    integer(ip)                            :: ipoin,isize,ielem
    integer(ip)                            :: pgaus,itime

    do isize = 1,nsize

       if( ADR(isize) % kfl_time_integration == 1 ) then

          if( ADR(isize) % kfl_time_scheme == TRAPEZOIDAL .and. ADR(isize) % kfl_time_order == 2 ) then
             !
             ! Crank-Nicolson: (:,1) <= 2*(:,1)-(:,3)
             ! 
             if( present(unkno) ) then
                !
                ! Unknown
                !
                do ipoin = 1,ADR(isize) % nunkn
                   unkno(isize,ipoin,1) = &
                      &   2.0_rp * unkno(isize,ipoin,1) &
                      &          - unkno(isize,ipoin,3)
                end do
             end if
             if( ADR(isize) % kfl_time_bubble /= 0 ) then
                !
                ! Bubble
                !             
                do ielem = 1,nelem 
                   ADR(isize) % bubble(ielem,1) = &
                      &  2.0_rp * ADR(isize) % bubble(ielem,1) &
                      &         - ADR(isize) % bubble(ielem,3)
                end do
             else if( ADR(isize) % kfl_time_sgs /= 0 ) then
                !
                ! SGS
                ! 
                if( ADR(isize) % kfl_first_order_sgs == 1 ) then
                   continue
                else
                   do ielem = 1,nelem 
                      pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
                      ADR(isize) % sgs(ielem) % a(1,1:pgaus,1) = &
                         &  2.0_rp * ADR(isize) % sgs(ielem) % a(1,1:pgaus,1) &
                         &         - ADR(isize) % sgs(ielem) % a(1,1:pgaus,3)
                   end do
                end if
             end if

          else if( ADR(isize) % kfl_time_scheme == BDF ) then
             !
             ! BDF scheme: (:,5) <= (:,4)
             !             (:,4) <= (:,3)
             !
             if( present(unkno) ) then
                !
                ! Unknown
                !
                do ipoin = 1,ADR(isize) % nunkn
                   do itime = 2 + ADR(isize) % kfl_time_order,4,-1
                      unkno(isize,ipoin,itime) = unkno(isize,ipoin,itime-1)
                   end do
                end do
             end if
             if( ADR(isize) % kfl_time_bubble /= 0 ) then
                !
                ! Bubble
                !
                do ielem = 1,nelem 
                   do itime = 2 + ADR(isize) % kfl_time_order,4,-1
                      ADR(isize) % bubble(ielem,itime) = ADR(isize) % bubble(ielem,itime-1) 
                   end do
                end do
             else if( ADR(isize) % kfl_time_sgs /= 0 ) then
                !
                ! SGS
                ! 
                if( ADR(isize) % kfl_first_order_sgs == 1 ) then
                   continue
                else
                   do ielem = 1,nelem 
                      pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
                      ADR(isize) % sgs(ielem) % a(1,1:pgaus,itime) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,itime-1) 
                   end do
                end if
             end if

          end if
          !
          ! New time old values (:,3) <= (:,1)
          !
          if( present(unkno) ) then
             !
             ! Unknown
             !
             do ipoin = 1,ADR(isize) % nunkn
                unkno(isize,ipoin,3) = unkno(isize,ipoin,1)
             end do
          end if
          if( ADR(isize) % kfl_time_bubble /= 0 ) then
             !
             ! Bubble
             !
             do ielem = 1,nelem 
                ADR(isize) % bubble(ielem,3) = ADR(isize) % bubble(ielem,1) 
             end do
          else if( ADR(isize) % kfl_time_sgs /= 0 ) then
             !
             ! SGS
             ! 
             do ielem = 1,nelem 
                pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
                ADR(isize) % sgs(ielem) % a(1,1:pgaus,3) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,1)
             end do
          end if
       end if
    end do

  end subroutine ADR_end_time_step_1

  subroutine ADR_end_time_step_s(ADR,unkno)

    type(ADR_typ), intent(inout)           :: ADR
    real(rp),      intent(inout), optional :: unkno(:,:)
    integer(ip)                            :: ipoin,ielem
    integer(ip)                            :: pgaus,itime

    if( ADR % kfl_time_integration == 1 ) then

       if( ADR % kfl_time_scheme == TRAPEZOIDAL .and. ADR % kfl_time_order == 2 ) then
          !
          ! Crank-Nicolson: (:,1) <= 2*(:,1)-(:,3)
          ! 
          if( present(unkno) ) then
             !
             ! Unknown
             !
             do ipoin = 1,ADR % nunkn
                unkno(ipoin,1) = &
                   &   2.0_rp * unkno(ipoin,1) &
                   &          - unkno(ipoin,3)
             end do
          end if
          if( ADR % kfl_time_bubble /= 0 ) then
             !
             ! Bubble
             !             
             do ielem = 1,nelem 
                ADR % bubble(ielem,1) = &
                   &  2.0_rp * ADR % bubble(ielem,1) &
                   &         - ADR % bubble(ielem,3)
             end do
          else if( ADR % kfl_time_sgs /= 0 ) then
             !
             ! SGS
             ! 
             if( ADR % kfl_first_order_sgs == 1 ) then
                continue
             else
                do ielem = 1,nelem 
                   pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
                   ADR % sgs(ielem) % a(1,1:pgaus,1) = &
                      &  2.0_rp * ADR % sgs(ielem) % a(1,1:pgaus,1) &
                      &         - ADR % sgs(ielem) % a(1,1:pgaus,3)
                end do
             end if
          end if

       else if( ADR % kfl_time_scheme == BDF ) then
          !
          ! BDF scheme: (:,5) <= (:,4)
          !             (:,4) <= (:,3)
          !
          if( present(unkno) ) then
             !
             ! Unknown
             !
             do itime = 2 + ADR % kfl_time_order_save,4,-1
                do ipoin = 1,ADR % nunkn
                   unkno(ipoin,itime) = unkno(ipoin,itime-1)
                end do
             end do
          end if

          if( ADR % kfl_time_bubble /= 0 ) then
             !
             ! Bubble
             !
             do ielem = 1,nelem 
                do itime = 2 + ADR % kfl_time_order_save,4,-1
                   ADR % bubble(ielem,itime) = ADR % bubble(ielem,itime-1) 
                end do
             end do
          else if( ADR % kfl_time_sgs /= 0 ) then
             !
             ! SGS
             ! 
             if( ADR % kfl_first_order_sgs == 1 ) then
                continue
             else             
                do ielem = 1,nelem 
                   pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
                   ADR % sgs(ielem) % a(1,1:pgaus,itime) = ADR % sgs(ielem) % a(1,1:pgaus,itime-1) 
                end do
             end if
          end if

       end if
       !
       ! New time old values (:,3) <= (:,1)
       !
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR % nunkn
             unkno(ipoin,3) = unkno(ipoin,1)
          end do
       end if
       if( ADR % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR % bubble(ielem,3) = ADR % bubble(ielem,1) 
          end do
       else if( ADR % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
             ADR % sgs(ielem) % a(1,1:pgaus,3) = ADR % sgs(ielem) % a(1,1:pgaus,1)
          end do
       end if

    end if

  end subroutine ADR_end_time_step_s

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Begin of an outer iteration 
  !> @details Begin of an outer iteration (:,2) <= (:,3)
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_begin_time_step_1(nsize,ADR,unkno)

    integer(ip),   intent(in)              :: nsize
    type(ADR_typ), intent(inout)           :: ADR(nsize)
    real(rp),      intent(inout), optional :: unkno(:,:,:)
    integer(ip)                            :: ipoin,isize,ielem
    integer(ip)                            :: pgaus

    do isize = 1,nsize
       
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR(isize) % nunkn
             unkno(isize,ipoin,2) = unkno(isize,ipoin,1)
          end do
       end if
       if( ADR(isize) % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR(isize) % bubble(ielem,2) = ADR(isize) % bubble(ielem,1) 
          end do
       else if( ADR(isize) % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
             ADR(isize) % sgs(ielem) % a(1,1:pgaus,2) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,1)
          end do
       end if
    end do

  end subroutine ADR_begin_time_step_1

  subroutine ADR_begin_time_step_s(ADR,unkno)

    type(ADR_typ), intent(inout)           :: ADR
    real(rp),      intent(inout), optional :: unkno(:,:)
    integer(ip)                            :: ipoin,ielem
    integer(ip)                            :: pgaus

    if( INOTMASTER ) then

       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR % nunkn
             unkno(ipoin,2) = unkno(ipoin,1)
          end do
       end if
       if( ADR % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR % bubble(ielem,2) = ADR % bubble(ielem,1) 
          end do
       else if( ADR % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
             ADR % sgs(ielem) % a(1,1:pgaus,2) = ADR % sgs(ielem) % a(1,1:pgaus,1)
          end do
       end if

    end if

  end subroutine ADR_begin_time_step_s

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Begin an inner iteration
  !> @details Begin an inner iteration (:,1) <= (:,2)
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_begin_inner_iteration_1(nsize,ADR,unkno)

    integer(ip),   intent(in)              :: nsize
    type(ADR_typ), intent(inout)           :: ADR(nsize)
    real(rp),      intent(inout), optional :: unkno(:,:,:)
    integer(ip)                            :: ipoin,isize,ielem,pgaus
 
    if( INOTMASTER ) then
    do isize = 1,nsize 
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR(isize) % nunkn
             unkno(isize,ipoin,1) = unkno(isize,ipoin,2)
          end do
       end if
       if( ADR(isize) % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR(isize) % bubble(ielem,1) = ADR(isize) % bubble(ielem,2) 
          end do
       else if( ADR(isize) % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
             ADR(isize) % sgs(ielem) % a(1,1:pgaus,1) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,2)
          end do
       end if
    end do
 end if

  end subroutine ADR_begin_inner_iteration_1

  subroutine ADR_begin_inner_iteration_s(ADR,unkno)

     type(ADR_typ), intent(inout)           :: ADR
     real(rp),      intent(inout), optional :: unkno(:,:)
     integer(ip)                            :: ipoin,ielem,pgaus

    if( INOTMASTER ) then
    if( present(unkno) ) then
       !
       ! Unknown
       !
       do ipoin = 1,ADR % nunkn
          unkno(ipoin,1) = unkno(ipoin,2)
       end do
    end if
    if( ADR % kfl_time_bubble /= 0 ) then
       !
       ! Bubble
       !
       do ielem = 1,nelem 
          ADR % bubble(ielem,1) = ADR % bubble(ielem,2) 
       end do
    else if( ADR % kfl_time_sgs /= 0 ) then
       !
       ! SGS
       ! 
       do ielem = 1,nelem 
          pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
          ADR % sgs(ielem) % a(1,1:pgaus,1) = ADR % sgs(ielem) % a(1,1:pgaus,2)
       end do
    end if
 end if

  end subroutine ADR_begin_inner_iteration_s

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   End an inner iteration
  !> @details End an inner iteration (:,2) <= (:,1)
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_end_inner_iteration_1(nsize,ADR,unkno)

    integer(ip),   intent(in)              :: nsize
    type(ADR_typ), intent(inout)           :: ADR(nsize)
    real(rp),      intent(inout), optional :: unkno(:,:,:)
    integer(ip)                            :: ipoin,isize,ielem
    integer(ip)                            :: pgaus

    if( INOTMASTER ) then
    do isize = 1,nsize
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR(isize) % nunkn
             unkno(isize,ipoin,2) = unkno(isize,ipoin,1)
          end do
       end if
       if( ADR(isize) % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR(isize) % bubble(ielem,2) = ADR(isize) % bubble(ielem,1) 
          end do
       else if( ADR(isize) % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
             ADR(isize) % sgs(ielem) % a(1,1:pgaus,2) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,1)
          end do
       end if
    end do
 end if

  end subroutine ADR_end_inner_iteration_1

  subroutine ADR_end_inner_iteration_s(ADR,unkno)

    type(ADR_typ), intent(inout)           :: ADR
    real(rp),      intent(inout), optional :: unkno(:,:)
    integer(ip)                            :: ipoin,ielem
    integer(ip)                            :: pgaus

    if( INOTMASTER ) then
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR % nunkn
             unkno(ipoin,2) = unkno(ipoin,1)
          end do
       end if
       if( ADR % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR % bubble(ielem,2) = ADR % bubble(ielem,1) 
          end do
       else if( ADR % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
             ADR % sgs(ielem) % a(1,1:pgaus,2) = ADR % sgs(ielem) % a(1,1:pgaus,1)
          end do
       end if
    end if

  end subroutine ADR_end_inner_iteration_s

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   After restart
  !> @details After restart (:,1) <= (:,3)
  !>
  !-----------------------------------------------------------------------

  subroutine ADR_after_restart_1(nsize,ADR,unkno)

    integer(ip),   intent(in)              :: nsize
    type(ADR_typ), intent(inout)           :: ADR(nsize)
    real(rp),      intent(inout), optional :: unkno(:,:,:)
    integer(ip)                            :: ipoin,isize,ielem
    integer(ip)                            :: pgaus,icomp

    if( INOTMASTER ) then

       do isize = 1,nsize
          if( ADR(isize) % kfl_time_integration == 0 ) then
             icomp = 2
          else
             icomp = 3
          end if

          if( present(unkno) ) then
             !
             ! Unknown
             !
             do ipoin = 1,ADR(isize) % nunkn
                unkno(isize,ipoin,1) = unkno(isize,ipoin,icomp)
             end do
          end if
          if( ADR(isize) % kfl_time_bubble /= 0 ) then
             !
             ! Bubble
             !
             do ielem = 1,nelem 
                ADR(isize) % bubble(ielem,1) = ADR(isize) % bubble(ielem,icomp) 
             end do
          else if( ADR(isize) % kfl_time_sgs /= 0 ) then
             !
             ! SGS
             ! 
             do ielem = 1,nelem 
                pgaus = size(ADR(isize) % sgs(ielem) % a,2,KIND=ip)
                ADR(isize) % sgs(ielem) % a(1,1:pgaus,1) = ADR(isize) % sgs(ielem) % a(1,1:pgaus,icomp)
             end do
          end if
       end do

    end if

   end subroutine ADR_after_restart_1

   subroutine ADR_after_restart_s(ADR,unkno)

      type(ADR_typ), intent(inout)           :: ADR
      real(rp),      intent(inout), optional :: unkno(:,:)
      integer(ip)                            :: ipoin,ielem
      integer(ip)                            :: pgaus,icomp

      if( ADR % kfl_time_integration == 0 ) then
         icomp = 2
      else
         icomp = 3
      end if

    if( INOTMASTER ) then
       if( present(unkno) ) then
          !
          ! Unknown
          !
          do ipoin = 1,ADR % nunkn
             unkno(ipoin,1) = unkno(ipoin,icomp)
          end do
       end if
       if( ADR % kfl_time_bubble /= 0 ) then
          !
          ! Bubble
          !
          do ielem = 1,nelem 
             ADR % bubble(ielem,1) = ADR % bubble(ielem,icomp) 
          end do
       else if( ADR % kfl_time_sgs /= 0 ) then
          !
          ! SGS
          ! 
          do ielem = 1,nelem 
             pgaus = size(ADR % sgs(ielem) % a,2,KIND=ip)
             ADR % sgs(ielem) % a(1,1:pgaus,1) = ADR % sgs(ielem) % a(1,1:pgaus,icomp)
          end do
       end if
    end if

   end subroutine ADR_after_restart_s

   !-----------------------------------------------------------------------
   !>
   !> @author  Guillaume Houzeaux
   !> @brief   Compute time step
   !> @details Compute time step
   !>
   !----------------------------------------------------------------------

   subroutine ADR_time_strategy(ittim,dtinv,dtinv_old,ADR)
      integer(ip),   intent(in)    :: ittim         !< Time step number
      real(rp),      intent(in)    :: dtinv         !< Kernel step size
      real(rp),      intent(in)    :: dtinv_old(*)  !< Kernel old step size
      type(ADR_typ), intent(inout) :: ADR           !< ADR type

      if( ADR % kfl_time_integration == 0 ) then
         ADR % dtinv           = 0.0_rp  
         ADR % dtinv_old       = 0.0_rp  
      else
         ADR % dtinv           = dtinv
         ADR % dtinv_old(1:10) = dtinv_old(1:10)
      end if

      if( ADR % kfl_time_scheme == TRAPEZOIDAL ) then
         !
         ! Trapezoidal rule: Euler iterations
         !
         if( ittim <= ADR % number_euler_steps ) then
            ADR % kfl_time_order = 1
         else
            ADR % kfl_time_order = ADR % kfl_time_order_save 
         end if
         ADR % ntime = 2

      else if( ADR % kfl_time_scheme == BDF ) then
         !
         ! BDF scheme: increase integration order at each time step
         !     
         if( ittim <= ADR % number_euler_steps ) then
            ADR % kfl_time_order = 1
         else
            ADR % kfl_time_order = min(ADR % kfl_time_order_save,ittim)
         end if
         ADR % ntime = 1 + ADR % kfl_time_order
      end if
      !
      ! Time coefficients
      !
      call ADR_time_coefficients(&
         ADR % kfl_time_scheme,ADR % kfl_time_order,ADR % time_parameters,dtinv,dtinv_old)

   end subroutine ADR_time_strategy

end module mod_ADR
!> @}
