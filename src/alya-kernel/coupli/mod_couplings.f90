!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling functions
!> @file    mod_couplings.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for coupli
!> @details ToolBox for coupli
!>          To create a coupling:
!>
!>          1. Initialize the structure:
!>          - call COU_INITIALIZE_COUPLING_STRUCTURE(COUPLING)
!>          2. Compute the following:
!>          - COUPLING % GEOME % NUMBER_WET_POINTS ....................... Number of wet points
!>          - COUPLING % GEOME % NPOIN_WET ............................... Number of wet nodes
!>          - COUPLING % GEOME % COORD_WET(:,1:NUMBER_WET_POINTS) ........ Coordinates of wet points
!>          - COUPLING % GEOME % LPOIN_WET(1:NPOIN_WET) .................. List of wet nodes
!>          - COUPLING % ITYPE ........................................... Vector projection
!>          - COUPLING % KIND ............................................ BETWEEN_SUBDOMAINS/BETWEEN_ZONES
!>          - COUPLING % COLOR_TARGET .................................... Target color
!>          - COUPLING % COLOR_SOURCE .................................... Source color
!>          3. Initialize the coupling:
!>          - call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_type(icoup) % geome % coord_wet,color_target,color_source,COUPLING)
!>          4. Generate transmission matrices:
!>          - call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling,what_to_do)
!>          - call COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES(coupling,what_to_do)
!>          - call COU_PARALLELIZE_TRANSMISSION_MATRICES(coupling)
!>
!>          Trick, if the NUMBER_WET_POINTS wet points do not have anything to do with the
!>          mesh, do the following:
!>
!>          - KIND = BETWEEN_ZONES
!>          - NPOIN_WET = NUMBER_WET_POINTS
!>          - LPOIN_WET(1:NPOIN_WET) = 1:NPOIN_WET
!>
!>          Implicit coupling
!>          -----------------
!>
!>          * Mass matrix
!>
!>          Let x = [1,1,1...]^t
!>          The lumped mass matrix Ml can be obtained from the mass matrix M as
!>          Ml = M.x
!>          Thus it can be computed just like the classical SpMV y = Ax
!>          which is valid as far as xd=Td.xn. This is our case a 1d=Td.1n due
!>          to the constant conservation of the Dirichlet transmission matrix.
!>
!> @{
!------------------------------------------------------------------------

module mod_couplings
  use def_kintyp_basic,   only : ip,rp,lg,r1p,r2p,i1p,i2p
  use def_kintyp_comm,    only : comm_data_par
  use def_spmat,          only : spmat
  use def_master,         only : ISEQUEN
  use def_master,         only : intost
  use def_master,         only : ittim
  use def_master,         only : current_code
  use def_master,         only : current_zone
  use def_master,         only : INOTMASTER
  use def_master,         only : IMASTER,kfl_paral,lninv_loc
  use def_master,         only : zeror,lzone,namda
  use def_master,         only : I_AM_IN_SUBD
  use def_master,         only : AT_BEGINNING
  use def_master,         only : THIS_NODE_IS_MINE
  use def_master,         only : modul
  use def_domain,         only : lnods
  use def_domain,         only : lesub,nelem,lnnod
  use def_domain,         only : mnodb,mnode,ndime
  use def_domain,         only : nbono,coord,ltopo,ltype
  use def_domain,         only : lbono,npoin,lnoch
  use def_domain,         only : nboun,lnodb,lnnob
  use def_domain,         only : lelch,meshe,lbset
  use def_domain,         only : nnode,kfl_codbo
  use def_domain,         only : xfiel
  use def_kermod,         only : ielse,relse,ndivi
  use def_elmtyp,         only : NOFRI
  use def_elmtyp,         only : ELHOL
  use mod_elmgeo,         only : elmgeo_natural_coordinates
  use mod_elmgeo,         only : elmgeo_natural_coordinates_on_boundaries
  use mod_elsest,         only : elsest_host_element
  use mod_parall,         only : PAR_COMM_COLOR_PERM
  use mod_parall,         only : par_part_in_color
  use mod_parall,         only : par_code_zone_subd_to_color
  use mod_parall,         only : color_target
  use mod_parall,         only : color_source
  use mod_parall,         only : par_bin_comin
  use mod_parall,         only : par_bin_comax
  use mod_parall,         only : par_bin_part
  use mod_parall,         only : par_bin_boxes
  use mod_parall,         only : par_bin_size
  use mod_parall,         only : PAR_COMM_CURRENT
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : PAR_COMM_COLOR_ARRAY
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_parall,         only : par_part_comin
  use mod_parall,         only : par_part_comax
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_alloca_min
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_resize
  use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_kdtree,         only : typ_kdtree
  use mod_kdtree,         only : kdtree_nearest_boundary
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_BARRIER
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_START_NON_BLOCKING_COMM
  use mod_communications, only : PAR_END_NON_BLOCKING_COMM
  use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_ALLGATHERV
  use mod_communications, only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_coupli,         only : scala_cou
  use def_coupli,         only : kdtree_typ
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_type
  use def_coupli,         only : UNKNOWN
  use def_coupli,         only : RESIDUAL
  use def_coupli,         only : DIRICHLET_IMPLICIT
  use def_coupli,         only : DIRICHLET_EXPLICIT
  use def_coupli,         only : RELAXATION_SCHEME
  use def_coupli,         only : UNKNOWN
  use def_coupli,         only : TRANSPOSE_MIRROR
  use def_coupli,         only : ON_WHOLE_MESH
  use def_coupli,         only : ON_IMMERSED_MESH
  use def_coupli,         only : ON_FLOATING_POINTS
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : typ_coupling_wet
  use def_coupli,         only : memor_cou
  use def_coupli,         only : AITKEN_SCHEME
  use def_coupli,         only : BROYDEN_SCHEME
  use def_coupli,         only : IQNLS_SCHEME
  use def_coupli,         only : STRESS_PROJECTION
  use def_coupli,         only : PROJECTION
  use def_coupli,         only : BETWEEN_SUBDOMAINS
  use def_coupli,         only : BETWEEN_ZONES
  use def_coupli,         only : ON_CHIMERA_MESH
  use def_coupli,         only : FIXED_UNKNOWN
  use def_coupli,         only : INTERFACE_MASS
  use def_coupli,         only : GLOBAL_MASS
  use def_coupli,         only : ncoup_implicit
  use def_coupli,         only : coupling_driver_couplings
  use def_coupli,         only : coupling_driver_iteration
  use def_coupli,         only : coupling_driver_max_iteration
  use def_coupli,         only : coupling_driver_number_couplings
  use def_coupli,         only : coupling_driver_tolerance
  use def_coupli,         only : nboun_cou
  use def_coupli,         only : lnodb_cou
  use def_coupli,         only : ltypb_cou
  use def_coupli,         only : lboch_cou
  use def_coupli,         only : lnnob_cou
  use def_coupli,         only : ON_SET,ON_FIELD,ON_CODE
  use def_coupli,         only : kfl_timco_cou
  use def_coupli,         only : kfl_efect
  use mod_matrix,         only : matrix_COO_spgemm
  use mod_matrix,         only : matrix_COO_aggregate
  use mod_matrix,         only : nullify_spmat
  use mod_iofile,         only : iofile_open_unit
  use mod_iofile,         only : iofile_close_unit
  use mod_iofile,         only : iofile_available_unit
  use mod_iofile,         only : iofile_flush_unit
  use mod_messages,       only : messages_live
  use mod_maths,          only : maths_matrix_vector_multiplication
  use mod_maths,          only : maths_outer_product
  use mod_projec,         only : projec_mass_conservation
  use mod_couplings_setup,   only : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_optional_argument, only : optional_argument
  use mod_std
  implicit none
  private
  
  integer(ip),   parameter :: my_huge = huge(1_ip)
  character(13), parameter :: vacal='mod_couplings'

  interface COU_INTERPOLATE_NODAL_VALUES
     module procedure COU_INTERPOLATE_NODAL_VALUES_11,&
          &           COU_INTERPOLATE_NODAL_VALUES_12,&
          &           COU_INTERPOLATE_NODAL_VALUES_22,&
          &           COU_INTERPOLATE_NODAL_VALUES_32,&
          &           COU_INTERPOLATE_NODAL_VALUES_33
  end interface COU_INTERPOLATE_NODAL_VALUES

  interface COU_PUT_VALUE_ON_TARGET
     module procedure COU_PUT_VALUE_ON_TARGET_IP_1,&
          &           COU_PUT_VALUE_ON_TARGET_IP_2,&
          &           COU_PUT_VALUE_ON_TARGET_IP_12
  end interface COU_PUT_VALUE_ON_TARGET

  public :: COU_INTERPOLATE_NODAL_VALUES
  public :: COU_INTERPOLATE_NODAL_VALUES_go
  public :: COU_RESIDUAL_FORCE
  public :: I_AM_IN_COUPLING
  public :: I_AM_INVOLVED_IN_A_COUPLING_TYPE
  public :: COU_INIT_INTERPOLATE_POINTS_VALUES                ! Initialize color coupling
  public :: COU_PRESCRIBE_DIRICHLET_IN_MATRIX
  public :: COU_LIST_SOURCE_NODES
  public :: COU_CHECK_CONVERGENCE
  public :: THERE_EXISTS_A_ZONE_COUPLING
  public :: I_HAVE_A_FRINGE_ELEMENT
  public :: MIRROR_COUPLING
  public :: I_HAVE_A_CHIMERA_COUPLING
  public :: COU_PUT_VALUE_ON_TARGET
  public :: COU_SET_FIXITY_ON_TARGET
  public :: COU_TEMPORAL_PREDICTOR
  public :: couplings_initialize_solver
  public :: couplings_impose_dirichlet
  public :: couplings_check_dirichlet                         ! Check if an array satisfies the Dirichlet coupling condition
  public :: couplings_time_step                               ! Compute time step according to strategy
  public :: couplings_initialization
  public :: couplings_select_boundaries
  public :: couplings_elements_graph
  public :: couplings_are_exhaustive
  public :: cou_wet_points_from_mirror
  public :: cou_residual_scalar
  
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-26
  !> @brief   Couplingsinitialization
  !> @details Initializaiton of coupling variables
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_initialization()

    use def_coupli
    !
    ! Memory
    !
    memor_cou              = 0
    !
    ! Read variables: initialization must be here as sometimes
    ! coupling is used without data file
    !
    mcoup                  =  0
    kfl_absolute_cou       = -1                      ! Basolute tolerance based on average element size
    toler_absolute_cou     = -1.00_rp                ! Element average size for partition bounding box
    toler_relative_cou     =  0.01_rp                ! Relative tolerance for partition bounding box
    kfl_lost_wet_point_cou = STOP_ALYA_WITH_WARNINGS ! Stop Alya with warnings
    number_of_holes        = 0                       ! No hole
    kfl_timco_cou          = 2                       ! Local time step (each Alya has its own time step)
    !
    ! Driver
    !
    coupling_driver_couplings     = 0
    coupling_driver_max_iteration = 1
    coupling_driver_iteration     = 0
    coupling_driver_tolerance     = 1.0e-12_rp
    !
    ! Others
    !
    nullify(coupling_type)
    !
    ! Implicit couplings
    !
    ncoup_implicit   = 0
    ncoup_implicit_d = 0
    ncoup_implicit_n = 0
    nullify(lcoup_implicit_d)
    nullify(lcoup_implicit_n)
    nullify(mask_cou)

  end subroutine couplings_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Find mirror coupling
  !> @details Obtain the mirror coupling MIRROR_COUPLING of ICOUP
  !>          It returns zero if no mirror has been found
  !>
  !----------------------------------------------------------------------

  function I_HAVE_A_CHIMERA_COUPLING()
    integer(ip) :: icoup
    logical(lg) :: I_HAVE_A_CHIMERA_COUPLING

    I_HAVE_A_CHIMERA_COUPLING = .false.
    do icoup = 1,mcoup
       if(    I_AM_IN_COUPLING(icoup) .and. &
            & coupling_type(icoup) % where_type == ON_CHIMERA_MESH ) then
          I_HAVE_A_CHIMERA_COUPLING = .true.
          return
       end if
    end do

  end function I_HAVE_A_CHIMERA_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Find mirror coupling
  !> @details Obtain the mirror coupling MIRROR_COUPLING of ICOUP
  !>          It returns zero if no mirror has been found
  !
  !----------------------------------------------------------------------

  function MIRROR_COUPLING(icoup)
    integer(ip), intent(in) :: icoup        !< Coupling
    integer(ip)             :: MIRROR_COUPLING
    logical(lg)             :: notfound

    MIRROR_COUPLING = 0
    notfound = .true.
    do while( notfound .and. MIRROR_COUPLING < mcoup )
       MIRROR_COUPLING = MIRROR_COUPLING + 1
       if(  &
            & coupling_type(icoup) % color_target           == coupling_type(MIRROR_COUPLING) % color_source .and. &
            & coupling_type(MIRROR_COUPLING) % color_target == coupling_type(icoup) % color_source ) then
          coupling_type(icoup) % mirror_coupling = MIRROR_COUPLING
          notfound = .false.
       end if
    end do

  end function MIRROR_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If I am in a coupling
  !> @details Check if I am involved in cupling icoup
  !
  !----------------------------------------------------------------------

  function I_AM_IN_COUPLING(icoup)
    integer(ip), intent(in) :: icoup        !< Coupling
    integer(ip)             :: icolo_source
    integer(ip)             :: icolo_target
    logical(lg)             :: I_AM_IN_COUPLING

    icolo_source = coupling_type(icoup) % color_source
    icolo_target = coupling_type(icoup) % color_target
    I_AM_IN_COUPLING = I_AM_IN_COLOR(icolo_source) .or. I_AM_IN_COLOR(icolo_target)

  end function I_AM_IN_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If I am in a coupling of a certain type
  !> @details If I am in a coupling of a certain type: RESIDUAL, UNKNOWN
  !>          DIRICHLET_EXPLICIT or DIRICHLET_IMPLICIT
  !
  !----------------------------------------------------------------------

  function I_AM_INVOLVED_IN_A_COUPLING_TYPE(ikind,iwhat,itype,micou_type)
     integer(ip), intent(in)           :: ikind    !< Coupling kind (between zones, between subdomain)
     integer(ip), intent(in)           :: iwhat    !< Coupling what (residual, unknown, dirichlet)
     integer(ip), intent(in), optional :: itype    !< Coupling type (projection, interpolation)
     integer(ip), intent(in), optional :: micou_type !< Mirror coupling type (projection, interpolation)
     integer(ip)                       :: icoup
     logical(lg)                       :: I_AM_INVOLVED_IN_A_COUPLING_TYPE
     integer(ip)                       :: micou

     I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.

     if( ikind == BETWEEN_SUBDOMAINS .and. ncoup_implicit == 0 ) then
        continue
     else
        do icoup = 1,mcoup
           if(    I_AM_IN_COUPLING(icoup) .and. &
              & coupling_type(icoup) % kind == ikind .and. &
              & coupling_type(icoup) % what == iwhat ) then
              I_AM_INVOLVED_IN_A_COUPLING_TYPE = .true.
              if(present(itype) .and. (.not.(coupling_type(icoup) % itype == itype))) &
                 & I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.
              if(present(micou_type)) then
                 micou = coupling_type(icoup) % mirror_coupling
                 if(.not.(coupling_type(micou) % itype == micou_type)) &
                    & I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.
              end if
           end if
        end do
     end if

  end function I_AM_INVOLVED_IN_A_COUPLING_TYPE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   Check if I have a fringe element
  !> @details A fringe element is a hole element with at least one if
  !>          fringe node
  !
  !----------------------------------------------------------------------

  function I_HAVE_A_FRINGE_ELEMENT()
    integer(ip) :: ipoin,ielem,inode
    logical(lg) :: I_HAVE_A_FRINGE_ELEMENT

    I_HAVE_A_FRINGE_ELEMENT = .false.
    do ielem = 1,nelem
       if( lelch(ielem) == ELHOL ) then
          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( lnoch(ipoin) == NOFRI ) then
                I_HAVE_A_FRINGE_ELEMENT = .true.
                return
             end if
          end do
       end if
    end do

  end function I_HAVE_A_FRINGE_ELEMENT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   Check if I have a hole element
  !> @details A hole element has only hole nodes
  !
  !----------------------------------------------------------------------

  function I_HAVE_A_HOLE_ELEMENT()
    integer(ip) :: ielem
    logical(lg) :: I_HAVE_A_HOLE_ELEMENT

    I_HAVE_A_HOLE_ELEMENT = .false.
    do ielem = 1,nelem
       if( lelch(ielem) == ELHOL ) then
          I_HAVE_A_HOLE_ELEMENT = .true.
          return
       end if
    end do

  end function I_HAVE_A_HOLE_ELEMENT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    30/09/2014
  !> @brief   Check convergence
  !> @details Check convergence of a coupling
  !
  !----------------------------------------------------------------------

  subroutine COU_CHECK_CONVERGENCE(iblok,kfl_gozon)
    integer(ip), intent(in)  :: iblok
    integer(ip), intent(out) :: kfl_gozon
    integer(ip)              :: icoup,kcoup
    integer(ip)              :: idime,idofn,ndofn,ipoin,kpoin
    integer(ip)              :: ntime_con,itime

    kfl_gozon = 0

    if( coupling_driver_iteration(iblok) >= coupling_driver_max_iteration(iblok) ) then
       !
       ! Not converged, max numb. of iteration overpassed
       !
       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)
          !
          ! Update values for subcycling coupling (frequency of exchanges different of one)
          !
          if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip )then
             !
             ! Only source code save exchanging values
             !
             if( current_code == coupling_type(icoup) % code_source )then

                do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source
                   ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
                   do idime = 1_ip, ndime
                      coupling_type(icoup) % values_frequ(idime,kpoin,1_ip) = coupling_type(icoup) % values_frequ(idime,kpoin,2_ip)
                   end do
                end do
             end if
          end if
          !
          ! Save last exchanged values for temporal predictor. Only the target code will make predictions
          !
          if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. current_code == coupling_type(icoup) % code_target .and. coupling_type(icoup) % temporal_predictor == 1_ip )then
             ndofn = size(coupling_type(icoup) % values_converged,1_ip)
             ntime_con = size(coupling_type(icoup) % values_converged,3)
             do itime = ntime_con,2,-1

                do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                   do idofn = 1,ndofn
                      coupling_type(icoup) % values_converged(idofn,ipoin,itime) = coupling_type(icoup) % values_converged(idofn,ipoin,itime-1)
                   end do
                end do

             end do
             do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                do idofn = 1,ndofn
                   coupling_type(icoup) % values_converged(idofn,ipoin,1_ip) = coupling_type(icoup) % values(idofn,ipoin,1_ip)
                end do
             end do

          end if

       end do
       return
    else
       
       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)
          !
          ! Only check if the current time step is a multiple of the frequency defined
          !
          if( mod( ittim,coupling_type(icoup) % frequ_send ) == 0_ip .and. current_code == coupling_type(icoup) % code_source )then
             if( coupling_type(icoup) % resid(1) > coupling_driver_tolerance(iblok) )then
                kfl_gozon = 1
             end if
          end if
          if( mod( ittim,coupling_type(icoup) % frequ_recv ) == 0_ip .and. current_code == coupling_type(icoup) % code_target )then
             if( coupling_type(icoup) % resid(1) > coupling_driver_tolerance(iblok) )then
                kfl_gozon = 1
             end if
          end if
       end do
       !
       ! Update values exchanged when convergence is reached and frequency of exchanges is larger than one
       !
       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)

          if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip )then
             ! do kcoup = 1,coupling_driver_number_couplings(iblok)
             !    icoup = coupling_driver_couplings(kcoup,iblok)
             if( mod( ittim,coupling_type(icoup) % frequ_send) == 0_ip .and. kfl_gozon==0 .and. current_code == coupling_type(icoup) % code_source )then
                do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source
                   ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
                   do idime = 1_ip, ndime
                      coupling_type(icoup) % values_frequ(idime,kpoin,1_ip) = coupling_type(icoup) % values_frequ(idime,kpoin,2_ip)
                   end do
                end do
             end if
             ! end do
          end if

          !
          ! Save converged values for temporal predictor. Only the target code will make predictions
          !

          if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. current_code == coupling_type(icoup) % code_target .and. kfl_gozon==0 .and. &
               coupling_type(icoup) % temporal_predictor == 1_ip )then

             ndofn = size(coupling_type(icoup) % values_converged,1_ip)
             ntime_con = size(coupling_type(icoup) % values_converged,3)

             ! do kcoup = 1,coupling_driver_number_couplings(iblok)
             !    icoup = coupling_driver_couplings(kcoup,iblok)

             do itime = ntime_con,2,-1
                do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                   do idofn = 1,ndime
                      coupling_type(icoup) % values_converged(idofn,ipoin,itime) = coupling_type(icoup) % values_converged(idofn,ipoin,itime-1)
                   end do
                end do
             end do

             do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                do idofn = 1,ndime
                   coupling_type(icoup) % values_converged(idofn,ipoin,1_ip) = coupling_type(icoup) % values(idofn,ipoin,1_ip)
                end do
             end do

             ! end do
          end if  ! coupling kind
       end do     ! loop over couplings
    end if

  end subroutine COU_CHECK_CONVERGENCE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   Actualize unknown according to scheme
  !> @details Actualize the unknown and save previous values according
  !>          to the scheme (relaxation, Aitken, etc.)
  !
  !----------------------------------------------------------------------
  
  subroutine COU_UPDATE_POINTS_VALUES(xxnew,coupling,xresi,mask)
    use mod_iqnls

    implicit none

    real(rp),      pointer,   intent(inout)        :: xxnew(:,:)
    type(typ_color_coupling), intent(inout)        :: coupling
    real(rp),                 intent(out)          :: xresi(2)
    integer(ip),   pointer,   intent(in), optional :: mask(:,:)
    integer(ip)                                    :: ndofn,ipoin,npoin_wet,ntime_wet
    integer(ip)                                    :: idofn,itime,kpoin,nsize
    real(rp)                                       :: relax,rela1,rip1,rip2,xfact
    real(rp)                                       :: numer,denom

    integer(ip)                                    :: current_in_iter, i_iter, i_col, shifted_columns
    integer(ip)                                    :: filled_columns_history_iqnls, filled_columns_increment_iqnls
    integer(ip)                                    :: ndofs, idofs
    integer(ip)                                    :: mindex
    integer(ip), save                              :: max_columns_increment_iqnls, filled_columns_inner_iqnls

    real(rp), allocatable                          :: alpha(:)
    real(rp)                                       :: vecaux(coupling % wet % npoin_wet * size(xxnew,1) )
    integer(ip), save                              :: inner_iter_counter
    integer(ip), allocatable,  save                :: history_tracking_iqnls(:)

    coupling % itera = coupling % itera + 1
    ntime_wet = 0

    if( INOTMASTER ) then
       !
       ! Degrees of freedom
       !
       ndofn     = size(xxnew,1)
       npoin_wet = coupling % wet % npoin_wet
       ndofs      = npoin_wet * ndofn
       !
       ! Allocate memory if required and initializate some values
       !
       if( .not. associated( coupling % values_converged ) .and. npoin_wet > 0 .and. coupling  % temporal_predictor == 1_ip )then
          !
          ! Three time steps values will be stored to make temporal predictions the allocation is done here because it is not
          ! know a priori the ndofn of the exchanged quantity
          !
          call memory_alloca(memor_cou,'values_converged',vacal,coupling % values_converged,ndofn,npoin_wet,3_ip)
          coupling % values_converged = 0.0_rp
       end if

       if( coupling % scheme == RELAXATION_SCHEME )then

          if( .not. associated(coupling % values) )then
             call memory_alloca(memor_cou,'values',vacal,coupling % values,ndofn,npoin_wet,2_ip)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == AITKEN_SCHEME )then

          if( .not. associated(coupling % values) )then
             call memory_alloca(memor_cou,'values',          vacal,coupling % values,ndofn,npoin_wet,3_ip)
             call memory_alloca(memor_cou,'values_predicted',vacal,coupling % values_predicted,ndofn,npoin_wet)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == BROYDEN_SCHEME )then

          if( .not. associated(coupling % values) )then
             !             call memory_alloca(memor_cou,vacal,'values',coupling % values,ndofn,npoin_wet,3_ip)
             !             call memory_alloca(memor_cou,vacal,'values_predicted',coupling % values_predicted,ndofn,npoin_wet)
             !             call memory_alloca(memor_cou,vacal,'jacobian_inverse',coupling % jacobian_inverse,ndofn,npoin_wet,npoin_wet)
             !             call memory_alloca(memor_cou,vacal,'dincr_predicted',coupling % dincr_predicted,npoin_wet)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == IQNLS_SCHEME )then

          if( .not. associated(coupling % relaxed_iqnls) )then
             !if(ndofs>0_ip) then

                !if(PAR_MY_CODE_RANK.eq.1)THEN
                !write(6,*) '|------------------------------------------------|'
                !write(6,*) 'ndof:', ndofn, '| npoin_wet', npoin_wet
                !write(6,*) 'ranku:', coupling % ranku_iqnls, 'history:', coupling % history_iqnls
                !write(6,*) '|------------------------------------------------|'
             !ENDIF
             nsize = max(1_ip,ndofn * npoin_wet)
                call memory_alloca(memor_cou,'values'                 ,vacal,coupling % values,ndofn,npoin_wet,1_ip)
                call memory_alloca(memor_cou,'relaxed'                ,vacal,coupling % relaxed_iqnls,           nsize, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'unrelaxed_iqnls'        ,vacal,coupling % unrelaxed_iqnls,         nsize, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'residues_iqnls'         ,vacal,coupling % residues_iqnls,          nsize, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'valincr_iqnls'          ,vacal,coupling % valincr_iqnls,           nsize, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'residincr_iqnls'        ,vacal,coupling % residincr_iqnls,         nsize, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'V_current_history_iqnls',vacal,coupling % V_current_history_iqnls, nsize, coupling % ranku_iqnls * (coupling % history_iqnls+1))
                call memory_alloca(memor_cou,'W_current_history_iqnls',vacal,coupling % W_current_history_iqnls, nsize, coupling % ranku_iqnls * (coupling % history_iqnls+1))

                if(coupling % history_iqnls == 0_ip)then
                  call memory_alloca_min(memor_cou,'valincr_history_iqnls'  ,vacal,coupling % valincr_history_iqnls )
                  call memory_alloca_min(memor_cou,'residincr_history_iqnls',vacal,coupling % residincr_history_iqnls )
                  call memory_alloca_min(memor_cou,'history_tracking_iqnls' ,vacal,coupling % history_tracking_iqnls )
                else
                  call memory_alloca(memor_cou,'valincr_history_iqnls'   ,vacal, coupling % valincr_history_iqnls,   nsize, coupling % ranku_iqnls, coupling % history_iqnls)
                  call memory_alloca(memor_cou,'residincr_history_iqnls' ,vacal, coupling % residincr_history_iqnls, nsize, coupling % ranku_iqnls, coupling % history_iqnls)
                  call memory_alloca(memor_cou,'history_tracking_iqnls'  ,vacal, coupling % history_tracking_iqnls, coupling % history_iqnls)

                endif
          endif

          if( associated(coupling % relaxed_iqnls) ) then
            ntime_wet = size(coupling % relaxed_iqnls,2)
          endif

       end if

    else
       ndofn     = 0_ip
       npoin_wet = 0_ip
       ntime_wet = 0_ip
       ndofs     = 0_ip

       if( .not. associated(coupling % values) )then
         call memory_alloca_min(memor_cou,'values'                 ,vacal, coupling % values )
         call memory_alloca_min(memor_cou,'relaxed'                ,vacal, coupling % relaxed_iqnls )
         call memory_alloca_min(memor_cou,'unrelaxed_iqnls'        ,vacal, coupling % unrelaxed_iqnls )
         call memory_alloca(memor_cou,'residues_iqnls'             ,vacal,coupling % residues_iqnls, 1_ip, coupling % ranku_iqnls)
         call memory_alloca(memor_cou,'valincr_iqnls'              ,vacal,coupling % valincr_iqnls,  1_ip, coupling % ranku_iqnls)
         call memory_alloca(memor_cou,'residincr_iqnls'            ,vacal,coupling % residincr_iqnls, 1_ip, coupling % ranku_iqnls)
         call memory_alloca_min(memor_cou,'v_current_history_iqnls',vacal, coupling % v_current_history_iqnls )
         call memory_alloca_min(memor_cou,'w_current_history_iqnls',vacal, coupling % w_current_history_iqnls )
         call memory_alloca(memor_cou,'valincr_history_iqnls'   ,vacal, coupling % valincr_history_iqnls,   1_ip, coupling % ranku_iqnls, coupling % history_iqnls)
         call memory_alloca(memor_cou,'residincr_history_iqnls' ,vacal, coupling % residincr_history_iqnls, 1_ip, coupling % ranku_iqnls, coupling % history_iqnls)

         call memory_alloca_min(memor_cou,'history_tracking_iqnls' ,vacal, coupling % history_tracking_iqnls )
         call memory_alloca_min(memor_cou,'values_predicted'       ,vacal, coupling % values_predicted )
         call memory_alloca_min(memor_cou,'values_converged'       ,vacal, coupling % values_converged )
       endif
       !!       elseif( (coupling % scheme == AITKEN_SCHEME) .or. (coupling % scheme == RELAXATION_SCHEME) )then
       !!         call memory_alloca_min( coupling % values )
       !!         call memory_alloca_min( coupling % values_converged )
       !!         call memory_alloca_min( coupling % values_predicted )
       !!       endif

    end if

    !
    ! Save old relaxed values
    !
    ! relaxed(1)=x^k+1         <- to be calculated
    ! relaxed(2)=x^k
    ! relaxed(3)=x^k-1
    ! relaxed(4)=x^k-2
    !        .
    !        .
    !        .
    !
    !
    ! And save old residues, if IQNLS
    ! residues_iqnls(1)=r^k    <- to be written
    ! residues_iqnls(2)=r^k-1
    ! residues_iqnls(3)=r^k-2
    !        .
    !        .
    !        .

!!!  if( coupling % scheme .ne. BROYDEN_SCHEME) then

    if ( (coupling % scheme == RELAXATION_SCHEME) .or. (coupling % scheme == AITKEN_SCHEME)) then
       !
       ! Rest of the schemes but broyden
       !
       do itime = ntime_wet,2,-1
          do ipoin = 1,npoin_wet
             do idofn = 1,ndofn
                coupling % values(idofn,ipoin,itime) = coupling % values(idofn,ipoin,itime-1)
             end do
          end do
       end do

    endif
!!!!  endif

    !
    ! Aitken: copy predicted value
    !
    if( coupling % scheme == RELAXATION_SCHEME ) then
       !
       ! Avoid to have first null Dirichlet guess 
       !
       if( coupling % what == UNKNOWN .and. coupling % itera == 1 ) then
          relax = 1.0_rp
       else
          relax = coupling % relax
       end if

    else if( coupling % scheme == AITKEN_SCHEME ) then
       !
       ! First two iterations are performed with constant relaxation
       !
       if( coupling_driver_iteration(1_ip) < 3_ip ) then
          ! if( coupling % itera < 3_ip ) then

          coupling % aitken = coupling % relax
          relax             = coupling % aitken
          ! print*, "DEBUG: AITKEN inicial ", relax

       else
          !
          ! Scalar product for aitken relaxation factor
          !
          numer = 0.0_rp
          denom = 0.0_rp
          do ipoin = 1,npoin_wet
             do idofn = 1,ndofn
                !
                ! rip1 = d_{i}   - d_{i+1}'
                !
                rip1  = coupling % values_predicted(idofn,ipoin) - coupling % values(idofn,ipoin,3_ip)
                ! rip1 = ( coupling % values(idofn,ipoin,2_ip) -  coupling % values(idofn,ipoin,3_ip) ) * xxnew(idofn,ipoin)
                !
                ! rip2 = d_{i+1} - d_{i+2}'
                !
                rip2  =  xxnew(idofn,ipoin) - coupling % values(idofn,ipoin,2_ip)
                ! rip2 = xxnew(idofn,ipoin) - coupling % values_predicted(idofn,ipoin)
                ! rip3 = xxnew(idofn,ipoin) - coupling % values(idofn,ipoin,2_ip)
                numer = numer + rip1 * (rip2-rip1) * coupling % wet % weight_wet(ipoin)
                denom = denom + (rip2-rip1)*(rip2-rip1) * coupling % wet % weight_wet(ipoin)
                ! numer = numer + rip1
                ! denom = denom + rip2 * rip3
             end do
          end do

          call PAR_SUM(numer,'IN CURRENT COLOR')
          call PAR_SUM(denom,'IN CURRENT COLOR')
          relax =-coupling % aitken * numer / (denom+zeror)
          ! if( relax < 0_rp .or. relax > 0.6 ) relax = 0.1_rp
          ! relax = - numer / (denom + zeror)
          ! relax = min(relax, 1._rp)
          ! relax = max(relax,-1._rp)
          coupling % aitken = relax
          ! if( relax < 0.01_rp ) relax = 0.02_rp
          ! if( relax > 1.00_rp ) relax = 0.03_rp
          ! aux = dabs(relax)
          ! print*, "DEBUG: AITKEN ", relax
       end if

    else if( coupling % scheme == BROYDEN_SCHEME ) then

       call COU_BROYDEN_BAD(xxnew,coupling)

    else if( coupling % scheme == IQNLS_SCHEME ) then
       ! Obtain the current inner iteration
       !
       current_in_iter = coupling_driver_iteration(1_ip)

       ! Initialize values in first iteration and if there's no history being saved
       !
       if( (ittim.eq.1).and.(current_in_iter .eq. 1_ip)) then
          coupling % relaxed_iqnls=0.0_rp
          coupling % unrelaxed_iqnls=0.0_rp
          coupling % residues_iqnls=0.0_rp
          coupling % valincr_iqnls=0.0_rp
          coupling % residincr_iqnls=0.0_rp
          coupling % residincr_history_iqnls=0.0_rp
          allocate(history_tracking_iqnls(coupling % history_iqnls))
          history_tracking_iqnls=0_ip
          inner_iter_counter=0_ip
       else

          ! Shifting of time columns
          !
          if ( coupling % history_iqnls > 0 .and. current_in_iter == 1 )then

            !Compute how many time columns to shift
            if (ittim .le. coupling % history_iqnls) then
              shifted_columns=ittim
            else
              shifted_columns=coupling % history_iqnls
            endif

            !IF(PAR_MY_CODE_RANK.eq.1) then
            !write(6,*) '|------------------------------------------------|'
            !write(6,*) '|---------------COLUMN SHIFTING ALGORITHM--------|'
            !write(6,*) '|------------------------------------------------|'
            !write(6,*) 'about to shift: ', shifted_columns, 'time steps'
            !write(6,*) 'tracking and history pre shifting'
            !write(6,*) coupling % history_tracking_iqnls(:)
            !write(6,*) history_tracking_iqnls(:)
            !IF(PAR_MY_CODE_RANK.eq.1)write(6,*) coupling % residincr_history_iqnls(5, : , :)
            !ENDIF

            ! Shift previous times
            do i_col = shifted_columns,2,-1

            ! IF(PAR_MY_CODE_RANK.eq.0)WRITE(6,*) 'i_col', i_col
            ! IF(PAR_MY_CODE_RANK.eq.0)write(6,*) 'shifting', coupling % history_tracking_iqnls(i_col-1), 'iterations from column', i_col-1, 'into', i_col
             if(INOTMASTER) then
               coupling % residincr_history_iqnls( : , :, i_col )  = 0.0_rp
               coupling % valincr_history_iqnls( : , : , i_col )  = 0.0_rp
               do i_iter=1, history_tracking_iqnls(i_col-1)
                  !IF(PAR_MY_CODE_RANK.eq.1)WRITE(6,*) '       moving iteration', i_iter
                  do idofs=1, ndofs
                      coupling % residincr_history_iqnls(idofs, i_iter, i_col )  = coupling % residincr_history_iqnls(idofs, i_iter, i_col-1)
                      coupling % valincr_history_iqnls(idofs, i_iter, i_col )  = coupling % valincr_history_iqnls(idofs, i_iter, i_col-1)
                  enddo
               enddo
             endif
             history_tracking_iqnls(i_col) = history_tracking_iqnls(i_col-1)
            enddo

            ! Save last time step in the newest column
            !IF(PAR_MY_CODE_RANK.eq.1) then
            !write(6,*) '|------------------------------------------------|'
            !WRITE(6,*) 'COPYING', inner_iter_counter, 'ITERATIONS FROM THE PREVIOUS TIME STEP'
            !ENDIF

            do i_iter=1, inner_iter_counter
               do idofs=1, ndofs
                   coupling % residincr_history_iqnls(idofs, i_iter,1)  = coupling % residincr_iqnls(idofs, i_iter)
                   coupling % valincr_history_iqnls(idofs, i_iter,1)  = coupling % valincr_iqnls(idofs, i_iter)
               enddo
            enddo

            history_tracking_iqnls(1) = filled_columns_inner_iqnls
            filled_columns_inner_iqnls=0
            inner_iter_counter=0_ip

!            IF(PAR_MY_CODE_RANK.eq.0) then
!            write(6,*) '|------------------------------------------------|'
!            write(6,*) 'tracking and history post algorithm'
!            write(6,*) coupling % history_tracking_iqnls(:)
!            write(6,*) history_tracking_iqnls(:)
!            IF(PAR_MY_CODE_RANK.eq.1)write(6,*) coupling % residincr_history_iqnls(5, : , :)
!            ENDIF

            coupling % relaxed_iqnls=0.0_rp
            coupling % unrelaxed_iqnls=0.0_rp
            coupling % residues_iqnls=0.0_rp
            coupling % valincr_iqnls=0.0_rp
            coupling % residincr_iqnls=0.0_rp

          endif

       ! Compute the maximum number of columns
       ! and how many of them are filled
       !
        filled_columns_inner_iqnls = current_in_iter
        if(filled_columns_inner_iqnls >= coupling % ranku_iqnls) filled_columns_inner_iqnls = coupling % ranku_iqnls

        filled_columns_history_iqnls=ittim-1
        if(filled_columns_history_iqnls >= coupling % history_iqnls) filled_columns_history_iqnls = coupling % history_iqnls

        max_columns_increment_iqnls = coupling % ranku_iqnls-1

       filled_columns_increment_iqnls=filled_columns_inner_iqnls-1

       if(coupling % history_iqnls > 0_ip) then
         max_columns_increment_iqnls = coupling % ranku_iqnls * ( coupling % history_iqnls + 1)-coupling % history_iqnls
         filled_columns_increment_iqnls= sum(history_tracking_iqnls(:))-filled_columns_history_iqnls+(filled_columns_inner_iqnls-1)
         if(filled_columns_increment_iqnls >= max_columns_increment_iqnls) filled_columns_increment_iqnls = max_columns_increment_iqnls
       endif

!IF(PAR_MY_CODE_RANK.eq.1) then
!write(6,*) ' '
!write(6,*) '|------------------------------------------------|'
!write(6,*) 'current_time_iter:', ittim,  ' | current_in_iter:', current_in_iter
!write(6,*)  'filled_inner:', filled_columns_inner_iqnls, '| filled_columns_history_iqnls: ', filled_columns_history_iqnls
!write(6,*) 'max_columns increment total (current+histories):', max_columns_increment_iqnls
!write(6,*) 'filled_increment', filled_columns_increment_iqnls
!write(6,*) '|------------------------------------------------|'
!! The last column is always zero
!write(6,*) 'residincr_pre_algorithm'
!write(6,*) coupling % residincr_iqnls(5, :)
!ENDIF
!write(6,*) 'IM RANK', PAR_MY_CODE_RANK, 'filled_increment', filled_columns_increment_iqnls

          !
          ! save residues and unrelaxed values
          ! Column 1 has the newest values
          !
          do i_col = coupling % ranku_iqnls,2,-1
                  ! write(6,*) '|------------------------------------------------|'
                  ! write(6,*) '|--- MOVING COLUMN', i_col-1, 'into', i_col
             do idofs=1, ndofs
                coupling % relaxed_iqnls(idofs,i_col)    = coupling % relaxed_iqnls(idofs,i_col-1)
                coupling % unrelaxed_iqnls(idofs, i_col) = coupling % unrelaxed_iqnls(idofs, i_col-1)
                coupling % residues_iqnls(idofs, i_col)  = coupling % residues_iqnls(idofs, i_col-1)
             enddo
          end do
       endif

       inner_iter_counter = inner_iter_counter + 1_ip

       if(inner_iter_counter > coupling % ranku_iqnls) then
          inner_iter_counter = coupling % ranku_iqnls
       endif
       !IF(PAR_MY_CODE_RANK.eq.1) THEN
       ! write(6,*) '|------------------------------------------------|'
       ! WRITE(6,*) 'COUPLING TRACKING: ', inner_iter_counter
       !ENDIF

       !!
       !! Save new values and compute residues
       !!

       mindex=0_ip
       do ipoin=1, npoin_wet
          do idofn=1, ndofn
             mindex = ipoin * ndofn + idofn - ndofn

             !
             ! In the first iteration of the current time step save the last
             ! converged value in the actual relaxed vector to compute
             ! the residue more precisely
             !
             if (current_in_iter == 1_ip .and. associated(coupling % values_converged) ) then
                coupling % relaxed_iqnls(mindex,2) = coupling % values_converged(idofn,ipoin,1_ip)
             endif

             ! Save new result, for notation consistency
             coupling % unrelaxed_iqnls(mindex,1) = xxnew(idofn, ipoin) * coupling % scaling_iqnls
             ! Compute new residue
             coupling % residues_iqnls(mindex,1) = coupling % unrelaxed_iqnls(mindex,1) - coupling % relaxed_iqnls (mindex, 2)
          enddo
       enddo

       !
       ! First two iterations are performed with constant relaxation
       !
        if( current_in_iter < 3_ip ) then
          !IF(PAR_MY_CODE_RANK.eq.1) then
          !write(6,*) '|------------------------------------------------|'
          !write(6,*) 'COMPUTING FIXED RELAXATION'
          !ENDIF

          relax = coupling % relax

          do idofs=1,ndofs
             coupling % relaxed_iqnls (idofs, 1) = relax * coupling % unrelaxed_iqnls(idofs, 1) + ( 1.0_rp -relax ) * coupling % relaxed_iqnls(idofs,2)
          enddo

       else
          !-----------------------------------------------------
          ! From the iteration 3, the actual IQN-LS
          ! will be computed
          !
          !-----------------------------------------------------

!IF(PAR_MY_CODE_RANK.eq.1) then
!write(6,*) '|------------------------------------------------|'
!write(6,*) ' GETTING INTO IQNLS'
!ENDIF
          ! Build increment vector for this iteration
          !
          ! Build vector V_i= r_{i+1} - r_i
          ! and vector   W_i = x*_{i+1} - x*_i
          !
          do i_iter=1,filled_columns_inner_iqnls-1
             do idofs=1, ndofs
                coupling % residincr_iqnls(idofs, i_iter) = coupling % residues_iqnls(idofs,i_iter+1)  - coupling % residues_iqnls(idofs,i_iter)
                coupling % valincr_iqnls(idofs, i_iter)   = coupling % unrelaxed_iqnls(idofs,i_iter+1) -  coupling % unrelaxed_iqnls(idofs,i_iter)
             enddo
          enddo

!IF(PAR_MY_CODE_RANK.eq.1) then
!write(6,*) '|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|'
!write(6,*) '|%%%%CALLING FILTERING ALGORITHM%%%%%%%%%%%%%%|'
!write(6,*) '|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|'
!write(6,*) 'epsi', coupling % efilter_iqnls
          !ENDIF
          
            call QRfiltering(coupling % residincr_iqnls,          &
                             coupling % residincr_history_iqnls,  &
                             coupling % valincr_iqnls,            &
                             coupling % valincr_history_iqnls,    &
                             filled_columns_inner_iqnls,          &
                             history_tracking_iqnls,              &
                             coupling % efilter_iqnls,            &
                             coupling % V_current_history_iqnls,  &
                             coupling % W_current_history_iqnls)

            ! Recompute filled increments
            !
            filled_columns_increment_iqnls=filled_columns_inner_iqnls-1
            if(coupling % history_iqnls > 0_ip) then
              filled_columns_increment_iqnls= sum(history_tracking_iqnls(:))-filled_columns_history_iqnls+(filled_columns_inner_iqnls-1)
              if(filled_columns_increment_iqnls >= max_columns_increment_iqnls) filled_columns_increment_iqnls = max_columns_increment_iqnls
            endif

!IF(PAR_MY_CODE_RANK.eq.1) then
!write(6,*) 'A MATRIX POST FILTERING:'
!write(6,*) 'filled_increment', filled_columns_increment_iqnls
!write(6,*) 'V_current_history:'
!write(6,*) coupling % V_current_history_iqnls(5,:)
!write(6,*) 'residincr:', filled_columns_inner_iqnls
!write(6,*) coupling % residincr_iqnls(5,:)
!write(6,*) 'residincr_history:', history_tracking_iqnls(:)
!write(6,*) coupling % residincr_history_iqnls(5,:,:)
!write(6,*) '|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|'
!ENDIF

          allocate(alpha(coupling % ranku_iqnls*(coupling % history_iqnls +1) ))
          alpha = 0.0_rp

          call compute_alpha( coupling % V_current_history_iqnls , &
               coupling % residues_iqnls(:,1)  , &
               ndofs, filled_columns_increment_iqnls, alpha )

!IF(PAR_MY_CODE_RANK.eq.1) then
!write(6,*) '|------------------------------------------------|'
!write(6,*) ' COMPUTED_ALPHA'
!write(6,*) alpha(:)
!ENDIF

          ! Matrix vector product W*alpha
          ! used for the new prediction
          !
          vecaux=0.0_rp
          !call maths_matrix_vector_multiplication( coupling % valincr_iqnls(:,:), alpha(:), vecaux(:), filled_columns_inner_iqnls - 1 )
          call maths_matrix_vector_multiplication( coupling % W_current_history_iqnls(:,:), alpha(:), vecaux(:), filled_columns_increment_iqnls )

          ! New prediction
          ! x^k+1 = x^k + W*alpha -r^k
          !
          do idofs=1,ndofs
             coupling % relaxed_iqnls (idofs, 1) = coupling % relaxed_iqnls (idofs, 2) + coupling % residues_iqnls(idofs, 1) + vecaux(idofs)
          enddo

          !! Here finishes the differentiation between iter<3 and iter>3
       endif

       !! We must save again the result in xxnew
       !!
       do ipoin=1, npoin_wet
          do idofn=1, ndofn
             mindex = ipoin * ndofn + idofn - ndofn
             xxnew(idofn,ipoin) = coupling % relaxed_iqnls (mindex, 1) / coupling % scaling_iqnls
             coupling % values(idofn,ipoin,1) = xxnew(idofn,ipoin)
          enddo
       enddo

       !! Here  finishes the SCHEME_IQNLS
    end if  ! RELAXATION SCHEME

    !
    ! Save unrelaxed results for next aitken calculation (must be performed in all iterations)
    !
    if(( coupling % scheme == AITKEN_SCHEME)) then
       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             coupling % values_predicted(idofn,ipoin) = xxnew(idofn,ipoin)
          end do
       end do
    end if
    !
    ! Relaxation of the solution
    !
    if( coupling % scheme == AITKEN_SCHEME .and. coupling_driver_iteration(1_ip) > 2_ip ) then
       rela1    = 1.0_rp - relax
       do kpoin = 1,npoin_wet
          do idofn = 1,ndofn
             xxnew(idofn,kpoin)               = rela1 * coupling % values(idofn,kpoin,2_ip) + relax * xxnew(idofn,kpoin)
             coupling % values(idofn,kpoin,1) = xxnew(idofn,kpoin)
          end do
       end do
    else if( coupling % scheme == BROYDEN_SCHEME ) then

    else if ( (coupling % scheme == RELAXATION_SCHEME) .or. &
         (coupling % scheme == AITKEN_SCHEME .and. coupling_driver_iteration(1_ip) <= 2_ip ) ) then
       rela1    = 1.0_rp - relax
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          do idofn = 1,ndofn
             xxnew(idofn,kpoin)               = relax * xxnew(idofn,kpoin) + rela1 * coupling % values(idofn,kpoin,2)
             coupling % values(idofn,kpoin,1) = xxnew(idofn,kpoin)
          end do
       end do
    end if
    !
    ! Compute residual
    !
    xresi = 0.0_rp
    if( coupling % scheme == IQNLS_SCHEME ) then
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          if( THIS_NODE_IS_MINE(ipoin) ) then
             do idofn = 1,ndofn
                mindex = kpoin * ndofn + idofn - ndofn
                xresi(2) = xresi(2) +  coupling % relaxed_iqnls(mindex,1)**2
                xresi(1) = xresi(1) + (coupling % relaxed_iqnls(mindex,1)-coupling % relaxed_iqnls(mindex,2))**2
             end do
          end if
       end do
    else
       xfact = 1.0_rp
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          if( THIS_NODE_IS_MINE(ipoin) ) then
             do idofn = 1,ndofn
                if( present(mask) ) then
                   xfact = 1.0_rp
                   if(       coupling % what == RESIDUAL ) then
                      if( mask(idofn,ipoin) > 0 ) xfact = 0.0_rp
                   else if( coupling %  what == UNKNOWN ) then
                      if( mask(idofn,ipoin) > 0 .and. mask(idofn,ipoin) /= FIXED_UNKNOWN )  xfact = 0.0_rp
                  end if
                end if
                xresi(2) = xresi(2) + xfact *  coupling % values(idofn,kpoin,1)**2
                xresi(1) = xresi(1) + xfact * (coupling % values(idofn,kpoin,1)-coupling % values(idofn,kpoin,2))**2
             end do
          end if
       end do
    end if

  end subroutine COU_UPDATE_POINTS_VALUES

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Interpolate and modify nodal arrays
  !> @details Do the following
  !>                               Interpolate
  !>          XTARGET(NDOFN,NPOIN)     <=       XSOURCE(NDOFN,NPOIN)
  !>
  !>          1. Allocate XINTERP(NDOFN,NPOIN_WET)
  !>
  !>          2. Interpolate XINTERP(NDOFN,NPOIN_WET) FROM XSOURCE:
  !>
  !>                     COU_GET_INTERPOLATE_POINTS_VALUES
  !>             XINTERP              <=                   XSOURCE
  !>
  !>
  !>          3. Scatter solution:
  !>                     LPOIN_WET
  !>             XTARGET    <=    XINTERP
  !>
  !----------------------------------------------------------------------

  subroutine COU_INTERPOLATE_NODAL_VALUES_11(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                       :: icoup
    integer(ip), intent(in)                       :: ndofn
    real(rp),    intent(inout), pointer           :: xtarget(:)
    real(rp),    intent(in),    pointer, optional :: xsource(:)
    integer(ip), intent(in),    pointer, optional :: mask(:,:)
    type(typ_color_coupling),            optional :: my_coupling
    real(rp)                                      :: xtarget_tmp(2)
    real(rp)                                      :: xsource_tmp(2)
    logical(lg)                                   :: ltarget,lsource

    if( memory_size(xtarget)>0 ) then
       ltarget = .true.
    else
       ltarget = .false.
    end if
    if( present(xsource) ) then
       if( memory_size(xsource) > 0 ) then
          lsource = .true.
       else
          lsource = .false.
       end if
    else
       lsource = .false.
    end if
    
    if( ltarget .and. lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)
    else if( ltarget ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource_tmp,mask,my_coupling)
    else if( lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource,mask,my_coupling)
    else
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource_tmp,mask,my_coupling)
    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES_11

  subroutine COU_INTERPOLATE_NODAL_VALUES_12(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                       :: icoup
    integer(ip), intent(in)                       :: ndofn
    real(rp),    intent(inout), pointer           :: xtarget(:)
    real(rp),    intent(in),    pointer           :: xsource(:,:)
    integer(ip), intent(in),    pointer, optional :: mask(:,:)
    type(typ_color_coupling),            optional :: my_coupling
    real(rp)                                      :: xtarget_tmp(2)
    real(rp)                                      :: xsource_tmp(2)
    logical(lg)                                   :: ltarget,lsource

    if( memory_size(xtarget)>0 ) then
       ltarget = .true.
    else
       ltarget = .false.
    end if
    if( memory_size(xsource) > 0 ) then
       lsource = .true.
    else
       lsource = .false.
    end if
    
    if( ltarget .and. lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)
    else if( ltarget ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource_tmp,mask,my_coupling)
    else if( lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource,mask,my_coupling)
    else
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource_tmp,mask,my_coupling)
    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES_12

  subroutine COU_INTERPOLATE_NODAL_VALUES_22(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                       :: icoup
    integer(ip), intent(in)                       :: ndofn
    real(rp),    intent(inout), pointer           :: xtarget(:,:)
    real(rp),    intent(in),    pointer, optional :: xsource(:,:)
    integer(ip), intent(in),    pointer, optional :: mask(:,:)
    type(typ_color_coupling),            optional :: my_coupling
    real(rp)                                      :: xtarget_tmp(2)
    real(rp)                                      :: xsource_tmp(2)
    logical(lg)                                   :: ltarget,lsource

    if( memory_size(xtarget)>0 ) then
       ltarget = .true.
    else
       ltarget = .false.
    end if
    if( present(xsource) ) then
       if( memory_size(xsource) > 0 ) then
          lsource = .true.
       else
          lsource = .false.
       end if
    else
       lsource = .false.
    end if
    
    if( ltarget .and. lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)
    else if( ltarget ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource_tmp,mask,my_coupling)
    else if( lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource,mask,my_coupling)
    else
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource_tmp,mask,my_coupling)
    end if
    
  end subroutine COU_INTERPOLATE_NODAL_VALUES_22

  subroutine COU_INTERPOLATE_NODAL_VALUES_32(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                       :: icoup
    integer(ip), intent(in)                       :: ndofn
    real(rp),    intent(inout), pointer           :: xtarget(:,:,:)
    real(rp),    intent(in),    pointer, optional :: xsource(:,:)
    integer(ip), intent(in),    pointer, optional :: mask(:,:)
    type(typ_color_coupling),            optional :: my_coupling
    real(rp)                                      :: xtarget_tmp(2)
    real(rp)                                      :: xsource_tmp(2)
    logical(lg)                                   :: ltarget,lsource

    if( memory_size(xtarget)>0 ) then
       ltarget = .true.
    else
       ltarget = .false.
    end if
    if( present(xsource) ) then
       if( memory_size(xsource) > 0 ) then
          lsource = .true.
       else
          lsource = .false.
       end if
    else
       lsource = .false.
    end if
    
    if( ltarget .and. lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)
    else if( ltarget ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource_tmp,mask,my_coupling)
    else if( lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource,mask,my_coupling)
    else
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource_tmp,mask,my_coupling)
    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES_32

  subroutine COU_INTERPOLATE_NODAL_VALUES_33(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                       :: icoup
    integer(ip), intent(in)                       :: ndofn
    real(rp),    intent(inout), pointer           :: xtarget(:,:,:)
    real(rp),    intent(in),    pointer           :: xsource(:,:,:)
    integer(ip), intent(in),    pointer, optional :: mask(:,:)
    type(typ_color_coupling),            optional :: my_coupling
    real(rp)                                      :: xtarget_tmp(2)
    real(rp)                                      :: xsource_tmp(2)
    logical(lg)                                   :: ltarget,lsource

    if( memory_size(xtarget)>0 ) then
       ltarget = .true.
    else
       ltarget = .false.
    end if
    if( memory_size(xsource) > 0 ) then
       lsource = .true.
    else
       lsource = .false.
    end if
    
    if( ltarget .and. lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)
    else if( ltarget ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource_tmp,mask,my_coupling)
    else if( lsource ) then
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource,mask,my_coupling)
    else
       call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget_tmp,xsource_tmp,mask,my_coupling)
    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES_33
  
  subroutine COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,xtarget,xsource,mask,my_coupling)

    integer(ip), intent(in)                    :: icoup 
    integer(ip), intent(in)                    :: ndofn
    real(rp),    intent(inout)                 :: xtarget(ndofn,*)
    real(rp),    intent(in),          optional :: xsource(ndofn,*)
    integer(ip), intent(in), pointer, optional :: mask(:,:)
    type(typ_color_coupling),         optional :: my_coupling
    real(rp),                pointer           :: xinterp(:,:)
    real(rp)                                   :: xresi(2)
    integer(ip)                                :: ipoin,kpoin,idofn,icolo,npoin_wet
    integer(ip)                                :: subdomain_target
    integer(ip)                                :: subdomain_source
    integer(ip)                                :: zone_target
    integer(ip)                                :: code_target
    integer(ip)                                :: whatis
    logical(lg)                                :: subdo_coupling
    real(rp)                                   :: weight
    real(rp),    pointer                       :: xsource_loc(:,:)
    logical(lg)                                :: if_mask

    if( present(my_coupling) ) then
       subdomain_target = my_coupling % subdomain_target
       subdomain_source = my_coupling % subdomain_source
       color_target     = my_coupling % color_target
       color_source     = my_coupling % color_source
       code_target      = my_coupling % code_target
       zone_target      = my_coupling % zone_target
       npoin_wet        = my_coupling % wet % npoin_wet
       whatis           = my_coupling % what
    else
       subdomain_target = coupling_type(icoup) % subdomain_target
       subdomain_source = coupling_type(icoup) % subdomain_source
       color_target     = coupling_type(icoup) % color_target
       color_source     = coupling_type(icoup) % color_source
       code_target      = coupling_type(icoup) % code_target
       zone_target      = coupling_type(icoup) % zone_target
       npoin_wet        = coupling_type(icoup) % wet % npoin_wet
       whatis           = coupling_type(icoup) % what
    end if
    !
    ! Initialization
    !
    subdo_coupling = .false.
    if( zone_target == 0 ) then
       if( I_AM_IN_SUBD( subdomain_target) ) subdo_coupling = .true.
       if( I_AM_IN_SUBD( subdomain_source) ) subdo_coupling = .true.
       icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    else
       icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
    end if

    xresi        = 0.0_rp
    if( IMASTER ) npoin_wet = 0
    nullify(xinterp)
    !
    ! Options
    !
    if_mask = .false.
    if( present(mask) ) then
       if( associated(mask) ) if_mask = .true.
    end if

    if( icolo == color_target .or. icolo == color_source .or. subdo_coupling ) then
       !
       ! Allocate values
       !
       if( current_code == code_target .and. INOTMASTER ) then
          call memory_alloca(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp,ndofn,max(1_ip,npoin_wet))
       else
          call memory_alloca(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp,1_ip,1_ip)
       end if
       !
       ! Interpolate values from XSOURCE
       !
       if( present(xsource) ) then
          if( present(my_coupling) ) then
             call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource,xinterp,my_coupling)
          else
             call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource,xinterp,coupling_type(icoup))
          end if
       else
          nullify(xsource_loc)
          call memory_alloca(memor_cou,'XSOURCE_LOC','cou_interpolate_nodal_values',xsource_loc,ndofn,max(1_ip,npoin))
          if( INOTMASTER ) xsource_loc(1:ndofn,1:npoin) = xtarget(1:ndofn,1:npoin)
          if( present(my_coupling) ) then
             call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource_loc,xinterp,my_coupling)
          else
             call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource_loc,xinterp,coupling_type(icoup))
          end if
          call memory_deallo(memor_cou,'XSOURCE_LOC','cou_interpolate_nodal_values',xsource_loc)
       end if
       !
       ! Permute value to target XTARGET
       !
       if( current_code == code_target ) then
          !
          ! Update solution according to scheme (relaxation, Aitken, etc.)
          !
          if( .not. subdo_coupling ) call COU_UPDATE_POINTS_VALUES(xinterp,coupling_type(icoup),xresi,mask)

          if(    whatis == UNKNOWN            .or. &
               & whatis == DIRICHLET_IMPLICIT .or. &
               & whatis == DIRICHLET_EXPLICIT ) then

             if( if_mask ) then
                !
                ! Unknown with mask
                !
                do kpoin = 1,npoin_wet
                   if( present(my_coupling) ) then
                      ipoin = my_coupling % wet % lpoin_wet(kpoin)
                   else
                      ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   end if
                   do idofn = 1,ndofn
                      if( mask(idofn,ipoin) <= 0 .or. mask(idofn,ipoin) == FIXED_UNKNOWN ) then
                         xtarget(idofn,ipoin) = xinterp(idofn,kpoin)
                      end if
                   end do
                end do
             else
                !
                ! Unknown without mask
                !
                do kpoin = 1,npoin_wet
                   if( present(my_coupling) ) then
                      ipoin = my_coupling % wet % lpoin_wet(kpoin)
                   else
                      ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   end if
                   do idofn = 1,ndofn
                      xtarget(idofn,ipoin) = xinterp(idofn,kpoin)
                   end do
                end do
             end if

          else if( whatis == RESIDUAL ) then

             if( if_mask ) then
                !
                ! Force with mask
                !
                do kpoin = 1,npoin_wet
                   if( present(my_coupling) ) then
                      ipoin  = my_coupling % wet % lpoin_wet(kpoin)
                      weight = my_coupling % wet % weight_wet(kpoin)
                   else
                      ipoin  = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                      weight = coupling_type(icoup) % wet % weight_wet(kpoin)
                   end if
                   do idofn = 1,ndofn
                      if( mask(idofn,ipoin) <= 0 ) then
                         xtarget(idofn,ipoin) = xtarget(idofn,ipoin) + weight * xinterp(idofn,kpoin)
                      end if
                   end do
                end do
             else
                !
                ! Force without mask
                !
                do kpoin = 1,npoin_wet
                   if( present(my_coupling) ) then
                      ipoin  = my_coupling % wet % lpoin_wet(kpoin)
                      weight = my_coupling % wet % weight_wet(kpoin)
                   else
                      ipoin  = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                      weight = coupling_type(icoup) % wet % weight_wet(kpoin)
                   end if                   
                   do idofn = 1,ndofn
                      xtarget(idofn,ipoin) = xtarget(idofn,ipoin) + weight * xinterp(idofn,kpoin)
                   end do
                end do
             end if
          else
             call runend('WRONG TAG')
          end if
          !
          ! Conservation
          !
          if( if_mask ) then
             if( present(my_coupling) ) then
                call COU_CONSERVATION(my_coupling,ndofn,xtarget,mask)
             else
                call COU_CONSERVATION(coupling_type(icoup),ndofn,xtarget,mask)
             end if
          end if
       end if
       call memory_deallo(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp)
       !
       ! It should be in color icolo and jcolo!
       !
       if( .not. subdo_coupling ) then
          call PAR_SUM(2_ip,xresi,'IN CURRENT COUPLING')
          if( present(my_coupling) ) then
             if(xresi(2)<1.0e-10_rp ) then
                xresi(1) = 0.0_rp
             end if
             my_coupling % resid(2) = sqrt(xresi(2))
             my_coupling % resid(1) = sqrt(xresi(1)) / ( sqrt(xresi(2)) + zeror )
          else
             coupling_type(icoup) % resid(2) = sqrt(xresi(2))
             coupling_type(icoup) % resid(1) = sqrt(xresi(1)) / ( sqrt(xresi(2)) + zeror )
          end if
       end if

    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES_GO
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Compute scalar coupling residual
  !> @details Compute scalar coupling residual
  !>
  !----------------------------------------------------------------------

  subroutine cou_residual_scalar(my_coupling)

    type(typ_color_coupling), intent(inout) :: my_coupling
    integer(ip)                             :: icoup

    icoup = my_coupling % number
    my_coupling % resid(1) = abs(scala_cou(icoup,1)-scala_cou(icoup,2))
    scala_cou(icoup,2)     = scala_cou(icoup,1)

  end subroutine cou_residual_scalar
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Conservation
  !> @details Apply a conservation algorithm
  !>
  !----------------------------------------------------------------------

  subroutine COU_CONSERVATION(coupling,ndofn,xtarget,mask)

    type(typ_color_coupling), intent(in)    :: coupling
    integer(ip),              intent(in)    :: ndofn
    integer(ip),              intent(in)    :: mask(ndofn,*)
    real(rp),                 intent(inout) :: xtarget(ndofn,*)
    integer(ip)                             :: kboun,iboun,inodb,ipoin
    integer(ip)                             :: idofn
    logical(lg)                             :: mark_node
    logical(lg),              pointer       :: gboun(:)
    logical(lg),              pointer       :: gpoin(:)

    if( coupling % conservation /= 0 ) then

       if( INOTMASTER ) then
          allocate( gboun(nboun) )
          allocate( gpoin(npoin) )
          gboun = .false.
          gpoin = .false.
          do kboun = 1,coupling % wet % nboun_wet
             iboun = coupling % wet % lboun_wet(kboun)
             gboun(iboun) = .true.
             do inodb = 1,lnnob(iboun)
                ipoin = lnodb(inodb,iboun)
                mark_node = .true.
                loop_idofn: do idofn = 1,ndofn
                   if( mask(idofn,ipoin) /= FIXED_UNKNOWN ) then
                      mark_node = .false.
                      exit loop_idofn
                   end if
                end do loop_idofn
                if( mark_node ) gpoin(ipoin) = .true.
             end do
          end do
       else
          allocate( gboun(1) )
          allocate( gpoin(1) )
       end if
       if( coupling % conservation == INTERFACE_MASS ) then
          !call projec_mass_conservation(xtarget,gboun,gpoin,'LOCAL MASS')
       else if( coupling % conservation == GLOBAL_MASS ) then
          !call projec_mass_conservation(xtarget,gboun,gpoin,'GLOBAL MASS')
       else
          call runend('COU_CONSERVATION: CONSERVATION NOT CODED')
       end if
       deallocate( gboun )
       deallocate( gpoin )

    end if

  end subroutine COU_CONSERVATION

  subroutine COU_RESIDUAL_FORCE(ndofn,ia,ja,amatr,rhsid,unkno,force_rhs)
    use def_kintyp, only        :  ip,rp
    use def_master, only        :  INOTMASTER
    use def_domain, only        :  npoin
    implicit none
    integer(ip),    intent(in)  :: ndofn
    integer(ip),    intent(in)  :: ia(*)
    integer(ip),    intent(in)  :: ja(*)
    real(rp),       intent(in)  :: amatr(ndofn,ndofn,*)
    real(rp),       intent(in)  :: rhsid(ndofn,*)
    real(rp),       intent(in)  :: unkno(ndofn,*)
    real(rp),       intent(out) :: force_rhs(ndofn,*)
    integer(ip)                 :: ipoin,iz,jpoin,idofn,jdofn
    !
    ! Marcar unicamente los nodos que tocan los elementos
    !
    if( INOTMASTER ) then
       do ipoin = 1,npoin
          force_rhs(1:ndofn,ipoin) = rhsid(1:ndofn,ipoin)
          do iz = ia(ipoin),ia(ipoin+1)-1
             jpoin = ja(iz)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   force_rhs(idofn,ipoin) = &
                        force_rhs(idofn,ipoin) - amatr(jdofn,idofn,iz)*unkno(jdofn,jpoin)
                end do
             end do
          end do
       end do
       call rhsmod(ndofn,force_rhs)
    end if

  end subroutine COU_RESIDUAL_FORCE

  subroutine COU_PRESCRIBE_DIRICHLET_IN_MATRIX(nbvar,npopo,ia,ja,an,bb,xx)
    integer(ip), intent(in)    :: nbvar
    integer(ip), intent(in)    :: npopo
    integer(ip), intent(in)    :: ia(*)
    integer(ip), intent(in)    :: ja(*)
    real(rp),    intent(inout) :: an(nbvar,nbvar,*)
    real(rp),    intent(inout) :: xx(nbvar,*)
    real(rp),    intent(inout) :: bb(nbvar,*)
    integer(ip)                :: nrows,icoup,ii,jj,kk,nn
    integer(ip)                :: izdom
    real(rp),    pointer       :: xx_tmp(:,:)

    if(    INOTMASTER .and. ( &
         & I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_IMPLICIT) .or. &
         & I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_EXPLICIT) ) ) then

       nrows = nbvar * npopo
       if( nbvar > 1 ) call runend('COU_PRESCRIBE_DIRICHLET_IN_MATRIX: NOT CODED')
       nullify( xx_tmp )
       allocate( xx_tmp(nbvar,npopo) )
       do ii = 1,npopo
          do nn = 1,nbvar
             xx_tmp(nn,ii) = xx(nn,ii)
          end do
       end do
       do icoup = 1,mcoup
          if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
             if( coupling_type(icoup) % what == DIRICHLET_EXPLICIT .or. coupling_type(icoup) % what == DIRICHLET_IMPLICIT ) then
                call COU_INTERPOLATE_NODAL_VALUES_go(icoup,1_ip,xx,xx_tmp)
             end if
          end if
       end do
       deallocate( xx_tmp )
       do icoup = 1,mcoup
          if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
             if( coupling_type(icoup) % what == DIRICHLET_EXPLICIT .or. coupling_type(icoup) % what == DIRICHLET_IMPLICIT ) then
                do kk = 1,coupling_type(icoup) % wet % npoin_wet
                   ii  = coupling_type(icoup) % wet % lpoin_wet(kk)
                   do izdom = ia(ii),ia(ii+1)-1
                      jj = ja(izdom)
                      if( ii == jj ) then
                         an(1,1,izdom) = 1.0_rp
                      else
                         an(1,1,izdom) = 0.0_rp
                      end if
                   end do
                   bb(1,ii) = xx(1,ii)
                end do
             end if
          end if
       end do
    end if

  end subroutine COU_PRESCRIBE_DIRICHLET_IN_MATRIX

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   List of wet nodes
  !> @details Give the list of source node for a specific coupling
  !>          of for all couplings
  !
  !----------------------------------------------------------------------

  subroutine COU_LIST_SOURCE_NODES(list_source_nodes,kcoup,what_to_do)
    logical(lg),  intent(inout), pointer  :: list_source_nodes(:)
    integer(ip),  intent(in),    optional :: kcoup
    character(*), intent(in),    optional :: what_to_do
    integer(ip)                           :: icoup_ini,icoup_end,icoup
    integer(ip)                           :: kpoin,ipoin
    !
    ! Allocate of necessary
    !
    if( associated(list_source_nodes) ) then
       if( present(what_to_do) ) then
          if( trim(what_to_do) == 'INITIALIZE' ) then
             do ipoin = 1,size(list_source_nodes)
                list_source_nodes(ipoin) = .false.
             end do
          end if
       end if
    else
       call memory_alloca(memor_cou,'LIST_SOURCE_NODES','cou_list_source_nodes',list_source_nodes,npoin)
    end if
    !
    ! Bounds
    !
    if( present(kcoup) ) then
       if( kcoup /= 0 ) then
          icoup_ini = kcoup
          icoup_end = kcoup
       else
          icoup_ini = 1
          icoup_end = mcoup
       end if
    else
       icoup_ini = 1
       icoup_end = mcoup
    end if
    !
    ! Activate node
    !
    do icoup = icoup_ini,icoup_end
       do kpoin = 1,coupling_type(icoup) % geome % npoin_source
          ipoin = coupling_type(icoup)  % geome % lpoin_source(kpoin)
          list_source_nodes(ipoin) = .true.
       end do
    end do

  end subroutine COU_LIST_SOURCE_NODES

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If there is at least one zone coupling
  !> @details If there is at least one zone coupling
  !
  !----------------------------------------------------------------------

  function THERE_EXISTS_A_ZONE_COUPLING()
    integer(ip) :: icoup
    logical(lg) :: THERE_EXISTS_A_ZONE_COUPLING

    THERE_EXISTS_A_ZONE_COUPLING = .false.
    do icoup = 1,mcoup
       if( coupling_type(icoup) % zone_target + coupling_type(icoup) % zone_source /= 0 ) then
          THERE_EXISTS_A_ZONE_COUPLING = .true.
          return
       end if
    end do
  end function THERE_EXISTS_A_ZONE_COUPLING

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   Initialize a value on target nodes
  !> @details Initialize a value on target nodes only if the coupling
  !>          is not on whole mesh
  !>
  !----------------------------------------------------------------------

  subroutine COU_PUT_VALUE_ON_TARGET_IP_1(value_in,xarray,kcoup,FORCE)
    
    integer(ip), intent(in)              :: value_in
    integer(ip), intent(inout), pointer  :: xarray(:)
    integer(ip), intent(in),    optional :: kcoup
    logical(lg), intent(in),    optional :: FORCE
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: kpoin,ipoin,icoup

    if( mcoup > 0 ) then
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       if( optional_argument(.false.,FORCE) ) then
          do icoup = icoup_ini,icoup_fin
             if(    (( coupling_type(icoup) % where_type /= ON_WHOLE_MESH      ) &
                  .and. ( coupling_type(icoup) % where_type /= ON_FLOATING_POINTS )) &
                  .or. ( coupling_type(icoup) % where_type == ON_IMMERSED_MESH ) )then
                do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                   ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   if( ipoin > 0 ) xarray(ipoin) = value_in
                end do
             end if
          end do 

       else
          do icoup = icoup_ini,icoup_fin
             if(  coupling_type(icoup) %  module_target == modul .or. coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                if(    (( coupling_type(icoup) % where_type /= ON_WHOLE_MESH      ) &
                     .and. ( coupling_type(icoup) % where_type /= ON_FLOATING_POINTS )) &
                     .or. ( coupling_type(icoup) % where_type == ON_IMMERSED_MESH ) )then
                   do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                      ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                      if( ipoin > 0 ) xarray(ipoin) = value_in
                   end do
                end if
             end if
          end do
       end if
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_1

  subroutine COU_PUT_VALUE_ON_TARGET_IP_2(value_in,xarray,kcoup)
    integer(ip), intent(in)              :: value_in(*)
    integer(ip), intent(inout), pointer  :: xarray(:,:)
    integer(ip), intent(in),    optional :: kcoup
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: ndofn,kpoin,ipoin,icoup

    if( mcoup > 0 ) then
       ndofn = size(xarray,1)
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       do icoup = icoup_ini,icoup_fin
          if(  coupling_type(icoup) %  module_target == modul .or. coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
             if(    ( coupling_type(icoup) % where_type /= ON_WHOLE_MESH    ) &
                  .or. ( coupling_type(icoup) % where_type == ON_IMMERSED_MESH ) )then
                do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                   ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   if( ipoin > 0 ) xarray(1:ndofn,ipoin) = value_in(1:ndofn)
                end do
             end if
          end if
       end do
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_2

  subroutine COU_PUT_VALUE_ON_TARGET_IP_12(value_in,xarray,kcoup)
    integer(ip), intent(in)              :: value_in
    integer(ip), intent(inout), pointer  :: xarray(:,:)
    integer(ip), intent(in),    optional :: kcoup
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: ndofn,kpoin,ipoin,icoup

    if( mcoup > 0 ) then
       ndofn = size(xarray,1)
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       do icoup = icoup_ini,icoup_fin
          if(  coupling_type(icoup) %  module_target == modul ) then
             if(    ( coupling_type(icoup) % where_type /= ON_WHOLE_MESH    ) &
                  .or. ( coupling_type(icoup) % where_type == ON_IMMERSED_MESH ) )then
                do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                   ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   if( ipoin > 0 ) xarray(1:ndofn,ipoin) = value_in
                end do
             end if
          end if
       end do
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_12

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    14/10/2014
  !> @brief   Change fixity array
  !> @details Modify fixity array on target:
  !>          - Coupling between zones:
  !>            UNKNOWN type: set fixity to FIXED_UNKNOWN whenener it
  !>            is different from 0. Can be forced by using "FORCE"
  !>          - Coupling between subdomains:
  !>            RESIDUAL type: free the nodes if FORCE SOBDOMAIN
  !             is present
  !>
  !----------------------------------------------------------------------

  subroutine COU_SET_FIXITY_ON_TARGET(variable,imodu,kfl_fixno,what)
    character(*), intent(in)              :: variable
    integer(ip),  intent(in)              :: imodu
    integer(ip),  intent(inout), pointer  :: kfl_fixno(:,:)
    character(*), intent(in), optional    :: what
    integer(ip)                           :: ipoin,kpoin,npoin_wet
    integer(ip)                           :: ndofn,kdofn,kpoin_fixed
    integer(ip)                           :: icoup,iffix_max,idofn
    integer(ip)                           :: code_target,zone_target
    logical(lg)                           :: force_subdomain

    if( mcoup == 0 ) return
    !
    ! What to do
    !
    iffix_max       = 0
    ndofn           = size(kfl_fixno,1)
    force_subdomain = .false.
    if( present(what) ) then
       if( trim(what) == 'FORCE ZONE' ) then
          iffix_max = huge(1_ip)
       else if( trim(what) == 'FREE BETWEEN SUBDOMAINS' ) then
          force_subdomain = .true.
       else if( trim(what) == 'FREE FIXITY' ) then
          iffix_max = -1_ip 
       end  if
    end if

    do icoup = 1,mcoup

       if(  coupling_type(icoup) %  module_target == imodu .or. coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then

          if(    ( coupling_type(icoup) % where_type /= ON_WHOLE_MESH    ) &
               .or. ( coupling_type(icoup) % where_type == ON_IMMERSED_MESH ) )then

             code_target  = coupling_type(icoup) % code_target
             zone_target  = coupling_type(icoup) % zone_target
             color_target = coupling_type(icoup) % color_target
             npoin_wet    = coupling_type(icoup) % wet % npoin_wet

             if( I_AM_IN_COLOR(color_target) ) then

                if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                   !
                   ! Between subdomains: just check not all dofs are prescribed
                   !
                   if( force_subdomain ) then
                      kpoin_fixed = 0
                      do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                         ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                         kfl_fixno(1:ndofn,ipoin) = 0
                      end do
                   else
                      kpoin_fixed = 0
                      do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                         ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                         if( minval(kfl_fixno(1:ndofn,ipoin)) > 0 ) kpoin_fixed = kpoin_fixed + 1
                      end do
                      if( kpoin_fixed == npoin_wet ) then
                         kpoin = 1
                      else
                         kpoin = 0
                      end if
                      call PAR_MIN(kpoin,'IN CURRENT TARGET COLOR')
                      if( kpoin == 1 ) then
                         call runend('ALL DOFS ARE PRESCIBED ON INTERFACE')
                      end if
                   end if

                else if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. lzone(imodu) == zone_target  ) then
                   !
                   ! Between zones
                   !
                   ! Unknown type:  force iffix
                   ! Residual type: just check not all dofs are prescribed
                   !
                   if( coupling_type(icoup) % variable == variable(1:5) ) then

                      if( coupling_type(icoup) % what == UNKNOWN  ) then

                         do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                            ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                            do idofn = 1,ndofn
                               if( kfl_fixno(idofn,ipoin) <= iffix_max ) kfl_fixno(idofn,ipoin) = FIXED_UNKNOWN
                            end do
                         end do

                      else if( coupling_type(icoup) % what == RESIDUAL ) then

                         kpoin_fixed = 0
                         do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                            ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                            kdofn = 0
                            if( minval(kfl_fixno(1:ndofn,ipoin)) > 0 ) kpoin_fixed = kpoin_fixed + 1
                         end do
                         if( kpoin_fixed == npoin_wet ) then
                            kpoin = 1
                         else
                            kpoin = 0
                         end if
                         call PAR_MIN(kpoin,'IN CURRENT TARGET COLOR')
                         if( kpoin == 1 ) then
                            print*,'popopo=',current_code,' ',variable
                            call runend('ALL DOFS ARE PRESCIBED ON INTERFACE')
                         end if

                      end if

                   end if

                end if

             end if

          end if

       end if

    end do

  end subroutine COU_SET_FIXITY_ON_TARGET

  !----------------------------------------------------------------------
  !>
  !> @author  J.C. Cajas
  !> @date    02/06/2016
  !> @brief   Predict interface values
  !> @details Predict coupling values when advancing in time
  !>
  !>
  !----------------------------------------------------------------------
  subroutine COU_TEMPORAL_PREDICTOR(icoup)

    integer(ip), intent(in) :: icoup
    integer(ip)             :: idime, kpoin, npoin_wet

    if( coupling_type(icoup) % temporal_predictor /= 1_ip ) return

    npoin_wet = coupling_type(icoup) % wet % npoin_wet

    if( .not. associated( coupling_type(icoup) % values_converged ) .and. npoin_wet > 0 ) &
         &call runend(' Sure you made two time steps before making a prediction? see cou_update_points_values, runend from cou_temporal_predictor ')
    !
    ! Call madame Sazu and make your prediction, it will be stored in
    ! coupling_type(icoup) % values_converged(:,:,1) in order to reuse the array
    !
    ! prediction(k+1) = 5/2 * values_conv(k) - 2 * values_conv(k-1) + 1/2 * values_conv(k-2)
    !

    if( coupling_type(icoup) % temporal_predictor_order == 0_ip ) then
      return ! Value is the one of the last time step. Saved in convergence subroutine

    elseif( coupling_type(icoup) % temporal_predictor_order == 1_ip ) then
      call runend('COU_TEMPORAL_PREDICTOR: first order predictor not coded')

    elseif( coupling_type(icoup) % temporal_predictor_order == 2_ip ) then
      do kpoin = 1_ip, npoin_wet
        do idime = 1_ip, ndime
          coupling_type(icoup) % values_converged(idime,kpoin,1_ip) = 2.5_rp * coupling_type(icoup) % values_converged(idime,kpoin,1_ip) - &
               & 2.0_rp * coupling_type(icoup) % values_converged(idime,kpoin,2_ip) + 0.5_rp * coupling_type(icoup) % values_converged(idime,kpoin,3_ip)
        end do
      end do

    else
      call runend('MOD_COUPLINGS: temporal predictor order not programmed')

    endif

  end subroutine COU_TEMPORAL_PREDICTOR

      !----------------------------------------------------------------------
      !>
      !> @author  J.C. Cajas
      !> @date    02/06/2016
      !> @brief   Q-N Broyden ('bad') approximation
      !> @details Broyden (bad) Q-N approximation
      !>
      !>
      !----------------------------------------------------------------------
  subroutine COU_BROYDEN_BAD(xxnew_o,coupli_o)

    real(rp),      pointer,   intent(inout) :: xxnew_o(:,:)
    type(typ_color_coupling), intent(inout) :: coupli_o
    integer(ip)                             :: ndofn,ipoin,npoin_wet,npoin_wet_total
    integer(ip)                             :: idofn,jpoin
    integer(ip)                             :: initial_index, itime
    integer(ip)                             :: ntime_wet,ipart
    integer(ip)                             :: target_rank, target_size

    ! Broyden auxiliary variables
    real(rp), pointer                       :: deltaf(:), deltag(:), deltav(:), auxjac(:,:)
    real(rp), pointer                       :: my_dex(:), my_def(:), my_dev(:)
    real(rp)                                :: normv, relax

    integer(ip), pointer                    :: npoin_wet_all(:)

    nullify(my_dex)
    nullify(my_def)
    nullify(my_dev)

    nullify(deltaf)
    nullify(deltav)

    nullify(deltag)
    nullify(auxjac)

    nullify(npoin_wet_all)

    if( INOTMASTER ) then
       !
       ! Degrees of freedom
       !
       ndofn     = size(xxnew_o,1)
       npoin_wet = coupli_o % wet % npoin_wet
       ntime_wet = 3

    else

       ndofn     = 0
       npoin_wet = 0
       ntime_wet = 0

    end if
    !
    ! Size of the whole wet surface
    !
    if( .not. associated(npoin_wet_all) )then

       call PAR_COMM_RANK_AND_SIZE(target_rank,target_size,'IN CURRENT TARGET COLOR')
       allocate(npoin_wet_all(0:target_size-1))
       !call memory_alloca(memor_cou,vacal,' npoin_all ',npoin_wet_all,target_size)

       call PAR_ALLGATHER(npoin_wet,npoin_wet_all,1_4,'IN CURRENT TARGET COLOR')

       npoin_wet_total = 0_rp
       do ipart = 0_ip, target_size-1_ip
          if( ipart == target_rank )initial_index = npoin_wet_total
          npoin_wet_total = npoin_wet_total + npoin_wet_all(ipart)
       end do

       npoin_wet_all = ndime * npoin_wet_all

    end if
    ! print*, "DEBUG: npoin_wet ", npoin_wet_all, target_rank, npoin_wet_total
    !
    ! Memory allocation
    !
    if( INOTMASTER )then

       if( .not. associated(coupli_o % values) )          &
            call memory_alloca(memor_cou,'values',vacal,coupli_o % values,ndofn,npoin_wet,3_ip)
       if( .not. associated(coupli_o % values_predicted) )&
            call memory_alloca(memor_cou,'values_predicted',vacal,coupli_o % values_predicted,ndofn,npoin_wet)

       if( .not. associated(auxjac) )                     &
            call memory_alloca(memor_cou,'auxjac',vacal,auxjac ,ndofn * npoin_wet, ndofn * npoin_wet_total)
       if( .not. associated(coupli_o % jacobian_inverse) )&
            call memory_alloca(memor_cou,'jacobian_inverse',vacal,coupli_o % jacobian_inverse,ndofn,npoin_wet,npoin_wet_total)

       if( .not. associated(my_dex) )                     &
            call memory_alloca(memor_cou,'my_dex',vacal, my_dex ,ndofn * npoin_wet)
       if( .not. associated(my_def) )     &
            call memory_alloca(memor_cou,'my_def',vacal, my_def ,ndofn * npoin_wet)
       if( .not. associated(deltaf) )     &
            call memory_alloca(memor_cou,'deltaf',vacal, deltaf ,ndofn * npoin_wet_total)

       if( .not. associated(my_dev) )     &
            call memory_alloca(memor_cou,'my_dev',vacal, my_dev ,ndofn * npoin_wet)
       if( .not. associated(deltav) )     &
            call memory_alloca(memor_cou,'deltav',vacal, deltav ,ndofn * npoin_wet_total)

       if( .not. associated(deltag) )     &
            call memory_alloca(memor_cou,'deltag',vacal, deltag ,ndofn * npoin_wet)

    else

       if( .not. associated(coupli_o % values) )             call memory_alloca_min(memor_cou,'values'           ,vacal, coupli_o % values )
       if( .not. associated(coupli_o % values_predicted) )   call memory_alloca_min(memor_cou,'values_predicted' ,vacal, coupli_o % values_predicted )
       if( .not. associated(coupli_o % jacobian_inverse) )   call memory_alloca_min(memor_cou,'jacobian_inverse' ,vacal, coupli_o % jacobian_inverse )

       call memory_alloca_min(memor_cou,'MY_DEX' ,vacal           ,  my_dex )
       call memory_alloca_min(memor_cou,'MUY_DEF',vacal           ,  my_def )
       call memory_alloca_min(memor_cou,'MY_DEV' ,vacal           ,  my_dev )

       if( .not. associated(deltaf) )                     &
            call memory_alloca(memor_cou,'deltaf',vacal,deltaf ,ndime * npoin_wet_total)
       if( .not. associated(deltav) )                     &
            call memory_alloca(memor_cou,'deltav',vacal,deltav ,ndime * npoin_wet_total)

       call memory_alloca_min(memor_cou,'deltag',vacal, deltag)
       call memory_alloca_min(memor_cou,'auxjac',vacal, auxjac)

       npoin_wet_total = 0_rp

    end if
    !
    ! Save old relaxed values
    !
    do itime = ntime_wet,2,-1
       do ipoin = 1,npoin_wet
          do idofn = 1,ndofn
             coupli_o % values(idofn,ipoin,itime) = coupli_o % values(idofn,ipoin,itime-1)
          end do
       end do
    end do
    !
    ! Search the roots of f(x)-x
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          xxnew_o(idofn,ipoin)                 = xxnew_o(idofn,ipoin) - coupli_o % values(idofn,ipoin,2)
          my_dev( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin)
       end do
    end do
    !
    ! First iteration is performed with an initial guess of the Jacobian inverse
    !
    if( coupling_driver_iteration(1_ip) < 2_ip ) then

       relax = coupli_o % relax

       if( coupling_driver_iteration(1_ip) == 1_ip ) then

          do ipoin = 1_ip,npoin_wet
             do jpoin = 1_ip, npoin_wet_total
                do idofn = 1_ip,ndofn
                   coupli_o % jacobian_inverse(idofn,ipoin,jpoin) = 0_rp
                end do
             end do
          end do

          do ipoin = 1_ip,npoin_wet
             do idofn = 1_ip,ndofn
                coupli_o % jacobian_inverse(idofn,ipoin,initial_index+ipoin) =-relax
             end do
          end do

       end if

       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn

             my_dex( (ipoin-1_ip)*ndofn + idofn ) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,3_ip)
             my_def( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin) - coupli_o % values_predicted(idofn,ipoin)

          end do
       end do

       call PAR_ALLGATHERV(my_dev,deltav,npoin_wet_all,'IN CURRENT TARGET COLOR')

    else

       ! print*, "DEBUG: mas de dos iter ", PAR_MY_CODE_RANK
       !
       ! Auxiliary vectors to compute the correction of the inverse jacobian
       !
       normv  = 0_rp
       do ipoin = 1_ip, npoin_wet_total
          deltaf(ipoin) = 0_rp
          deltav(ipoin) = 0_rp
       end do
       do ipoin = 1_ip, npoin_wet
          deltag(ipoin) = 0_rp
       end do

       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             my_dex( (ipoin-1_ip)*ndofn + idofn ) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,3_ip)
             my_def( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin) - coupli_o % values_predicted(idofn,ipoin)

             ! deltaf( (initial_index+ipoin-1_ip)*ndofn + idofn ) =  my_def( (ipoin-1_ip)*ndofn + idofn )
             ! deltav( (initial_index+ipoin-1_ip)*ndofn + idofn ) =  my_dev( (ipoin-1_ip)*ndofn + idofn )

          end do
       end do

       call PAR_ALLGATHERV(my_def,deltaf,npoin_wet_all,'IN CURRENT TARGET COLOR')
       call PAR_ALLGATHERV(my_dev,deltav,npoin_wet_all,'IN CURRENT TARGET COLOR')

       do ipoin = 1_ip, npoin_wet_total
          do idofn = 1_ip, ndofn
             normv = normv + deltaf( (ipoin-1_ip)*ndofn + idofn ) * deltaf( (ipoin-1_ip)*ndofn + idofn )
          end do
       end do

       normv = 1_rp / ( normv + zeror )

       do ipoin = 1_ip, npoin_wet
          do jpoin = 1_ip, npoin_wet_total
             do idofn = 1_ip, ndofn
                deltag( (ipoin-1_ip)*ndofn + idofn ) = deltag( (ipoin-1_ip)*ndofn + idofn ) + &
                     & coupli_o % jacobian_inverse(idofn,ipoin,jpoin) * deltaf((jpoin-1_ip)*ndofn + idofn )
             end do
          end do
       end do
       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             deltag((ipoin-1_ip)*ndofn + idofn ) = ( my_dex( (ipoin-1_ip)*ndofn + idofn ) - &
                  &deltag( (ipoin-1_ip)*ndofn + idofn ) ) * normv
          end do
       end do
       auxjac = 0_rp
       if( INOTMASTER .and. npoin_wet /= 0 )then
          call maths_outer_product(deltag,deltaf,auxjac)
       end if
       !
       ! Update the Jacobian approximation
       !
       do jpoin = 1_ip, npoin_wet_total
          do ipoin = 1_ip, npoin_wet
             do idofn = 1_ip, ndofn

                coupli_o % jacobian_inverse(idofn,ipoin,jpoin) = coupli_o % jacobian_inverse(idofn,ipoin,jpoin) + &
                     &auxjac( (ipoin-1_ip)*ndofn+idofn ,(jpoin-1_ip)*ndofn+idofn  )

             end do
             ! print*, "DEBUG: jacobian_inverse ", ipoin, jpoin, coupli_o % jacobian_inverse(:,ipoin,jpoin)
          end do
       end do

    end if      ! Iterations greater thar 1
    !
    ! Save unmodified results
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          coupli_o % values_predicted(idofn,ipoin) = xxnew_o(idofn,ipoin)
       end do
    end do
    !
    ! Update of the solution
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn

          coupli_o % values(idofn,ipoin,1_ip) = 0_rp

       end do
    end do
    do ipoin = 1_ip,npoin_wet
       do jpoin = 1_ip,npoin_wet_total
          do idofn = 1_ip,ndofn
             coupli_o % values(idofn,ipoin,1_ip) = coupli_o % values(idofn,ipoin,1_ip)+&
                  coupli_o % jacobian_inverse(idofn,ipoin,jpoin) * deltav( (jpoin-1_ip)*ndofn+idofn )
          end do
       end do
    end do

    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          coupli_o % values(idofn,ipoin,1_ip) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,1_ip)
          xxnew_o(idofn,ipoin) = coupli_o % values(idofn,ipoin,1_ip)
       end do

    end do
    ! print*, "DEBUG: salgo Broyden ", PAR_MY_CODE_RANK

    if( associated(my_def) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',my_def )
    if( associated(my_dev) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltag ',my_dev )

    ! if( associated(deltax) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltax ',deltax )
    if( associated(deltaf) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',deltaf )
    if( associated(deltav) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',deltav )

    if( associated(deltag) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltag ',deltag )
    if( associated(auxjac) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' auxjac ',auxjac )

    if( associated(npoin_wet_all) ) call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' npoin_wet_all ',npoin_wet_all )

  end subroutine COU_BROYDEN_BAD

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize mask for dot product
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_initialize_solver_mask(solve)

    use def_kintyp_solvers, only : soltyp
    use def_coupli, only : ncoup_implicit_n
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : mask_cou
    type(soltyp), intent(inout) :: solve          !< Solver structure

    if( ncoup_implicit_n + ncoup_implicit_d > 0 ) then
       solve % kfl_mask =  1
       solve % mask     => mask_cou
    end if

  end subroutine couplings_initialize_solver_mask

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize mask for dot product
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_initialize_solver_dirichlet_condition(solve)

    use def_kintyp_solvers, only : soltyp
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : lcoup_implicit_d
    type(soltyp), optional, intent(inout) :: solve          !< Solver structure
    integer(ip)                           :: kcoup,icoup

    if( INOTMASTER .and. mcoup > 0 ) then
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup)
        !!!   solve % kfl_dirichlet = 2
        end do
     end if

   end subroutine couplings_initialize_solver_dirichlet_condition

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Detect mnodes where reaction is required
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_initialize_solver_reaction(solve)

    use def_kintyp, only : ip
    use def_kintyp_solvers, only : soltyp
    use def_master, only : modul
    use def_master, only : mem_modul
    use def_domain, only : npoin
    use mod_parall, only : color_target
    use mod_parall, only : color_source

    type(soltyp), intent(inout), pointer :: solve(:)          !< Solver structure
    integer(ip)                          :: ireaction
    integer(ip)                          :: num_blocks
    integer(ip)                          :: ndofn,icoup
    integer(ip)                          :: kblok
    integer(ip)                          :: ndofn_block
    logical(lg)                          :: kfl_residual

    ireaction =  0
    ndofn      = solve(1) % ndofn
    num_blocks = solve(1) % num_blocks

    if( solve(1) % block_num == 1 ) then
       !
       ! SOLVE(1) % LPOIN_REACTION(1:NPOIN): Mark the nodes where reaction is required
       ! They are the source nodes
       !
       do icoup = 1,mcoup
          kfl_residual = coupling_type(icoup) % what       == RESIDUAL      .and. &
                       & coupling_type(icoup) % kind       == BETWEEN_ZONES
          if ( kfl_residual .or. kfl_efect ) then
             color_target = coupling_type(icoup) % color_target
             color_source = coupling_type(icoup) % color_source

             if( ireaction == 0 ) then
                ! TODO: Check if commenting this if-statement doesn't fuck up other cases!
                ! if( I_AM_IN_COLOR(color_source) ) then
                   solve(1) % kfl_react = max(solve(1) % kfl_react,1_ip)
                   call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % LPOIN_REACTION','inivar',solve(1) % lpoin_reaction,npoin)
                ! end if

                if( I_AM_IN_COLOR(color_target) ) then
                   solve(1) % kfl_bvnat = 1
                   call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % BVNAT','inivar',solve(1) % block_array(1) % bvnat,ndofn,npoin)
                   solve(1) % bvnat => solve(1) % block_array(1) % bvnat
                   do kblok = 2,num_blocks
                      ndofn_block = solve(1) % block_dimensions(kblok)
                      call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % BVNAT','inivar',solve(1) % block_array(kblok) % bvnat,ndofn_block,npoin)
                   end do
                end if

             end if

             ireaction = ireaction + 1

             if( INOTMASTER ) call COU_LIST_SOURCE_NODES(solve(1) % lpoin_reaction,icoup)
          end if
       end do
    end if

  end subroutine couplings_initialize_solver_reaction

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize solver according to coupling options
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_initialize_solver()

    use def_kintyp, only : ip
    use def_kintyp_solvers, only : soltyp
    use def_master, only : modul,momod,mmodu
    use def_master, only : current_zone
    use def_master, only : I_AM_IN_ZONE
    use def_master, only : lzone
    use def_master, only : kfl_modul
    use def_solver, only : solve_sol

    integer(ip) :: ivari,jvari

    do modul = 1,mmodu
       current_zone = lzone(modul)
       if( kfl_modul(modul) == 1 .and. associated(momod(modul) % solve) .and. I_AM_IN_ZONE(current_zone) ) then
          do ivari = 1,size(momod(modul) % solve)
             solve_sol => momod(modul) % solve(ivari:)
             call couplings_initialize_solver_mask(solve_sol(1))
             call couplings_initialize_solver_dirichlet_condition(solve_sol(1))
          end do
          jvari = 1
          do while( jvari <= size(momod(modul) % solve) )
             solve_sol => momod(modul) % solve(jvari:)
             call couplings_initialize_solver_reaction(solve_sol)
             jvari = jvari + momod(modul) % solve(jvari) % block_num
          end do
       end if
    end do

  end subroutine couplings_initialize_solver

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-05-02
  !> @brief   Impose the Dirichlet condition
  !> @details Impose the Dirichlet condition on implicit couplings
  !>
  !-----------------------------------------------------------------------

  subroutine couplings_impose_dirichlet(solve,unkno)

    use def_kintyp_solvers, only : soltyp
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : lcoup_implicit_d

    type(soltyp), intent(in)    :: solve
    real(rp),     intent(inout) :: unkno(*)
    integer(ip)                 :: icoup,kcoup,ndofn

    ndofn = solve % ndofn

    if( INOTMASTER .and. mcoup > 0 ) then
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup)
           if( solve % kfl_iffix /= 0 ) then
              call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,unkno,mask=solve % kfl_fixno)
           else
              call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,unkno)
           end if
        end do
     end if

   end subroutine couplings_impose_dirichlet

   !-----------------------------------------------------------------------
   !>
   !> @author  houzeaux
   !> @date    2018-05-02
   !> @brief   Check the Dirichlet condition
   !> @details Check if a variable satisfies the Dirichlet condition
   !>          on implicit couplings
   !>
   !-----------------------------------------------------------------------

   subroutine couplings_check_dirichlet(ndofn,unkno,solve)

     use def_kintyp_solvers, only : soltyp
     use def_master, only : lninv_loc
     use def_master, only : intost
     use def_domain, only : npoin
     use def_coupli, only : ncoup_implicit_d
     use def_coupli, only : lcoup_implicit_d

     integer(ip),  intent(in)             :: ndofn
     real(rp),     intent(inout)          :: unkno(*)
     type(soltyp), intent(in),   optional :: solve
     integer(ip)                          :: icoup,kcoup,ipoin,kpoin,idofn,itotn
     real(rp),     pointer                :: unkno_cpy(:)
     real(rp)                             :: eps
     logical(lg)                          :: mask

     if( INOTMASTER .and. mcoup > 0 ) then
        allocate(unkno_cpy(npoin*ndofn))
        unkno_cpy(1:ndofn*npoin) = unkno(1:ndofn*npoin)
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup)
           call COU_INTERPOLATE_NODAL_VALUES_go(icoup,ndofn,unkno)
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              do idofn = 1,ndofn
                 itotn = (ipoin-1)*ndofn+idofn
                 eps   = abs(unkno_cpy(itotn)-unkno(itotn))
                 mask  = .false.
                 if( present(solve) ) then
                    if( solve % kfl_iffix /= 0 ) then
                       if( solve % kfl_fixno(idofn,ipoin) == 1 ) mask = .true.
                    end if
                 end if
                 if( eps > 1.0e-6_rp .and. .not. mask ) then
                    print*,'WE DO NOT SATISFY THE DIRICHLET CONDITION AT NODE '//intost(lninv_loc(ipoin))//', VALUE =',eps
                 end if
              end do
           end do
        end do
        deallocate(unkno_cpy)
     end if

   end subroutine couplings_check_dirichlet

   subroutine couplings_test_transmission_conditions()
     use def_master
     use def_domain
     use def_coupli
     implicit none
     real(rp), pointer :: xx(:)
     integer(ip)       :: ii,icoup,kcoup,kpoin

     allocate(xx(max(1_ip,npoin)))
     do ii = 1,npoin
        xx(ii) = 2.0_rp*coord(2,ii)+3.0_rp
     end do

     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)
        do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
           ii = coupling_type(icoup) % wet % lpoin_wet(kpoin)
           xx(ii) = 0.0_rp
        end do
     end do

     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)
        call COU_INTERPOLATE_NODAL_VALUES_go(icoup,1_ip,xx)
        do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
           ii = coupling_type(icoup) % wet % lpoin_wet(kpoin)
           print*,lninv_loc(ii),xx(ii),abs(xx(ii)-(2.0_rp*coord(2,ii)+3.0_rp))
        end do
     end do

   end subroutine couplings_test_transmission_conditions

   subroutine couplings_select_boundaries(coupling,boundary_mask,where_type,where_number)

     type(typ_color_coupling),           intent(inout) :: coupling
     logical(lg),              pointer,  intent(inout) :: boundary_mask(:)
     integer(ip),              optional, intent(in)    :: where_type
     integer(ip),              optional, intent(in)    :: where_number
     integer(ip)                                       :: where_type_loc
     integer(ip)                                       :: where_number_loc
     integer(ip)                                       :: bcode_target
     integer(ip)                                       :: set_target
     integer(ip)                                       :: field_target
     integer(ip)                                       :: iboun
     
     if( .not. associated(boundary_mask) ) then
        call memory_alloca(memor_cou,'boundary_masj','MOD_COUPLINGS',boundary_mask,nboun)       
     end if
     
     if( present(where_type) ) then
        where_type_loc = where_type
     else
        where_type_loc = coupling % where_type 
     end if
     if( present(where_number) ) then
        where_number_loc = where_number
     else
        where_number_loc = coupling % where_number 
     end if

     if( where_type_loc == ON_SET ) then
        !
        ! Coupling on a set
        !
        set_target  = where_number_loc
        do iboun = 1,nboun
           if( lbset(iboun) == set_target ) then
              boundary_mask(iboun) = .true.
           end if
        end do
        
     else if( where_type_loc == ON_FIELD ) then
        !
        ! Coupling on a field
        !
        field_target = where_number_loc 
        do iboun = 1,nboun
           if( int(xfiel(field_target) % a(1,iboun,1),ip) /= 0 ) then
              boundary_mask(iboun) = .true.
           end if
        end do

     else if( where_type_loc == ON_CODE ) then
        !
        ! Coupling on a boundary code
        !
        bcode_target = where_number_loc 
        do iboun = 1,nboun
           if( kfl_codbo(iboun) == bcode_target ) then
              boundary_mask(iboun) = .true.
           end if
        end do

     end if

   end subroutine couplings_select_boundaries
   
   !-----------------------------------------------------------------------
   !>
   !> @author  houzeaux
   !> @date    2019-05-02
   !> @brief   Time step
   !> @details Time step strategy for coupling
   !>
   !-----------------------------------------------------------------------

   subroutine couplings_time_step(dtinv,MESSAGE)

     real(rp),              intent(inout) :: dtinv
     logical(lg), optional, intent(in)    :: MESSAGE
     logical(lg)                          :: if_message

     if( mcoup > 0 .and. kfl_timco_cou == 0 ) then
        if_message = .false.
        if( present(MESSAGE) ) if_message = MESSAGE
        if( if_message) call messages_live('COUPLI: TIME STEP SYNCHRONIZATION')
        call PAR_MAX(dtinv,'IN THE WORLD')
     end if

   end subroutine couplings_time_step
   
   !-----------------------------------------------------------------------
   !>
   !> @author  borrell
   !> @date    2020-03-10
   !> @brief   Check if all the copulings are exaustive
   !> @details Check if all the couplings are exaustive, i.e if all wet
   !>          points for all copulings receive a contribution
   !>
   !-----------------------------------------------------------------------
   function  couplings_are_exhaustive() result(istat)
      
      use def_domain,         only :  npoin
      use def_master,         only :  INOTMASTER
      use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE

      implicit none

      logical(lg)                 :: istat
      real(rp),       pointer     :: xsource(:)
      real(rp),       pointer     :: xtarget(:)
      integer(ip)                 :: icoup, ipoin
      integer(ip)                 :: color_target_aux
      character(100), PARAMETER   :: vacal = "coupling_toolbox_is_exhaustive"

      nullify(xsource,xtarget)
      
      color_target_aux = color_target
      call memory_alloca(memor_cou,'XSOURCE',vacal,xsource,max(1_ip,npoin))
      call memory_alloca(memor_cou,'XTARGET',vacal,xtarget,max(1_ip,npoin))

      istat = .true.
      do icoup = 1,mcoup
        
         color_target = coupling_type(icoup) % color_target 
         xsource(:) = 1.0_rp
         xtarget(:) = 0.0_rp
         if( I_AM_IN_COUPLING(icoup) .and. INOTMASTER ) then

            call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,xtarget,xsource)

            if (coupling_type(icoup) % kfl_par_transm == 0) then
               call PAR_INTERFACE_NODE_EXCHANGE(1_ip,xsource,'SUM','IN CURRENT TARGET COLOR')
            endif

            do ipoin = 1,coupling_type(icoup) % wet % npoin_wet
               if(xtarget(ipoin) == 0_ip ) then
                  istat = .false.
               endif
            end do
         end if

      end do
      call memory_deallo(memor_cou,'XSOURCE',vacal,xsource)
      call memory_deallo(memor_cou,'XTARGET',vacal,xtarget)
      color_target = color_target_aux

      return

   end function couplings_are_exhaustive

   !----------------------------------------------------------------------
   !>
   !> @author  Ricard Borrell
   !> @date    13/06/2017
   !> @brief   Element graph of coupled problelms 
   !> @details Generates the graph (base on nodal connectivity) of the coupled 
   !>          problem. Two elements of different codes are coupled if any of 
   !>          their nodes are connected through the transmission matrices.  
   !>          Actually this tool does not work properly if combined with mesh 
   !>          multiplication, because the "internal" graph (lelel, pelel)  is 
   !>          evaluated before carring out the multiplicaiton, and the "external"
   !>          couplings after
   !>
   !>   OUTPUT (saved in case.cou.dat_graph) 
   !>
   !>   NELEM_COU:  #elmements of the graph representing overall coupled problem
   !>   NEDGE_COU:  #edges of the graph
   !>   PELEL_COU:  Pointers to the LELEL_COU array (of size NELEM_COU+1)
   !>   PELEL_COU:  Adjacencies (edges) of each element
   !>   LELEW_COU:  Edges weights
   !>
   !----------------------------------------------------------------------
   subroutine couplings_elements_graph()

      use def_master,         only :  namda
      use def_master,         only :  IMASTER
      use def_master,         only :  current_code
      use def_master,         only :  leinv_loc
      use def_domain,         only :  npoin
      use def_domain,         only :  nelem
      use def_domain,         only :  pelpo
      use def_domain,         only :  lelpo
      use def_domain,         only :  pelel
      use def_domain,         only :  lelel
      use def_domain,         only :  nedge
      use def_coupli,         only :  memor_cou
      use mod_memory,         only :  memory_alloca
      use mod_memory,         only :  memory_deallo
      use mod_parall,         only :  color_target
      use mod_parall,         only :  color_source
      use mod_parall,         only :  PAR_MY_WORLD_RANK 
      use mod_parall,         only :  PAR_WORLD_SIZE
      use mod_communications, only :  PAR_SUM
      use mod_communications, only :  PAR_MAX
      use mod_communications, only :  PAR_GATHER
      use mod_communications, only :  PAR_GATHERV
      use mod_communications, only :  PAR_ALLGATHER
      use mod_communications, only :  PAR_SEND_RECEIVE_TO_ALL
      use mod_maths,          only :  maths_heap_sort

      implicit none

      logical(lg)                 :: ICOUMASTER
      logical(lg)                 :: auxl
      character(150)              :: fname       !output file name
      integer(ip)                 :: u=189962    !output file unit     
      integer(ip)                 :: icoup,ipoin,irank,icode
      integer(ip)                 :: ielms,inewe,ineig,itmat
      integer(ip)                 :: icont,jcont,ia,ja
      integer(ip)                 :: ncode,nnewe,nneig
      integer(ip)                 :: nelem_world, nelem_code
      integer(ip)                 :: nedge_world, nedge_code
      integer(ip)                 :: npoin_recv,npelel,nlelel
      integer(ip)                 :: melpo
      integer(ip)                 :: cuelm
      integer(ip)                 :: sumnew
      integer(ip)                 :: auxi,auxj
      integer(ip)                 :: iaux,kaux
      integer(ip)                 :: auxs(2)
      integer(ip),pointer         :: nelco(:)          ! #elements per core
      integer(ip),pointer         :: nselc(:)          ! sum of #elements in "previous" codes
      integer(ip),pointer         :: lcora(:)          ! code associated to each rank
      real(rp),pointer            :: lnels(:,:)        ! list number of elements per point in source
      real(rp),pointer            :: lpels(:,:)        ! list elements indices per point source
      type(r2p),pointer           :: lnelt(:)          ! list number of elements per point in target
      type(r2p),pointer           :: lpelt(:)          ! list elements indices per point target
      integer(ip),pointer         :: lnewe(:,:)
      integer(4),pointer          :: lnnewe(:)         ! number of new edges generated on each rank
      integer(ip),pointer         :: lnewem(:,:)       ! list of new edged stored in COUMASTER
      type(i1p),pointer           :: extco(:)          ! external couplings
      integer(4),pointer          :: lnpel(:)          ! list of pelel sizes
      integer(4),pointer          :: lnlel(:)          ! list of lelel sizes
      integer(ip),pointer         :: pelel_buf(:)      ! buffer used by COUMASTER to receive pelels
      integer(ip),pointer         :: lelel_buf(:)      ! buffer used by COUMASTER to receive lelels  
      integer(ip),pointer         :: lelel_aux(:)
      type(typ_color_coupling)    :: coupling
      character(100), PARAMETER :: vacal = "cou_elements_graph "
      !
      !   Outputs (written to disc)
      !
      integer(ip)                 :: nelem_cou 
      integer(ip)                 :: nedge_cou
      integer(ip),pointer         :: pelel_cou(:)
      integer(ip),pointer         :: lelel_cou(:)
      integer(ip),pointer         :: lelew_cou(:)


      !
      !   Initializations
      !
      ICOUMASTER = .false.
      if(PAR_MY_WORLD_RANK == 0_ip) ICOUMASTER = .true.
      nullify(nelco,nselc,lcora,lnnewe,lnpel,lnlel,lnewe,lnewem,extco,pelel_cou)
      nullify(lelel_buf,pelel_buf,lelel_aux,lpels,lelel_cou,lnels,lnelt,lelew_cou)
      nullify(lpelt)
      sumnew = 0_ip
      !
      ! Evaluate: nedge_world and nelem_world
      !
      auxs(1)=nelem
      auxs(2)=nedge
      call PAR_SUM(2_ip, auxs,"IN THE WORLD")
      nelem_world=auxs(1)
      nedge_world=auxs(2)
      !
      ! Evaluate: nedge_code and nelem_code
      !
      auxs(1)=nelem
      auxs(2)=nedge
      call PAR_SUM(2_ip, auxs,"IN MY CODE")
      nelem_code=auxs(1)
      nedge_code=auxs(2)
      !
      ! Evaluate number of codes (ncode) elements per code (nelco(:))
      ! sum of elements in previous codes (nselc(:)) and the 
      ! code to which each rank is associated (lcora(:))
      !
      ncode = current_code
      call PAR_MAX(ncode,"IN THE WORLD")

      call memory_alloca(memor_cou,'nelco',vacal,nelco,ncode)
      call memory_alloca(memor_cou,'nselc',vacal,nselc,ncode)
      call memory_alloca(memor_cou,'lcora',vacal,lcora,PAR_WORLD_SIZE)

      if(IMASTER) then
         nelco(current_code) = nelem_code
      endif
      call PAR_SUM(nelco,"IN THE WORLD")
      auxi = 0
      do icode = 1,ncode
         nselc(icode) = auxi
         auxi = auxi + nelco(icode)
      enddo
      if(auxi /= nelem_world) then
         call runend("cou_initialize_coupling: something wrong at summing up code elements")
      endif
      call PAR_ALLGATHER(current_code,lcora,1_4,"IN THE WORLD")
      !
      ! Send through the coupling commuincaiton, #elements (lnels(:)) 
      ! and elements containing each source point (lpels(:)). Those are 
      ! received in lnelt(:) and lpelt(:), respectively. Used melpo, maximum
      ! elements containing any point of the coupling
      !
      if(.not.IMASTER) then
         melpo = 0_ip 
         call memory_alloca(memor_cou,'lnels',vacal,lnels,1_ip,npoin)
         call memory_alloca(memor_cou,'lnelt',vacal,lnelt,mcoup)
         do ipoin = 1,npoin
            lnels(1_ip,ipoin) = real(pelpo(ipoin+1)-pelpo(ipoin),rp)
         enddo
         do icoup = 1,mcoup

            if( I_AM_IN_COUPLING(icoup) ) then

               coupling         = coupling_type(icoup)
               color_source     = coupling % color_source
               color_target     = coupling  % color_target

               npoin_recv       = coupling % commd % lrecv_dim
               nullify(lnelt(icoup) % a)
               call memory_alloca(memor_cou,'lnelt % a',vacal,lnelt(icoup) % a,1_ip,max(1_ip,npoin_recv))

               call PAR_SEND_RECEIVE_TO_ALL(1_ip,lnels,lnelt(icoup) % a,coupling  % commd&
                  ,'ASYNCHRONOUS',coupling % commd % lsend_perm)
               melpo = max(melpo,int(maxval(lnelt(icoup) % a(1_ip,:)),ip))

            endif
         enddo
      endif

      call PAR_MAX(melpo,"IN THE WORLD")

      if(.not.IMASTER) then
         call memory_alloca(memor_cou,'lpels',vacal,lpels,melpo,npoin)
         call memory_alloca(memor_cou,'lpelt',vacal,lpelt,mcoup)
         lpels = -1_rp
         do ipoin = 1,npoin
            icont = 1_ip
            do ielms = pelpo(ipoin),pelpo(ipoin+1)-1_ip
               if(icont <= melpo) then
                  lpels(icont,ipoin) = real(leinv_loc(lelpo(ielms)),rp) + real(nselc(lcora(PAR_MY_WORLD_RANK+1)),rp)
                  icont = icont + 1_ip
               endif
            enddo
         enddo
         do icoup = 1,mcoup

            if( I_AM_IN_COUPLING(icoup) ) then

               coupling         = coupling_type(icoup)
               color_source     = coupling % color_source
               color_target     = coupling % color_target

               npoin_recv       = coupling % commd % lrecv_dim
               nullify(lpelt(icoup) % a)
               call memory_alloca(memor_cou,'lpelt % a',vacal,lpelt(icoup) % a,melpo,max(1_ip,npoin_recv))

               call PAR_SEND_RECEIVE_TO_ALL(mcoup,lpels,lpelt(icoup) % a,coupling % commd&
                  ,'ASYNCHRONOUS',coupling % commd % lsend_perm)

            endif
         enddo
      endif
      !
      ! Count new edges (nnewe) and store them in lnewe(:,:)
      ! Note: for each new edge (i,j), the oposite edge (j,i) is also included
      !
      nnewe = 0_ip
      do icoup = 1,mcoup
         coupling         = coupling_type(icoup)
         if(associated(coupling % ltransmat_target)) then
            do ineig = 1,coupling % commd % nneig
               if(associated(coupling % ltransmat_target(ineig) % iA)) then

                  do itmat = 1, size(coupling % ltransmat_target(ineig) % jA)

                     ja = coupling % ltransmat_target(ineig) % jA(itmat)
                     ia = coupling % wet % lpoin_wet( coupling % ltransmat_target(ineig) % iA(itmat))
                     nnewe = nnewe + int(lnelt(icoup) % a(1, ja) * lnels(1,ia),ip)

                  enddo
               end if
            end do
         end if
      enddo

      call memory_alloca(memor_cou,'lnewe',vacal,lnewe,2_ip,2_ip*nnewe)
      inewe = 1_ip
      do icoup = 1,mcoup
         coupling         = coupling_type(icoup)
         if(associated(coupling % ltransmat_target)) then
            do ineig = 1,coupling % commd % nneig
               if(associated(coupling % ltransmat_target(ineig) % iA)) then

                  do itmat = 1, size(coupling % ltransmat_target(ineig) % jA)

                     ja = coupling % ltransmat_target(ineig) % jA(itmat)
                     ia = coupling % wet % lpoin_wet( coupling % ltransmat_target(ineig) % iA(itmat))

                     do icont= 1,int(lnels(1,ia),ip)
                        do jcont = 1,int(lnelt(icoup) % a(1_ip,ja),ip) 
                           lnewe(1_ip,inewe) = int(lpels(icont,ia),ip)
                           lnewe(2_ip,inewe) = int(lpelt(icoup) % a(jcont, ja),ip)
                           inewe = inewe + 1_ip
                           lnewe(1_ip,inewe) = int(lpelt(icoup) % a(jcont, ja),ip)
                           lnewe(2_ip,inewe) = int(lpels(icont,ia),ip)
                           inewe = inewe + 1_ip
                        enddo
                     enddo

                  enddo
               end if
            end do
         end if
      enddo
      if(inewe-1_ip /= nnewe*2_ip) call runend("Error in contruction of new edges of coupled graph")

      !
      ! Send all "external" couplings edges to COUMASTER
      !

      if(ICOUMASTER) call memory_alloca(memor_cou,'lnnewe',vacal,lnnewe,PAR_WORLD_SIZE)
      call PAR_GATHER(int(2*nnewe,4),lnnewe,"IN THE WORLD")

      if(ICOUMASTER)then
         nnewe = 0_ip
         nnewe = nnewe + sum(lnnewe)
         lnnewe(:) = 2_ip*lnnewe(:)
         call memory_alloca(memor_cou,'lnewem',vacal,lnewem,2_ip,nnewe)
      endif
      call PAR_GATHERV(lnewe,lnewem,lnnewe,"IN THE WORLD")

      !
      ! Store "external" couplings in extco double pointer:
      !
      !   extco[i] -> allocated in case there are couplings for element i
      !   extco[i][1] -> #edges
      !   extco[i][2] ... extco[i][exco[i][1]+1] -> edges with element i as initial point
      !

      if(ICOUMASTER) then
         !
         ! Order new edges
         !
         call maths_heap_sort(2_ip,nnewe,ivin=lnewem(1,:),ivo1=lnewem(2,:))
         !
         ! Allocate structure to estore graph of external couplings (count and allocate)  
         !
         call memory_alloca(memor_cou,'extco',vacal,extco,nelem_world)
         nneig = 0_ip
         cuelm = lnewem(1,1)
         do inewe=1,nnewe
            if(lnewem(1_ip,inewe)==cuelm) then
               nneig = nneig + 1_ip
            else
               call memory_alloca(memor_cou,'extco(cuelm) % l',vacal,extco(cuelm) % l,nneig+1)
               extco(cuelm) % l(1_ip) = 0_ip 
               cuelm = lnewem(1_ip,inewe)
               nneig = 1_ip
            endif
         enddo
         call memory_alloca(memor_cou,'extco(cuelm) % l',vacal,extco(cuelm) % l,nneig+1)
         extco(cuelm) % l(1_ip) = 0_ip 
         !
         ! Store couplings in allocated structure   
         !
         sumnew = 0_ip
         do inewe=1,nnewe
            cuelm = lnewem(1_ip,inewe)
            auxl = .false.
            do icont = 1,extco(cuelm) % l(1)
               if(extco(cuelm) % l(icont+1) == lnewem(2_ip,inewe)) then
                  auxl = .true.
               endif
            enddo
            if(.not.auxl) then
               extco(cuelm) % l(1) = extco(cuelm) % l(1) + 1_ip
               extco(cuelm) % l(extco(cuelm) % l(1)+1) = lnewem(2_ip,inewe)
               sumnew = sumnew + 1_ip
            endif
         enddo

      endif
      !
      ! Send all "pelels" to COUMASTER (this could be optimized, 
      !                                 with less comms)
      !
      if(ICOUMASTER) call memory_alloca(memor_cou,'lnpel',vacal,lnpel,PAR_WORLD_SIZE)
      call PAR_GATHER(int(size(pelel),4),lnpel,"IN THE WORLD")

      if(ICOUMASTER) then
         jcont = 1_ip
         do icont = 1,PAR_WORLD_SIZE
            if(lnpel(icont) == 1_ip) then
               lnpel(icont) = 0_ip
            else
               if(jcont /= lcora(icont)) then
                  call runend("Codes disorderd across ranks. Code not ready for this")
               endif
               jcont = jcont + 1_ip
            endif
         enddo
         npelel = 0_ip
         npelel = npelel + sum(lnpel)
         call memory_alloca(memor_cou,'pelel_buf',vacal,pelel_buf,npelel)
      endif
      call PAR_GATHERV(pelel,pelel_buf,lnpel,"IN THE WORLD")
      !
      ! Send all "lelels" to COUMASTER 
      !
      if(ICOUMASTER) then
         call memory_alloca(memor_cou,'lnlel',vacal,lnlel,PAR_WORLD_SIZE)
         lnlel = 0_ip
         auxi = 0_ip
         do icont = 1,PAR_WORLD_SIZE
            if(lnpel(icont) > 0_ip) then
               auxi = auxi + lnpel(icont)
               lnlel(icont) = pelel_buf(auxi)-1_ip
            endif
         enddo
         nlelel = 0_ip
         nlelel = nlelel + sum(lnlel)
         call memory_alloca(memor_cou,'lelel_buf',vacal,lelel_buf,nlelel)
      endif
      if(size(lelel) > 1_ip) then
         call memory_alloca(memor_cou,'lele_aux',vacal,lelel_aux,int(size(lelel),ip))
         lelel_aux(1:size(lelel)) = lelel(1:size(lelel)) + nselc(lcora(PAR_MY_WORLD_RANK+1))
      endif
      call PAR_GATHERV(lelel_aux,lelel_buf,lnlel,"IN THE WORLD")
      !
      ! Generate pelel_cou and lelel_cou
      !
      if(ICOUMASTER) then

         call memory_alloca(memor_cou,'pelel_cou',vacal,pelel_cou,npelel-ncode+1)
         call memory_alloca(memor_cou,'lelel_cou',vacal,lelel_cou,size(lelel_buf)+sumnew)
         call memory_alloca(memor_cou,'lelew_cou',vacal,lelew_cou,size(lelel_buf)+sumnew)

         pelel_cou(1) = 1_ip
         iaux = 1_ip
         kaux = 2_ip
         auxi = 0_ip
         auxj = 1_ip

         do irank = 1_ip,PAR_WORLD_SIZE
            do icont = auxi+2_ip,auxi+lnpel(irank)
               pelel_cou(kaux) = pelel_cou(kaux-1)+pelel_buf(icont)-pelel_buf(icont-1)
               do jcont = pelel_cou(kaux-1),pelel_cou(kaux)-1_ip
                  lelel_cou(auxj) = lelel_buf(iaux)
                  lelew_cou(auxj) = 1_ip
                  iaux = iaux + 1_ip
                  auxj = auxj + 1_ip
               enddo
               if(associated(extco(kaux-1) % l)) then
                  pelel_cou(kaux) = pelel_cou(kaux) + extco(kaux-1_ip) % l(1_ip)
                  do jcont = 2_ip,extco(kaux-1) % l(1) + 1_ip
                     lelel_cou(auxj) = extco(kaux-1) % l(jcont)
                     lelew_cou(auxj) = 2_ip
                     auxj = auxj + 1_ip
                  enddo
               endif
               kaux = kaux + 1_ip              
            enddo
            auxi = auxi + lnpel(irank)
         enddo

         nelem_cou = nelem_world
         nedge_cou = size(lelel_cou)
         !
         ! Writte to disc
         !
         fname = adjustl(trim(namda))//".cou.dat_graph"
         open(u,file=fname,action="write")
         write(u,"(I6,I6,A6)") nelem_cou,nedge_cou,"001"
         do icont = 1,nelem_cou
            do jcont = pelel_cou(icont),pelel_cou(icont+1)-1_ip
               write(u,"(I6,I6)",advance='no') lelel_cou(jcont),lelew_cou(jcont)
            enddo
            write(u,*)
         enddo
         close(u)


      endif
      !
      ! Deallo memory
      !   
      if(associated(nelco)) call memory_deallo(memor_cou,'nelco',vacal,nelco)
      if(associated(nselc)) call memory_deallo(memor_cou,'nselc',vacal,nselc)
      if(associated(lcora)) call memory_deallo(memor_cou,'lcora',vacal,lcora)
      if(associated(lnels)) call memory_deallo(memor_cou,'lnels',vacal,lnels)
      if(associated(lpels)) call memory_deallo(memor_cou,'lpels',vacal,lpels)
      if(associated(lnelt)) call memory_deallo(memor_cou,'lnelt',vacal,lnelt)
      if(associated(lpelt)) call memory_deallo(memor_cou,'lpelt',vacal,lpelt)
      if(associated(lnewe)) call memory_deallo(memor_cou,'lnewe',vacal,lnewe)
      if(associated(lnnewe)) call memory_deallo(memor_cou,'lnnewe',vacal,lnnewe)
      if(associated(lnewem)) call memory_deallo(memor_cou,'lnewem',vacal,lnewem)
      if(associated(extco)) call memory_deallo(memor_cou,'extco',vacal,extco)
      if(associated(lnpel)) call memory_deallo(memor_cou,'lnpel',vacal,lnpel)
      if(associated(lnlel)) call memory_deallo(memor_cou,'lnlel',vacal,lnlel)
      if(associated(pelel_buf)) call memory_deallo(memor_cou,'pelel_buf',vacal,pelel_buf)
      if(associated(lelel_buf)) call memory_deallo(memor_cou,'lelel_buf',vacal,lelel_buf)
      if(associated(lelel_aux)) call memory_deallo(memor_cou,'lelel_aux',vacal,lelel_aux)
      if(associated(pelel_cou)) call memory_deallo(memor_cou,'pelel_cou',vacal,pelel_cou)
      if(associated(lelel_cou)) call memory_deallo(memor_cou,'lelel_cou',vacal,lelel_cou)
      if(associated(lelew_cou)) call memory_deallo(memor_cou,'lelew_cou',vacal,lelew_cou)

   end subroutine couplings_elements_graph

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2020-06-30
   !> @brief   Compute wet nodes
   !> @details Define wet nodes for a mirror coupling
   !> 
   !-----------------------------------------------------------------------

   subroutine cou_wet_points_from_mirror(coupling)

     type(typ_color_coupling), intent(inout) :: coupling
     integer(ip)                             :: jcoup,kpoin,ipoin
     integer(ip),              pointer       :: kweight_wet(:)
     integer(ip),              pointer       :: list_wet_nodes(:)

     nullify(kweight_wet)
     nullify(list_wet_nodes)
     if( INOTMASTER ) then
     !
     ! Compute weight points
     ! 
     call memory_alloca(memor_cou,'LIST_WET_NODES',vacal,list_wet_nodes,npoin)
     jcoup                              = coupling % mirror_coupling
     coupling % wet % number_wet_points = coupling_type(jcoup) % geome % npoin_source
     do kpoin = 1,coupling_type(jcoup) % geome % npoin_source
        ipoin                 = coupling_type(jcoup) % geome % lpoin_source(kpoin)       
        list_wet_nodes(ipoin) = 1
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(list_wet_nodes,'MAX','IN CURRENT TARGET COLOR')
     coupling % wet % number_wet_points = 0     
     do ipoin = 1,npoin
        coupling % wet % number_wet_points = coupling % wet % number_wet_points + list_wet_nodes(ipoin)
     end do
     coupling % wet % npoin_wet = coupling % wet % number_wet_points
     
     call memory_deallo(memor_cou,'COUPLING % WET % COORD_WET'     ,vacal,coupling % wet % coord_wet)
     call memory_deallo(memor_cou,'COUPLING % WET % LPOIN_WET'     ,vacal,coupling % wet % lpoin_wet)
     call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET'    ,vacal,coupling % wet % weight_wet)
     call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET_IMP',vacal,coupling % wet % weight_wet_imp)
     call memory_deallo(memor_cou,'KWEIGHT_WET'                    ,vacal,kweight_wet)
     
     call memory_alloca(memor_cou,'COUPLING % WET % COORD_WET'     ,vacal,coupling % wet % coord_wet,ndime,coupling % wet % number_wet_points)
     call memory_alloca(memor_cou,'COUPLING % WET % LPOIN_WET'     ,vacal,coupling % wet % lpoin_wet,      coupling % wet % number_wet_points)
     call memory_alloca(memor_cou,'COUPLING % WET % WEIGHT_WET'    ,vacal,coupling % wet % weight_wet,     coupling % wet % number_wet_points)
     call memory_alloca(memor_cou,'COUPLING % WET % WEIGHT_WET_IMP',vacal,coupling % wet % weight_wet_imp, coupling % wet % number_wet_points)
     call memory_alloca(memor_cou,'KWEIGHT_WET'                    ,vacal,kweight_wet,npoin)

     kpoin = 0
     do ipoin = 1,npoin
        if( list_wet_nodes(ipoin) == 1 ) then
           kpoin                                     = kpoin + 1
           coupling % wet % coord_wet(1:ndime,kpoin) = coord(1:ndime,ipoin)      
           coupling % wet % lpoin_wet(kpoin)         = ipoin
           kweight_wet(ipoin)                        = 1
        end if
     end do
     !
     ! Define weights of wet nodes             
     !
     call PAR_INTERFACE_NODE_EXCHANGE(kweight_wet,'SUM','IN CURRENT TARGET COLOR')
     if( coupling % kind == BETWEEN_SUBDOMAINS ) then
        !
        ! Subdomain coupling: weight is one, because coupling is carried out AFTER exchange
        !
        do kpoin = 1,coupling % wet % npoin_wet
           coupling % wet % weight_wet(kpoin) = 1.0_rp
        end do
     else
        !
        ! Weight depends on the number of neighbors, because coupling is carried out BEFORE exchange
        !   
        do kpoin = 1,coupling % wet % npoin_wet
           ipoin = coupling % wet % lpoin_wet(kpoin)
           coupling % wet % weight_wet(kpoin) = 1.0_rp / real(kweight_wet(ipoin),rp)
        end do
     end if
     !
     ! Weight used when transmision matrices are not parallelized
     !
     do kpoin = 1,coupling % wet % npoin_wet
        ipoin = coupling % wet % lpoin_wet(kpoin)
        coupling % wet % weight_wet_imp(kpoin) = 1.0_rp / real(kweight_wet(ipoin),rp)
     end do
     !
     ! Compute fringe wet nodes
     !
     if ( coupling % wet % kfl_get_fringe ) then
        ! Allocate fringe wet nodes
        call memory_alloca(memor_cou,'COUPLING % WET % KFL_FRINGE_WETNODES','mod_couplings', &
           coupling % wet % kfl_fringe_wetnodes, npoin)
        ! Get fringe wet nodes
        call cou_get_fringe_wetnodes(list_wet_nodes, coupling)
     end if
     !
     ! Deallocate
     !
     call memory_deallo(memor_cou,'KWEIGHT_WET'   ,vacal,kweight_wet)
     call memory_deallo(memor_cou,'LIST_WET_NODES',vacal,list_wet_nodes)

     end if

   end subroutine cou_wet_points_from_mirror

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Get fringe wet nodes of background
   !> @details Get fringe wet nodes of background for multi-code
   !>          immersed-type problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_get_fringe_wetnodes(list_wet_nodes, coupling)
     ! ----------------------------------------------- !
     ! Identify fringe nodes as:
     !
     ! kfl_fringe = 0: dry node
     !              1: fringe wet node
     !             -1: interior wet node
     !
     ! ----------------------------------------------- !

     use def_domain, only : c_dom, r_dom

     implicit none

     integer(ip),              intent(in)    :: list_wet_nodes(*)
     type(typ_color_coupling), intent(inout) :: coupling
     integer(ip)                             :: ipoin, jpoin, kpoin, izdom

     ! Set to 0 as default for all nodes
     coupling % wet % kfl_fringe_wetnodes(1:npoin) = 0_ip

     ! Loop over nodes
     do kpoin = 1,coupling % wet % npoin_wet
        ipoin = coupling % wet % lpoin_wet(kpoin)

        ! Set to -1 as default for wet nodes
        coupling % wet % kfl_fringe_wetnodes(ipoin) = -1_ip

        ! Loop over neighbors
        do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
           jpoin = c_dom(izdom)

           ! If neighbor is not wet then ipoin is a fringe node
           if (list_wet_nodes(jpoin) == 0_ip) then
              coupling % wet % kfl_fringe_wetnodes(ipoin) = 1_ip
           end if
           !
        end do
        !
     end do
     !
     ! Exchange fringe nodes
     !
     call PAR_INTERFACE_NODE_EXCHANGE(coupling % wet % kfl_fringe_wetnodes,'MAX','IN CURRENT TARGET COLOR')

   end subroutine cou_get_fringe_wetnodes

end module mod_couplings
