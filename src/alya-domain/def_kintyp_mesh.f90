!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_mesh_type.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Mesh type
!> @details Mesh type definitions
!>
!>          Basic type
!>          ----------
!>          This type is intended to contain all element
!>          types. Thus boundaries are not represented specifically.
!>          This type is similar to gmsh structure.
!>
!>          |-> Basic procedures are defined here. More complex ones
!>              are in the associated module: mod_mesh_type_basic.f90
!>          !-> Postprocess is in mod_postpr_mesh.f90
!>
!>          Extended type
!>          -------------
!>          Lot more info about the mesh. For example, boundaries are
!>          explicitly declared.
!>          This is the native mesh type of Alya.
!>
!>   
!-----------------------------------------------------------------------

module def_kintyp_mesh

  use def_kintyp_basic,         only : ip,rp,i1p,i1pp,r3p,lg
  use def_kintyp_mesh_basic,    only : mesh_type_basic
  use def_kintyp_mesh_basic,    only : mesh_name_default
  use def_kintyp_domain,        only : mfiel
  use mod_memory_basic,         only : memory_alloca
  use mod_memory_basic,         only : memory_deallo 
  use mod_memory_basic,         only : memory_copy
  use mod_elmgeo,               only : element_type
  use mod_elmgeo,               only : elmgeo_element_distance
  use def_master,               only : IPARALL
  use mod_optional_argument,    only : optional_argument
  implicit none
  private

  !----------------------------------------------------------------------
  !
  ! Extended mesh type
  !
  !----------------------------------------------------------------------

  type, extends(mesh_type_basic) :: mesh_type
     !
     ! GEOMETRY 
     !
     integer(ip)               :: ntens                     ! Number Hessian components
     integer(ip)               :: nboun                     ! Number of boundaries
     integer(ip)               :: nfiel                     ! Number of fields
     integer(ip)               :: kfl_field(7,mfiel)        ! Fields dimensions
     integer(ip)               :: nedge                     ! Number of edges
     integer(ip)               :: mnodb                     ! Max number of node per boundary
     integer(ip)               :: mgaus                     ! Max number of Gauss points
     integer(ip)               :: medge                     ! Max number of edge per element
     integer(ip)               :: nbopo                     ! Number of boundary nodes
     integer(ip)               :: npoi1                     ! Number of interior node
     integer(ip)               :: npoi2                     ! First own boundary node
     integer(ip)               :: npoi3                     ! Last own boundary node
     integer(ip)               :: nedg1                     ! Number of interior node
     integer(ip)               :: nedg2                     ! First own boundary node
     integer(ip)               :: nedg3                     ! Last own boundary node
     integer(ip)               :: npoin_2                   ! Number of nodes including halo nodes
     integer(ip)               :: nelem_2                   ! Number of elements including halo nodes
     integer(ip)               :: nboun_2                   ! Number of boundaries including halo nodes
     integer(ip)               :: npoin_own                 ! Own nodes = npoi3
     integer(ip)               :: npoin_halo                ! Number of nodes up to halo nodes
     integer(ip)               :: kfl_ngrou                 ! How groups are constructed
     integer(ip),      pointer :: npoin_par(:)
     integer(ip),      pointer :: nelem_par(:)
     integer(ip),      pointer :: nboun_par(:)
     integer(ip)               :: npoin_origi
     integer(ip)               :: npoin_total
     integer(ip)               :: nelem_total
     integer(ip)               :: nboun_total
     integer(ip),      pointer :: lelch(:)                  ! NELEM
     integer(ip),      pointer :: lnnod(:)                  ! NELEM
     integer(ip),      pointer :: lesub(:)                  ! NELEM
     integer(ip),      pointer :: lmate(:)                  ! NELEM
     integer(ip),      pointer :: lgaus(:)                  ! NELEM
     integer(ip),      pointer :: lbinv_loc(:)              ! NBOUN
     integer(ip),      pointer :: lnodb(:,:)                ! NBOUN
     integer(ip),      pointer :: lboel(:,:)                ! NBOUN
     integer(ip),      pointer :: lelbo(:)                  ! NBOUN
     integer(ip),      pointer :: ltypb(:)                  ! NBOUN
     integer(ip),      pointer :: lboch(:)                  ! NBOUN
     integer(ip),      pointer :: lnnob(:)                  ! NBOUN
     integer(ip),      pointer :: lnoch(:)                  ! NPOIN: List of node characteristics
     integer(ip),      pointer :: lmast(:)                  ! NPOIN: List of master nodes
     integer(ip),      pointer :: lpoty(:)                  ! NPOIN: list of boudnary nodes
     type(r3p)                 :: xfiel(mfiel)              ! FIELDS
     ! 
     ! Edge structure
     !
     integer(ip),      pointer :: edge_to_node(:,:)         ! NEDGE
     integer(ip),      pointer :: ledgs(:,:)                ! NELEM: element edge connectivity  (see LNODS)
     integer(ip),      pointer :: ledgb(:,:)                ! NBOUN: boundary edge connectivity (see LNODB)
     integer(ip),      pointer :: lnned(:)                  ! NELEM: Number of edge by element  (see LNNOD)
     integer(ip),      pointer :: lnneb(:)                  ! NBOUN: Number of edge by element  (see LNNOB)
     !
     ! SETS
     !
     integer(ip),      pointer :: leset(:)                  ! NELEM
     integer(ip),      pointer :: lbset(:)                  ! NBOUN
     integer(ip),      pointer :: lnset(:)                  ! NPOIN
     !
     ! BOUNDARY CONDITIONS
     !
     integer(ip),      pointer :: kfl_codno(:,:)            ! NPOIN
     integer(ip),      pointer :: kfl_codbo(:)              ! NBOUN
     !
     ! GROUPS
     !
     integer(ip),      pointer :: lgrou_dom(:)              ! NPOIN
     !
     ! INTERPOLATION
     !
     type(i1pp),       pointer :: linno(:)                  ! NPOIN
     !
     ! Graphs
     !
     integer(ip)               :: nzdom                     ! Size of node graph
     integer(ip)               :: nzsym                     ! Size of node symmetric graph
     integer(ip)               :: nzedg                     ! Size of edge graph
     integer(ip)               :: nzelm_2                   ! Size of element graph (including halos)
     integer(ip)               :: nzdom_own                 ! Size of node graph of own nodes
     integer(ip)               :: nzdom_ell                 ! Size of ell graph
     integer(ip),      pointer :: c_dom(:)                  ! Node CSR graph
     integer(ip),      pointer :: r_dom(:)                  ! Node CSR graph
     integer(ip),      pointer :: coo_rows(:)               ! Node COO graph
     integer(ip),      pointer :: coo_cols(:)               ! Node COO graph
     integer(ip),      pointer :: ell_cols(:,:)             ! Node ELL graph
     integer(ip),      pointer :: c_sym(:)                  ! Symmetric node graph
     integer(ip),      pointer :: r_sym(:)                  ! Symmetric node graph
     integer(ip),      pointer :: c_edg(:)                  ! Edge graph
     integer(ip),      pointer :: r_edg(:)                  ! Edge graph
     integer(ip),      pointer :: c_elm_2(:)                ! Element graph including halos
     integer(ip),      pointer :: r_elm_2(:)                ! Element graph including halos
     integer(ip),      pointer :: c_dom_2(:)                ! Node CSR graph including halos
     integer(ip),      pointer :: r_dom_2(:)                ! Node CSR graph including halos
     integer(ip),      pointer :: r_dom_own(:)              ! Full row graph
     integer(ip),      pointer :: r_dom_end(:)              ! Full row graph
     integer(ip),      pointer :: r_dom_ini(:)              ! Full row graph
     integer(ip),      pointer :: c_dom_own(:)              ! Full row graph
     !
     ! Geometrical arrays
     !
     real(rp),         pointer :: exnor(:,:,:)              ! NDIME,NDIME,NBOPO: Local basis
     real(rp),         pointer :: vmass(:)                  ! NPOIN: Lumped mass matrix
     real(rp),         pointer :: vmasc(:)                  ! NPOIN: Mass matrix with closed rule
     
   contains
     
     procedure,        pass    :: init                      ! Initialize mesh
     procedure,        pass    :: deallo                    ! Deallocate
     procedure,        pass    :: boundary_bb               ! Boundary bounding box
     procedure,        pass    :: associate_boundary        ! Allocate and associated boundary

  end type mesh_type
  
  real(rp), parameter :: epsil = epsilon(1.0_rp)

  public :: mesh_type_basic
  public :: mesh_type

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Deallocate mesh type
  !> @details Deallocate mesh type
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(mesh,MEMORY_COUNTER,MESH_NAME)

    class(mesh_type)                           :: mesh
    integer(8),       optional,  intent(inout) :: MEMORY_COUNTER(2)
    character(len=*), optional,  intent(in)    :: MESH_NAME
    character(len=:), allocatable              :: my_mesh_name
    integer(ip)                                :: ifiel
    integer(8)                                 :: memor_loc(2)

    my_mesh_name = optional_argument(mesh % name,MESH_NAME)
    memor_loc    = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % NPOIN_PAR','deallo',mesh % npoin_par)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % NELEM_PAR','deallo',mesh % nelem_par)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % NBOUN_PAR','deallo',mesh % nboun_par)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LEINV_LOC','deallo',mesh % leinv_loc)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNODS','deallo',mesh % lnods)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LTYPE','deallo',mesh % ltype)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LELCH','deallo',mesh % lelch)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNNOD','deallo',mesh % lnnod)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LESUB','deallo',mesh % lesub)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LMATE','deallo',mesh % lmate)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LGAUS','deallo',mesh % lgaus)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERME','deallo',mesh % perme)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LBINV_LOC','deallo',mesh % lbinv_loc)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNODB','deallo',mesh % lnodb)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LBOEL','deallo',mesh % lboel)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LELBO','deallo',mesh % lelbo)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LTYPB','deallo',mesh % ltypb)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LBOCH','deallo',mesh % lboch)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNNOB','deallo',mesh % lnnob)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNINV_LOC','deallo',mesh % lninv_loc)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COORD','deallo',mesh % coord)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNOCH','deallo',mesh % lnoch)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LMAST','deallo',mesh % lmast)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LPOTY','deallo',mesh % lpoty)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERMN','deallo',mesh % permn)

    do ifiel = 1,mesh % nfiel
       call memory_deallo(memor_loc,trim(my_mesh_name)//' % XFIEL % A','deallo',mesh % xfiel(ifiel) % a)
    end do

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % edge_to_node','deallo',mesh % edge_to_node)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % ledgs','deallo',mesh % ledgs)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % ledgb','deallo',mesh % ledgb)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % lnned','deallo',mesh % lnned)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % lnneb','deallo',mesh % lnneb)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % leset','deallo',mesh % leset)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % lbset','deallo',mesh % lbset)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % lnset','deallo',mesh % lnset)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % kfl_codno','deallo',mesh % kfl_codno)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % kfl_codbo','deallo',mesh % kfl_codbo)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LGROU_DOM','deallo',mesh % lgrou_dom)

    !call memory_deallo(memor_loc,trim(my_mesh_name)//' % LINNO','deallo',mesh % linno)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_dom   ' ,'deallo',mesh % c_dom)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_dom   ' ,'deallo',mesh % r_dom)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % coo_rows' ,'deallo',mesh % coo_rows)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % coo_cols' ,'deallo',mesh % coo_cols)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % ell_cols' ,'deallo',mesh % ell_cols)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_sym   ' ,'deallo',mesh % c_sym)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_sym   ' ,'deallo',mesh % r_sym)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_edg   ' ,'deallo',mesh % c_edg)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_edg   ' ,'deallo',mesh % r_edg)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_elm_2 ' ,'deallo',mesh % c_elm_2)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_elm_2 ' ,'deallo',mesh % r_elm_2)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_dom_2 ' ,'deallo',mesh % c_dom_2)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_dom_2 ' ,'deallo',mesh % r_dom_2)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % c_dom_own','deallo',mesh % c_dom_own)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_dom_own','deallo',mesh % r_dom_own)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_dom_end','deallo',mesh % r_dom_end)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % r_dom_ini','deallo',mesh % r_dom_ini)

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % EXNOR'    ,'deallo',mesh % exnor)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % VMASS'    ,'deallo',mesh % vmass)
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % VMASC'    ,'deallo',mesh % vmasc)

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Initialization
  !> @details Initialization. Shoudl be called prior to any operations
  !>          on a mesh.
  !> 
  !-----------------------------------------------------------------------

  subroutine init(mesh,wname)

    class(mesh_type),                intent(inout) :: mesh
    character(len=*),      optional, intent(in)    :: wname
    integer(ip)                                    :: ifiel

    if( present(wname) ) then
       mesh % name = trim(wname)
    else
       mesh % name = mesh_name_default
    end if

    mesh % id         = 0
    mesh % ndime      = 0
    mesh % ntens      = 0
    mesh % npoin      = 0
    mesh % nelem      = 0
    mesh % nboun      = 0
    mesh % ntags      = 0
    mesh % nfiel      = 0
    mesh % kfl_field  = 0
    mesh % nedge      = 0
    mesh % mnode      = 0
    mesh % mnodb      = 0
    mesh % mgaus      = 0
    mesh % medge      = 0
    mesh % nbopo      = 0
    mesh % npoi1      = 0
    mesh % npoi2      = 0
    mesh % npoi3      = 0
    mesh % npoin_own  = 0
    mesh % npoin_halo = 0
    mesh % npoin_2    = 0
    mesh % nelem_2    = 0
    mesh % nboun_2    = 0
    mesh % nzdom      = 0
    mesh % nzsym      = 0
    mesh % nzedg      = 0
    mesh % nzelm_2    = 0
    mesh % nzdom_own  = 0
    mesh % nzdom_ell  = 0
    mesh % kfl_ngrou  = 0

    nullify(mesh % npoin_par)
    nullify(mesh % nelem_par)
    nullify(mesh % nboun_par)

    nullify(mesh % leinv_loc)
    nullify(mesh % lnods)
    nullify(mesh % ltype)
    nullify(mesh % lelch)
    nullify(mesh % lnnod)
    nullify(mesh % lesub)
    nullify(mesh % lmate)
    nullify(mesh % lgaus)
    nullify(mesh % perme)

    nullify(mesh % lbinv_loc)
    nullify(mesh % lnodb)
    nullify(mesh % lboel)
    nullify(mesh % lelbo)
    nullify(mesh % ltypb)
    nullify(mesh % lboch)
    nullify(mesh % lnnob)

    nullify(mesh % lninv_loc)
    nullify(mesh % coord)
    nullify(mesh % lnoch)
    nullify(mesh % lmast)
    nullify(mesh % lpoty)
    nullify(mesh % permn)
    nullify(mesh % tags )

    nullify(mesh % boundary)
    nullify(mesh % parent)
    
    do ifiel = 1,mesh % nfiel
       nullify(mesh % xfiel(ifiel) % a)
    end do

    nullify(mesh % edge_to_node)
    nullify(mesh % ledgs)
    nullify(mesh % ledgb)
    nullify(mesh % lnned)
    nullify(mesh % lnneb)

    nullify(mesh % leset)
    nullify(mesh % lbset)
    nullify(mesh % lnset)

    nullify(mesh % kfl_codno)
    nullify(mesh % kfl_codbo)

    nullify(mesh % lgrou_dom)

    nullify(mesh % linno)

    nullify(mesh % c_dom)
    nullify(mesh % r_dom)
    nullify(mesh % coo_rows)
    nullify(mesh % coo_cols)
    nullify(mesh % ell_cols)
    nullify(mesh % c_sym)
    nullify(mesh % r_sym)
    nullify(mesh % c_edg)
    nullify(mesh % r_edg)
    nullify(mesh % c_elm_2)
    nullify(mesh % r_elm_2)
    nullify(mesh % c_dom_2)
    nullify(mesh % r_dom_2)
    nullify(mesh % c_dom_own)
    nullify(mesh % r_dom_own)
    nullify(mesh % r_dom_end)
    nullify(mesh % r_dom_ini)

    nullify(mesh % exnor)
    nullify(mesh % vmass)
    nullify(mesh % vmasc)

    if( IPARALL ) call mesh % comm % init(COMM_NAME='MESH % COMM')

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Compute boundary bounding boxes
  !> @details Compute boundary bounding boxes. Check for possible
  !>          zero bounding boxes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine boundary_bb(mesh,bobox,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(mesh_type),                intent(in)    :: mesh
    real(rp),               pointer, intent(inout) :: bobox(:,:,:)   
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),  optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                    :: pelty,ielem,inode,ipoin
    real(rp)                                       :: comin(3),comax(3),rdime
    real(rp)                                       :: element_distance
    real(rp)                                       :: elcod(mesh % ndime,mesh % mnode)
    integer(8)                                     :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       call memory_deallo(memor_loc,'BOBOX','element_bb',bobox)
    else
       rdime     = 1.0_rp / real(mesh % ndime,rp)

       if( .not. associated(bobox) ) & 
            call memory_alloca(memor_loc,'BOBOX','boundary_bb',bobox,2_ip,mesh % ndime,mesh % nboun)
       
       do ielem = 1,mesh % nboun
          pelty =  mesh % ltypb(ielem)
          comin =  huge(1.0_rp)*0.1_rp
          comax = -huge(1.0_rp)*0.1_rp
          do inode = 1,element_type(pelty) % number_nodes
             ipoin                       = mesh % lnodb(inode,ielem)
             elcod(1:mesh % ndime,inode) = mesh % coord(1:mesh % ndime,ipoin)
             comin(1:mesh % ndime)       = min(comin(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
             comax(1:mesh % ndime)       = max(comax(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
          end do
          call elmgeo_element_distance(mesh % ndime,mesh % ltypb(ielem),elcod,element_distance)
          bobox(1,1:mesh % ndime,ielem) = comin(1:mesh % ndime) - 1.0e-16_rp * element_distance - epsil
          bobox(2,1:mesh % ndime,ielem) = comax(1:mesh % ndime) + 1.0e-16_rp * element_distance + epsil
       end do
       
    end if

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine boundary_bb

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Boundary mesh
  !> @details Associate boundary mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine associate_boundary(mesh,ONLY_DEALLOCATE)
    
    class(mesh_type),           intent(inout) :: mesh
    logical(lg),      optional, intent(in)    :: ONLY_DEALLOCATE

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       if( associated(mesh % boundary) ) deallocate(mesh % boundary)
    else
       if( associated(mesh % boundary) ) then
          write(*,'(a)') 'BOUNDARY MESH ALREADY ASSOCIATED'
          stop -1
       else
          allocate(mesh % boundary)
          call mesh % boundary % init()
          mesh % boundary % ndime =  mesh % ndime
          mesh % boundary % mnode =  mesh % mnodb
          mesh % boundary % nelem =  mesh % nboun
          mesh % boundary % ltype => mesh % ltypb
          mesh % boundary % lnods => mesh % lnodb
          mesh % boundary % coord => mesh % coord
          mesh % boundary % lboel => mesh % lboel
          mesh % boundary % lelbo => mesh % lelbo
       end if
    end if
    
  end subroutine associate_boundary
  
end module def_kintyp_mesh
!> @}
