!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Mesh_Type_Toolbox
!> @{
!> @name    ToolBox for mesh type
!> @file    mod_mesh_type.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for mesh type
!> @details ToolBox for mesh type
!
!-----------------------------------------------------------------------

module mod_mesh_type

  use def_kintyp_basic,          only : ip,rp,lg
  use def_kintyp_mesh,           only : mesh_type_basic
  use def_kermod,                only : ndivi
  use def_kermod,                only : kfl_posdi

  use def_domain,                only : npoin_origi
  use def_domain,                only : npoin_total
  use def_domain,                only : nelem_total
  use def_domain,                only : mesh_type
  use def_domain,                only : nboun_total
  use def_domain,                only : ndime,ntens,npoin,nfiel
  use def_domain,                only : mfiel,kfl_field
  use def_domain,                only : nelem,nboun,mnode,mnodb
  use def_domain,                only : mgaus,nbopo,npoin_own
  use def_domain,                only : npoin_2,nelem_2
  use def_domain,                only : nboun_2,mnodb
  use def_domain,                only : nelem,nboun
  use def_domain,                only : nedge,medge
  use def_domain,                only : npoin_halo,npoin
  use def_domain,                only : kfl_ngrou

  use def_master,                only : intost
  use def_master,                only : npoi1,npoi2,npoi3
  use def_master,                only : nedg1,nedg2,nedg3
  use def_master,                only : npoin_par,nelem_par
  use def_master,                only : nboun_par,npart
  use def_master,                only : INOTSLAVE
  use def_master,                only : IPARALL
  use def_master,                only : INOTMASTER
  use def_master,                only : IMASTER,ISEQUEN
  use def_master,                only : intost
  use def_master,                only : title
  use def_master,                only : leinv_loc
  use def_master,                only : lbinv_loc
  use def_master,                only : lninv_loc
  use def_master,                only : THIS_NODE_IS_MINE

  use def_domain,                only : memor_dom
  use def_domain,                only : meshe
  use def_domain,                only : lnods
  use def_domain,                only : ltype
  use def_domain,                only : lnnod
  use def_domain,                only : lelch
  use def_domain,                only : lesub
  use def_domain,                only : lmate
  use def_domain,                only : lgaus
  use def_domain,                only : lnodb
  use def_domain,                only : ltypb
  use def_domain,                only : lboch
  use def_domain,                only : lboel
  use def_domain,                only : lelbo
  use def_domain,                only : lnnob
  use def_domain,                only : coord
  use def_domain,                only : lnoch
  use def_domain,                only : lmast
  use def_domain,                only : lpoty
  use def_domain,                only : xfiel
  use def_domain,                only : nzdom
  use def_domain,                only : r_dom
  use def_domain,                only : c_dom
  use def_domain,                only : exnor
  use def_domain,                only : vmass
  use def_domain,                only : vmasc

  use def_domain,                only : edge_to_node
  use def_domain,                only : ledgs
  use def_domain,                only : ledgb
  use def_domain,                only : lnned
  use def_domain,                only : lnneb

  use mod_memory,                only : memory_copy
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo

  use mod_parall,                only : PAR_COMM_MY_CODE
  use mod_parall,                only : PAR_MY_CODE_RANK
  use mod_parall,                only : PAR_CODE_SIZE
  use mod_parall,                only : commd
  
  use mod_communications_global, only : PAR_SUM
  use mod_elmgeo,                only : element_type
  use mod_optional_argument,     only : optional_argument

  implicit none 

  private

  integer(ip), target :: int_min_1(1),int_min_2(1,1)
  real(rp),    target :: rea_min_2(1,1)
  character(100)      :: vacal='mod_mesh_type'

  public :: mesh_type_allocate_initialize
  public :: mesh_type_save_original_mesh
  public :: mesh_type_update_last_mesh
  public :: mesh_type_basic_to_complete
  public :: mesh_type_allocate_minimum
  public :: mesh_type_deallocate
  public :: mesh_type_initialize
  public :: mesh_type_copy
  public :: mesh_coherent_boundary_mesh
  public :: mesh_number_boundary_nodes
  
  interface mesh_type_initialize
     module procedure mesh_type_initialize_s,&
          &           mesh_type_initialize_1
  end interface mesh_type_initialize

  interface mesh_type_deallocate
     module procedure mesh_type_deallocate_all,&
          &           mesh_type_deallocate_s
  end interface mesh_type_deallocate

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-18
  !> @brief   Boundary nodes
  !> @details Compute number of boundary nodes, local and global
  !> 
  !-----------------------------------------------------------------------

  function mesh_number_boundary_nodes(mesh) result(nn)

    type(mesh_type),  intent(in) :: mesh
    integer(ip)                  :: nn(2)
    logical(lg),     pointer     :: lmask(:)
    integer(ip)                  :: inodb,iboun,ipoin,pblty
    
    nullify(lmask)
    call memory_alloca(memor_dom,'LMASK','mod_mesh_type',lmask,mesh % npoin)
    do iboun = 1,nboun
       pblty = abs(mesh % ltypb(iboun))
       do inodb = 1,element_type(pblty) % number_nodes
          ipoin = mesh % lnodb(inodb,iboun)
          if(THIS_NODE_IS_MINE(ipoin)) lmask(ipoin) = .true.
       end do
    end do
    if( mesh % npoin > 0 ) then
       nn(1) = count(lmask)
    else
       nn(1) = 0
    end if
    nn(2) = nn(1)
    call PAR_SUM(nn(2))
    call memory_deallo(memor_dom,'LMASK','mod_mesh_type',lmask)
    
  end function mesh_number_boundary_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-18
  !> @brief   Create a coherent boundary mesh
  !> @details Reorder nodes and connectivity to have stand alone coherent
  !>          boundary mesh of basic type. Halo boundaries can also
  !>          be included
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_coherent_boundary_mesh(mesh_vol,mesh_bou,INCLUDE_HALOS) 

    type(mesh_type),                 intent(in)    :: mesh_vol
    type(mesh_type_basic),           intent(inout) :: mesh_bou
    logical(lg),           optional, intent(in)    :: INCLUDE_HALOS 
    integer(ip),           pointer                 :: lperm(:)
    integer(ip)                                    :: kpoin,ipoin
    integer(ip)                                    :: inodb,iboun
    integer(ip)                                    :: pblty,nboun_loc
    integer(ip)                                    :: npoin_loc
    logical(lg)                                    :: if_halos

    if_halos = optional_argument(.false.,INCLUDE_HALOS)

    if( if_halos ) then
       npoin_loc = mesh_vol % npoin_2
       nboun_loc = mesh_vol % nboun_2
    else
       npoin_loc = mesh_vol % npoin
       nboun_loc = mesh_vol % nboun
    end if
    
    nullify(lperm)
    call mesh_bou % init(trim(mesh_vol % name)//'_BOUNDARY')
    call memory_alloca(memor_dom,'LPERM','mod_mesh_type',lperm,npoin_loc)

    kpoin = 0
    do iboun = 1,nboun_loc
       pblty = abs(mesh_vol % ltypb(iboun))
       do inodb = 1,element_type(pblty) % number_nodes
          ipoin = mesh_vol % lnodb(inodb,iboun)
          if( lperm(ipoin) == 0 ) then
             kpoin = kpoin + 1
             lperm(ipoin) = kpoin
          end if
       end do
    end do

    mesh_bou % ndime = mesh_vol % ndime
    mesh_bou % mnode = mesh_vol % mnode
    mesh_bou % nelem = nboun_loc
    mesh_bou % npoin = kpoin

    call mesh_bou % alloca()

    do ipoin = 1,npoin_loc
       kpoin = lperm(ipoin)
       if( kpoin /= 0 ) & 
            mesh_bou % coord(1:mesh_bou % ndime,kpoin) = mesh_vol % coord(1:mesh_bou % ndime,ipoin)
    end do
    
    do iboun = 1,nboun_loc
       pblty = abs(mesh_vol % ltypb(iboun))
       mesh_bou % ltype(iboun) = pblty
       do inodb = 1,element_type(pblty) % number_nodes
          ipoin = mesh_vol % lnodb(inodb,iboun)       
          kpoin = lperm(ipoin)
          mesh_bou % lnods(inodb,iboun) = kpoin
       end do
    end do
    
    call memory_deallo(memor_dom,'LPERM','mod_mesh_type',lperm)
    
  end subroutine mesh_coherent_boundary_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-15
  !> @brief   Copy an extended mesh type
  !> @details Copy an extended mesh typw 
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_copy(meshe_in,meshe_out) 

    type(mesh_type), intent(inout) :: meshe_in
    type(mesh_type), intent(inout) :: meshe_out

    meshe_out % name       = meshe_in % name
    meshe_out % id         = meshe_in % id
    meshe_out % ndime      = meshe_in % ndime
    meshe_out % ntens      = meshe_in % ntens
    meshe_out % npoin      = meshe_in % npoin
    meshe_out % nelem      = meshe_in % nelem
    meshe_out % nboun      = meshe_in % nboun
    meshe_out % nfiel      = meshe_in % nfiel
    meshe_out % kfl_field  = meshe_in % kfl_field
    meshe_out % nedge      = meshe_in % nedge
    meshe_out % mnode      = meshe_in % mnode
    meshe_out % mnodb      = meshe_in % mnodb
    meshe_out % mgaus      = meshe_in % mgaus
    meshe_out % medge      = meshe_in % medge
    meshe_out % nbopo      = meshe_in % nbopo
    meshe_out % npoi1      = meshe_in % npoi1
    meshe_out % npoi2      = meshe_in % npoi2
    meshe_out % npoi3      = meshe_in % npoi3
    meshe_out % nedg1      = meshe_in % nedg1
    meshe_out % nedg2      = meshe_in % nedg2
    meshe_out % nedg3      = meshe_in % nedg3
    meshe_out % npoin_own  = meshe_in % npoin_own
    meshe_out % npoin_halo = meshe_in % npoin_halo
    meshe_out % npoin_2    = meshe_in % npoin
    meshe_out % nelem_2    = meshe_in % nelem
    meshe_out % nboun_2    = meshe_in % nboun
    meshe_out % kfl_ngrou  = meshe_in % kfl_ngrou

    meshe_out % comm % RANK4          = meshe_in % comm % RANK4         
    meshe_out % comm % SIZE4          = meshe_in % comm % SIZE4         
    meshe_out % comm % PAR_COMM_WORLD = meshe_in % comm % PAR_COMM_WORLD

    call memory_copy(memor_dom,'MESH_OUT % LNODS'    ,'mesh_type_save_original_mesh',meshe_in % lnods,    meshe_out % lnods,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LTYPE'    ,'mesh_type_save_original_mesh',meshe_in % ltype,    meshe_out % ltype,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNNOD'    ,'mesh_type_save_original_mesh',meshe_in % lnnod,    meshe_out % lnnod,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LELCH'    ,'mesh_type_save_original_mesh',meshe_in % lelch,    meshe_out % lelch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LESUB'    ,'mesh_type_save_original_mesh',meshe_in % lesub,    meshe_out % lesub,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LMATE'    ,'mesh_type_save_original_mesh',meshe_in % lmate,    meshe_out % lmate,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LEINV_LOC','mesh_type_save_original_mesh',meshe_in % leinv_loc,meshe_out % leinv_loc,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % PERME'    ,'mesh_type_save_original_mesh',meshe_in % perme    ,meshe_out % perme     ,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_OUT % COORD'    ,'mesh_type_save_original_mesh',meshe_in % coord,    meshe_out % coord,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNOCH'    ,'mesh_type_save_original_mesh',meshe_in % lnoch,    meshe_out % lnoch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LMAST'    ,'mesh_type_save_original_mesh',meshe_in % lmast    ,meshe_out % lmast,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNINV_LOC','mesh_type_save_original_mesh',meshe_in % lninv_loc,meshe_out % lninv_loc,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % PERMN'    ,'mesh_type_save_original_mesh',meshe_in % permn    ,meshe_out % permn     ,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_OUT % LNODB'    ,'mesh_type_save_original_mesh',meshe_in % lnodb,    meshe_out % lnodb,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LTYPB'    ,'mesh_type_save_original_mesh',meshe_in % ltypb,    meshe_out % ltypb,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LBOCH'    ,'mesh_type_save_original_mesh',meshe_in % lboch,    meshe_out % lboch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LELBO'    ,'mesh_type_save_original_mesh',meshe_in % lelbo,    meshe_out % lelbo,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LBINV_LOC','mesh_type_save_original_mesh',meshe_in % lbinv_loc,meshe_out % lbinv_loc,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_OUT % EDGE_TO_NODE','mesh_type_save_original_mesh',meshe_in % edge_to_node,    meshe_out % edge_to_node,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LEDGS'    ,'mesh_type_save_original_mesh',meshe_in % ledgs,    meshe_out % ledgs,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LEDGB'    ,'mesh_type_save_original_mesh',meshe_in % ledgb,    meshe_out % ledgb,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNNED'    ,'mesh_type_save_original_mesh',meshe_in % lnned,    meshe_out % lnned,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNNEB'    ,'mesh_type_save_original_mesh',meshe_in % lnneb,    meshe_out % lnneb,    'DO_NOT_DEALLOCATE')

  end subroutine mesh_type_copy


  !-----------------------------------------------------------------------      
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate a minimum size
  !> @details Allocate a minimum size to mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_allocate_minimum(CURRENT_MESH)

    integer(ip), intent(in) :: CURRENT_MESH

    meshe(CURRENT_MESH) % lninv_loc => int_min_1
    meshe(CURRENT_MESH) % leinv_loc => int_min_1
    meshe(CURRENT_MESH) % lbinv_loc => int_min_1

    meshe(CURRENT_MESH) % lnods     => int_min_2
    meshe(CURRENT_MESH) % ltype     => int_min_1
    meshe(CURRENT_MESH) % lelch     => int_min_1
    meshe(CURRENT_MESH) % lnnod     => int_min_1
    meshe(CURRENT_MESH) % lesub     => int_min_1
    meshe(CURRENT_MESH) % lmate     => int_min_1

    meshe(CURRENT_MESH) % coord     => rea_min_2
    meshe(CURRENT_MESH) % lnoch     => int_min_1
    meshe(CURRENT_MESH) % lmast     => int_min_1

    meshe(CURRENT_MESH) % lnodb     => int_min_2
    meshe(CURRENT_MESH) % ltypb     => int_min_1
    meshe(CURRENT_MESH) % lboch     => int_min_1
    meshe(CURRENT_MESH) % lelbo     => int_min_1
    
    meshe(CURRENT_MESH) % edge_to_node => int_min_2
    meshe(CURRENT_MESH) % ledgs     => int_min_2
    meshe(CURRENT_MESH) % ledgb     => int_min_2
    meshe(CURRENT_MESH) % lnned     => int_min_1
    meshe(CURRENT_MESH) % lnneb     => int_min_1

  end subroutine mesh_type_allocate_minimum

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Deallocate mesh type
  !> @details Deallocate mesh type
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_deallocate_s(meshe_in,MESH_NAME)

    type(mesh_type),             intent(inout) :: meshe_in
    character(len=*),  optional, intent(in)    :: MESH_NAME
    character(20)                              :: my_mesh_name
    integer(ip)                                :: ifiel

    if( present(MESH_NAME) ) then
       my_mesh_name = trim(MESH_NAME)
    else
       my_mesh_name = meshe_in % name
    end if

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NPOIN_PAR',trim(vacal),meshe_in % npoin_par)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NELEM_PAR',trim(vacal),meshe_in % nelem_par)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NBOUN_PAR',trim(vacal),meshe_in % nboun_par)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LEINV_LOC',trim(vacal),meshe_in % leinv_loc)               
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnods'    ,trim(vacal),meshe_in % lnods)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ltype'    ,trim(vacal),meshe_in % ltype)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lelch'    ,trim(vacal),meshe_in % lelch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnnod'    ,trim(vacal),meshe_in % lnnod)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lesub'    ,trim(vacal),meshe_in % lesub)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lmate'    ,trim(vacal),meshe_in % lmate)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lgaus'    ,trim(vacal),meshe_in % lgaus)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % perme'    ,trim(vacal),meshe_in % perme)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LBINV_LOC',trim(vacal),meshe_in % lbinv_loc)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LNODB'    ,trim(vacal),meshe_in % lnodb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lboel'    ,trim(vacal),meshe_in % lboel)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lelbo'    ,trim(vacal),meshe_in % lelbo)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ltypb'    ,trim(vacal),meshe_in % ltypb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lboch'    ,trim(vacal),meshe_in % lboch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnnob'    ,trim(vacal),meshe_in % lnnob)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LNINV_LOC',trim(vacal),meshe_in % lninv_loc)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coord'    ,trim(vacal),meshe_in % coord)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnoch'    ,trim(vacal),meshe_in % lnoch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lmast'    ,trim(vacal),meshe_in % lmast)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lpoty'    ,trim(vacal),meshe_in % lpoty)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % permn'    ,trim(vacal),meshe_in % permn)

    do ifiel = 1,meshe_in % nfiel
       call memory_deallo(memor_dom,trim(my_mesh_name)//' % XFIEL % A',trim(vacal),meshe_in % xfiel(ifiel) % a)
    end do

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % edge_to_node',trim(vacal),meshe_in % edge_to_node)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ledgs',trim(vacal),meshe_in % ledgs)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ledgb',trim(vacal),meshe_in % ledgb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnned',trim(vacal),meshe_in % lnned)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnneb',trim(vacal),meshe_in % lnneb)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % leset',trim(vacal),meshe_in % leset)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lbset',trim(vacal),meshe_in % lbset)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnset',trim(vacal),meshe_in % lnset)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % kfl_codno',trim(vacal),meshe_in % kfl_codno)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % kfl_codbo',trim(vacal),meshe_in % kfl_codbo)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LGROU_DOM',trim(vacal),meshe_in % lgrou_dom)

    !call memory_deallo(memor_dom,trim(my_mesh_name)//' % LINNO',trim(vacal),meshe_in % linno)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom   ',trim(vacal),meshe_in % c_dom)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom   ',trim(vacal),meshe_in % r_dom)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coo_rows',trim(vacal),meshe_in % coo_rows)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coo_cols',trim(vacal),meshe_in % coo_cols)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ell_cols',trim(vacal),meshe_in % ell_cols)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_sym   ',trim(vacal),meshe_in % c_sym)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_sym   ',trim(vacal),meshe_in % r_sym)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_edg   ',trim(vacal),meshe_in % c_edg)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_edg   ',trim(vacal),meshe_in % r_edg)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_elm_2 ',trim(vacal),meshe_in % c_elm_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_elm_2 ',trim(vacal),meshe_in % r_elm_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom_2 ',trim(vacal),meshe_in % c_dom_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_2 ',trim(vacal),meshe_in % r_dom_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom_own',trim(vacal),meshe_in % c_dom_own)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_own',trim(vacal),meshe_in % r_dom_own)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_end',trim(vacal),meshe_in % r_dom_end)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_ini',trim(vacal),meshe_in % r_dom_ini)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % EXNOR',trim(vacal),meshe_in % exnor)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % VMASS',trim(vacal),meshe_in % vmass)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % VMASC',trim(vacal),meshe_in % vmasc)


  end subroutine mesh_type_deallocate_s

  subroutine mesh_type_deallocate_all(NUMBER_MESHES)

    integer(ip), optional, intent(in) :: NUMBER_MESHES
    integer(ip)                       :: idivi
    integer(ip)                       :: ndivi_loc

    if( present(NUMBER_MESHES) ) then
       ndivi_loc = NUMBER_MESHES
    else
       ndivi_loc = ndivi
    end if

    if( INOTSLAVE ) then
       do idivi = 0,ndivi_loc
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NPOIN_PAR','mesh_type_save_original_mesh',meshe(idivi) % npoin_par)
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NELEM_PAR','mesh_type_save_original_mesh',meshe(idivi) % nelem_par)
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NBOUN_PAR','mesh_type_save_original_mesh',meshe(idivi) % nboun_par)
       end do
    end if
    deallocate( meshe )

  end subroutine mesh_type_deallocate_all

  !-----------------------------------------------------------------------
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate and initialize
  !> @details Allocate and initiaize mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_initialize_1(meshe_in)

    type(mesh_type), pointer, intent(inout) :: meshe_in(:)
    integer(ip)                    :: ii

    if( associated(meshe_in) ) then
       do ii = lbound(meshe_in,1),ubound(meshe_in,1)
          call mesh_type_initialize_s(meshe_in(ii))
       end do
    end if

  end subroutine mesh_type_initialize_1

  subroutine mesh_type_initialize_s(meshe_in,wname)

    type(mesh_type),            intent(inout) :: meshe_in
    character(len=*), optional, intent(in)    :: wname

    call meshe_in % init(wname)

  end subroutine mesh_type_initialize_s

  !-----------------------------------------------------------------------
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate and initialize
  !> @details Allocate and initiaize mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_allocate_initialize(NUMBER_MESHES)

    integer(ip), optional, intent(in) :: NUMBER_MESHES
    integer(ip)                       :: idivi
    integer(ip)                       :: ndivi_loc

    if( present(NUMBER_MESHES) ) then
       ndivi_loc = NUMBER_MESHES
    else
       ndivi_loc = ndivi
    end if

    allocate( meshe(-1:ndivi_loc) )

    do idivi = -1,ndivi_loc
       call mesh_type_initialize(meshe(idivi))
       meshe(idivi) % name = trim(title)//'_MM'//trim(intost(idivi))
    end do

    if( INOTSLAVE ) then
       do idivi = 0,ndivi_loc
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NPOIN_PAR','mesh_type_save_original_mesh',meshe(idivi) % npoin_par,npart)
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NELEM_PAR','mesh_type_save_original_mesh',meshe(idivi) % nelem_par,npart)
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NBOUN_PAR','mesh_type_save_original_mesh',meshe(idivi) % nboun_par,npart) 
       end do
    end if

  end subroutine mesh_type_allocate_initialize

  !-----------------------------------------------------------------------
  !>
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Save original mesh
  !> @details Mesh has been divided and postprocess is on original mesh
  !>          Save the original mesh in MESHE(0)
  !>
  !-----------------------------------------------------------------------

  subroutine mesh_type_save_original_mesh()

    integer(ip) :: ipart,ipoin,ielem,iboun
    !
    ! Save LNINV_LOC of original mesh
    ! It is useful to refer to original mesh numbering
    !   - IPOIN = lninv_loc(lnlev(JPOIN))
    !     - IPOIN:  original numbering
    !     - JPOIN:  local numbering
    !     - Check that lnlev(JPOIN) /= 0
    ! ISEQUEN does it too to avoid many if's
    ! MASTER allocate minimum memory
    ! Save also LEINV_LOC of original mesh
    !
    if( INOTMASTER ) then

       !if( ISEQUEN ) then
       !   call memory_alloca(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',lninv_loc,npoin,'IDENTITY')
       !   call memory_alloca(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',leinv_loc,nelem,'IDENTITY')
       !   call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',lbinv_loc,nboun,'IDENTITY')
       !end if

       if( ndivi > 0 ) then
          call memory_alloca(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',meshe(0) % lninv_loc,npoin)
          do ipoin = 1,npoin
             meshe(0) % lninv_loc(ipoin) = lninv_loc(ipoin)
          end do
          call memory_alloca(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',meshe(0) % leinv_loc,nelem)
          do ielem = 1,nelem
             meshe(0) % leinv_loc(ielem) = leinv_loc(ielem)
          end do
          call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',meshe(0) % lbinv_loc,max(1_ip,nboun))
          do iboun = 1,nboun
             meshe(0) % lbinv_loc(iboun) = lbinv_loc(iboun)
          end do
       else
          meshe(0) % lninv_loc => lninv_loc
          meshe(0) % leinv_loc => leinv_loc
          meshe(0) % lbinv_loc => lbinv_loc
       end if

    else if( IMASTER ) then

       call mesh_type_allocate_minimum(0_ip)

    end if

    if( ndivi > 0 .and. kfl_posdi == 0 ) then
       !
       ! Mesh has been divided and postprocess is on original mesh
       !
       if( IMASTER ) then
          !
          ! Master and sequen save global geometry parameters
          !
          meshe(0) % npoin_origi = npoin_origi
          meshe(0) % npoin_total = npoin_total
          meshe(0) % nelem_total = nelem_total
          meshe(0) % nboun_total = nboun_total
          do ipart = 1,npart
             meshe(0) % npoin_par(ipart) = npoin_par(ipart)
             meshe(0) % nelem_par(ipart) = nelem_par(ipart)
             meshe(0) % nboun_par(ipart) = nboun_par(ipart)
          end do
       end if

       if( INOTMASTER ) then
          !
          ! Slaves and sequen save original geometry
          !
          meshe(0) % ndime      = ndime
          meshe(0) % ntens      = ntens
          meshe(0) % npoin      = npoin
          meshe(0) % nelem      = nelem
          meshe(0) % nboun      = nboun
          meshe(0) % nfiel      = nfiel
          meshe(0) % kfl_field  = kfl_field
          meshe(0) % nedge      = nedge
          meshe(0) % mnode      = mnode
          meshe(0) % mnodb      = mnodb
          meshe(0) % mgaus      = mgaus
          meshe(0) % medge      = medge
          meshe(0) % nbopo      = nbopo
          meshe(0) % npoi1      = npoi1
          meshe(0) % npoi2      = npoi2
          meshe(0) % npoi3      = npoi3
          meshe(0) % nedg1      = nedg1
          meshe(0) % nedg2      = nedg2
          meshe(0) % nedg3      = nedg3
          meshe(0) % npoin_own  = npoin_own
          meshe(0) % npoin_halo = npoin_halo
          meshe(0) % npoin_2    = npoin
          meshe(0) % nelem_2    = nelem
          meshe(0) % nboun_2    = nboun
          meshe(0) % kfl_ngrou  = kfl_ngrou

          call memory_copy(memor_dom,'LNODS','mesh_type_save_original_mesh',lnods,meshe(0) % lnods,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LTYPE','mesh_type_save_original_mesh',ltype,meshe(0) % ltype,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNNOD','mesh_type_save_original_mesh',lnnod,meshe(0) % lnnod,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LELCH','mesh_type_save_original_mesh',lelch,meshe(0) % lelch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LESUB','mesh_type_save_original_mesh',lesub,meshe(0) % lesub,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LMATE','mesh_type_save_original_mesh',lmate,meshe(0) % lmate,'DO_NOT_DEALLOCATE')

          call memory_copy(memor_dom,'COORD','mesh_type_save_original_mesh',coord,meshe(0) % coord,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNOCH','mesh_type_save_original_mesh',lnoch,meshe(0) % lnoch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LMAST','mesh_type_save_original_mesh',lmast,meshe(0) % lmast,'DO_NOT_DEALLOCATE')

          call memory_copy(memor_dom,'LNODB','mesh_type_save_original_mesh',lnodb,meshe(0) % lnodb,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LTYPB','mesh_type_save_original_mesh',ltypb,meshe(0) % ltypb,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LBOCH','mesh_type_save_original_mesh',lboch,meshe(0) % lboch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LELBO','mesh_type_save_original_mesh',lelbo,meshe(0) % lelbo,'DO_NOT_DEALLOCATE')

          call memory_copy(memor_dom,'EDGE_TO_NODE','mesh_type_save_original_mesh',edge_to_node,meshe(0) % edge_to_node,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LEDGS','mesh_type_save_original_mesh',ledgs,meshe(0) % ledgs,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LEDGB','mesh_type_save_original_mesh',ledgb,meshe(0) % ledgb,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNNED','mesh_type_save_original_mesh',lnned,meshe(0) % lnned,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNNEB','mesh_type_save_original_mesh',lnneb,meshe(0) % lnneb,'DO_NOT_DEALLOCATE')

       end if
       !
       ! Communication
       !
       meshe(0) % comm % RANK4          = int(PAR_MY_CODE_RANK,4)        
       meshe(0) % comm % SIZE4          = int(PAR_CODE_SIZE,4) 
       meshe(0) % comm % PAR_COMM_WORLD = PAR_COMM_MY_CODE
       !
       ! Name
       !
       meshe(0) % name = 'ORIGINAL'
       meshe(0) % id   = 0
    end if

  end subroutine mesh_type_save_original_mesh

  !-----------------------------------------------------------------------
  !>
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Update last mesh
  !> @details Update last mesh upon changes and deallocations
  !>          e.g. repoint to last mesh after computing halo geometry
  !>          Called by mod_ghost_geometry.f90
  !>
  !-----------------------------------------------------------------------

  subroutine mesh_type_update_last_mesh(CURRENT_MESH)

    integer(ip), optional, intent(in) :: CURRENT_MESH
    integer(ip)                       :: idivi,ipart,ifiel

    if( present(CURRENT_MESH) ) then
       idivi = CURRENT_MESH
    else
       idivi = ndivi
    end if

    if( associated(meshe) .and. size(meshe,kind=ip) >= idivi ) then       
       !
       ! Dimensions
       !
       meshe(idivi) % ndime       =  ndime
       meshe(idivi) % ntens       =  ntens
       meshe(idivi) % npoin       =  npoin
       meshe(idivi) % nelem       =  nelem
       meshe(idivi) % nboun       =  nboun
       meshe(idivi) % nfiel       =  nfiel
       meshe(idivi) % nedge       =  nedge
       meshe(idivi) % kfl_field   =  kfl_field
       meshe(idivi) % mnode       =  mnode
       meshe(idivi) % mnodb       =  mnodb
       meshe(idivi) % mgaus       =  mgaus
       meshe(idivi) % medge       =  medge
       meshe(idivi) % nbopo       =  nbopo     ! Computed in extnor
       meshe(idivi) % npoi1       =  npoi1
       meshe(idivi) % npoi2       =  npoi2
       meshe(idivi) % npoi3       =  npoi3
       meshe(idivi) % nedg1       =  nedg1
       meshe(idivi) % nedg2       =  nedg2
       meshe(idivi) % nedg3       =  nedg3
       meshe(idivi) % npoin_own   =  npoin_own
       meshe(idivi) % npoin_halo  =  npoin_halo
       meshe(idivi) % kfl_ngrou   =  kfl_ngrou

       meshe(idivi) % nelem_2     =  nelem_2
       meshe(idivi) % nboun_2     =  nboun_2
       meshe(idivi) % npoin_2     =  npoin_2
       !
       ! Dimensions of master
       !
       if( INOTSLAVE ) then
          meshe(ndivi) % npoin_origi = npoin_origi
          meshe(ndivi) % npoin_total = npoin_total
          meshe(ndivi) % nelem_total = nelem_total
          meshe(ndivi) % nboun_total = nboun_total
          do ipart = 1,npart
             meshe(ndivi) % npoin_par(ipart) = npoin_par(ipart)
             meshe(ndivi) % nelem_par(ipart) = nelem_par(ipart)
             meshe(ndivi) % nboun_par(ipart) = nboun_par(ipart)
          end do
       end if
       !
       ! Element arrays
       !
       meshe(idivi) % leinv_loc   => leinv_loc
       meshe(idivi) % lnods       => lnods
       meshe(idivi) % ltype       => ltype
       meshe(idivi) % lnnod       => lnnod
       meshe(idivi) % lelch       => lelch
       meshe(idivi) % lesub       => lesub
       meshe(idivi) % lmate       => lmate
       meshe(idivi) % lgaus       => lgaus
       !
       ! Boundary arrays
       !     
       meshe(idivi) % lbinv_loc   => lbinv_loc
       meshe(idivi) % lnodb       => lnodb
       meshe(idivi) % ltypb       => ltypb
       meshe(idivi) % lboch       => lboch
       meshe(idivi) % lboel       => lboel
       meshe(idivi) % lelbo       => lelbo
       meshe(idivi) % lnnob       => lnnob
       !
       ! Nodal arrays
       !     
       meshe(idivi) % lninv_loc   => lninv_loc
       meshe(idivi) % coord       => coord
       meshe(idivi) % lnoch       => lnoch
       meshe(idivi) % lmast       => lmast
       meshe(idivi) % lpoty       => lpoty
       !
       ! Edge arrays
       !     
       meshe(idivi) % edge_to_node   => edge_to_node
       meshe(idivi) % ledgs       => ledgs
       meshe(idivi) % ledgb       => ledgb
       meshe(idivi) % lnned       => lnned
       meshe(idivi) % lnneb       => lnneb
       !
       ! Others
       !
       do ifiel = 1,meshe(idivi) % nfiel
          meshe(idivi) % xfiel(ifiel) % a => xfiel(ifiel) % a
       end do
       !
       ! Graphs
       !
       meshe(idivi) % nzdom       =  nzdom
       meshe(idivi) % r_dom       => r_dom        ! Computed in domgra
       meshe(idivi) % c_dom       => c_dom        ! Computed in domgra
       !
       ! Geometrical arrays
       !
       meshe(idivi) % exnor       => exnor        ! Computed in extnor
       meshe(idivi) % vmass       => vmass        ! Computed in massma
       meshe(idivi) % vmasc       => vmasc        ! Computed in massmc
       !
       ! Communication
       !
       meshe(idivi) % comm % RANK4          =  int(PAR_MY_CODE_RANK,4)        
       meshe(idivi) % comm % SIZE4          =  int(PAR_CODE_SIZE,4)    
       meshe(idivi) % comm % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
       if( IPARALL ) then
          meshe(idivi) % comm % nneig          =  commd % nneig
          meshe(idivi) % comm % bound_dim      =  commd % bound_dim
          meshe(idivi) % comm % neights        => commd % neights
          meshe(idivi) % comm % bound_perm     => commd % bound_perm
          meshe(idivi) % comm % bound_size     => commd % bound_size
       end if
       !
       ! Boundary
       !
       if( .not. associated(meshe(idivi) % boundary) ) allocate(meshe(idivi) % boundary)
       meshe(idivi) % boundary % ndime        =  meshe(idivi) % ndime
       meshe(idivi) % boundary % mnode        =  meshe(idivi) % mnodb
       meshe(idivi) % boundary % nelem        =  meshe(idivi) % nboun
       meshe(idivi) % boundary % npoin        =  meshe(idivi) % npoin
       meshe(idivi) % boundary % lnods        => meshe(idivi) % lnodb
       meshe(idivi) % boundary % ltype        => meshe(idivi) % ltypb
       meshe(idivi) % boundary % leinv_loc    => meshe(idivi) % lbinv_loc
       meshe(idivi) % boundary % lninv_loc    => meshe(idivi) % lninv_loc
       meshe(idivi) % boundary % coord        => meshe(idivi) % coord
       meshe(idivi) % boundary % lelbo        => meshe(idivi) % lelbo
       meshe(idivi) % boundary % lboel        => meshe(idivi) % lboel
       meshe(idivi) % boundary % mesh         => meshe(idivi)
       !
       ! Name
       !
       meshe(idivi) % name = 'MESHE('//trim(intost(idivi))//')' !'LAST_MESH'
       meshe(idivi) % id   = 0
       
    end if

  end subroutine mesh_type_update_last_mesh

  subroutine mesh_type_basic_to_complete(mesh_basic,mesh_complete)

    type(mesh_type_basic), intent(inout) :: mesh_basic
    type(mesh_type),       intent(inout) :: mesh_complete

    mesh_complete % ndime                 =   mesh_basic % ndime
    mesh_complete % mnode                 =   mesh_basic % mnode
    mesh_complete % mnodb                 =   mnodb
    mesh_complete % npoin                 =   mesh_basic % npoin
    mesh_complete % nelem                 =   mesh_basic % nelem
    mesh_complete % lnods                 =>  mesh_basic % lnods
    mesh_complete % ltype                 =>  mesh_basic % ltype
    mesh_complete % leinv_loc             =>  mesh_basic % leinv_loc
    mesh_complete % lninv_loc             =>  mesh_basic % lninv_loc
    mesh_complete % coord                 =>  mesh_basic % coord

    mesh_complete % comm % RANK4          =   mesh_basic % comm % RANK4
    mesh_complete % comm % SIZE4          =   mesh_basic % comm % SIZE4
    mesh_complete % comm % PAR_COMM_WORLD =   mesh_basic % comm % PAR_COMM_WORLD
    mesh_complete % comm % bound_dim      =   mesh_basic % comm % bound_dim
    mesh_complete % comm % nneig          =   mesh_basic % comm % nneig
    mesh_complete % comm % neights        =>  mesh_basic % comm % neights
    mesh_complete % comm % bound_size     =>  mesh_basic % comm % bound_size
    mesh_complete % comm % bound_perm     =>  mesh_basic % comm % bound_perm

  end subroutine mesh_type_basic_to_complete

end module mod_mesh_type
!> @}
