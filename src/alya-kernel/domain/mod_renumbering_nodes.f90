!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_renumbering_nodes.f90
!> @author  houzeaux
!> @date    2019-05-07
!> @brief   Renumber nodes
!> @details Renumber nodes in this order, with npoi2 = npoi1+1.
!>          Interior nodes are renumbered using METIS.
!>
!>          1       -> npoi1: interior
!>          npoi2   -> npoi3: own boundary
!>          npoi3+1 -> npoin: oth boundary
!>
!>          Then, halos will also be renumbered later on. New and old
!>          notations are:
!>
!>          NEW RENUMBERING:
!>
!>          +---------------+------------+------------+----------+-------+
!>          |   interior    | own bound. | oth bound. |      geom halo   |
!>          +---------------+------------+------------+----------+-------+
!>
!>          ---------------->------------>------------>------------------>
!>          1          npoi1 npoi2  npoi3         npoin            npoin_2
!>
!>
!>          AFTER RENUMBERING:
!>
!>          +----------------------------+-----------------+-------------+
!>          |    Own nodes               |    comp halo    | rest halo   |
!>          +----------------------------+-----------------+-------------+
!>
!>          ----------------------------->----------------->------------->
!>          1                    npoin_own       npoin_halo        npoin_2
!>
!-----------------------------------------------------------------------

module mod_renumbering_nodes

  use def_kintyp,      only : ip,rp
  use def_master,      only : ISLAVE
  use def_master,      only : INOTSLAVE
  use def_master,      only : ISEQUEN
  use def_master,      only : IMASTER
  use def_master,      only : INOTMASTER
  use def_master,      only : npoi1
  use def_master,      only : npoi2
  use def_master,      only : npoi3
  use def_kermod,      only : kfl_renumbering_npoin
  use def_kermod,      only : nsfc_renumbering_npoin
  use def_domain,      only : npoin_own
  use def_domain,      only : npoin_halo
  use def_domain,      only : npoin
  use def_domain,      only : mnode
  use def_domain,      only : nelem
  use def_domain,      only : lnods
  use def_domain,      only : ltype
  use def_domain,      only : lnnod
  use def_domain,      only : coord
  use def_domain,      only : memor_dom
  use mod_parall,      only : par_memor
  use mod_parall,      only : PAR_COMM_MY_CODE_ARRAY
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_renumbering, only : renumbering_node_arrays
  use mod_renumbering, only : renumbering_sfc_recursive
  use mod_renumbering, only : renumbering_reverse_cuthill_mckee
  use mod_graphs,      only : graphs_poipoi
  use mod_mesh_type,   only : mesh_type_update_last_mesh
  use mod_elmgeo,      only : element_type
  use mod_messages,    only : messages_live
  use mod_alya2metis,  only : alya2metis_METIS_NodeND

  use def_master
  implicit none
  private

  integer(ip), pointer :: permR_nodes(:)

  public :: renumbering_nodes
  public :: renumbering_nodes_initialization
  public :: permR_nodes

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-07
  !> @brief   Renumber nodes
  !> @details Renumber nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine renumbering_nodes_initialization()

    nullify(permR_nodes)

  end subroutine renumbering_nodes_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-07
  !> @brief   Renumber nodes
  !> @details Renumber nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine renumbering_nodes()

    integer(ip)          :: ipoin,kpoin,ielem,pelty
    integer(ip), pointer :: lnnod_loc(:)
    character(15)        :: my_method
    !
    ! Interior graph
    !
    integer(ip), pointer :: ia(:)
    integer(ip), pointer :: ja(:)
    
    nullify(lnnod_loc)
    nullify(ia)
    nullify(ja)

    if( INOTSLAVE ) then
       if(      kfl_renumbering_npoin == 0 ) then
          my_method = 'NOTHING!...'
       else if( kfl_renumbering_npoin == 1 ) then
          my_method = 'METIS'
       else if( kfl_renumbering_npoin == 2 ) then
          my_method = 'SFC'
       else if( kfl_renumbering_npoin == 3 ) then
          my_method = 'CUTHILL-MCKEE'
       endif
       call messages_live('NODE RENUMBERING WITH '//trim(my_method))
    end if

    if( INOTMASTER ) then

       call memory_alloca(par_memor,'PERMR_NODES','par_renumber_nodes',permr_nodes,npoin)
       !
       ! PERMR_NODES: Identify boundary nodes. Others'=-1, Own=-2
       !         
       do ipoin = npoi1+1,npoi2-1
          permr_nodes(ipoin) = -1
       end do
       do ipoin = npoi2,npoi3
          permr_nodes(ipoin) = -2
       end do
       do ipoin = npoi3+1,npoin
          permr_nodes(ipoin) = -1
       end do
       !
       ! Renumber interior nodes
       !
       kpoin = npoi1
       do ipoin = 1,npoi1
          permr_nodes(ipoin) = ipoin
       end do
       !
       ! Renumber own boundary nodes
       !
       npoi1 = kpoin
       npoi2 = kpoin + 1
       do ipoin = 1,npoin
          if( permr_nodes(ipoin) == -2 ) then
             kpoin              = kpoin + 1
             permr_nodes(ipoin) = kpoin
          end if
       end do
       npoi3 = kpoin
       !
       ! Renumber others' boundary nodes
       !
       do ipoin = 1,npoin
          if( permr_nodes(ipoin) == -1 ) then
             kpoin              = kpoin + 1
             permr_nodes(ipoin) = kpoin
          end if
       end do

       if( npoi1 > 0 ) then
          !
          ! IA, JA: graph of interior nodes without diagonal
          !
          if( kfl_renumbering_npoin == 1 .or. kfl_renumbering_npoin == 3 ) then
             if( .not. associated(lnnod) ) then
                call memory_alloca(par_memor,'LNNOD_LOC','par_renumber_nodes',lnnod_loc,nelem)
                do ielem = 1,nelem
                   pelty = abs(ltype(ielem))
                   lnnod_loc(ielem) = element_type(pelty) % number_nodes
                end do
             else
                lnnod_loc => lnnod
             end if
             call graphs_poipoi(&
                  npoi1,nelem,mnode,lnods,lnnod_loc,ltype,ia,ja,'SQUARE REMOVE DIAGONAL',memor=memor_dom)        
          end if
          !
          ! Reorder interior nodes 
          !
          if(      kfl_renumbering_npoin == 1 ) then
             !
             ! METIS renumbering
             !
             if( associated(ja) ) then
                call alya2metis_METIS_NodeND(npoi1,ia,ja,permr_nodes)
             else
                do ipoin = 1,npoi1
                   permr_nodes(ipoin) = ipoin
                end do
             end if
             
          else if( kfl_renumbering_npoin == 2 ) then
             !
             ! SFC renumbering
             !
             call renumbering_sfc_recursive(nsfc_renumbering_npoin,coord,permr_nodes,npoi1)

          else if( kfl_renumbering_npoin == 3 ) then
             !
             ! Cuthill-McKee renumbering
             !
             if( associated(ja) ) then
                call renumbering_reverse_cuthill_mckee(npoi1,ia,ja,permr_nodes)
             else
                do ipoin = 1,npoi1
                   permr_nodes(ipoin) = ipoin
                end do
             end if
             
          else
             !
             ! No numbering
             !
          end if
          !
          ! Deallocate if necessary
          !
          if( .not. associated(lnnod) ) &
               call memory_deallo(par_memor,'LNNOD_LOC','par_renumber_nodes',lnnod_loc)
       end if
       !
       ! Permute nodal arrays
       !
       call renumbering_node_arrays(permr_nodes)
       !
       ! Deallocate memory
       !
       !call memory_deallo(par_memor,'PERMR_NODES','par_renumber_nodes',permr_nodes)
       call memory_deallo(par_memor,'IA'   ,'par_renumber_nodes',ia)
       call memory_deallo(par_memor,'JA'   ,'par_renumber_nodes',ja)

    end if

    if( ISLAVE ) then
       !
       ! Useful parameters
       !
       npoin_own  = npoi3
       npoin_halo = npoin_own
       !
       ! Update main communicator
       !
       PAR_COMM_MY_CODE_ARRAY(1) % npoi1 = npoi1
       PAR_COMM_MY_CODE_ARRAY(1) % npoi2 = npoi2
       PAR_COMM_MY_CODE_ARRAY(1) % npoi3 = npoi3

    else if( ISEQUEN ) then

       npoi3      = npoin
       npoin_own  = npoin
       npoin_halo = npoin

    else if( IMASTER ) then

       npoi3      = 0
       npoin_own  = 0
       npoin_halo = 0

    end if
    !
    ! Point to mesh structure => MESHE(NDIVI)
    !
    call mesh_type_update_last_mesh()

  end subroutine renumbering_nodes

end module mod_renumbering_nodes
!> @}
