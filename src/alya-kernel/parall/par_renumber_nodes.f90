!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 !-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_renumber_nodes.f90
!> @date    21/03/2017
!> @author  Guillaume Houzeaux
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
!> @}
!-----------------------------------------------------------------------

subroutine par_renumber_nodes()

  use def_kintyp,      only : ip,rp,lg
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
  use def_domain,      only : npoin
  use def_domain,      only : mnode
  use def_domain,      only : nelem
  use def_domain,      only : lnods
  use def_domain,      only : ltype
  use def_domain,      only : lnnod
  use def_domain,      only : coord
  use def_domain,      only : memor_dom
  use mod_parall,      only : par_memor
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_renumbering, only : renumbering_node_arrays
  use mod_renumbering, only : renumbering_nodes
  use mod_renumbering, only : renumbering_sfc_recursive
  use mod_renumbering, only : renumbering_reverse_cuthill_mckee
  use mod_renumbering, only : renumbering_update
  use mod_graphs,      only : graphs_poipoi
  use mod_graphs,      only : graphs_poipoi_deallocate
  use mod_mesh_type,   only : mesh_type_update_last_mesh
  use mod_messages,    only : messages_live
  use mod_alya2metis,  only : alya2metis_METIS_NodeND
 
  use def_master
  implicit none

  integer(ip)          :: ipoin,kpoin
  integer(ip), pointer :: permR(:)
  character(15)        :: my_method
  !
  ! Interior graph
  !
  integer(ip), pointer :: ia(:)
  integer(ip), pointer :: ja(:)

  nullify(permR)
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
   
  if( npoin > 0 ) then

     call memory_alloca(par_memor,'permR','par_renumber_nodes',permR,npoin)
     !
     ! PERMR: Identify boundary nodes. Others'=-1, Own=-2
     !
     do ipoin = npoi1+1,npoi2-1
        permR(ipoin) = -1
     end do
     do ipoin = npoi2,npoi3
        permR(ipoin) = -2
     end do
     do ipoin = npoi3+1,npoin
        permR(ipoin) = -1
     end do
     !
     ! Renumber interior nodes
     !
     kpoin = npoi1
     do ipoin = 1,npoi1
        permR(ipoin) = ipoin
     end do
     !
     ! Renumber own boundary nodes
     !
     npoi1 = kpoin
     npoi2 = kpoin + 1
     do ipoin = 1,npoin
        if( permR(ipoin) == -2 ) then
           kpoin        = kpoin + 1
           permR(ipoin) = kpoin
        end if
     end do
     npoi3 = kpoin
     !
     ! Renumber others' boundary nodes
     !
     do ipoin = 1,npoin
        if( permR(ipoin) == -1 ) then
           kpoin        = kpoin + 1
           permR(ipoin) = kpoin
        end if
     end do
     
     if( npoi1 >  0 ) then     
        !
        ! IA, JA: graph of interior nodes without diagonal
        !
        call graphs_poipoi(&
             npoi1,nelem,mnode,lnods,lnnod,ltype,ia,ja,'SQUARE REMOVE DIAGONAL',memor=memor_dom)        
        !
        ! Reorder interior nodes 
        !
        if(      kfl_renumbering_npoin == 1 .and. associated(ja) ) then
           !
           ! METIS renumbering
           !
           call alya2metis_METIS_NodeND(npoi1,ia,ja,permr)
           
        else if( kfl_renumbering_npoin == 2 ) then
           !
           ! SFC renumbering
           !
           call renumbering_sfc_recursive(nsfc_renumbering_npoin,coord,permr,npoi1)
           
        else if( kfl_renumbering_npoin == 3 ) then
           !
           ! Cuthill-McKee renumbering
           !
           call renumbering_reverse_cuthill_mckee(npoi1,ia,ja,permr)
           
        else
           !
           ! No numbering
           !
        end if
     end if
     !
     ! Permute nodal arrays
     !
     call renumbering_node_arrays(permR)
     !
     ! Deallocate memory
     !
     call graphs_poipoi_deallocate(ia,ja,memor_dom)
     call memory_deallo(par_memor,'PERMR','par_renumber_nodes',permR)

  end if
  !
  ! Update some mesh parameters, NPOIN_OWN, NPOIN_HALO, etc.  
  !
  call renumbering_update()  

end subroutine par_renumber_nodes
