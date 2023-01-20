!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    edge_data_structures.f90
!> @date    04/02/2016
!> @author  Guillaume Houzeaux
!> @brief   Edge data structures
!> @details Compute edge data structures in MESHE(NDIVI):
!>          \verbatim
!>          MEDGE ....................... Maximum number of edges per elements in mesh
!>          NEDGE ....................... Tolal number of edges in the mesh
!>          LNNED(1:NELEM) .............. Element number of edges
!>          LNNEB(1:NBOUN) .............. Boundary number of edges
!>          LEDGS(1:PEDGE,1:NELEM) ...... List of global edges for all elements
!>          LEDGB(1:PEDGE,1:NBOUN) ...... List of global edges for all boundaries
!>          R_EDG(:), C_EDG(:) .......... Edge graph
!>          EDGE_TO_NODE(1:2,1:NEDGE) ... Edge nodes
!>          \endverbatim
!>
!> @} 
!-----------------------------------------------------------------------

subroutine edge_data_structures()
  use def_kintyp,                        only : ip
  use def_kermod,                        only : kfl_edge_elements
  use def_kermod,                        only : kfl_edge_graph
  use def_kermod,                        only : ndivi
  use def_domain,                        only : meshe
  use def_domain,                        only : medge
  use def_domain,                        only : nedge
  use def_domain,                        only : edge_to_node
  use def_domain,                        only : ledgs
  use def_domain,                        only : lnned
  use def_domain,                        only : ledgb
  use def_domain,                        only : lnneb
  use def_domain,                        only : r_edg
  use def_domain,                        only : c_edg
  use def_domain,                        only : nedg3
  use def_domain,                        only : memor_dom
  use mod_domain,                        only : domain_memory_allocate
  use mod_graphs,                        only : graphs_edges
  use def_master,                        only : INOTMASTER
  use def_master,                        only : IPARALL
  use def_master,                        only : kfl_paral
  use def_master,                        only : npart
  use def_master,                        only : lginv_loc
  use def_master,                        only : lninv_loc
  use def_master,                        only : nedg1,nedg2,nedg3
  use def_master,                        only : namda
  use mod_parall,                        only : PAR_COMM_MY_CODE_ARRAY
  use mod_messages,                      only : messages_live
  use mod_memory,                        only : memory_alloca
  use mod_memory,                        only : memory_deallo
  use mod_communications_global,         only : PAR_ALLGATHER
  use mod_communications_point_to_point, only : PAR_INTERFACE_EDGE_EXCHANGE
  use mod_matrix_market,                 only : matrix_market_matrix
  use mod_matrix_market,                 only : matrix_market_vector
  use mod_strings,                       only : integer_to_string
  implicit none

  integer(ip)          :: nedge_offset,iedge
  integer(ip), pointer :: nedge_own_tot(:)

  if( kfl_edge_elements == 1 ) then
     
     call messages_live('COMPUTE EDGE DATA STRUCTURES')
     !
     ! Edge arrays
     !
     if( INOTMASTER ) call graphs_edges(&
          meshe(ndivi),medge,nedge,edge_to_node,&
          ledgs,lnned,ledgb,lnneb,r_edg,&
          c_edg,POINT_TO_TYPE=.true.,memor=memor_dom)
     !call check_edge_data_structures()
     !
     ! Communication arrays: edges are renumbered to account
     ! for interior, own and other boundary edges
     !
     if( IPARALL ) then
        call par_edge_communication_array(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY)
     else
        nedg1 = nedge
        nedg2 = nedge
        nedg3 = nedge
     end if
     !
     ! Create a global numbering     
     !
     nullify(nedge_own_tot)
     call domain_memory_allocate('LGINV_LOC')
     call memory_alloca(memor_dom,'nedge_own_tot','renumbering_lexical_order',nedge_own_tot,npart+1,'INITIALIZE',0_ip)
     call PAR_ALLGATHER(nedg3,nedge_own_tot)
     if( INOTMASTER ) then
        nedge_offset = sum(nedge_own_tot(0:kfl_paral-1))
        do iedge = 1,nedg3        
           lginv_loc(iedge) = iedge + nedge_offset
        end do
        call PAR_INTERFACE_EDGE_EXCHANGE(lginv_loc,'SUM','IN MY CODE')
     end if
     call memory_deallo(memor_dom,'nedge_own_tot','renumbering_lexical_order',nedge_own_tot)
     !
     ! Output edge data structure
     !
     if( kfl_edge_graph == 1 ) then
        call matrix_market_matrix(1_ip,1_ip,nedge,r_edg,c_edg,FILENAME=trim(namda)//'-edge-graph'               //integer_to_string(kfl_paral)//'.mtx')
        call matrix_market_vector(lginv_loc,                  FILENAME=trim(namda)//'edge-global-numbering'     //integer_to_string(kfl_paral)//'.mtx')
        call matrix_market_vector(meshe(ndivi) % edge_to_node,FILENAME=trim(namda)//'edge-node-connectivity-loc'//integer_to_string(kfl_paral)//'.mtx')
        call matrix_market_vector(meshe(ndivi) % edge_to_node,FILENAME=trim(namda)//'edge-node-connectivity-glo'//integer_to_string(kfl_paral)//'.mtx',PERMX=lninv_loc)
     end if
     
  end if

end subroutine edge_data_structures

subroutine check_edge_data_structures()
  use def_domain
  use def_master
  use def_kermod
  implicit none
  integer(ip) :: ielem,iedge,izdom

  print*,'MEDGE=',meshe(ndivi) % medge
  print*,'NEDGE=',meshe(ndivi) % nedge
  print*,'CONNECTIVITY:'
  do ielem = 1,nelem
     print*,ielem,meshe(ndivi) % ledgs(:,ielem)
  end do
  print*,'EDGES:'
  do iedge = 1,meshe(ndivi) % nedge
     print*,iedge,meshe(ndivi) % edge_to_node(1:2,iedge)
  end do 
  print*,'NUMBER:'
  do ielem = 1,meshe(ndivi) % nelem
     print*,iedge,meshe(ndivi) % lnned(ielem)
  end do 
  print*,'GRAPHS:'
  do iedge = 1,meshe(ndivi) % nedge
     do izdom = meshe(ndivi) % r_edg(iedge),meshe(ndivi) % r_edg(iedge+1)-1
        print*,iedge,meshe(ndivi) % c_edg(izdom)
     end do
  end do 

  call runend('O.K.!')

end subroutine check_edge_data_structures
