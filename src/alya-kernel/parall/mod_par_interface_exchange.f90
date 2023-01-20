!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_interface_exchange.f90
!> @author  houzeaux
!> @date    2018-05-30
!> @brief   Interface exchange
!> @details Set up node interface exchange communication arrays
!>          in parallel. It uses the coupling
!>
!-----------------------------------------------------------------------

module mod_par_interface_exchange

  use def_kintyp_basic,          only : ip,rp,lg
  use def_kintyp_comm,           only : comm_data_par
  use def_master,                only : IMASTER
  use def_master,                only : INOTMASTER
  use def_master,                only : kfl_paral
  use def_master,                only : zeror
  use def_master,                only : current_code
  use def_master,                only : npart
  use def_master,                only : npoin_par
  use def_master,                only : nelem_par
  use def_master,                only : nboun_par
  use def_master,                only : npoi1,npoi2,npoi3
  use def_master,                only : lninv_loc
  use def_domain,                only : nelem,npoin,lnods,ltype
  use def_domain,                only : nboun,nboun_total,ndime
  use def_domain,                only : coord,npoin_own
  use def_coupli,                only : typ_color_coupling
  use def_coupli,                only : GLOBAL_NUMBERING
  use def_coupli,                only : SAME_COORDINATE
  use mod_parall,                only : par_memor
  use mod_parall,                only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                only : par_code_zone_subd_to_color
  use mod_parall,                only : PAR_COMM_MY_CODE
  use def_parall,                only : nneig
  use mod_maths,                 only : maths_heap_sort
  use mod_maths,                 only : maths_geometrical_sort_using_coordinates
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_copy
  use mod_memory,                only : memory_deallo
  use mod_memory,                only : memory_size
  use mod_renumbering,           only : renumbering_node_arrays
  use mod_coupling_memory,       only : cou_initialization
  use mod_coupling_memory,       only : cou_deallocate
  use mod_couplings,             only : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_meshes,                only : meshes_list_boundary_nodes
  use mod_graphs,                only : graphs_number_to_linked_list
  use mod_graphs,                only : graphs_dual_graph
  use mod_graphs,                only : graphs_dual_graph_deallocate
  use mod_graphs,                only : graphs_color_graph
  use mod_graphs,                only : graphs_color_graph_deallocate
  use mod_communications,        only : PAR_ALLGATHER
  use mod_communications,        only : PAR_ALLGATHERV
  use mod_communications,        only : PAR_SUM,PAR_MAX
  use mod_communications,        only : PAR_GATHER
  use mod_communications,        only : PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications,        only : PAR_COMM_RANK_AND_SIZE
  use mod_messages,              only : messages_live
  use def_mpi
#include "def_mpi.inc"
  implicit none
  private

  public :: par_interface_exchange
  public :: par_interface_exchange_scheduling

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-05-30
  !> @brief   Interface exchange communication arrays
  !> @details Generate communication data structure for independent
  !>          meshes, based on the the GLOBAL_NODE_NUMBERING option
  !>          of coupling module
  !>
  !>
  !-----------------------------------------------------------------------

  subroutine par_interface_exchange(comm,TYPE_OF_COUPLING,COMM_NAME)

    type(comm_data_par),     intent(inout), target   :: comm
    integer(ip),             intent(in),    optional :: TYPE_OF_COUPLING
    character(len=*),        intent(in),    optional :: COMM_NAME
    character(20)                                    :: my_comm_name
    type(typ_color_coupling)                         :: my_coupling
    integer(ip)                                      :: ii,kk,ineig,lsize,dom_i
    integer(ip)                                      :: kpoin,ipoin,nsubd
    integer(ip)                                      :: color_target
    integer(ip)                                      :: color_source
    integer(ip)                                      :: PAR_CURRENT_RANK
    integer(ip)                                      :: PAR_CURRENT_SIZE
    integer(ip), pointer                             :: node_type(:)
    integer(ip), pointer                             :: node_mult(:)
    integer(ip), pointer                             :: size_inte(:)
    integer(ip), pointer                             :: permr(:)
    logical(lg), pointer                             :: source_nodes(:)
    integer(ip), pointer                             :: npoin_par0(:)
    integer(ip), pointer                             :: nelem_par0(:)
    integer(ip), pointer                             :: nboun_par0(:)
    real(rp),    pointer                             :: xcoor(:,:)
    integer(ip), pointer                             :: lninv_sav(:)

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if

    call PAR_COMM_RANK_AND_SIZE(comm % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
 
    nsubd = PAR_CURRENT_SIZE-1
    !
    ! Nullify pointers
    !
    nullify(source_nodes)
    nullify(node_type)
    nullify(node_mult)
    nullify(size_inte)
    nullify(permr)
    nullify(npoin_par0)
    nullify(nelem_par0)
    nullify(nboun_par0)

    nullify(xcoor)
    nullify(lninv_sav)

    !--------------------------------------------------------------------
    !
    ! Periodicity
    !
    !--------------------------------------------------------------------
    !
    ! Put a negative number to enable coupling with myself
    !

    !--------------------------------------------------------------------
    !
    ! Create a coupling
    !
    !--------------------------------------------------------------------

    call cou_initialization(my_coupling)
    if( present(TYPE_OF_COUPLING) ) then
       my_coupling % itype               =  TYPE_OF_COUPLING
    else
       my_coupling % itype               =  GLOBAL_NUMBERING           ! Coupling based on global numbering
    end if
    my_coupling % kfl_multi_source       = -1                          ! Enable multiple source nodes
    my_coupling % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
    my_coupling % subdomain_target       =  0
    my_coupling % subdomain_source       =  0
    my_coupling % number                 =  1000
    my_coupling % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    my_coupling % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    my_coupling % kfl_symmetry           =  1                          ! Coupling is symmetric my wet nodes are also target
    my_coupling % commd % PAR_COMM_WORLD =  comm % PAR_COMM_WORLD
    color_target                         =  my_coupling % color_target
    color_source                         =  my_coupling % color_source
    !
    ! Define wet geometry: select boundary nodes only
    !
    call messages_live('COMPUTE LIST OF CANDIDATE INTERFACE NODES')
    if( INOTMASTER ) then
       call meshes_list_boundary_nodes(&
            nelem,npoin,lnods,ltype,&
            my_coupling % wet % npoin_wet,&
            my_coupling % wet % lpoin_wet,&
            LIST_BOUNDARY_NODE_NAME='COUPLING % WET % LPOIN_WET')
       call memory_alloca(par_memor,'COUPLING % WET % COORD_WET' ,'par_interface_exchange',my_coupling % wet % coord_wet,ndime,my_coupling % wet % npoin_wet)
       call memory_alloca(par_memor,'SOURCE_NODES'               ,'par_interface_exchange',source_nodes,npoin)
       !
       ! Restrict source nodes as target nodes!
       !
       do ii = 1,my_coupling % wet % npoin_wet
          ipoin = my_coupling % wet % lpoin_wet(ii)
          source_nodes(ipoin) = .true.
          my_coupling % wet % coord_wet(1:ndime,ii) = coord(1:ndime,ipoin)
       end do
    end if

!!$    block
!!$      use def_kintyp_mesh_basic, only: mesh_type_basic
!!$      use def_master
!!$      use def_domain
!!$      type(mesh_type_basic)           :: mesh
!!$      real(rp),             pointer   :: res(:)
!!$      character(len=5)                :: names
!!$      if(INOTMASTER) then 
!!$         nullify(res)
!!$         allocate(res(npoin))
!!$         res = 0.0_rp
!!$         names = 'WET  '
!!$         do ii = 1,my_coupling % wet % npoin_wet
!!$            ipoin = my_coupling % wet % lpoin_wet(ii)
!!$            res(ipoin) = 1.0_rp
!!$         end do
!!$         call mesh % init()
!!$         mesh % ndime = ndime
!!$         mesh % npoin = npoin
!!$         mesh % nelem = nelem
!!$         mesh % mnode = mnode
!!$         mesh % coord => coord
!!$         mesh % lnods => lnods
!!$         mesh % ltype => ltype
!!$         call mesh % output(FILENAME='test-'//trim(intost(kfl_paral)))
!!$         call mesh % results(res,names,FILENAME='test-'//trim(intost(kfl_paral)),where='ON NODES')
!!$      end if
!!$    end block
    !
    ! Coupling, it's better to take a minimum relative tolerance (coupling data has not been read!)
    !
    call messages_live('COMPUTE INTERFACE NODES WITH COUPLING')
    call COU_INIT_INTERPOLATE_POINTS_VALUES(&
         my_coupling,PAR_COMM_MY_CODE,CANDIDATE_SOURCE_NODES=source_nodes)
    call memory_deallo(par_memor,'SOURCE_NODES','par_interface_exchange',source_nodes)
    !
    ! We do not have LNINV_LOC... use distance to order nodes
    !
    if( my_coupling % itype == SAME_COORDINATE ) then
       call memory_alloca(par_memor,'XCOOR','par_interface_exchange',xcoor,ndime,my_coupling % commd % lsend_dim)
    end if

    call messages_live('COMPUTE INTERFACE NODE PROPERTY')
    if( INOTMASTER ) then
       !
       ! Allocate memory
       ! NODE_TYPE(II) = 0 ... interior
       ! NODE_TYPE(II) = 1 ... own boundary
       ! NODE_TYPE(II) = 2 ... boundary
       !
       call memory_alloca(par_memor,'NODE_TYPE','par_interface_exchange',node_type,npoin)
       call memory_alloca(par_memor,'NODE_MULT','par_interface_exchange',node_mult,npoin)
       call memory_alloca(par_memor,'SIZE_INTE','par_interface_exchange',size_inte,my_coupling % commd % nneig)
       !
       ! Order nodes on both sides of the interface
       !
       do ineig = 1,my_coupling % commd % nneig
          lsize            = my_coupling % commd % lsend_size(ineig+1)-my_coupling % commd % lsend_size(ineig)
          size_inte(ineig) = lsize
          kk               = my_coupling % commd % lsend_size(ineig)
          if( lsize > 0 ) then
             if( my_coupling % itype == SAME_COORDINATE ) then
                do ii =  1,lsize
                   ipoin = my_coupling % geome % lpoin_source(kk+ii-1)
                   xcoor(1:ndime,ii) = coord(1:ndime,ipoin)
                end do
                call maths_geometrical_sort_using_coordinates(2_ip,ndime,lsize,xcoor,my_coupling % geome % lpoin_source(kk:))
                do ii =  1,lsize
                   ipoin = my_coupling % geome % lpoin_source(kk+ii-1)
                end do
             else
                do ii =  1,lsize
                   ipoin = my_coupling % geome % lpoin_source(kk+ii-1)
                end do
                call maths_heap_sort(2_ip,lsize,my_coupling % geome % lpoin_source(kk:),PERMUTATION=lninv_loc)
             end if
             do kk = my_coupling % commd % lsend_size(ineig),my_coupling % commd % lsend_size(ineig+1)-1
                ipoin = my_coupling % geome % lpoin_source(kk)
                node_type(ipoin) = 2
                node_mult(ipoin) = node_mult(ipoin) + 1
             end do
          end if
       end do
       if( my_coupling % itype == SAME_COORDINATE ) then
          call memory_deallo(par_memor,'XCOOR','par_interface_exchange',xcoor)
       end if
       !
       ! Size of interface with neighbors, removing nodes with multiplicity > 2
       !
       do ineig = 1,my_coupling % commd % nneig
          do kk = my_coupling % commd % lsend_size(ineig),my_coupling % commd % lsend_size(ineig+1)-1
             ipoin = my_coupling % geome % lpoin_source(kk)
             if( node_mult(ipoin) > 1 ) size_inte(ineig) = size_inte(ineig) - 1
          end do
       end do
       !
       ! Split boundary between neighbors
       !
       do ineig = 1,my_coupling % commd % nneig
          dom_i = my_coupling % commd % neights(ineig)
          if( size_inte(ineig) > 0 ) then
             lsize = size_inte(ineig)/2
             if( PAR_CURRENT_RANK < dom_i ) then
                !
                ! Take the first chunk
                !
                do kk = my_coupling % commd % lsend_size(ineig),my_coupling % commd % lsend_size(ineig)+lsize
                   ipoin = my_coupling % geome % lpoin_source(kk)
                   if( node_mult(ipoin) <= 1 ) node_type(ipoin) = 1
                end do
             else
                !
                ! Take the second chunk
                !
                do kk = my_coupling % commd % lsend_size(ineig)+lsize+1,my_coupling % commd % lsend_size(ineig+1)-1
                   ipoin = my_coupling % geome % lpoin_source(kk)
                   if( node_mult(ipoin) <= 1 ) node_type(ipoin) = 1
                end do
             end if
          end if
       end do
       !
       ! Mark nodes with multiplicity > 2 (see node * in next example). It will be
       ! given to the lowest rank MPI process
       !
       !       o
       !       |
       !   1   o
       !       |    3
       !  -----*
       !       |
       !   2   o
       !
       do ipoin = 1,npoin
          if( node_mult(ipoin) > 1 ) then
             node_mult(ipoin) = PAR_CURRENT_RANK
          else
             node_mult(ipoin) = 0
          end if
       end do
       do ineig = 1,my_coupling % commd % nneig
          dom_i = my_coupling % commd % neights(ineig)
          do kk = my_coupling % commd % lsend_size(ineig),my_coupling % commd % lsend_size(ineig+1)-1
             ipoin = my_coupling % geome % lpoin_source(kk)
             if( node_mult(ipoin) > 0 ) node_mult(ipoin) = min(node_mult(ipoin),dom_i)
          end do
       end do
       do ipoin = 1,npoin
          if( node_mult(ipoin) > 0 ) then
             if( PAR_CURRENT_RANK == node_mult(ipoin) ) then
                node_type(ipoin) = 1
             else
                node_type(ipoin) = 2
             end if
          end if
       end do
       !
       ! Compute permutation to irder interior and boundary nodes
       !
       call memory_alloca(par_memor,'PERMR1','par_interface_exchange',permr,npoin)
       kpoin = 0
       do ipoin = 1,npoin
          if( node_type(ipoin) == 0 ) then
             kpoin = kpoin + 1
             permr(ipoin) = kpoin
          end if
       end do
       npoi1 = kpoin
       npoi2 = kpoin+1
       do ipoin = 1,npoin
          if( node_type(ipoin) == 1 ) then
             kpoin = kpoin + 1
             permr(ipoin) = kpoin
          end if
       end do
       npoi3     = kpoin
       npoin_own = npoi3
       do ipoin = 1,npoin
          if( node_type(ipoin) == 2 ) then
             kpoin = kpoin + 1
             permr(ipoin) = kpoin
          end if
       end do
       call renumbering_node_arrays(permR,PERMUTE_LMAST=.false.)
    end if
    !
    ! Permute remaining array... LPOIN_SOURCE
    !
    do kk = 1,memory_size(my_coupling % geome % lpoin_source)
       ipoin = my_coupling % geome % lpoin_source(kk)
       my_coupling % geome % lpoin_source(kk) = permR(ipoin)
    end do
    !
    ! Copy communication arrays
    !
    call messages_live('COMPUTE COMMUNICATION STRATEGY')
    call comm % deallo(MEMORY_COUNTER=par_memor,COMM_NAME='COMMD')


    if( INOTMASTER ) then

       comm % nneig     = my_coupling % commd % nneig
       comm % bound_dim = my_coupling % commd % lsend_size(my_coupling % commd % nneig+1)-1

       call comm % alloca(par_memor)
       
       ii = 0 
       do ineig = 1,my_coupling % commd % nneig
          dom_i = my_coupling % commd % neights(ineig)
          comm % bound_size(ineig) = my_coupling % commd % lsend_size(ineig)
          comm % neights(ineig)    = dom_i
          do kk = my_coupling % commd % lsend_size(ineig),my_coupling % commd % lsend_size(ineig+1)-1
             ii    = ii+1
             ipoin = my_coupling % geome % lpoin_source(kk)
             comm % bound_perm(ii) = ipoin
             comm % bound_invp(ii) = ipoin
          end do
       end do
       comm % bound_size(my_coupling % commd % nneig+1)=my_coupling % commd % lsend_size(my_coupling % commd % nneig+1)

       call memory_deallo(par_memor,'PERMR1'   ,'par_interface_exchange',permr    )
       call memory_deallo(par_memor,'NODE_TYPE','par_interface_exchange',node_type)
       call memory_deallo(par_memor,'NODE_MULT','par_interface_exchange',node_mult)
       call memory_deallo(par_memor,'SIZE_INTE','par_interface_exchange',size_inte)
    end if
    !
    ! Compute scheduling
    !
    call par_interface_exchange_scheduling(comm)
    !
    ! Deallocate
    !
    call cou_deallocate(my_coupling)

  end subroutine par_interface_exchange

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-05-30
  !> @brief   Compute scheduling
  !> @details Comptue communication scheduling for node interface exchange
  !>
  !-----------------------------------------------------------------------

  subroutine par_interface_exchange_scheduling(comm)

    type(comm_data_par), intent(inout)        :: comm
    integer(ip)                               :: ipart,ii,iz,nn
    integer(ip)                               :: total_neighbors
    integer(ip)                               :: nnDual
    integer(ip)                               :: nncolors,lsize,kneig
    integer(ip)                               :: i1_new,i2_new,i1,i2
    integer(ip)                               :: dom_i,ineig,ineig_new
    integer(ip)                               :: PAR_CURRENT_RANK
    integer(ip)                               :: PAR_CURRENT_SIZE
    MY_MPI_COMM                               :: PAR_COMM4
    integer(ip),         pointer              :: number_neighbors(:)
    integer(ip),         pointer              :: list_neighbors(:)
    integer(ip),         pointer              :: my_list_neighbors(:)
    integer(ip),         pointer              :: iaDom(:)
    integer(ip),         pointer              :: jaDom(:)
    integer(ip),         pointer              :: iaDomDual(:)
    integer(ip),         pointer              :: jaDomDual(:)
    integer(ip),         pointer              :: translDual(:)
    integer(ip),         pointer              :: colors(:)
    integer(ip),         pointer              :: lcomm(:)
    integer(ip),         pointer              :: bound_size(:)
    integer(ip),         pointer              :: neights(:)
    integer(ip),         pointer              :: bound_perm(:)
    integer(ip),         pointer              :: bound_invp(:)

    nullify(number_neighbors)
    nullify(list_neighbors)
    nullify(my_list_neighbors)
    nullify(iaDom)
    nullify(jaDom)
    nullify(iaDomDual)
    nullify(jaDomDual)
    nullify(translDual)
    nullify(colors)
    nullify(lcomm)
    nullify(bound_size)
    nullify(neights)
    nullify(bound_perm)
    nullify(bound_invp)
    !
    ! Communicator
    !
    !call PAR_COMM_RANK_AND_SIZE(comm % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    PAR_CURRENT_RANK = comm % RANK4
    PAR_CURRENT_SIZE = comm % SIZE4
    PAR_COMM4 = comm % PAR_COMM_WORLD
    nn        = PAR_CURRENT_SIZE-1
    kneig     = comm % nneig

    !--------------------------------------------------------------------
    !
    ! Neighbor graph: (IADOM,JADOM)
    !
    !--------------------------------------------------------------------
    !
    ! Gather neigboring information
    !
    call memory_alloca(par_memor,'NUMBER_NEIGHBORS','par_interface_exchange_scheduling',number_neighbors,nn+1_ip,lboun=0_ip)
    call PAR_ALLGATHER(comm % nneig,number_neighbors,PAR_COMM_IN=PAR_COMM4)
    total_neighbors = sum(number_neighbors)
    call memory_alloca(par_memor,'LIST_NEIGHBORS','par_interface_exchange_scheduling',list_neighbors,total_neighbors)
    if( INOTMASTER ) my_list_neighbors => comm % neights
    call PAR_ALLGATHERV(my_list_neighbors,list_neighbors,number_neighbors,PAR_COMM_IN=PAR_COMM4)

    if( INOTMASTER ) then
       !
       ! Compute subdomains graph (iaDom,jaDom)
       !
       call memory_alloca(par_memor,'IADOM','par_interface_exchange_scheduling',iaDom,nn+1_ip)
       call memory_alloca(par_memor,'JADOM','par_interface_exchange_scheduling',jaDom,total_neighbors)
       do ipart = 1,nn
          iaDom(ipart) = number_neighbors(ipart)
       end do
       call graphs_number_to_linked_list(nn,iaDom)
       ii = 0
       do ipart = 1,nn
          do iz = iaDom(ipart),iaDom(ipart+1)-1
             ii = ii + 1
             jaDom(ii) = list_neighbors(ii)
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Dual graph: (DUAL_IADOM,DUAL_JADOM)
       !
       !--------------------------------------------------------------------

       call graphs_dual_graph(nn,iaDom,jaDom,nnDual,iaDomDual,jaDomDual,translDual,par_memor)

       !--------------------------------------------------------------------
       !
       ! Color dual graph: NNCOLORS, COLORS
       !
       !--------------------------------------------------------------------

       call graphs_color_graph(nnDual,iaDomDual,jaDomDual,nncolors,colors,par_memor)

       !--------------------------------------------------------------------
       !
       ! Communication order
       !
       !--------------------------------------------------------------------

       call memory_alloca(par_memor,'LCOMM','par_interface_exchange_scheduling',lcomm,nncolors)
       call par_interface_exchange_order_communications(nncolors,nn,PAR_CURRENT_RANK,nnDual,iaDom,jaDom,translDual,colors,lcomm)

       !--------------------------------------------------------------------
       !
       ! Rearrange communication arrays to account for scheduling
       !
       !--------------------------------------------------------------------

       allocate( bound_size(kneig+1) )
       allocate( neights(kneig)      )
       allocate( bound_perm(comm % bound_dim))
       allocate( bound_invp(comm % bound_dim))
       !
       ! Compute new arrays
       !
       ineig_new = 0
       bound_size(1) = 1

       do ii = 1,nncolors
          dom_i = lcomm(ii)
          if( dom_i /= 0 ) then
             ineig = 1
             do while( comm % neights(ineig) /= dom_i )
                ineig = ineig + 1
             end do
             ineig_new                 = ineig_new + 1
             neights(ineig_new)        = dom_i
             lsize                     = comm % bound_size(ineig+1)-comm % bound_size(ineig)
             bound_size(ineig_new+1)   = bound_size(ineig_new) + lsize
             i1                        = comm % bound_size(ineig)
             i2                        = comm % bound_size(ineig+1)-1
             i1_new                    = bound_size(ineig_new)
             i2_new                    = bound_size(ineig_new+1)-1
             bound_perm(i1_new:i2_new) = comm % bound_perm(i1:i2)
             bound_invp(i1_new:i2_new) = comm % bound_invp(i1:i2)
          end if
       end do
       !
       ! Copy
       !
       do ineig = 1,comm % nneig
          comm % neights(ineig)    = neights(ineig)
          comm % bound_size(ineig) = bound_size(ineig)
       end do
       comm % bound_size(kneig+1) = bound_size(kneig+1)
       do ii = 1,comm % bound_dim
          comm % bound_perm(ii) = bound_perm(ii)
          comm % bound_invp(ii) = bound_invp(ii)
       end do
       !
       ! Deallocate
       !
       deallocate( bound_size )
       deallocate( neights    )
       deallocate( bound_perm )
       deallocate( bound_invp )

    end if
    
    call memory_deallo(par_memor,'NUMBER_NEIGHBORS' ,'par_interface_exchange_scheduling',number_neighbors)
    call memory_deallo(par_memor,'LIST_NEIGHBORS'   ,'par_interface_exchange_scheduling',list_neighbors)
    call memory_deallo(par_memor,'IADOM'            ,'par_interface_exchange_scheduling',iaDom)
    call memory_deallo(par_memor,'JADOM'            ,'par_interface_exchange_scheduling',jaDom)
    call memory_deallo(par_memor,'LCOMM'            ,'par_interface_exchange_scheduling',lcomm)

    call graphs_color_graph_deallocate(colors,MEMOR=par_memor)
    call graphs_dual_graph_deallocate (iaDomDual,JaDomDual,translDual,MEMOR=par_memor)

  end subroutine par_interface_exchange_scheduling

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    2018-06-01
  !> @brief   Compute the communication strategy
  !> @details Compute the communication strategy. Some variables
  !>          have been computed in previous subroutines:
  !>
  !>          LCOMM: Communication scheduling, organized by color
  !>          LCOMM(COLOR,SUBDOMAIN)= ADJACENT SUBDOMAIN
  !>
  !>          LCOMM( 1, 1 ) = 2  <= During the first communication stage, color 1,
  !>          LCOMM( 1, 2 ) = 1     we have: 1 <=> 2 and 3 <=> 4
  !>          LCOMM( 1, 3 ) = 4
  !>          LCOMM( 1, 4 ) = 3
  !>          LCOMM( 2, 1 ) = 4  <= During the 2nd communcation statge, color 3,
  !>          LCOMM( 2, 2 ) = 3     we have: 1 <=> 4 and 2 <=> 3
  !>          LCOMM( 2, 3 ) = 2
  !>          LCOMM( 2, 4 ) = 1
  !>          LCOMM( 3, 1 ) = 0  <= subdomain 1 is not involved in color 3
  !>          LCOMM( 3, 2 ) = 4  <= Stage 4, only 2 and 4 communicate
  !>          LCOMM( 3, 3 ) = 0  <= subdomain 3 is not involved in color 3
  !>          LCOMM( 3, 4 ) = 2
  !>
  !>          \endverbatim
  !>
  !------------------------------------------------------------------------

  subroutine par_interface_exchange_order_communications(nncolors,nsubd,isubd,nnDual,ia,ja,translDual,colors,lcomm)

    integer(ip),    intent(in)             :: nncolors            !< Number of colors
    integer(ip),    intent(in)             :: nsubd               !< Number of nodes of the dual graph
    integer(ip),    intent(in)             :: isubd               !< Number of nodes of the dual graph
    integer(ip),    intent(in)             :: nnDual              !< Number of nodes of the dual graph
    integer(ip),    intent(in),    pointer :: ia(:)               !< IA
    integer(ip),    intent(in),    pointer :: ja(:)               !< JA
    integer(ip),    intent(in),    pointer :: translDual(:)
    integer(ip),    intent(in),    pointer :: colors(:)
    integer(ip),    intent(inout), pointer :: lcomm(:)
    integer(ip)                            :: adjncy,ii,jj,vv,ww
    integer(ip),    pointer                :: invTransl(:,:)

    nullify(invTransl)
    if( .not. associated(lcomm) ) call memory_alloca(par_memor,'LCOMM'     ,'par_interface_exchange_order_communications',lcomm,nncolors)
    call memory_alloca(par_memor,'invTranslS','par_interface_exchange_order_communications',invTransl,2_ip,nnDual)
    !
    ! Construct inverse of TRANSL
    ! TRANSDUAL gives the dual graph node JJ (edge connecting adjacent edges)
    ! INVTRANSL(JJ) % NODE1/NODE2 gives for a dual graph node JJ the two adjacent
    ! subdomains NODE1 and NODE2 forming this communication edge
    !
    do vv = 1,nsubd
       do ii = ia(vv),ia(vv+1)-1
          ww = ja(ii)
          if( vv < ww ) then
             jj              = translDual(ii)
             invTransl(1,jj) = vv
             invTransl(2,jj) = ww
          end if
       end do
    end do
    !
    ! Construct communication array LCOMM
    ! Nodes of the dual graph of the same color can communicate at the same time
    ! so the strategy is the following:
    ! for each color, make adjacent subdomains of this color communicate. Each
    ! color consists of a stage of the communcation scheduling.
    ! Communications are automatically ordered.
    !
    do ii = 1,nncolors
       do adjncy = 1,nnDual
          if( colors(adjncy) == ii ) then
             vv = invTransl(1,adjncy)
             ww = invTransl(2,adjncy)
             !lcomm(ii,vv) = ww
             !lcomm(ii,ww) = vv
             if( vv == isubd ) lcomm(ii) = ww
             if( ww == isubd ) lcomm(ii) = vv
          endif
       end do
    end do

    call memory_deallo(par_memor,'invTranslS','par_interface_exchange_order_communications',invTransl)

  end subroutine par_interface_exchange_order_communications

end module mod_par_interface_exchange
!> @}
