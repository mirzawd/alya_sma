!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Mesh_Multiplication_Toolbox
!> @{
!> @name    ToolBox for mesh multiplication
!> @file    mod_mesh_multiplication.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for mesh multiplication
!> @details Idem
!
!-----------------------------------------------------------------------

module mod_mesh_multiplication

  use def_kintyp
  use def_domain
  use def_master
  use def_kermod
  use def_master,            only : lfacg,nfacg
  use def_elmtyp,            only : element_max
  use def_elmgeo,            only : element_type
  use def_elmgeo,            only : list_edges_HEX27
  use def_elmgeo,            only : list_edges_TET10
  use def_elmgeo,            only : QUA09_TO_QUA04
  use def_elmtyp,            only : BAR02,BAR03,QUA09,TRI06
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_SUM
  use mod_messages,          only : messages_live
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_outfor,            only : outfor
  use mod_messages,          only : messages_live
  use mod_messages,          only : livinf
  use mod_elmgeo,            only : elmgeo_mm_elements
  use mod_elmgeo,            only : elmgeo_mm_interior_nodes
  use mod_elmgeo,            only : elmgeo_mm_face_nodes
  use mod_elmgeo,            only : elmgeo_number_nodes
  use mod_elmgeo,            only : elmgeo_edge
  use mod_elmgeo,            only : elmgeo_segment
  use mod_elmgeo,            only : elmgeo_shapf_deriv_heslo
  use mod_domain,            only : domain_memory_allocate
  use mod_domain,            only : domain_memory_deallocate
  use mod_domain,            only : domain_memory_reallocate
  use mod_parall,            only : par_memor
  use def_elmtyp,            only : element_end
  use def_elmtyp,            only : element_max
  use mod_strings,           only : integer_to_string
  use mod_maths_arrays,      only : maths_findloc
  use mod_maths_sort,        only : maths_heap_sort

#ifdef __PGI
#define MEMPGI )
#else
#define MEMPGI ,memor=memor_dom)
#endif
  implicit none

  private

  integer(ip)            :: howdi(-element_max:element_max)   ! How to divide elements according to type
  integer(ip)            :: newno(-element_max:element_max)   ! New nodes according to type  
  real(rp)               :: shapf(64)
  integer(ip)            :: nfaco(30)
  integer(ip)            :: mnodm
  integer(ip)            :: face_node(21,6)
  integer(ip), pointer   :: pedgp(:)                          ! List of boundary faces (mesh subdivision)
  logical(lg)            :: if_high_order
  type(i1p),   pointer   :: face_divide(:,:)
  integer(ip), parameter :: subQUA09(4,8) = reshape ( [ 1,5,9,8, 5,9,0,0, 5,2,6,9, 8,9,0,0, &
       &                                                9,6,0,0, 8,9,7,4, 9,7,0,0, 9,6,3,7  ], [4,8] )
  integer(ip), parameter :: numQUA09(8)   =           [ 4,2,4,2,2,4,2,4 ]
  integer(ip), parameter :: subTRI06(2,3) = reshape ( [ 6,4, 4,5, 6,5 ], [2,3] )
  integer(ip), parameter :: numTRI06(3)   =           [ 2,2,2 ]
  
  public      :: mesh_multiplication
  public      :: mesh_multiplication_node_codes

contains

  subroutine mesh_multiplication()

    !-----------------------------------------------------------------------
    !****f* domain/submsh
    ! NAME
    !    submsh
    ! DESCRIPTION
    !    This subroutine recursively subdivides the mesh
    ! OUTPUT
    ! USED BY
    !    submsh
    !***
    !-----------------------------------------------------------------------

    integer(ip) :: idivi,kelem,ielem,ipoin,kpoin,iboun,kboun
    real(rp)    :: time1,time2,time3,time4

    if( ndivi > 0 ) then

       if( multiply_with_curvature == 0 ) then
          call messages_live('MESH MULTIPLICATION','START SECTION')
       else
          call messages_live('CURVED MESH MULTIPLICATION','START SECTION')
       end if
       !
       ! Subdivide recursively the mesh
       !
       call par_memory(18_ip)                               ! Parall: Allocate some memory
       call mesh_multiplication_memory(1_ip)                ! Allocate memory

       do idivi = 1,ndivi
          call messages_live('LEVEL '//trim(intost(idivi)),'START SECTION')
          call cputim(time1)
          call mesh_multiplication_list_edges()             ! Edge table
          call cputim(time2)
          call mesh_multiplication_list_faces()             ! Face table
          call cputim(time3)
          call mesh_multiplication_divide_elements()        ! Subdivide elements
          call cputim(time4)
          if( kfl_mmpar == 0 .and. kfl_elint == 0 .and. kfl_elndi == 0 ) then
             call par_submsh()
          else
             call mesh_multiplication_parallel_interfaces()
          end if
          call mesh_multiplication_output_info()            ! New mesh dimensions
          call mesh_multiplication_memory(2_ip)             ! Deallocate memory
          call domvar(2_ip)                                 ! LNUTY...
          call par_checks()                                 ! Parall: Perform some checkings

          call messages_live('LEVEL '//trim(intost(idivi)),'END SECTION')
          cpu_other(1) = cpu_other(1) + time2 - time1
          cpu_other(2) = cpu_other(2) + time3 - time2
          cpu_other(3) = cpu_other(3) + time4 - time3
       end do

       call messages_live('MESH MULTIPLICATION','END SECTION')

    else

       call mesh_multiplication_memory(1_ip)           ! Allocate memory

    end if
    !
    ! Save levels
    !
    if( INOTMASTER .and. ndivi > 0 .and. kfl_posdi == 0 ) then
       call memory_alloca(memor_dom,'LPMSH','mesh_multiplication',lpmsh,meshe(0) % npoin)
       do ipoin = 1,npoin
          kpoin = lnlev(ipoin)
          if( kpoin > 0 ) lpmsh(kpoin) = ipoin
       end do
       call memory_alloca(memor_dom,'LEMSH','mesh_multiplication',lemsh,meshe(0) % nelem)
       do ielem = 1,nelem
          kelem = lelev(ielem)
          if( kelem > 0 ) lemsh(kelem) = ielem
       end do
       call memory_alloca(memor_dom,'LBMSH','mesh_multiplication',lbmsh,meshe(0) % nboun)
       do iboun = 1,nboun
          kboun = lblev(iboun)
          if( kboun > 0 ) lbmsh(kboun) = iboun
       end do
    end if
    !
    ! Output CPU info
    !
    call mesh_multiplication_cpu()

  end subroutine mesh_multiplication

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-06-06
  !> @brief   Reconstruct interface
  !> @details Reconstruct interface using coupling
  !>          WARNING: This technique is incompatible with coupling, as
  !>          if nodes coincide at interface, Alya will try to glue the
  !>          nodes!
  !>
  !-----------------------------------------------------------------------

  subroutine mesh_multiplication_parallel_interfaces()

    use def_coupli,                 only : SAME_COORDINATE
    use mod_par_interface_exchange, only : par_interface_exchange
    use mod_par_global_numbering,   only : par_global_numbering_nodes
    use mod_par_global_numbering,   only : par_global_numbering_elements_boundaries
    use def_domain,                 only : npoin,nelem,nboun
    use def_master,                 only : lninv_loc,leinv_loc,lbinv_loc
    use mod_parall,                 only : PAR_COMM_MY_CODE_ARRAY
    use mod_par_additional_arrays,  only : par_global_variables_arrays
    use mod_par_additional_arrays,  only : par_ordered_exchange_update

    if( IPARALL ) then

       call PAR_COMM_MY_CODE_ARRAY(1) % deallo()
       if(      kfl_mmpar == 2 ) then
          call par_interface_exchange(PAR_COMM_MY_CODE_ARRAY(1))
       else if( kfl_mmpar == 1 ) then
          call par_interface_exchange(PAR_COMM_MY_CODE_ARRAY(1),TYPE_OF_COUPLING=SAME_COORDINATE)
       end if
       !
       ! Global numbering
       !
       call par_global_numbering_nodes              (npoin,lninv_loc)
       call par_global_numbering_elements_boundaries(nelem,leinv_loc)
       call par_global_numbering_elements_boundaries(nboun,lbinv_loc)
       !
       ! Parallelization variables
       !
       call par_global_variables_arrays()
       !
       ! Mesh dimensions
       !
       call mesh_multiplication_dimensions()

    end if

  end subroutine mesh_multiplication_parallel_interfaces

  !-----------------------------------------------------------------------
  !
  !> @brief   Remove nodes from a graph
  !> @details Remove nodes from a graph using a mask
  !>
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine mesh_multiplication_cpu()

    real(rp) :: cpu_multi

    if( ndivi > 0 .and. .not. PART_AND_WRITE() ) then

       call PAR_MAX(9_ip,cpu_other,'IN MY CODE')

       cpu_multi  =   cpu_other(1) + cpu_other(2) + cpu_other(3) &
            &       + cpu_other(4) + cpu_other(5) + cpu_other(6) &
            &       + cpu_other(7) + cpu_other(8) + cpu_other(9)
       cpu_multi  = max(zeror,cpu_multi)
       routp( 1)  = cpu_multi
       routp( 2)  = cpu_other(1)
       routp( 3)  = 100.0_rp*routp( 2)/cpu_multi
       routp( 4)  = cpu_other(2)
       routp( 5)  = 100.0_rp*routp( 4)/cpu_multi
       routp( 6)  = cpu_other(3)
       routp( 7)  = 100.0_rp*routp( 6)/cpu_multi
       routp( 8)  = cpu_other(4)
       routp( 9)  = 100.0_rp*routp( 8)/cpu_multi
       routp(10)  = cpu_other(5)
       routp(11)  = 100.0_rp*routp(10)/cpu_multi
       routp(12)  = cpu_other(6)
       routp(13)  = 100.0_rp*routp(12)/cpu_multi
       routp(14)  = cpu_other(7)
       routp(15)  = 100.0_rp*routp(14)/cpu_multi
       routp(16)  = cpu_other(8)
       routp(17)  = 100.0_rp*routp(16)/cpu_multi
       routp(18)  = cpu_other(9)
       routp(19)  = 100.0_rp*routp(18)/cpu_multi
       call outfor(49_ip,lun_outpu,' ')
    end if

  end subroutine mesh_multiplication_cpu

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-08-03
  !> @brief   Create edge table
  !> @details Create edge table
  !>          NNEDG ... Number of edges
  !>          LEDGG ... Edge table
  !>          LEDGB ... Boundary edge table (when Parall is on)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_multiplication_list_edges()

    use def_kintyp
    use def_parame
    use def_domain
    use def_master
    use mod_graphs, only : graphs_elepoi
    use mod_graphs, only : graphs_elepoi_deallocate
    use mod_memory, only : memory_alloca
    use mod_memory, only : memory_deallo
    use mod_graphs, only : graphs_number_to_linked_list
    implicit none
    integer(ip)          :: ipoin,iedgg,ielem
    integer(ip)          :: ilisn,jpoin,dummi,pelty,pblty,lnode(2)
    integer(ip)          :: iedge,ielty,jedgg,kedgg,jlisn,isegm
    integer(ip), pointer :: nedpo(:)
    logical(lg)          :: found
    integer(ip), pointer :: lnnod_loc(:)
    integer(ip), pointer :: ledgp(:,:) 
  
    nullify(lnnod_loc)
    nullify(ledgp)
    nullify(pedgp)

    call messages_live('EDGE TABLE')

    if( INOTMASTER ) then
       !
       ! Deallocate node-element graph if needed
       !
       call graphs_elepoi_deallocate(pelpo,lelpo,memor=memor_dom)
       !
       ! Compute node-element graph
       !
       if( .not. associated(lnnod) ) then
          call memory_alloca(memor_dom,'LNNOD_LOC','par_renumber_nodes',lnnod_loc,nelem)
          do ielem = 1,nelem
             pelty = abs(ltype(ielem))
             if( pelty < element_end) & 
                  lnnod_loc(ielem) = element_type(pelty) % number_nodes
          end do
       else
          lnnod_loc => lnnod
       end if
       call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod_loc,dummi,pelpo,lelpo MEMPGI
       !
       !
       nullify(nedpo)
       call memory_alloca(memor_dom,'PEDGP','ledges',pedgp,npoin+1)
       call memory_alloca(memor_dom,'NEDPO','ledges',nedpo,npoin)
       do ielem = 1,nelem
          ielty = ltype(ielem)
          do iedge = 1,element_type(ielty) % number_edges
             pblty = element_type(ielty) % type_edges(iedge)

             select case ( pblty )
                
             case ( 0_ip )
                
                ipoin = lnods(element_type(ielty) % list_edges(2,iedge),ielem)
                jpoin = lnods(element_type(ielty) % list_edges(1,iedge),ielem)
                if( ipoin < jpoin ) then
                   nedpo(ipoin) = nedpo(ipoin) + 1
                else if( ipoin > jpoin ) then
                   nedpo(jpoin) = nedpo(jpoin) + 1
                end if

             case default
                
                do isegm = 1,element_type(pblty) % order
                   lnode = elmgeo_segment(ielty,iedge,isegm)
                   ipoin = lnods(lnode(1),ielem)
                   jpoin = lnods(lnode(2),ielem)
                   if( ipoin < jpoin ) then
                      nedpo(ipoin) = nedpo(ipoin) + 1
                   else if( ipoin > jpoin ) then
                      nedpo(jpoin) = nedpo(jpoin) + 1
                   end if
                end do
                
             end select
             
          end do
       end do
       !
       ! Allocate memory for LELPO and construct the list
       !
       pedgp(1) = 1
       do ipoin = 1,npoin
          pedgp(ipoin+1) = pedgp(ipoin) + nedpo(ipoin)
       end do
       nedgg = pedgp(npoin+1)-1
       call memory_alloca(memor_dom,'ledgp','elepoi',ledgp,2_ip,nedgg)
       do ielem = 1,nelem
          ielty = ltype(ielem)
          do iedge = 1,element_type(ielty) % number_edges
             pblty = element_type(ielty) % type_edges(iedge)
             
             select case ( pblty )
                
             case ( 0_ip )

                ipoin = lnods(element_type(ielty) % list_edges(2,iedge),ielem)
                jpoin = lnods(element_type(ielty) % list_edges(1,iedge),ielem)
                if( ipoin < jpoin ) then
                   ledgp(1,pedgp(ipoin)) = jpoin
                   ledgp(2,pedgp(ipoin)) = ielem
                   pedgp(ipoin)          = pedgp(ipoin)+1
                else if( ipoin > jpoin ) then
                   ledgp(1,pedgp(jpoin)) = ipoin
                   ledgp(2,pedgp(jpoin)) = ielem
                   pedgp(jpoin)          = pedgp(jpoin)+1
                end if
                
             case default
                
                do isegm = 1,element_type(pblty) % order
                   lnode = elmgeo_segment(ielty,iedge,isegm)
                   ipoin = lnods(lnode(1),ielem)
                   jpoin = lnods(lnode(2),ielem)
                   if( ipoin < jpoin ) then
                      ledgp(1,pedgp(ipoin)) = jpoin
                      ledgp(2,pedgp(ipoin)) = ielem
                      pedgp(ipoin)          = pedgp(ipoin)+1
                   else if( ipoin > jpoin ) then
                      ledgp(1,pedgp(jpoin)) = ipoin
                      ledgp(2,pedgp(jpoin)) = ielem
                      pedgp(jpoin)          = pedgp(jpoin)+1
                   end if
                end do
                
             end select
          end do
       end do
       !
       ! Fill in edge table:
       ! LEDGG(1,IEDGE) = EDGE NODE 1
       ! LEDGG(2,IEDGE) = EDGE NODE 2
       ! LEDGG(3,IEDGE) = NEW NODE
       ! LEDGG(4,IEDGE) = 0/1 interface edge in parallel
       ! LEDGG(5,IEDGE) = IELEM
       !
       call memory_alloca(memor_dom,'LEDGG','ledges',ledgg,5_ip,nedgg)
       
       iedgg = 0
       nedgg = 0
       do ipoin = 1,npoin
          pedgp(ipoin) = 0
       end do

       do ipoin = 1,npoin 
          kedgg = iedgg
          do ilisn = 1,nedpo(ipoin)
             iedgg = iedgg + 1
             jpoin = ledgp(1,iedgg)     
             jlisn = 0
             jedgg = kedgg
             found = .false.
             edge: do while( jlisn < ilisn-1 )
                jlisn = jlisn + 1
                jedgg = jedgg + 1
                if( jpoin == ledgp(1,jedgg) ) then
                   found = .true.
                   exit edge
                end if
             end do edge
             if( .not. found ) then
                nedgg          = nedgg + 1
                ledgg(1,nedgg) = ipoin
                ledgg(2,nedgg) = jpoin                
                ledgg(5,nedgg) = ledgp(2,iedgg)                         
                pedgp(ipoin)   = pedgp(ipoin)+1
             end if
          end do
       end do
       call graphs_number_to_linked_list(npoin,pedgp)
       call memory_deallo(memor_dom,'NEDPO','ledges',nedpo)
       !
       ! Deallocate
       !
       if( .not. associated(lnnod) ) then
          call memory_deallo(memor_dom,'LNNOD_LOC','par_renumber_nodes',lnnod_loc)
       end if
       !
       ! Deallocate memory
       ! PEDGP, LELPO, PELPO are deallocated at the end of each multiplicaiton
       ! in subelm
       !
       call memory_deallo(memor_dom,'LEDGP','ledges',ledgp)

    end if

  end subroutine mesh_multiplication_list_edges

  subroutine mesh_multiplication_divide_elements()
    !-----------------------------------------------------------------------
    !****f* domain/subelm
    ! NAME
    !    domain
    ! DESCRIPTION
    !    Subdivide elements
    !    Nodes are numbered as follows:
    !    1 ................... NPOIN_OLD:       Old nodes
    !    NPOIN_OLD+1 ......... NPOIN_OLD+NEDDG: Edge nodes
    !    NPOIN_OLD+NEDDG+1 ... NPOIN:           Central nodes
    ! OUTPUT
    ! USED BY
    !    submsh
    !***
    !-----------------------------------------------------------------------
    use def_kintyp
    use def_parame
    use def_elmtyp
    use def_domain
    use def_master
    use def_kermod
    use mod_memory
    use mod_mpio_config,           only : mpio_config
    use mod_curved_multiplication, only : curve_subdivided_element
    implicit none
    integer(ip)          :: ielem,inode,ipoin,idime,pnode
    integer(ip)          :: ipoi1,ipoi2,ipoi3,ipoi4,istep
    integer(ip)          :: node1,node2,node3,node4,neldi,nbodi,npoif,iface
    integer(ip)          :: iedgg,kelem,pelty,ipoif,ipoic,pblty
    integer(ip)          :: iboun,inodb,kboun,ii,npoic,jpoin,pflty,ifacl
    integer(ip)          :: lnodx(100),igrou,jgrou,ifacg,ifiel,iedge,isegm
    integer(ip)          :: code1,code2,code3,code4,mcod1
    integer(ip)          :: knode,knodb,pnodb,kpoin
    integer(ip)          :: lninv(4),lnode(20),kk
    real(rp)             :: dummr
    !
    ! PERMUTATIONS
    !
    integer(ip), pointer :: perm_elem(:)       ! Element permutation
    integer(ip), pointer :: perm_boun(:)       ! Boundary permutation
    !
    ! GEOMETRY
    !
    integer(ip), pointer :: lnods_old(:,:)     ! NELEM
    integer(ip), pointer :: ltype_old(:)       ! NELEM
    integer(ip), pointer :: lelch_old(:)       ! NELEM
    integer(ip), pointer :: lnnod_old(:)       ! NELEM
    integer(ip), pointer :: lesub_old(:)       ! NELEM
    integer(ip), pointer :: lnodb_old(:,:)     ! NBOUN
    integer(ip), pointer :: lelbo_old(:)       ! NBOUN
    integer(ip), pointer :: ltypb_old(:)       ! NBOUN
    integer(ip), pointer :: lboch_old(:)       ! NBOUN
    real(rp),    pointer :: coord_old(:,:)     ! NPOIN
    integer(ip), pointer :: lgrou_dom_old(:)   ! NPOIN
    integer(ip), pointer :: lnoch_old(:)       ! NPOIN
    integer(ip), pointer :: lmast_old(:)       ! NPOIN
    !
    ! MATERIALS
    !
    integer(ip), pointer :: lmate_old(:)       ! NELEM
    !
    ! SETS
    !
    integer(ip), pointer :: leset_old(:)       ! NELEM
    integer(ip), pointer :: lbset_old(:)       ! NBOUN
    integer(ip), pointer :: lnset_old(:)       ! NPOIN
    !
    ! BOUNDARY CONDITIONS
    !
    integer(ip), pointer :: kfl_codno_old(:,:) ! NPOIN
    integer(ip), pointer :: kfl_codbo_old(:)   ! NBOUN
    !
    ! FIELDS
    !
    type(r3p),   pointer :: xfiel_old(:)       ! NPOIN/NBOUN
    !
    ! INTERPOLATION
    !
    type(i1pp),  pointer :: linno(:)           ! NPOIN
    !
    ! Permutation with original mesh
    !
    integer(ip), pointer :: lnlev_old(:)       ! NPOIN
    integer(ip), pointer :: lelev_old(:)       ! NELEM
    integer(ip), pointer :: lblev_old(:)       ! NBOUN
    logical(lg)          :: if_lnnod
    !
    ! Edge/face division
    !
    logical(lg), pointer :: divide_edge(:)     ! Flag for edge division
    logical(lg), pointer :: divide_face(:)     ! Flag for face division
    logical(lg), pointer :: divide_boun(:)     ! Flag for boundary division
    !
    ! Characteristic elements
    !
    integer(ip)          :: howdi_ELINT        ! Number elements division for ELINT
    integer(ip)          :: howbo_ELINT        ! Number boundaries division for ELINT
    integer(ip)          :: nelem_ELINT        ! Number of ELINT element types
    integer(ip)          :: nboun_ELINT        ! Number of boundaries belonging for ELINT element types
    integer(ip)          :: nedgf_ELINT        ! Number of edges from ELINT elements without division
    integer(ip)          :: nfacf_ELINT        ! Number of faces from ELINT elemnts without division
    !
    ! Boundary ordering for QUA04
    !
    integer(ip), parameter :: list_inodb_QUA04(4,4) = reshape((/1,2,3,4,2,3,4,1,3,4,1,2,4,1,2,3/),(/4,4/)) ! Boundary ordering cases
    !
    !-|CURVED| Mesh division with curvature variables
    !   real(rp), pointer :: curved_geometry(:,:,:) ! Elemental field with the curved geometry corresponding to each element
    !   real(rp), pointer :: subdivided_curved_geometry(:,:,:) ! Elemental field with the curved geometry corresponding to the subdivided mesh
    real(rp), allocatable :: curved_geometry(:,:,:) ! Elemental field with the curved geometry corresponding to each element
    real(rp), allocatable :: subdivided_curved_geometry(:,:,:)  ! Elemental field with the curved geometry corresponding to the subdivided mesh
    real(rp), allocatable :: subdivided_coordinates(:,:)        ! Curved coordinates of the subdivided elements
    real(rp), allocatable :: subdivided_element_geometry(:,:,:) ! Curved coordinates of the subdivided elements
    integer(ip)           :: geometry_order(nelem)              ! Interpolation order of the geometry for each element
    integer(ip)           :: geometry_numNodes(nelem)           ! Number of geometry nodes for each element
    integer(ip)           :: geometry_dimension(nelem)          ! Spatial dimension of the geometry for each element
    integer(ip)           :: geometry_fieldComponents(nelem)    ! Number of elemental fields for the geometry for each element
    integer(ip)           :: maxNumNodesGeometry
    integer(ip)           :: maxSubdividedCoordinates
    integer(ip)           :: maxSubdividedElements
    integer(ip)           :: beginDimCurvature
    integer(ip)           :: endDimCurvature
    integer(ip)           :: numFieldCurvaturePerDimension
    integer(ip)           :: numNewCoordinates, numNewElements
    integer(ip)           :: maxTotalNumberOfDividedElements
    integer(ip)           :: nsteps

    !--------------------------------------------------------------------
    !
    ! Nullify pointers
    !
    !--------------------------------------------------------------------

    if( associated(lnnod) ) then
       if_lnnod = .true.
    else
       if_lnnod = .false.
    end if

    nullify(perm_elem)
    nullify(perm_boun)

    nullify(lnods_old)
    nullify(ltype_old)
    nullify(lelch_old)
    nullify(lnnod_old)
    nullify(lesub_old)
    nullify(lnodb_old)
    nullify(lelbo_old)
    nullify(ltypb_old)
    nullify(lboch_old)
    nullify(coord_old)
    nullify(lgrou_dom_old)
    nullify(lnoch_old)
    nullify(lmast_old)
    nullify(xfiel_old)
    nullify(lmate_old)
    nullify(leset_old)
    nullify(lbset_old)
    nullify(lnset_old)
    nullify(kfl_codno_old)
    nullify(kfl_codbo_old)
    nullify(linno)
    nullify(lnlev_old)
    nullify(lblev_old)
    nullify(lelev_old)

    nullify(divide_edge)
    nullify(divide_face)
    nullify(divide_boun)

    !--------------------------------------------------------------------
    !
    ! Check possible errors
    !
    !--------------------------------------------------------------------

    if( lexis(PYR05) > 0 .and. lexis(TET04) == 0 ) then
       call runend('MESH MULTIPLICATION: YOU HAVE PYRAMIDS, AND HAVE '//&
            'TO EXPLICITLY DECLARE THE PRESENCE OF TET04 ELEMENTS IN DIMENSIONS FIELD')
    end if

    if( (kfl_elint == 1 .and. kfl_mmpar == 0 .and. kfl_elndi == 0) .or. &
         (kfl_elint == 1 .and. kfl_mmpar == 0 .and. kfl_elndi == 1) .or. &
         (kfl_elint == 0 .and. kfl_mmpar == 0 .and. kfl_elndi == 1) ) then
       call runend('MESH MULTIPLICATION: WRONG PARALLELIZATION METHOD FOR INTERFACE ELEMENTS'//&
            'OR FOR ELEMENT NORMAL DIVISION')
    end if

    !--------------------------------------------------------------------
    !
    ! Edge/face division
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then

       call memory_alloca(memor_dom,'DIVIDE_EDGE','mesh_multiplication_divide_elements',divide_edge,nedgg)
       call memory_alloca(memor_dom,'DIVIDE_BOUN','mesh_multiplication_divide_elements',divide_boun,nboun)
       call memory_alloca(memor_dom,'DIVIDE_FACE','mesh_multiplication_divide_elements',divide_face,nfacg)
       if( nfacg > 0 ) divide_face(1:nfacg) = .true.
       if( nedgg > 0 ) divide_edge(1:nedgg) = .true.
       if( nboun > 0 ) divide_boun(1:nboun) = .true.

    end if

    !--------------------------------------------------------------------
    !
    ! How elements are divided
    !
    !--------------------------------------------------------------------

    do pelty = -element_max,element_max 
       howdi(pelty) = elmgeo_mm_elements(pelty)
    end do
    do pelty = -element_max,element_max 
       newno(pelty) = elmgeo_mm_interior_nodes(pelty)
    end do

    !--------------------------------------------------------------------
    !
    ! How characteristic elements are divided
    !
    !--------------------------------------------------------------------

    howdi_ELINT  = 0
    howbo_ELINT  = 0
    if( kfl_elint == 1 .or. kfl_elndi == 1 ) then
       if( ndime == 2 ) then
          howdi_ELINT  = 2
          howbo_ELINT  = 1
       else
          howdi_ELINT  = 4
          howbo_ELINT  = 2
       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Special treatment for characteristic elements / element normal division
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then

       if( kfl_elint == 1 .or. kfl_elndi == 1 ) then
          !
          ! Interface elements
          !
          call mesh_multiplication_divide_ELINT(divide_edge,divide_face,divide_boun)

       end if

    end if

    !--------------------------------------------------------------------
    !
    ! Copy old mesh data
    !
    !--------------------------------------------------------------------

    call messages_live('ADD NEW NODES AND ELEMENTS')

    if( INOTMASTER ) then

       nelem_old = nelem
       nboun_old = nboun
       npoin_old = npoin

       call memory_alloca(memor_dom,'LNODS_OLD','mesh_multiplication_divide_elements',lnods_old,mnode,nelem,  'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LTYPE_OLD','mesh_multiplication_divide_elements',ltype_old,nelem,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LELCH_OLD','mesh_multiplication_divide_elements',lelch_old,nelem,        'DO_NOT_INITIALIZE')
       if( if_lnnod ) call memory_alloca(memor_dom,'LNNOD_OLD','mesh_multiplication_divide_elements',lnnod_old,nelem,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LESUB_OLD','mesh_multiplication_divide_elements',lesub_old,nelem,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LMATE_OLD','mesh_multiplication_divide_elements',lmate_old,nelem,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LNODB_OLD','mesh_multiplication_divide_elements',lnodb_old,mnodb,nboun,  'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LELBO_OLD','mesh_multiplication_divide_elements',lelbo_old,nboun,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LTYPB_OLD','mesh_multiplication_divide_elements',ltypb_old,nboun,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LBOCH_OLD','mesh_multiplication_divide_elements',lboch_old,nboun,        'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'COORD_OLD','mesh_multiplication_divide_elements',coord_old,ndime,npoin,  'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LNOCH_OLD','mesh_multiplication_divide_elements',lnoch_old,npoin)
       call memory_alloca(memor_dom,'LMAST_OLD','mesh_multiplication_divide_elements',lmast_old,npoin)

       do ielem = 1,nelem
          ltype_old(ielem) = ltype(ielem)
          lelch_old(ielem) = lelch(ielem)
          lesub_old(ielem) = lesub(ielem)
          lmate_old(ielem) = lmate(ielem)
          do inode = 1,mnode
             lnods_old(inode,ielem) = lnods(inode,ielem)
          end do
       end do
       if( if_lnnod ) then
          do ielem = 1,nelem
             lnnod_old(ielem) = lnnod(ielem)
          end do
       end if
       do iboun = 1,nboun
          do inodb = 1,mnodb
             lnodb_old(inodb,iboun) = lnodb(inodb,iboun)
          end do
          lelbo_old(iboun)         = lelbo(iboun)
          ltypb_old(iboun)         = ltypb(iboun)
          lboch_old(iboun)         = lboch(iboun)
       end do
       do ipoin = 1,npoin
          lnoch_old(ipoin) = lnoch(ipoin)
          lmast_old(ipoin) = lmast(ipoin)
          do idime = 1,ndime
             coord_old(idime,ipoin) = coord(idime,ipoin)
          end do
       end do
       !
       ! Characteristic elements counters
       !
       nelem_ELINT = 0
       nedgf_ELINT = 0
       nfacf_ELINT = 0
       nboun_ELINT = 0
       if( kfl_elint == 1 .and. kfl_elndi == 0 ) then
          nelem_ELINT = count(lelch_old(1:nelem_old)==ELINT)
          nedgf_ELINT = count(divide_edge(1:nedgg) .eqv. .false.)
          if( associated(divide_face) ) nfacf_ELINT = count(divide_face(1:nfacg) .eqv. .false.)
          do iboun = 1,nboun
             ielem = lelbo_old(iboun)
             if( lelch_old(ielem) == ELINT ) then
                nboun_ELINT = nboun_ELINT + 1
             end if
          end do
       else if( kfl_elndi == 1 ) then
          nelem_ELINT = nelem_old
          nedgf_ELINT = count(divide_edge(1:nedgg) .eqv. .false.)
          if( associated(divide_face) ) nfacf_ELINT = count(divide_face(1:nfacg) .eqv. .false.)
          if( associated(divide_boun) ) nboun_ELINT = count(divide_boun(1:nboun) .eqv. .false.)
       end if
       !
       ! Allocate memory fo new mesh
       !
       call domain_memory_deallocate('GEOMETRY')
       if( if_lnnod ) call domain_memory_deallocate('LNNOD')

       select case ( ndime )

       case ( 1_ip ) ! 1D

          nelem = ( lnuty(BAR02) + lnuty(BAR03) + lnuty(BAR04) ) * 2
          nboun = nboun_old
          npoin = npoin_old + lnuty(BAR02)*newno(BAR02) + lnuty(BAR03)*newno(BAR03) + lnuty(BAR04)*newno(BAR04)  

       case ( 2_ip ) ! 2D

          nelem = ( lnuty(TRI03) + lnuty(QUA04) + lnuty(TRI06) + lnuty(QUA09) + lnuty(QUA16) ) * 4 
          nboun = nboun_old * 2
          npoin = npoin_old + nedgg &
               + elmgeo_mm_interior_nodes(QUA04) * lnuty(QUA04) &
               + elmgeo_mm_interior_nodes(TRI06) * lnuty(TRI06) &
               + elmgeo_mm_interior_nodes(QUA09) * lnuty(QUA09) &
               + elmgeo_mm_interior_nodes(QUA16) * lnuty(QUA16)

       case ( 3_ip ) ! 3D

          nelem = ( lnuty(TET04) + lnuty(PEN06) + lnuty(HEX08) + lnuty(HEX27) + lnuty(TET10) ) * 8 &
               + elmgeo_mm_elements(PYR05) * lnuty(PYR05) &
               + elmgeo_mm_elements(BAR3D) * lnuty(BAR3D) &
               + lnuty(-DDDNE)
          nboun = ( lnuty(TRI03) + lnuty(QUA04) + lnuty(TRI06) + lnuty(QUA09) ) * 4 + lnuty(POINT) 
          npoin = npoin_old + nedgg &
               + elmgeo_mm_face_nodes    (TRI03) * nfaco(TRI03) &
               + elmgeo_mm_face_nodes    (TRI06) * nfaco(TRI06) &
               + elmgeo_mm_face_nodes    (QUA04) * nfaco(QUA04) &
               + elmgeo_mm_face_nodes    (QUA09) * nfaco(QUA09) &
               + elmgeo_mm_interior_nodes(HEX08) * lnuty(HEX08) &
               + elmgeo_mm_interior_nodes(HEX27) * lnuty(HEX27) &
               + elmgeo_mm_interior_nodes(TET10) * lnuty(TET10) 

       end select

       if( kfl_elint == 1 .or. kfl_elndi == 1_ip ) then
          if( ndime == 2 ) then
             nelem = nelem - nelem_ELINT * 2
             nboun = nboun - nboun_ELINT
             npoin = npoin - nelem_ELINT - nedgf_ELINT
          else
             nelem = nelem - nelem_ELINT * 4
             nboun = nboun - nboun_ELINT * 2
             npoin = npoin - nelem_ELINT - nedgf_ELINT - nfacf_ELINT
          end if
       end if

       !--------------------------------------------------------------------
       !
       ! Local permutation arrays NEW <= OLD
       !
       !--------------------------------------------------------------------

       call memory_alloca(memor_dom,'PERM_ELEM','mesh_multiplication_divide_elements',perm_elem,nelem)
       call memory_alloca(memor_dom,'PERM_BOUN','mesh_multiplication_divide_elements',perm_boun,nboun)
       kelem = 0
       do ielem = 1,nelem_old
          neldi = howdi(abs(ltype_old(ielem)))
          if( lelch_old(ielem) == ELINT .or. kfl_elndi == 1_ip ) neldi = howdi_ELINT
          do ii = 1,neldi
             kelem = kelem + 1
             perm_elem(kelem) = ielem
          end do
       end do
       kboun = 0
       do iboun = 1,nboun_old
          nbodi = howdi(abs(ltypb_old(iboun)))
          if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
          do ii = 1,nbodi
             kboun = kboun + 1
             perm_boun(kboun) = iboun
          end do
       end do

       !--------------------------------------------------------------------
       !
       ! Permutation arrays
       !
       !--------------------------------------------------------------------

       call memory_alloca(memor_dom,'LNLEV_OLD','mesh_multiplication_divide_elements',lnlev_old,npoin,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LELEV_OLD','mesh_multiplication_divide_elements',lelev_old,nelem,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'LBLEV_OLD','mesh_multiplication_divide_elements',lblev_old,nboun,'DO_NOT_INITIALIZE')

       do ipoin = 1,npoin_old
          lnlev_old(ipoin) = lnlev(ipoin)
       end do
       do ielem = 1,nelem_old
          lelev_old(ielem) = lelev(ielem)
       end do
       do iboun = 1,nboun_old
          lblev_old(iboun) = lblev(iboun)
       end do

       call mesh_multiplication_memory(-1_ip)
       call mesh_multiplication_memory( 1_ip)

       do ipoin = 1,npoin_old
          lnlev(ipoin) = lnlev_old(ipoin)
       end do

       kelem = 0
       do ielem = 1,nelem_old
          neldi = howdi(abs(ltype_old(ielem)))
          if( lelch_old(ielem) == ELINT .or. kfl_elndi == 1 ) neldi = howdi_ELINT
          do ii = 1,neldi
             kelem = kelem + 1
             if( ii == 1 ) then
                lelev(kelem) = lelev_old(ielem)
             else
                lelev(kelem) = 0
             end if
          end do
       end do

       kboun = 0
       do iboun = 1,nboun_old
          nbodi = howdi(abs(ltypb_old(iboun)))
          if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
          do ii = 1, nbodi
             kboun = kboun + 1
             if( ii == 1 ) then
                lblev(kboun) = lblev_old(iboun)
             else
                lblev(kboun) = 0
             end if
          end do
       end do

       if( ISEQUEN ) then
          !
          ! OJO: should be taken off when renumbering sequential run
          !
          do ipoin = npoin_old+1,npoin
             lnlev(ipoin) = 0
          end do
          do ielem = nelem_old+1,nelem
             lelev(ielem) = 0
          end do
          do iboun = nboun_old+1,nboun
             lblev(iboun) = 0
          end do
       end if

       call memory_deallo(memor_dom,'LBLEV_OLD','mesh_multiplication_divide_elements',lblev_old)
       call memory_deallo(memor_dom,'LELEV_OLD','mesh_multiplication_divide_elements',lelev_old)
       call memory_deallo(memor_dom,'LNLEV_OLD','mesh_multiplication_divide_elements',lnlev_old)

       !--------------------------------------------------------------------
       !
       !-|CURVED| Read curved geometry field (if exists)
       !
       !--------------------------------------------------------------------

       if(multiply_with_curvature.eq.1) then
          geometry_order           = int(xfiel(curvatureDataField)%a(1,:,1),ip)
          geometry_numNodes        = int(xfiel(curvatureDataField)%a(2,:,1),ip)
          geometry_dimension       = int(xfiel(curvatureDataField)%a(3,:,1),ip)
          geometry_fieldComponents = geometry_numNodes*geometry_dimension

          maxNumNodesGeometry = maxval(geometry_numNodes,size(geometry_numNodes,1,kind=ip))

          if(lexis(TRI03).eq.1_ip) then
             maxSubdividedCoordinates = 3_ip
             maxSubdividedElements = 4_ip
          end if
          if(lexis(QUA04).eq.1_ip) then
             maxSubdividedCoordinates = 5_ip
             maxSubdividedElements = 4_ip
          end if
          if(lexis(TET04).eq.1_ip) then
             maxSubdividedCoordinates = 6_ip
             maxSubdividedElements = 8_ip
          end if
          if(lexis(HEX08).eq.1_ip) then
             maxSubdividedCoordinates = 19_ip
             maxSubdividedElements = 8_ip
          end if

          !        if(ndime.ne.maxval(geometry_dimension,nelem)) then !if((ndime.ne.maxval(geometry_dimension,nelem)).or.(ndime.ne.minval(geometry_dimension,nelem))) then
          !          call runend('SUBMSH: Wrong curve geometry dimension')
          !        end if

          maxTotalNumberOfDividedElements = nelem!=> this is exact!!!! not as ...nelem_old*maxSubdividedElements
          allocate(curved_geometry(maxNumNodesGeometry,ndime,nelem_old))
          allocate(subdivided_curved_geometry(maxNumNodesGeometry,ndime,maxTotalNumberOfDividedElements))
          !        call memory_alloca(memor_dom,'curved_geometry','curved_geometry',maxNumNodesGeometry,ndime,nelem,'DO_NOT_INITIALIZE')
          !        call memory_alloca(memor_dom,'subdivided_curved_geometry','subdivided_curved_geometry',maxNumNodesGeometry,ndime,nelem*maxSubdividedElements,'DO_NOT_INITIALIZE')

          beginDimCurvature = 1_ip
          numFieldCurvaturePerDimension = maxNumNodesGeometry ! size(xfiel(curvatureField)%a,1,KIND=ip)/ndime
          do idime = 1_ip,ndime
             endDimCurvature = beginDimCurvature+numFieldCurvaturePerDimension-1_ip
             curved_geometry(:,idime,:) = xfiel(curvatureField)%a(beginDimCurvature:endDimCurvature,:,1)
             beginDimCurvature = endDimCurvature+1_ip
          end do

       end if

       !--------------------------------------------------------------------
       !
       ! New mesh: do not know what to do with LMAST and LNOCH
       ! Nodes are ordered with LNINV_LOC to have the same results on
       ! parallel interfaces (to avoid round-off errors)
       !
       !--------------------------------------------------------------------

       call domain_memory_allocate('GEOMETRY')
       if( if_lnnod ) call domain_memory_allocate('LNNOD')
       !
       ! NPOIN arrays
       !
       ! Old nodes
       do ipoin = 1,npoin_old
          lnoch(ipoin) = lnoch_old(ipoin)
          lmast(ipoin) = lmast_old(ipoin)
          do idime = 1,ndime
             coord(idime,ipoin) = coord_old(idime,ipoin)
          end do
       end do
       !
       ! Edge nodes
       !
       ipoin = npoin_old
       do iedgg = 1,nedgg
          if( divide_edge(iedgg) ) then
             ipoin          = ipoin + 1
             node1          = ledgg(1,iedgg)
             node2          = ledgg(2,iedgg)
             ledgg(3,iedgg) = ipoin
             lnoch(ipoin)   = NOFEM
             lmast(ipoin)   = 0

             if( if_high_order ) then
                ielem = ledgg(5,iedgg)
                pelty = abs(ltype_old(ielem))
                pnode = element_type(pelty) % number_nodes
                node1 = maths_findloc(pnode,lnods_old(:,ielem),node1) 
                node2 = maths_findloc(pnode,lnods_old(:,ielem),node2)
                call coord_new_point([node1,node2],ielem,pelty,lnods_old,coord_old,coord(1:ndime,ipoin))                                   
             else
                if( lninv_loc(node1) < lninv_loc(node2) ) then
                   node3 = node1
                   node1 = node2
                   node2 = node3
                end if
                coord(1:ndime,ipoin) = 0.5_rp * ( coord(1:ndime,node1) + coord(1:ndime,node2) )
             end if
          end if
       end do
       !
       ! Face nodes
       !
       if( if_high_order ) then
          call memory_alloca(memor_dom,'FACE_DIVIDE','renpoi',face_divide,8_ip,nfacg)
          do ifacg = 1,nfacg
             do ii = 1,8
                call memory_alloca(memor_dom,'FACE_DIVIDE % L','renpoi',face_divide(ii,ifacg) % l,4_ip)
             end do
          end do
       end if

       npoif = ipoin
       do ifacg = 1,nfacg
          if( divide_face(ifacg) ) then
             pflty = lfacg(mnodm+2,ifacg)
             ielem = lfacg(mnodm+3,ifacg)
             ifacl = lfacg(mnodm+4,ifacg)
             pelty = abs(ltype_old(ielem))

             select case ( pflty ) 

             case ( QUA04 )
                !
                ! QUA04
                !
                ipoin                = ipoin + 1
                lfacg(mnodm+1,ifacg) = ipoin             
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0           
                lnode(1:4)           = (/lfacg(1,ifacg),lfacg(2,ifacg),lfacg(3,ifacg),lfacg(4,ifacg)/)
                lninv(1:4)           = (/lninv_loc(lnode(1)),lninv_loc(lnode(2)),lninv_loc(lnode(3)),lninv_loc(lnode(4))/)
                do ii = 1,3
                   kk = minloc(lninv(ii:),1) + ii - 1
                   if( lninv(ii) > lninv(kk) ) then
                      node1     = lninv(ii)
                      lninv(ii) = lninv(kk)
                      lninv(kk) = node1
                      node1     = lnode(ii)
                      lnode(ii) = lnode(kk)
                      lnode(kk) = node1   
                   end if
                end do
                coord(1:ndime,ipoin) = 0.25_rp * ( coord(1:ndime,lnode(1)) + coord(1:ndime,lnode(2)) &
                     &                           + coord(1:ndime,lnode(3)) + coord(1:ndime,lnode(4)) )

             case ( TRI06 )
                !
                ! TRI06
                !
                do ii = 1,6
                   lnode(ii)         = element_type(pelty) % list_faces(ii,ifacl)
                end do

                lfacg(mnodm+1,ifacg) = ipoin + 1

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(1,ifacg) % l(1:2) = lnods_old([lnode(6),lnode(4)],ielem)
                call coord_new_point([lnode(6),lnode(4)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(2,ifacg) % l(1:2) = lnods_old([lnode(4),lnode(5)],ielem)
                call coord_new_point([lnode(4),lnode(5)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(3,ifacg) % l(1:2) = lnods_old([lnode(6),lnode(5)],ielem)
                call coord_new_point([lnode(6),lnode(5)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                !kk = 0
                !print*,'----------------',ifacl,ifacg
                !do kpoin = ipoin-2,ipoin
                !   kk = kk + 1
                !   write(*,'(a,i1,a,i1,a,3(1x,f4.1),a,i4)') 'face_node=(',kk,',',ifacl,')',coord(1:ndime,kpoin),': ',kpoin
                !end do

             case ( QUA09 )
                !
                ! QUA09: Keep only the first point for each face... the other ones are consecutively numbered
                !
                do ii = 1,9
                   lnode(ii)         = element_type(pelty) % list_faces(ii,ifacl)
                end do

                lfacg(mnodm+1,ifacg) = ipoin + 1

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(1,ifacg) % l(1:4) = lnods_old([lnode(1),lnode(5),lnode(9),lnode(8)],ielem)
                call coord_new_point([lnode(1),lnode(5),lnode(9),lnode(8)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(2,ifacg) % l(1:2) = lnods_old([lnode(5),lnode(9)],ielem)
                call coord_new_point([lnode(5),lnode(9)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(3,ifacg) % l(1:4) = lnods_old([lnode(5),lnode(2),lnode(6),lnode(9)],ielem)
                call coord_new_point([lnode(5),lnode(2),lnode(6),lnode(9)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(4,ifacg) % l(1:2) = lnods_old([lnode(8),lnode(9)],ielem)
                call coord_new_point([lnode(8),lnode(9)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(5,ifacg) % l(1:2) = lnods_old([lnode(9),lnode(6)],ielem)
                call coord_new_point([lnode(9),lnode(6)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0               
                face_divide(6,ifacg) % l(1:4) = lnods_old([lnode(8),lnode(9),lnode(7),lnode(4)],ielem)                
                call coord_new_point([lnode(8),lnode(9),lnode(7),lnode(4)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(7,ifacg) % l(1:2) = lnods_old([lnode(7),lnode(9)],ielem)                
                call coord_new_point([lnode(7),lnode(9)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))

                ipoin                = ipoin + 1
                lnoch(ipoin)         = NOFEM
                lmast(ipoin)         = 0
                face_divide(8,ifacg) % l(1:4) = lnods_old([lnode(9),lnode(6),lnode(3),lnode(7)],ielem)                
                call coord_new_point([lnode(9),lnode(6),lnode(3),lnode(7)],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))            

             end select
          end if
       end do
       ! stop
       !
       ! Center nodes: these are nodes located inside the elements, that do not belong to faces netiher edges
       !
       npoic = ipoin
       do ielem = 1,nelem_old

          pelty = abs(ltype_old(ielem))
          if( lelch_old(ielem) /= ELINT .and. kfl_elndi == 0 ) then

             select case ( pelty )

             case ( BAR02 )
                !
                ! BAR02
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,2_ip],ielem,pelty,lnods_old,coord_old,coord(1:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( BAR03 )
                !
                ! BAR03
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,3_ip],ielem,pelty,lnods_old,coord_old,coord(:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2
                ipoin            = ipoin + 1
                call coord_new_point([3_ip,2_ip],ielem,pelty,lnods_old,coord_old,coord(:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( BAR04 )
                !
                ! BAR04
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,3_ip],ielem,pelty,lnods_old,coord_old,coord(:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2
                ipoin            = ipoin + 1
                call coord_new_point([3_ip,4_ip],ielem,pelty,lnods_old,coord_old,coord(:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 3
                ipoin            = ipoin + 1
                call coord_new_point([4_ip,2_ip],ielem,pelty,lnods_old,coord_old,coord(:,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( TRI06 )
                !
                ! TRI06: 3 new nodes
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([6_ip,4_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2
                ipoin            = ipoin + 1
                call coord_new_point([4_ip,5_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 3
                ipoin            = ipoin + 1
                call coord_new_point([5_ip,6_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( QUA04 )
                !
                ! QUA04
                !
                ipoin = ipoin + 1

                coord(1,ipoin) = 0.0_rp
                coord(2,ipoin) = 0.0_rp
                do inode = 1,4
                   jpoin = lnods_old(inode,ielem)
                   coord(1,ipoin) = coord(1,ipoin) + coord_old(1,jpoin)
                   coord(2,ipoin) = coord(2,ipoin) + coord_old(2,jpoin)
                end do
                coord(1,ipoin) = 0.25_rp * coord(1,ipoin)
                coord(2,ipoin) = 0.25_rp * coord(2,ipoin)

                lnoch(ipoin) = NOFEM
                lmast(ipoin) = 0

             case ( QUA09 )
                !
                ! QUA09: 8 new nodes
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,5_ip,9_ip,8_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))            
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2
                ipoin            = ipoin + 1
                call coord_new_point([2_ip,5_ip,6_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin)) 
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 3
                ipoin            = ipoin + 1
                call coord_new_point([9_ip,6_ip,3_ip,7_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin)) 
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 4
                ipoin            = ipoin + 1
                call coord_new_point([8_ip,9_ip,7_ip,4_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin)) 
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 5
                ipoin            = ipoin + 1
                call coord_new_point([5_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 6
                ipoin            = ipoin + 1
                call coord_new_point([9_ip,6_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin)) 
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 7
                ipoin            = ipoin + 1
                call coord_new_point([7_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin)) 
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 8
                ipoin            = ipoin + 1
                call coord_new_point([8_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( QUA16 )
                !
                ! QUA16: 21 new nodes
                !
                ! 1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,5_ip,13_ip,12_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2
                ipoin            = ipoin + 1
                call coord_new_point([5_ip,13_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 3
                ipoin            = ipoin + 1
                call coord_new_point([5_ip,6_ip,14_ip,13_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 4
                ipoin            = ipoin + 1
                call coord_new_point([6_ip,14_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 5
                ipoin            = ipoin + 1
                call coord_new_point([6_ip,2_ip,7_ip,14_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 6
                ipoin            = ipoin + 1
                call coord_new_point([12_ip,13_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 7
                ipoin            = ipoin + 1
                call coord_new_point([13_ip,14_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 8
                ipoin            = ipoin + 1
                call coord_new_point([14_ip,7_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 9
                ipoin            = ipoin + 1
                call coord_new_point([12_ip,13_ip,16_ip,11_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 10
                ipoin            = ipoin + 1
                call coord_new_point([13_ip,16_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                lmast(ipoin)     = 0
                ! 11
                ipoin            = ipoin + 1
                call coord_new_point([13_ip,14_ip,15_ip,16_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 12
                ipoin            = ipoin + 1
                call coord_new_point([14_ip,15_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 13
                ipoin            = ipoin + 1
                call coord_new_point([14_ip,7_ip,8_ip,15_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 14
                ipoin            = ipoin + 1
                call coord_new_point([11_ip,16_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 15
                ipoin            = ipoin + 1
                call coord_new_point([16_ip,15_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 16
                ipoin            = ipoin + 1
                call coord_new_point([15_ip,8_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 17
                ipoin            = ipoin + 1
                call coord_new_point([11_ip,16_ip,10_ip,4_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 18
                ipoin            = ipoin + 1
                call coord_new_point([16_ip,10_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 19
                ipoin            = ipoin + 1
                call coord_new_point([16_ip,15_ip,9_ip,10_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 20
                ipoin            = ipoin + 1
                call coord_new_point([15_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 21
                ipoin            = ipoin + 1
                call coord_new_point([15_ip,8_ip,3_ip,9_ip],ielem,pelty,lnods_old,coord_old,coord(1:2,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( TET10 )
                !
                ! HEX27: 24 new nodes
                !
                ! 6-8
                ipoin            = ipoin + 1
                call coord_new_point([6_ip,8_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             case ( HEX08 )
                !
                ! HEX08
                !
                ipoin = ipoin + 1
                coord(1,ipoin) = 0.0_rp
                coord(2,ipoin) = 0.0_rp
                coord(3,ipoin) = 0.0_rp
                do inode = 1,8
                   jpoin          = lnods_old(inode,ielem)
                   coord(1,ipoin) = coord(1,ipoin) + coord_old(1,jpoin)
                   coord(2,ipoin) = coord(2,ipoin) + coord_old(2,jpoin)
                   coord(3,ipoin) = coord(3,ipoin) + coord_old(3,jpoin)
                end do
                coord(1,ipoin) = 0.125_rp * coord(1,ipoin)
                coord(2,ipoin) = 0.125_rp * coord(2,ipoin)
                coord(3,ipoin) = 0.125_rp * coord(3,ipoin)

                lnoch(ipoin) = NOFEM
                lmast(ipoin) = 0

             case ( HEX27 )
                !
                ! HEX27: 24 new nodes
                !
                ! 1, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([9_ip,21_ip,27_ip,22_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 2, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([21_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 3, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([22_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 4, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([12_ip,21_ip,27_ip,25_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 5, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,25_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 6, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([13_ip,22_ip,27_ip,25_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 7, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([1_ip,9_ip,22_ip,13_ip,12_ip,21_ip,27_ip,25_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 8, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,23_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 9 y=-1
                ipoin            = ipoin + 1
                call coord_new_point([21_ip,10_ip,23_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 10, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([22_ip,14_ip,23_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 11, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([9_ip,2_ip,10_ip,21_ip,22_ip,14_ip,23_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 12, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,24_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 13, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,23_ip,15_ip,24_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 14, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([21_ip,10_ip,3_ip,11_ip,27_ip,23_ip,15_ip,24_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 15, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([25_ip,27_ip,24_ip,16_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 16, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([25_ip,27_ip,24_ip,16_ip,12_ip,21_ip,11_ip,4_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 17, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([26_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 18, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([22_ip,27_ip,26_ip,17_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 19, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([25_ip,27_ip,26_ip,20_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 20, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([13_ip,22_ip,27_ip,25_ip,5_ip,17_ip,26_ip,20_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 21, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,23_ip,18_ip,26_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 22, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([22_ip,14_ip,23_ip,27_ip,17_ip,6_ip,18_ip,26_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 23, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,24_ip,19_ip,26_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                lmast(ipoin)     = 0
                ! 24, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([25_ip,27_ip,24_ip,16_ip,20_ip,26_ip,19_ip,8_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 25, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([27_ip,23_ip,15_ip,24_ip,26_ip,18_ip,7_ip,19_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0
                ! 26, y=-1
                ipoin            = ipoin + 1
                call coord_new_point([21_ip,11_ip,24_ip,27_ip],ielem,pelty,lnods_old,coord_old,coord(1:3,ipoin))
                lnoch(ipoin)     = NOFEM
                lmast(ipoin)     = 0

             end select
          end if
       end do
       !
       ! NELEM arrays
       !
       kelem = 0
       ipoic = npoic
       ipoif = npoif

       if( ndime == 1 ) then
          do ielem = 1,nelem_old

             pelty = abs(ltype_old(ielem))
             select case ( pelty )
                
             case ( BAR02 ) 
                !
                ! 1           2 
                ! o-----*-----o
                !       3
                !

                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = ipoic + 1
                ltype(kelem)   = BAR02
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = ipoic + 1
                lnods(2,kelem) = lnods_old(2,ielem)
                ltype(kelem)   = BAR02
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + newno(BAR02)

             case ( BAR03 ) 
                !
                ! 1     3     2
                ! o--*--o--*--o
                !    4     5
                !
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnods_old(3,ielem)
                lnods(3,kelem) = ipoic + 1
                ltype(kelem)   = BAR03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(3,ielem)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = ipoic + 2
                ltype(kelem)   = BAR03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + newno(BAR03)

             case ( BAR04 )
                !
                ! 1     3     4     2
                ! o--*--o--*--o--*--o
                !    5     6     7
                !
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = ipoic +2
                lnods(3,kelem) = ipoic + 1
                lnods(4,kelem) = lnods_old(3,ielem)
                ltype(kelem)   = BAR04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = ipoic + 2
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnods_old(4,ielem)
                lnods(4,kelem) = ipoic + 3
                ltype(kelem)   = BAR04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + newno(BAR04)

             end select
          end do

       else if( ndime == 2 ) then
          do ielem = 1,nelem_old

             pelty =  abs(ltype_old(ielem))
             select case ( pelty )
                
             case ( TRI03 ) 
                !
                !        3
                !        o
                !  (3)  /4\  (2)
                !      *---*
                !    / 1\2/3 \
                !   o----*----o
                !   1   (1)   2
                !
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx(2)) ! Edge 2-3
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(3,ielem),lnodx(3)) ! Edge 1-3

                !-|CURVED| Curved mesh division for TRI03
                if( multiply_with_curvature.eq.1_ip ) then
                   if(geometry_order(ielem).gt.1_ip) then
                      numNewCoordinates = 3_ip
                      numNewElements    = 4_ip

                      allocate( subdivided_coordinates(numNewCoordinates,ndime) )
                      allocate( subdivided_element_geometry(geometry_numNodes(ielem),ndime,numNewElements) )
                      call curve_subdivided_element(&
                           TRI03,&
                           numNewCoordinates,&
                           numNewElements,&
                           geometry_order(ielem),&
                           geometry_numNodes(ielem),&
                           ndime,&
                           curved_geometry(1:geometry_numNodes(ielem),1:ndime,ielem),&
                           subdivided_coordinates,&
                           subdivided_element_geometry)

                      coord(1:ndime,lnodx(1:3)) = transpose(subdivided_coordinates)
                      subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+4)) = subdivided_element_geometry

                      deallocate(subdivided_coordinates)
                      deallocate(subdivided_element_geometry)
                   else
                      subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+4)) = 0.0_rp
                   end if
                end if
                !
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnodx(1)
                lnods(3,kelem) = lnodx(3)
                ltype(kelem)   = TRI03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnodx(2)
                lnods(3,kelem) = lnodx(3)
                ltype(kelem)   = TRI03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnodx(2)
                ltype(kelem)   = TRI03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(3)
                lnods(2,kelem) = lnodx(2)
                lnods(3,kelem) = lnods_old(3,ielem)
                ltype(kelem)   = TRI03
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                if( lelch_old(ielem) == ELEXT ) then
                   lelch(kelem)   = ELEXT
                   lelch(kelem-1) = ELFEM
                   lelch(kelem-2) = ELFEM
                   lelch(kelem-2) = ELFEM
                end if

             case ( TRI06 ) 
                !
                ! TRI06
                !
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(4,ielem),lnodx( 1)) ! Edge 1-4
                call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(2,ielem),lnodx( 2)) ! Edge 4-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(5,ielem),lnodx( 3)) ! Edge 2-5
                call mesh_multiplication_edge_node(lnods_old(5,ielem),lnods_old(3,ielem),lnodx( 4)) ! Edge 5-3
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(6,ielem),lnodx( 5)) ! Edge 3-6
                call mesh_multiplication_edge_node(lnods_old(6,ielem),lnods_old(1,ielem),lnodx( 6)) ! Edge 6-1
                ! 1
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnods_old(4,ielem)
                lnods(3,kelem) = lnods_old(6,ielem)
                lnods(4,kelem) = lnodx( 1)
                lnods(5,kelem) = ipoic + 1
                lnods(6,kelem) = lnodx( 6)

                ltype(kelem)   = TRI06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(4,ielem)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnods_old(5,ielem)
                lnods(4,kelem) = lnodx( 2)
                lnods(5,kelem) = lnodx( 3)
                lnods(6,kelem) = ipoic + 2

                ltype(kelem)   = TRI06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(6,ielem)
                lnods(2,kelem) = lnods_old(5,ielem)
                lnods(3,kelem) = lnods_old(3,ielem)
                lnods(4,kelem) = ipoic + 3
                lnods(5,kelem) = lnodx( 4)
                lnods(6,kelem) = lnodx( 5)

                ltype(kelem)   = TRI06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(4,ielem)
                lnods(2,kelem) = lnods_old(5,ielem)
                lnods(3,kelem) = lnods_old(6,ielem)
                lnods(4,kelem) = ipoic + 2
                lnods(5,kelem) = ipoic + 3
                lnods(6,kelem) = ipoic + 1

                ltype(kelem)   = TRI06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + 3

             case ( QUA04 ) 
                !
                ! QUA04
                !
                if( (lelch_old(ielem) == ELFEM .or. lelch_old(ielem) == ELEXT) .and. kfl_elndi == 0  ) then
                   !
                   !   4           3
                   !   o-----*-----o
                   !   |     |     |
                   !   |     |     |
                   !   *-----*-----*
                   !   |     |     |
                   !   |     |     |
                   !   o-----*-----o
                   !   1           2
                   !
                   call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx(2)) ! Edge 2-3
                   call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx(3)) ! Edge 3-4
                   call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(1,ielem),lnodx(4)) ! Edge 4-1

                   ipoic          = ipoic + 1
                   !-|CURVED| Curved mesh division for QUA04
                   if( multiply_with_curvature.eq.1_ip ) then
                      if(geometry_order(ielem).gt.1_ip) then

                         numNewCoordinates = 5_ip
                         numNewElements    = 4_ip

                         allocate( subdivided_coordinates(numNewCoordinates,ndime) )
                         allocate( subdivided_element_geometry(geometry_numNodes(ielem),ndime,numNewElements) )

                         call curve_subdivided_element(&
                              QUA04,&
                              numNewCoordinates,&
                              numNewElements,&
                              geometry_order(ielem),&
                              geometry_numNodes(ielem),&
                              ndime,&
                              curved_geometry(1:geometry_numNodes(ielem),1:ndime,ielem),&
                              subdivided_coordinates,&
                              subdivided_element_geometry)
                         coord(1:ndime,(/lnodx(1:4),ipoic/)) = transpose(subdivided_coordinates)
                         subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+4)) = subdivided_element_geometry

                         deallocate(subdivided_coordinates)
                         deallocate(subdivided_element_geometry)
                      else
                         subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+4)) = 0.0_rp
                      end if
                   end if
                   ! 1
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnods_old(1,ielem)
                   lnods(2,kelem) = lnodx(1)
                   lnods(3,kelem) = ipoic
                   lnods(4,kelem) = lnodx(4)
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 2
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(1)
                   lnods(2,kelem) = lnods_old(2,ielem)
                   lnods(3,kelem) = lnodx(2)
                   lnods(4,kelem) = ipoic
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 3
                   kelem          = kelem + 1
                   lnods(1,kelem) = ipoic
                   lnods(2,kelem) = lnodx(2)
                   lnods(3,kelem) = lnods_old(3,ielem)
                   lnods(4,kelem) = lnodx(3)
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 4
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(4)
                   lnods(2,kelem) = ipoic
                   lnods(3,kelem) = lnodx(3)
                   lnods(4,kelem) = lnods_old(4,ielem)
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)

                   if( lelch_old(ielem) == ELEXT ) then
                      lelch(kelem)   = ELHOL
                      lelch(kelem-1) = ELHOL
                      lelch(kelem-2) = ELHOL
                   end if

                else if ( lelch_old(ielem) == ELINT .or. (lelch_old(ielem) == ELFEM .and. kfl_elndi == 1) ) then
                   !
                   !   4           3
                   !   o-----*-----o
                   !   |     |     |
                   !   |     |     |
                   !   |     |     |
                   !   |     |     |
                   !   |     |     |
                   !   o-----*-----o
                   !   1           2
                   !
                   call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx(2)) ! Edge 3-4
                   !-|CURVED| Curved mesh division for QUA04 - ELINT
                   if ( multiply_with_curvature .eq. 1_ip ) call runend('SUBMSH: not implemented curved division for QUA04-ELINT')

                   ! 1
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnods_old(1,ielem)
                   lnods(2,kelem) = lnodx(1)
                   lnods(3,kelem) = lnodx(2)
                   lnods(4,kelem) = lnods_old(4,ielem)
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 2
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(1)
                   lnods(2,kelem) = lnods_old(2,ielem)
                   lnods(3,kelem) = lnods_old(3,ielem)
                   lnods(4,kelem) = lnodx(2)
                   ltype(kelem)   = QUA04
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                end if

             case ( QUA09 ) 
                !
                ! QUA09
                !
                !   4     7     3
                !   o-----*-----o
                !   |     |     |
                !   |     |     |
                ! 8 *-----9-----* 6
                !   |     |     |
                !   |     |     |
                !   o-----*-----o
                !   1     5     2
                !
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(5,ielem),lnodx( 1)) ! Edge 1-5
                call mesh_multiplication_edge_node(lnods_old(5,ielem),lnods_old(2,ielem),lnodx( 2)) ! Edge 5-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(6,ielem),lnodx( 3)) ! Edge 2-6
                call mesh_multiplication_edge_node(lnods_old(6,ielem),lnods_old(3,ielem),lnodx( 4)) ! Edge 6-3
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(7,ielem),lnodx( 5)) ! Edge 3-7
                call mesh_multiplication_edge_node(lnods_old(7,ielem),lnods_old(4,ielem),lnodx( 6)) ! Edge 7-4
                call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(8,ielem),lnodx( 7)) ! Edge 4-8
                call mesh_multiplication_edge_node(lnods_old(8,ielem),lnods_old(1,ielem),lnodx( 8)) ! Edge 8-1
                ! 1
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnods_old(5,ielem)
                lnods(3,kelem) = lnods_old(9,ielem)
                lnods(4,kelem) = lnods_old(8,ielem)
                lnods(5,kelem) = lnodx( 1)
                lnods(6,kelem) = ipoic + 5
                lnods(7,kelem) = ipoic + 8
                lnods(8,kelem) = lnodx( 8)
                lnods(9,kelem) = ipoic + 1

                ltype(kelem)   = QUA09
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(5,ielem)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnods_old(6,ielem)
                lnods(4,kelem) = lnods_old(9,ielem)
                lnods(5,kelem) = lnodx( 2)
                lnods(6,kelem) = lnodx( 3)
                lnods(7,kelem) = ipoic + 6
                lnods(8,kelem) = ipoic + 5
                lnods(9,kelem) = ipoic + 2

                ltype(kelem)   = QUA09
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(9,ielem)
                lnods(2,kelem) = lnods_old(6,ielem)
                lnods(3,kelem) = lnods_old(3,ielem)
                lnods(4,kelem) = lnods_old(7,ielem)
                lnods(5,kelem) = ipoic + 6
                lnods(6,kelem) = lnodx( 4)
                lnods(7,kelem) = lnodx( 5)
                lnods(8,kelem) = ipoic + 7
                lnods(9,kelem) = ipoic + 3

                ltype(kelem)   = QUA09
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(8,ielem)
                lnods(2,kelem) = lnods_old(9,ielem)
                lnods(3,kelem) = lnods_old(7,ielem)
                lnods(4,kelem) = lnods_old(4,ielem)
                lnods(5,kelem) = ipoic + 8
                lnods(6,kelem) = ipoic + 7
                lnods(7,kelem) = lnodx( 6)
                lnods(8,kelem) = lnodx( 7)
                lnods(9,kelem) = ipoic + 4

                ltype(kelem)   = QUA09
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem) 

                if( lelch_old(ielem) == ELEXT ) then
                   lelch(kelem)   = ELHOL
                   lelch(kelem-1) = ELHOL
                   lelch(kelem-2) = ELHOL
                end if
                ipoic = ipoic + 8

             case ( QUA16 ) 
                !
                ! QUA16
                !
                call mesh_multiplication_edge_node(lnods_old( 1,ielem),lnods_old( 5,ielem),lnodx( 1)) ! Edge 1-5
                call mesh_multiplication_edge_node(lnods_old( 5,ielem),lnods_old( 6,ielem),lnodx( 2)) ! Edge 5-6
                call mesh_multiplication_edge_node(lnods_old( 6,ielem),lnods_old( 2,ielem),lnodx( 3)) ! Edge 6-2
                call mesh_multiplication_edge_node(lnods_old( 2,ielem),lnods_old( 7,ielem),lnodx( 4)) ! Edge 2-7
                call mesh_multiplication_edge_node(lnods_old( 7,ielem),lnods_old( 8,ielem),lnodx( 5)) ! Edge 7-8
                call mesh_multiplication_edge_node(lnods_old( 8,ielem),lnods_old( 3,ielem),lnodx( 6)) ! Edge 8-3
                call mesh_multiplication_edge_node(lnods_old( 3,ielem),lnods_old( 9,ielem),lnodx( 7)) ! Edge 3-9
                call mesh_multiplication_edge_node(lnods_old( 9,ielem),lnods_old(10,ielem),lnodx( 8)) ! Edge 9-10
                call mesh_multiplication_edge_node(lnods_old(10,ielem),lnods_old( 4,ielem),lnodx( 9)) ! Edge 10-4
                call mesh_multiplication_edge_node(lnods_old( 4,ielem),lnods_old(11,ielem),lnodx(10)) ! Edge 4-11
                call mesh_multiplication_edge_node(lnods_old(11,ielem),lnods_old(12,ielem),lnodx(11)) ! Edge 11-12
                call mesh_multiplication_edge_node(lnods_old(12,ielem),lnods_old( 1,ielem),lnodx(12)) ! Edge 12-1
                ! 1
                kelem           = kelem + 1

                lnods( 1,kelem) = lnods_old( 1,ielem) 
                lnods( 2,kelem) = lnodx( 2)  
                lnods( 3,kelem) = ipoic + 11        
                lnods( 4,kelem) = lnodx(11)  

                lnods( 5,kelem) = lnodx( 1)
                lnods( 6,kelem) = lnods_old( 5,ielem)
                lnods( 7,kelem) = ipoic +  3
                lnods( 8,kelem) = ipoic +  7
                lnods( 9,kelem) = ipoic + 10
                lnods(10,kelem) = ipoic +  9
                lnods(11,kelem) = lnods_old(12,ielem)
                lnods(12,kelem) = lnodx(12)

                lnods(13,kelem) = ipoic + 1
                lnods(14,kelem) = ipoic + 2
                lnods(15,kelem) = lnods_old(13,ielem)
                lnods(16,kelem) = ipoic + 6

                ltype(kelem)   = QUA16
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem           = kelem + 1

                lnods( 1,kelem) = lnodx( 2) 
                lnods( 2,kelem) = lnods_old( 2,ielem)
                lnods( 3,kelem) = lnodx( 5)
                lnods( 4,kelem) = ipoic + 11

                lnods( 5,kelem) = lnods_old( 6,ielem)
                lnods( 6,kelem) = lnodx( 3)
                lnods( 7,kelem) = lnodx( 4)
                lnods( 8,kelem) = lnods_old( 7,ielem)
                lnods( 9,kelem) = ipoic + 13
                lnods(10,kelem) = ipoic + 12
                lnods(11,kelem) = ipoic +  7
                lnods(12,kelem) = ipoic +  3

                lnods(13,kelem) = ipoic +  4
                lnods(14,kelem) = ipoic +  5
                lnods(15,kelem) = ipoic +  8
                lnods(16,kelem) = lnods_old(14,ielem)

                ltype(kelem)   = QUA16
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem           = kelem + 1

                lnods( 1,kelem) = lnodx(11)
                lnods( 2,kelem) = ipoic + 11
                lnods( 3,kelem) = lnodx( 8)
                lnods( 4,kelem) = lnods_old( 4,ielem)

                lnods( 5,kelem) = ipoic +  9
                lnods( 6,kelem) = ipoic + 10
                lnods( 7,kelem) = ipoic + 15
                lnods( 8,kelem) = ipoic + 19
                lnods( 9,kelem) = lnods_old(10,ielem)
                lnods(10,kelem) = lnodx( 9)
                lnods(11,kelem) = lnodx(10)
                lnods(12,kelem) = lnods_old(11,ielem)

                lnods(13,kelem) = ipoic + 14
                lnods(14,kelem) = lnods_old(16,ielem)
                lnods(15,kelem) = ipoic + 18
                lnods(16,kelem) = ipoic + 17

                ltype(kelem)   = QUA16
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem           = kelem + 1

                lnods( 1,kelem) = ipoic + 11
                lnods( 2,kelem) = lnodx( 5)
                lnods( 3,kelem) = lnods_old( 3,ielem)
                lnods( 4,kelem) = lnodx( 8)

                lnods( 5,kelem) = ipoic + 12
                lnods( 6,kelem) = ipoic + 13
                lnods( 7,kelem) = lnods_old( 8,ielem)
                lnods( 8,kelem) = lnodx( 6)

                lnods( 9,kelem) = lnodx( 7)
                lnods(10,kelem) = lnods_old( 9,ielem)
                lnods(11,kelem) = ipoic + 19
                lnods(12,kelem) = ipoic + 15

                lnods(13,kelem) = lnods_old(15,ielem)
                lnods(14,kelem) = ipoic + 16
                lnods(15,kelem) = ipoic + 21
                lnods(16,kelem) = ipoic + 20

                ltype(kelem)   = QUA16
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + 21

             case default

                call runend('SUBELM: ELEMENT NOT CODED: '//integer_to_string(ltype_old(ielem)))

             end select

          end do

       else

          do ielem = 1,nelem_old

             pelty = abs(ltype_old(ielem))

             select case ( pelty )
                
             case ( TET04 ) 

                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx(2)) ! Edge 2-3
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(3,ielem),lnodx(3)) ! Edge 1-3
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(4,ielem),lnodx(4)) ! Edge 1-4
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(4,ielem),lnodx(5)) ! Edge 2-4
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx(6)) ! Edge 3-4
                !
                !-|CURVED| Curved mesh division for TET04
                !-|CURVED| ******HERE IS WHERE I HAVE TO CODE for TET04
                if( multiply_with_curvature.eq.1_ip ) then
                   call runend('SUBMSH: not implemented curved division for TET04')
                   !                 if(geometry_order(ielem).gt.1) then
                   !                 end if
                end if
                ! 1
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnodx(1)
                lnods(3,kelem) = lnodx(3)
                lnods(4,kelem) = lnodx(4)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(5)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(2)
                lnods(2,kelem) = lnods_old(3,ielem)
                lnods(3,kelem) = lnodx(3)
                lnods(4,kelem) = lnodx(6)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnodx(5)
                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(4)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 5
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(5)
                lnods(2,kelem) = lnodx(6)
                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(4)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 6
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(3)
                lnods(2,kelem) = lnodx(2)
                lnods(3,kelem) = lnodx(6)
                lnods(4,kelem) = lnodx(4)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 7
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(3)
                lnods(2,kelem) = lnodx(1)

                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(4)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 8
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(6)
                lnods(2,kelem) = lnodx(4)
                lnods(3,kelem) = lnodx(5)
                lnods(4,kelem) = lnods_old(4,ielem)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

             case ( TET10 ) 
                !
                ! TET10
                !
                ! Edge nodes
                !
                !
                ii = 0
                do iedge = 1,size(list_edges_TET10,2)
                   pblty = element_type(pelty) % type_edges(iedge)
                   do isegm = 1,element_type(pblty) % order
                      lnode(1:2) = elmgeo_segment(pelty,iedge,isegm)
                      ipoin      = lnods_old(lnode(1),ielem)
                      jpoin      = lnods_old(lnode(2),ielem)
                      ii         = ii + 1
                      call mesh_multiplication_edge_node(ipoin,jpoin,lnodx(ii))
                   end do
                end do
                !
                ! Face nodes
                !
                do iface = 1,element_type(pelty) % number_faces
                   pblty = element_type(pelty) % type_faces(iface)
                   call mesh_multiplication_face_node(lnods_old(element_type(pelty) % list_faces(:,iface),ielem),ipoin,ifacg)
                   call mesh_multiplication_subface(iface,ifacg,pblty,lnods_old(element_type(pelty) % list_faces(:,iface),ielem))
                end do
                ! 1: 1-5-7-8
                kelem           = kelem + 1
                ii = kelem
                lnods( 1,kelem) = lnods_old(1,ielem)
                lnods( 2,kelem) = lnods_old(5,ielem)
                lnods( 3,kelem) = lnods_old(7,ielem)
                lnods( 4,kelem) = lnods_old(8,ielem)
                lnods( 5,kelem) = lnodx(1)
                lnods( 6,kelem) = face_node(1,4)
                lnods( 7,kelem) = lnodx(6)
                lnods( 8,kelem) = lnodx(7)
                lnods( 9,kelem) = face_node(1,1)
                lnods(10,kelem) = face_node(2,3)
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 2: 5-2-6-9
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(5,ielem)
                lnods( 2,kelem) = lnods_old(2,ielem)
                lnods( 3,kelem) = lnods_old(6,ielem)
                lnods( 4,kelem) = lnods_old(9,ielem)
                lnods( 5,kelem) = lnodx(2)           ! 5-2
                lnods( 6,kelem) = lnodx(3)           ! 2-6
                lnods( 7,kelem) = face_node(3,4)     ! 6-5
                lnods( 8,kelem) = face_node(2,1)     ! 5-9
                lnods( 9,kelem) = lnodx(9)           ! 2-9
                lnods(10,kelem) = face_node(1,2)     ! 6-9
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 3: 6-3-7-10
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 6,ielem)
                lnods( 2,kelem) = lnods_old( 3,ielem)
                lnods( 3,kelem) = lnods_old( 7,ielem)
                lnods( 4,kelem) = lnods_old(10,ielem)
                lnods( 5,kelem) = lnodx(4)            ! 6-3
                lnods( 6,kelem) = lnodx(5)            ! 3-7
                lnods( 7,kelem) = face_node(2,4)      ! 7-6
                lnods( 8,kelem) = face_node(2,2)      ! 6-10
                lnods( 9,kelem) = lnodx(11)           ! 3-10
                lnods(10,kelem) = face_node(1,3)      ! 7-10
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 4: 5-9-6-8
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 5,ielem)
                lnods( 2,kelem) = lnods_old( 9,ielem)
                lnods( 3,kelem) = lnods_old( 6,ielem)
                lnods( 4,kelem) = lnods_old( 8,ielem)
                lnods( 5,kelem) = face_node(2,1)
                lnods( 6,kelem) = face_node(1,2)
                lnods( 7,kelem) = face_node(3,4)
                lnods( 8,kelem) = face_node(1,1)
                lnods( 9,kelem) = face_node(3,1)
                lnods(10,kelem) = ipoic + 1
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 5:  9-10-6-8
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 9,ielem) 
                lnods( 2,kelem) = lnods_old(10,ielem)
                lnods( 3,kelem) = lnods_old( 6,ielem)
                lnods( 4,kelem) = lnods_old( 8,ielem)
                lnods( 5,kelem) = face_node(3,2) ! 9-10
                lnods( 6,kelem) = face_node(2,2) ! 10-6
                lnods( 7,kelem) = face_node(1,2) ! 6-9
                lnods( 8,kelem) = face_node(3,1) ! 9-8
                lnods( 9,kelem) = face_node(3,3) ! 10-8
                lnods(10,kelem) = ipoic + 1      ! 6-8
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 6: 7-6-10-8
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 7,ielem) 
                lnods( 2,kelem) = lnods_old( 6,ielem)
                lnods( 3,kelem) = lnods_old(10,ielem)
                lnods( 4,kelem) = lnods_old( 8,ielem)
                lnods( 5,kelem) = face_node(2,4) ! 7-6
                lnods( 6,kelem) = face_node(2,2) ! 6-10
                lnods( 7,kelem) = face_node(1,3) ! 10-7
                lnods( 8,kelem) = face_node(2,3) ! 7-8
                lnods( 9,kelem) = ipoic + 1      ! 6-8
                lnods(10,kelem) = face_node(3,3) ! 10-8
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 7
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 7,ielem) 
                lnods( 2,kelem) = lnods_old( 5,ielem)
                lnods( 3,kelem) = lnods_old( 6,ielem)
                lnods( 4,kelem) = lnods_old( 8,ielem)
                lnods( 5,kelem) = face_node(1,4) ! 7-5
                lnods( 6,kelem) = face_node(3,4) ! 5-6
                lnods( 7,kelem) = face_node(2,4) ! 6-7
                lnods( 8,kelem) = face_node(2,3) ! 7-8 ok
                lnods( 9,kelem) = face_node(1,1) ! 5-8
                lnods(10,kelem) = ipoic + 1      ! 6-8
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)
                ! 8
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(10,ielem) 
                lnods( 2,kelem) = lnods_old( 8,ielem)
                lnods( 3,kelem) = lnods_old( 9,ielem)
                lnods( 4,kelem) = lnods_old( 4,ielem)
                lnods( 5,kelem) = face_node(3,3) ! 10-8
                lnods( 6,kelem) = face_node(3,1) ! 8-9
                lnods( 7,kelem) = face_node(3,2) ! 9-10
                lnods( 8,kelem) = lnodx(12)      ! 10-4
                lnods( 9,kelem) = lnodx(8)       ! 8-4
                lnods(10,kelem) = lnodx(10)      ! 9-4
                ltype(kelem)    = TET10
                lelch(kelem)    = lelch_old(ielem)
                lesub(kelem)    = lesub_old(ielem)
                lmate(kelem)    = lmate_old(ielem)
                if( if_lnnod ) lnnod(kelem) = lnnod_old(ielem)

                do kk = kelem+1,ii+7
                   lnods(:,kk) = lnods(:,ii)
                   ltype(kk)   = ltype(ii)
                   lmate(kk)   = lmate(ii)
                   lelch(kk)   = lelch(ii)
                   lesub(kk)   = lesub(ii)
                end do

                ipoic = ipoic + newno(TET10)

             case ( HEX08 ) 

                if ( lelch_old(ielem) == ELFEM .and. kfl_elndi == 0 ) then
                   !
                   !       8           7
                   !       o-----*-----o
                   !      /     /     /|
                   !     *-----x-----* |
                   !   5/     /    6/| *
                   !   o-----*-----o |/|
                   !   |     |     | x |
                   !   |     |     |/| o
                   !   *-----x-----* |/3
                   !   |     |     | *
                   !   |     |     |/
                   !   o-----*-----o
                   !   1           2
                   !
                   call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx( 1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx( 2)) ! Edge 2-3
                   call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx( 3)) ! Edge 3-4
                   call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(1,ielem),lnodx( 4)) ! Edge 4-1
                   call mesh_multiplication_edge_node(lnods_old(5,ielem),lnods_old(6,ielem),lnodx( 5)) ! Edge 5-6
                   call mesh_multiplication_edge_node(lnods_old(6,ielem),lnods_old(7,ielem),lnodx( 6)) ! Edge 6-7
                   call mesh_multiplication_edge_node(lnods_old(7,ielem),lnods_old(8,ielem),lnodx( 7)) ! Edge 7-8
                   call mesh_multiplication_edge_node(lnods_old(8,ielem),lnods_old(5,ielem),lnodx( 8)) ! Edge 8-5
                   call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(5,ielem),lnodx( 9)) ! Edge 1-5
                   call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(6,ielem),lnodx(10)) ! Edge 2-6
                   call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(7,ielem),lnodx(11)) ! Edge 3-7
                   call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(8,ielem),lnodx(12)) ! Edge 4-8
                   !
                   ipoif = 12
                   do iface = 1,6
                      node1 = element_type(HEX08) % list_faces(1,iface) 
                      ipoi1 = lnods_old(node1,ielem)
                      node2 = element_type(HEX08) % list_faces(2,iface) 
                      ipoi2 = lnods_old(node2,ielem)
                      node3 = element_type(HEX08) % list_faces(3,iface) 
                      ipoi3 = lnods_old(node3,ielem)
                      node4 = element_type(HEX08) % list_faces(4,iface) 
                      ipoi4 = lnods_old(node4,ielem)
                      ipoif = ipoif + 1
                      ifacg = facel(1,iface,ielem)
                      lnodx(ipoif) = lfacg(mnodm+1,ifacg)
                      !call mesh_multiplication_face_node(ipoi1,ipoi2,ipoi3,ipoi4,lnodx(ipoif))         ! Face IFACE
                   end do
                   ipoic          = ipoic + 1
                   lnodx(19)      = ipoic

                   !-|CURVED| Curved mesh division for HEX08
                   if( multiply_with_curvature.eq.1_ip ) then
                      if(geometry_order(ielem).gt.1_ip) then
                         numNewCoordinates = 19_ip
                         numNewElements    = 8_ip

                         allocate( subdivided_coordinates(numNewCoordinates,ndime) )
                         allocate( subdivided_element_geometry(geometry_numNodes(ielem),ndime,numNewElements) )

                         call curve_subdivided_element(&
                              HEX08,&
                              numNewCoordinates,&
                              numNewElements,&
                              geometry_order(ielem),&
                              geometry_numNodes(ielem),&
                              ndime,&
                              curved_geometry(1:geometry_numNodes(ielem),1:ndime,ielem),&
                              subdivided_coordinates,&
                              subdivided_element_geometry)

                         coord(1:ndime,lnodx(1:19)) = transpose(subdivided_coordinates)
                         subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+8)) = subdivided_element_geometry

                         deallocate(subdivided_coordinates)
                         deallocate(subdivided_element_geometry)
                      else
                         subdivided_curved_geometry(1:geometry_numNodes(ielem),1:ndime,(kelem+1):(kelem+8)) = 0.0_rp
                      end if
                   end if

                   ! 1
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnods_old(1,ielem)
                   lnods(2,kelem) = lnodx( 1)
                   lnods(3,kelem) = lnodx(13)
                   lnods(4,kelem) = lnodx( 4)
                   lnods(5,kelem) = lnodx( 9)
                   lnods(6,kelem) = lnodx(17)
                   lnods(7,kelem) = lnodx(19)
                   lnods(8,kelem) = lnodx(16)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 2
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 1)
                   lnods(2,kelem) = lnods_old(2,ielem)
                   lnods(3,kelem) = lnodx( 2)
                   lnods(4,kelem) = lnodx(13)
                   lnods(5,kelem) = lnodx(17)
                   lnods(6,kelem) = lnodx(10)
                   lnods(7,kelem) = lnodx(14)
                   lnods(8,kelem) = lnodx(19)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 3
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(13)
                   lnods(2,kelem) = lnodx( 2)
                   lnods(3,kelem) = lnods_old(3,ielem)
                   lnods(4,kelem) = lnodx( 3)
                   lnods(5,kelem) = lnodx(19)
                   lnods(6,kelem) = lnodx(14)
                   lnods(7,kelem) = lnodx(11)
                   lnods(8,kelem) = lnodx(18)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 4
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 4)
                   lnods(2,kelem) = lnodx(13)
                   lnods(3,kelem) = lnodx( 3)
                   lnods(4,kelem) = lnods_old(4,ielem)
                   lnods(5,kelem) = lnodx(16)
                   lnods(6,kelem) = lnodx(19)
                   lnods(7,kelem) = lnodx(18)
                   lnods(8,kelem) = lnodx(12)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 5
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 9)
                   lnods(2,kelem) = lnodx(17)
                   lnods(3,kelem) = lnodx(19)
                   lnods(4,kelem) = lnodx(16)
                   lnods(5,kelem) = lnods_old(5,ielem)
                   lnods(6,kelem) = lnodx( 5)
                   lnods(7,kelem) = lnodx(15)
                   lnods(8,kelem) = lnodx( 8)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 6
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(17)
                   lnods(2,kelem) = lnodx(10)
                   lnods(3,kelem) = lnodx(14)
                   lnods(4,kelem) = lnodx(19)
                   lnods(5,kelem) = lnodx( 5)
                   lnods(6,kelem) = lnods_old(6,ielem)
                   lnods(7,kelem) = lnodx( 6)
                   lnods(8,kelem) = lnodx(15)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 7
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(19)
                   lnods(2,kelem) = lnodx(14)
                   lnods(3,kelem) = lnodx(11)
                   lnods(4,kelem) = lnodx(18)
                   lnods(5,kelem) = lnodx(15)
                   lnods(6,kelem) = lnodx( 6)
                   lnods(7,kelem) = lnods_old(7,ielem)
                   lnods(8,kelem) = lnodx( 7)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 8
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx(16)
                   lnods(2,kelem) = lnodx(19)
                   lnods(3,kelem) = lnodx(18)
                   lnods(4,kelem) = lnodx(12)
                   lnods(5,kelem) = lnodx( 8)
                   lnods(6,kelem) = lnodx(15)
                   lnods(7,kelem) = lnodx( 7)
                   lnods(8,kelem) = lnods_old(8,ielem)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)

                else if ( lelch_old(ielem) == ELINT .or. (lelch_old(ielem) == ELFEM .and. kfl_elndi == 1) ) then
                   !
                   !       8           7
                   !       o-----*-----o
                   !      /     /     /|
                   !     *-----x-----* |
                   !   5/     /    6/| |
                   !   o-----*-----o | |
                   !   |     |     | | |
                   !   |     |     | | o
                   !   |     |     | |/3
                   !   |     |     | *
                   !   |     |     |/
                   !   o-----*-----o
                   !   1           2
                   !
                   call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx(2)) ! Edge 2-3
                   call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx(3)) ! Edge 3-4
                   call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(1,ielem),lnodx(4)) ! Edge 4-1
                   call mesh_multiplication_edge_node(lnods_old(5,ielem),lnods_old(6,ielem),lnodx(5)) ! Edge 5-6
                   call mesh_multiplication_edge_node(lnods_old(6,ielem),lnods_old(7,ielem),lnodx(6)) ! Edge 6-7
                   call mesh_multiplication_edge_node(lnods_old(7,ielem),lnods_old(8,ielem),lnodx(7)) ! Edge 7-8
                   call mesh_multiplication_edge_node(lnods_old(8,ielem),lnods_old(5,ielem),lnodx(8)) ! Edge 8-5
                   !
                   ipoif = 8

                   ! Face 1
                   node1 = element_type(HEX08) % list_faces(1,1) 
                   ipoi1 = lnods_old(node1,ielem)
                   node2 = element_type(HEX08) % list_faces(2,1) 
                   ipoi2 = lnods_old(node2,ielem)
                   node3 = element_type(HEX08) % list_faces(3,1) 
                   ipoi3 = lnods_old(node3,ielem)
                   node4 = element_type(HEX08) % list_faces(4,1) 
                   ipoi4 = lnods_old(node4,ielem)
                   ipoif = ipoif + 1
                   ifacg = facel(1,1,ielem)
                   lnodx(ipoif) = lfacg(mnodm+1,ifacg)
                   ! Face 3
                   node1 = element_type(HEX08) % list_faces(1,3) 
                   ipoi1 = lnods_old(node1,ielem)
                   node2 = element_type(HEX08) % list_faces(2,3) 
                   ipoi2 = lnods_old(node2,ielem)
                   node3 = element_type(HEX08) % list_faces(3,3) 
                   ipoi3 = lnods_old(node3,ielem)
                   node4 = element_type(HEX08) % list_faces(4,3) 
                   ipoi4 = lnods_old(node4,ielem)
                   ipoif = ipoif + 1
                   ifacg = facel(1,3,ielem)
                   lnodx(ipoif) = lfacg(mnodm+1,ifacg)

                   !-|CURVED| Curved mesh division for HEX08
                   if ( multiply_with_curvature .eq. 1_ip ) call runend('SUBMSH: not implemented curved division for HEX08-ELINT')

                   ! 1
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnods_old(1,ielem)
                   lnods(2,kelem) = lnodx( 1)
                   lnods(3,kelem) = lnodx( 9)
                   lnods(4,kelem) = lnodx( 4)
                   lnods(5,kelem) = lnods_old(5,ielem)
                   lnods(6,kelem) = lnodx( 5)
                   lnods(7,kelem) = lnodx(10)
                   lnods(8,kelem) = lnodx( 8)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 2
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 1)
                   lnods(2,kelem) = lnods_old(2,ielem)
                   lnods(3,kelem) = lnodx( 2)
                   lnods(4,kelem) = lnodx( 9)
                   lnods(5,kelem) = lnodx( 5)
                   lnods(6,kelem) = lnods_old(6,ielem)
                   lnods(7,kelem) = lnodx( 6)
                   lnods(8,kelem) = lnodx(10)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 3
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 9)
                   lnods(2,kelem) = lnodx( 2)
                   lnods(3,kelem) = lnods_old(3,ielem)
                   lnods(4,kelem) = lnodx( 3)
                   lnods(5,kelem) = lnodx(10)
                   lnods(6,kelem) = lnodx( 6)
                   lnods(7,kelem) = lnods_old(7,ielem)
                   lnods(8,kelem) = lnodx( 7)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                   ! 4
                   kelem          = kelem + 1
                   lnods(1,kelem) = lnodx( 4)
                   lnods(2,kelem) = lnodx( 9)
                   lnods(3,kelem) = lnodx( 3)
                   lnods(4,kelem) = lnods_old(4,ielem)
                   lnods(5,kelem) = lnodx( 8)
                   lnods(6,kelem) = lnodx(10)
                   lnods(7,kelem) = lnodx( 7)
                   lnods(8,kelem) = lnods_old(8,ielem)
                   ltype(kelem)   = HEX08
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)

                end if

             case ( HEX27 ) 
                !
                ! HEX27
                !
                ! Edge nodes numbered first, then face nodes
                !
                ii = 0
                do iedge = 1,size(list_edges_HEX27,2)
                   pblty = element_type(pelty) % type_edges(iedge)
                   do isegm = 1,element_type(pblty) % order
                      lnode(1:2) = elmgeo_segment(pelty,iedge,isegm)
                      ipoin      = lnods_old(lnode(1),ielem)
                      jpoin      = lnods_old(lnode(2),ielem)
                      ii         = ii + 1
                      call mesh_multiplication_edge_node(ipoin,jpoin,lnodx(ii))
                   end do
                end do
                !
                ! Face nodes
                !
                do iface = 1,element_type(pelty) % number_faces
                   pblty = element_type(pelty) % type_faces(iface)
                   call mesh_multiplication_face_node(lnods_old(element_type(pelty) % list_faces(:,iface),ielem),ipoin,ifacg)
                   call mesh_multiplication_subface(iface,ifacg,pblty,lnods_old(element_type(pelty) % list_faces(:,iface),ielem))
                   !ipoin = face_node(1,iface)
                   !write(*,'(a,i1,a,i1,a,3(1x,f4.1),a,i4)') 'face_node=(',1_ip,',',iface,')',coord(1:ndime,ipoin),': ',ipoin
                   !do ii = 1,7
                   !   ipoin = face_node(ii+1,iface)
                   !   write(*,'(a,i1,a,i1,a,3(1x,f4.1),a,i4)') 'face_node=(',ii+1,',',iface,')',coord(1:ndime,ipoin+ii),': ',ipoin+ii
                   !end do
                end do

                ! 1                
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 1,ielem)
                lnods( 2,kelem) = lnods_old( 9,ielem)
                lnods( 3,kelem) = lnods_old(21,ielem)
                lnods( 4,kelem) = lnods_old(12,ielem)
                lnods( 5,kelem) = lnods_old(13,ielem)
                lnods( 6,kelem) = lnods_old(22,ielem)
                lnods( 7,kelem) = lnods_old(27,ielem)
                lnods( 8,kelem) = lnods_old(25,ielem)
                lnods( 9,kelem) = lnodx(1)
                lnods(10,kelem) = face_node(4,5) 
                lnods(11,kelem) = face_node(2,5)
                lnods(12,kelem) = lnodx(8)

                lnods(13,kelem) = lnodx(17)
                lnods(14,kelem) = face_node(2,1) 
                lnods(15,kelem) = ipoic + 2
                lnods(16,kelem) = face_node(2,4)

                lnods(17,kelem) = face_node(4,1)
                lnods(18,kelem) = ipoic + 3
                lnods(19,kelem) = ipoic + 5
                lnods(20,kelem) = face_node(5,4)

                lnods(21,kelem) = face_node(1,5)
                lnods(22,kelem) = face_node(1,1)
                lnods(23,kelem) = ipoic +  1
                lnods(24,kelem) = ipoic +  4
                lnods(25,kelem) = face_node(3,4)

                lnods(26,kelem) = ipoic + 6
                lnods(27,kelem) = ipoic + 7

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2               
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old( 9,ielem)
                lnods( 2,kelem) = lnods_old( 2,ielem)
                lnods( 3,kelem) = lnods_old(10,ielem)
                lnods( 4,kelem) = lnods_old(21,ielem)
                lnods( 5,kelem) = lnods_old(22,ielem)
                lnods( 6,kelem) = lnods_old(14,ielem)
                lnods( 7,kelem) = lnods_old(23,ielem)
                lnods( 8,kelem) = lnods_old(27,ielem)
                lnods( 9,kelem) = lnodx(2)
                lnods(10,kelem) = lnodx(3)
                lnods(11,kelem) = face_node(7,5)
                lnods(12,kelem) = face_node(4,5)
                lnods(13,kelem) = face_node(2,1)
                lnods(14,kelem) = lnodx(19)
                lnods(15,kelem) = face_node(2,2)
                lnods(16,kelem) = ipoic +  2

                lnods(17,kelem) = face_node(5,1)
                lnods(18,kelem) = face_node(4,2)
                lnods(19,kelem) = ipoic +  8
                lnods(20,kelem) = ipoic +  3

                lnods(21,kelem) = face_node(6,5)

                lnods(22,kelem) = face_node(3,1)
                lnods(23,kelem) = face_node(1,2)
                lnods(24,kelem) = ipoic +  9
                lnods(25,kelem) = ipoic +  1

                lnods(26,kelem) = ipoic + 10
                lnods(27,kelem) = ipoic + 11

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3           
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(12,ielem)
                lnods( 2,kelem) = lnods_old(21,ielem)
                lnods( 3,kelem) = lnods_old(11,ielem)
                lnods( 4,kelem) = lnods_old( 4,ielem)
                lnods( 5,kelem) = lnods_old(25,ielem)
                lnods( 6,kelem) = lnods_old(27,ielem)
                lnods( 7,kelem) = lnods_old(24,ielem)
                lnods( 8,kelem) = lnods_old(16,ielem)

                lnods( 9,kelem) = face_node(2,5)
                lnods(10,kelem) = face_node(5,5)
                lnods(11,kelem) = lnodx(6)
                lnods(12,kelem) = lnodx(7)

                lnods(13,kelem) = face_node(2,4) 
                lnods(14,kelem) = ipoic +  2
                lnods(15,kelem) = face_node(2,3)
                lnods(16,kelem) = lnodx(23)

                lnods(17,kelem) = ipoic +  5
                lnods(18,kelem) = ipoic + 12
                lnods(19,kelem) = face_node(5,3)
                lnods(20,kelem) = face_node(4,4) 

                lnods(21,kelem) = face_node(3,5)

                lnods(22,kelem) = ipoic +  4
                lnods(23,kelem) = ipoic + 26
                lnods(24,kelem) = face_node(3,3)
                lnods(25,kelem) = face_node(1,4)

                lnods(26,kelem) = ipoic + 15
                lnods(27,kelem) = ipoic + 16

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4           
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(21,ielem)
                lnods( 2,kelem) = lnods_old(10,ielem)
                lnods( 3,kelem) = lnods_old( 3,ielem)
                lnods( 4,kelem) = lnods_old(11,ielem)
                lnods( 5,kelem) = lnods_old(27,ielem)
                lnods( 6,kelem) = lnods_old(23,ielem)
                lnods( 7,kelem) = lnods_old(15,ielem)
                lnods( 8,kelem) = lnods_old(24,ielem)

                lnods( 9,kelem) = face_node(7,5)
                lnods(10,kelem) = lnodx(4)      
                lnods(11,kelem) = lnodx(5)              
                lnods(12,kelem) = face_node(5,5)

                lnods(13,kelem) = ipoic +  2
                lnods(14,kelem) = face_node(2,2)
                lnods(15,kelem) = lnodx(21)
                lnods(16,kelem) = face_node(2,3)

                lnods(17,kelem) = ipoic +  8
                lnods(18,kelem) = face_node(5,2)
                lnods(19,kelem) = face_node(4,3)
                lnods(20,kelem) = ipoic + 12

                lnods(21,kelem) = face_node(8,5)

                lnods(22,kelem) = ipoic +  9
                lnods(23,kelem) = face_node(3,2)
                lnods(24,kelem) = face_node(1,3)
                lnods(25,kelem) = ipoic + 26

                lnods(26,kelem) = ipoic + 13
                lnods(27,kelem) = ipoic + 14

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 5          
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(13,ielem)
                lnods( 2,kelem) = lnods_old(22,ielem)
                lnods( 3,kelem) = lnods_old(27,ielem)
                lnods( 4,kelem) = lnods_old(25,ielem)
                lnods( 5,kelem) = lnods_old( 5,ielem)
                lnods( 6,kelem) = lnods_old(17,ielem)
                lnods( 7,kelem) = lnods_old(26,ielem)
                lnods( 8,kelem) = lnods_old(20,ielem)

                lnods( 9,kelem) = face_node(4,1)
                lnods(10,kelem) = ipoic +  3
                lnods(11,kelem) = ipoic +  5
                lnods(12,kelem) = face_node(5,4)

                lnods(13,kelem) = lnodx(18)
                lnods(14,kelem) = face_node(7,1)
                lnods(15,kelem) = ipoic + 17
                lnods(16,kelem) = face_node(7,4)

                lnods(17,kelem) = lnodx(9)
                lnods(18,kelem) = face_node(2,6)
                lnods(19,kelem) = face_node(4,6)
                lnods(20,kelem) = lnodx(16)

                lnods(21,kelem) = ipoic +  6

                lnods(22,kelem) = face_node(6,1)
                lnods(23,kelem) = ipoic + 18
                lnods(24,kelem) = ipoic + 19
                lnods(25,kelem) = face_node(8,4)

                lnods(26,kelem) = face_node(1,6)
                lnods(27,kelem) = ipoic + 20

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 6         
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(22,ielem)
                lnods( 2,kelem) = lnods_old(14,ielem)
                lnods( 3,kelem) = lnods_old(23,ielem)
                lnods( 4,kelem) = lnods_old(27,ielem)
                lnods( 5,kelem) = lnods_old(17,ielem)
                lnods( 6,kelem) = lnods_old( 6,ielem)
                lnods( 7,kelem) = lnods_old(18,ielem)
                lnods( 8,kelem) = lnods_old(26,ielem)

                lnods( 9,kelem) = face_node(5,1)
                lnods(10,kelem) = face_node(4,2)
                lnods(11,kelem) = ipoic + 8
                lnods(12,kelem) = ipoic + 3

                lnods(13,kelem) = face_node(7,1)
                lnods(14,kelem) = lnodx(20)
                lnods(15,kelem) = face_node(7,2)
                lnods(16,kelem) = ipoic + 17

                lnods(17,kelem) = lnodx(10)
                lnods(18,kelem) = lnodx(11)
                lnods(19,kelem) = face_node(5,6)
                lnods(20,kelem) = face_node(2,6)

                lnods(21,kelem) = ipoic + 10

                lnods(22,kelem) = face_node(8,1)
                lnods(23,kelem) = face_node(6,2)
                lnods(24,kelem) = ipoic + 21
                lnods(25,kelem) = ipoic + 18

                lnods(26,kelem) = face_node(3,6)
                lnods(27,kelem) = ipoic + 22

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 7        
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(25,ielem)
                lnods( 2,kelem) = lnods_old(27,ielem)
                lnods( 3,kelem) = lnods_old(24,ielem)
                lnods( 4,kelem) = lnods_old(16,ielem)
                lnods( 5,kelem) = lnods_old(20,ielem)
                lnods( 6,kelem) = lnods_old(26,ielem)
                lnods( 7,kelem) = lnods_old(19,ielem)
                lnods( 8,kelem) = lnods_old( 8,ielem)

                lnods( 9,kelem) = ipoic +  5
                lnods(10,kelem) = ipoic + 12
                lnods(11,kelem) = face_node(5,3)
                lnods(12,kelem) = face_node(4,4)

                lnods(13,kelem) = face_node(7,4)
                lnods(14,kelem) = ipoic + 17
                lnods(15,kelem) = face_node(7,3)
                lnods(16,kelem) = lnodx(24)

                lnods(17,kelem) = face_node(4,6)
                lnods(18,kelem) = face_node(7,6)
                lnods(19,kelem) = lnodx(14)
                lnods(20,kelem) = lnodx(15)

                lnods(21,kelem) = ipoic + 15

                lnods(22,kelem) = ipoic + 19
                lnods(23,kelem) = ipoic + 23
                lnods(24,kelem) = face_node(8,3)
                lnods(25,kelem) = face_node(6,4)

                lnods(26,kelem) = face_node(6,6)
                lnods(27,kelem) = ipoic + 24

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 8        
                kelem           = kelem + 1
                lnods( 1,kelem) = lnods_old(27,ielem)
                lnods( 2,kelem) = lnods_old(23,ielem)
                lnods( 3,kelem) = lnods_old(15,ielem)
                lnods( 4,kelem) = lnods_old(24,ielem)
                lnods( 5,kelem) = lnods_old(26,ielem)
                lnods( 6,kelem) = lnods_old(18,ielem)
                lnods( 7,kelem) = lnods_old( 7,ielem)
                lnods( 8,kelem) = lnods_old(19,ielem)

                lnods( 9,kelem) = ipoic +  8
                lnods(10,kelem) = face_node(5,2)
                lnods(11,kelem) = face_node(4,3)
                lnods(12,kelem) = ipoic + 12

                lnods(13,kelem) = ipoic + 17
                lnods(14,kelem) = face_node(7,2)
                lnods(15,kelem) = lnodx(22)
                lnods(16,kelem) = face_node(7,3)

                lnods(17,kelem) = face_node(5,6)
                lnods(18,kelem) = lnodx(12)
                lnods(19,kelem) = lnodx(13)
                lnods(20,kelem) = face_node(7,6)

                lnods(21,kelem) = ipoic + 13

                lnods(22,kelem) = ipoic + 21
                lnods(23,kelem) = face_node(8,2)
                lnods(24,kelem) = face_node(6,3)
                lnods(25,kelem) = ipoic + 23

                lnods(26,kelem) = face_node(8,6)
                lnods(27,kelem) = ipoic + 25

                ltype(kelem)   = HEX27
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                ipoic = ipoic + newno(HEX27)

             case ( PEN06 ) 

                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx( 1)) ! Edge 1-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx( 2)) ! Edge 2-3
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(1,ielem),lnodx( 3)) ! Edge 3-1

                call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(5,ielem),lnodx( 4)) ! Edge 4-5
                call mesh_multiplication_edge_node(lnods_old(5,ielem),lnods_old(6,ielem),lnodx( 5)) ! Edge 5-6
                call mesh_multiplication_edge_node(lnods_old(6,ielem),lnods_old(4,ielem),lnodx( 6)) ! Edge 6-4

                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(4,ielem),lnodx( 7)) ! Edge 1-4
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(5,ielem),lnodx( 8)) ! Edge 2-5
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(6,ielem),lnodx( 9)) ! Edge 3-6

                !-|CURVED| Curved mesh division for PEN06
                !-|CURVED| ******HERE IS WHERE I HAVE TO CODE for PEN06
                if( multiply_with_curvature.eq.1_ip ) then
                   call runend('SUBMSH: not implemented curved division for PEN06')
                   !                 if(geometry_order(ielem).gt.1) then
                   !                 end if
                end if

                ipoif = 9
                !do iface = 1,3
                do iface = 3,5
                   node1 = element_type(PEN06) % list_faces(1,iface) 
                   ipoi1 = lnods_old(node1,ielem)
                   node2 = element_type(PEN06) % list_faces(2,iface) 
                   ipoi2 = lnods_old(node2,ielem)
                   node3 = element_type(PEN06) % list_faces(3,iface) 
                   ipoi3 = lnods_old(node3,ielem)
                   node4 = element_type(PEN06) % list_faces(4,iface) 
                   ipoi4 = lnods_old(node4,ielem)
                   ipoif = ipoif + 1
                   ifacg = facel(1,iface,ielem)
                   lnodx(ipoif) = lfacg(mnodm+1,ifacg)
                   !call mesh_multiplication_face_node(ipoi1,ipoi2,ipoi3,ipoi4,lnodx(ipoif))         ! Face IFACE
                end do
                ! 1
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnodx( 1)
                lnods(3,kelem) = lnodx( 3)
                lnods(4,kelem) = lnodx( 7)
                lnods(5,kelem) = lnodx(10)
                lnods(6,kelem) = lnodx(12)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx( 1)
                lnods(2,kelem) = lnodx( 2)
                lnods(3,kelem) = lnodx( 3)
                lnods(4,kelem) = lnodx(10)
                lnods(5,kelem) = lnodx(11)
                lnods(6,kelem) = lnodx(12)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx( 1)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnodx( 2)
                lnods(4,kelem) = lnodx(10)
                lnods(5,kelem) = lnodx( 8)
                lnods(6,kelem) = lnodx(11)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx( 3)
                lnods(2,kelem) = lnodx( 2)
                lnods(3,kelem) = lnods_old(3,ielem)
                lnods(4,kelem) = lnodx(12)
                lnods(5,kelem) = lnodx(11)
                lnods(6,kelem) = lnodx( 9)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 5
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx( 7)
                lnods(2,kelem) = lnodx(10)
                lnods(3,kelem) = lnodx(12)
                lnods(4,kelem) = lnods_old(4,ielem)
                lnods(5,kelem) = lnodx( 4)
                lnods(6,kelem) = lnodx( 6)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 6
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(10)
                lnods(2,kelem) = lnodx(11)
                lnods(3,kelem) = lnodx(12)
                lnods(4,kelem) = lnodx( 4)
                lnods(5,kelem) = lnodx( 5)
                lnods(6,kelem) = lnodx( 6)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 7
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(10)
                lnods(2,kelem) = lnodx( 8)
                lnods(3,kelem) = lnodx(11)
                lnods(4,kelem) = lnodx( 4)
                lnods(5,kelem) = lnods_old(5,ielem)
                lnods(6,kelem) = lnodx( 5)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 8
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(12)
                lnods(2,kelem) = lnodx(11)
                lnods(3,kelem) = lnodx( 9)
                lnods(4,kelem) = lnodx( 6)
                lnods(5,kelem) = lnodx( 5)
                lnods(6,kelem) = lnods_old(6,ielem)
                ltype(kelem)   = PEN06
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

             case ( PYR05 ) 

                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx( 1)) ! Edge 1-2
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(3,ielem),lnodx( 2)) ! Edge 2-3
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(4,ielem),lnodx( 3)) ! Edge 3-4
                call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(1,ielem),lnodx( 4)) ! Edge 4-1

                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(5,ielem),lnodx( 5)) ! Edge 1-5
                call mesh_multiplication_edge_node(lnods_old(2,ielem),lnods_old(5,ielem),lnodx( 6)) ! Edge 2-5
                call mesh_multiplication_edge_node(lnods_old(3,ielem),lnods_old(5,ielem),lnodx( 7)) ! Edge 3-5
                call mesh_multiplication_edge_node(lnods_old(4,ielem),lnods_old(5,ielem),lnodx( 8)) ! Edge 4-5
                !
                !-|CURVED| Curved mesh division for PYR05
                !-|CURVED| ******HERE IS WHERE I HAVE TO CODE for PYR05
                if( multiply_with_curvature.eq.1_ip ) then
                   call runend('SUBMSH: not implemented curved division for PYR05')
                   !                 if(geometry_order(ielem).gt.1) then
                   !                 end if
                end if
                !
                ifacg    = facel(1,1,ielem)
                lnodx(9) = lfacg(mnodm+1,ifacg)
                ! 1
                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnodx(1)
                lnods(3,kelem) = lnodx(9)
                lnods(4,kelem) = lnodx(4)
                lnods(5,kelem) = lnodx(5)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 2
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnods_old(2,ielem)
                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(9)
                lnods(5,kelem) = lnodx(6)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 3
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(9)
                lnods(2,kelem) = lnodx(2)
                lnods(3,kelem) = lnods_old(3,ielem)
                lnods(4,kelem) = lnodx(3)
                lnods(5,kelem) = lnodx(7)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 4
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(4)
                lnods(2,kelem) = lnodx(9)
                lnods(3,kelem) = lnodx(3)
                lnods(4,kelem) = lnods_old(4,ielem)
                lnods(5,kelem) = lnodx(8)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 5
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(8)
                lnods(2,kelem) = lnodx(7)
                lnods(3,kelem) = lnodx(6)
                lnods(4,kelem) = lnodx(5)
                lnods(5,kelem) = lnodx(9)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 6
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(5)
                lnods(2,kelem) = lnodx(6)
                lnods(3,kelem) = lnodx(7)
                lnods(4,kelem) = lnodx(8)
                lnods(5,kelem) = lnods_old(5,ielem)
                ltype(kelem)   = PYR05
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 7
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(5)
                lnods(2,kelem) = lnodx(6)
                lnods(3,kelem) = lnodx(1)
                lnods(4,kelem) = lnodx(9)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = 4
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 8
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(6)
                lnods(2,kelem) = lnodx(7)
                lnods(3,kelem) = lnodx(2)
                lnods(4,kelem) = lnodx(9)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = 4
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 9
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(7)
                lnods(2,kelem) = lnodx(8)
                lnods(3,kelem) = lnodx(3)
                lnods(4,kelem) = lnodx(9)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = 4
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)
                ! 10
                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(8)
                lnods(2,kelem) = lnodx(5)
                lnods(3,kelem) = lnodx(4)
                lnods(4,kelem) = lnodx(9)
                ltype(kelem)   = TET04
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = 4
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

             case ( BAR3D ) 
                !
                ! o-----*-----o
                !
                call mesh_multiplication_edge_node(lnods_old(1,ielem),lnods_old(2,ielem),lnodx(1)) ! Edge 1-2

                kelem          = kelem + 1
                lnods(1,kelem) = lnods_old(1,ielem)
                lnods(2,kelem) = lnodx(1)
                ltype(kelem)   = BAR3D
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

                kelem          = kelem + 1
                lnods(1,kelem) = lnodx(1)
                lnods(2,kelem) = lnods_old(2,ielem)
                ltype(kelem)   = BAR3D
                lelch(kelem)   = lelch_old(ielem)
                if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                lesub(kelem)   = lesub_old(ielem)
                lmate(kelem)   = lmate_old(ielem)

             case default 
                
                if( ltype_old(ielem) < 0 ) then
                   !
                   ! Virtual element
                   !
                   kelem          = kelem + 1
                   do inode = 1,elmgeo_number_nodes(ltype_old(ielem),lnods_old(:,ielem))
                      lnods(inode,kelem) = lnods_old(inode,ielem)
                   end do
                   ltype(kelem)   = ltype_old(ielem)
                   lelch(kelem)   = lelch_old(ielem)
                   if( if_lnnod ) lnnod(kelem)   = lnnod_old(ielem)
                   lesub(kelem)   = lesub_old(ielem)
                   lmate(kelem)   = lmate_old(ielem)
                else
                   call runend('NOT CODED')
                end if
                
             end select
          end do
       end if
       !
       ! Deallocate memory of FACEL
       !
       if( nfacg > 0 ) then
          call memory_deallo(memor_dom,'FACEL','mesh_multiplication_divide_elements',facel)
       end if
       !
       ! NBOUN arrays
       !
       kboun = 0
       if( ndime == 1 ) then
          do iboun = 1,nboun_old
             kboun          = kboun + 1
             lnodb(1,kboun) = lnodb_old(1,iboun)
             ltypb(kboun)   = ltypb_old(iboun)
             lboch(kboun)   = lboch_old(iboun)
          end do

       else if( ndime == 2 ) then

          do iboun = 1,nboun_old

             pblty = abs(ltypb_old(iboun))

             if( divide_boun(iboun) ) then

                select case ( pblty )

                case ( BAR02 ) 
                   !
                   ! BAR02
                   !
                   call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(2,iboun),lnodx(1)) ! Edge 1-2

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodb_old(1,iboun)
                   lnodb(2,kboun) = lnodx(1)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)
                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodx(1)
                   lnodb(2,kboun) = lnodb_old(2,iboun)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                case ( BAR03 ) 
                   !
                   ! BAR03
                   !
                   call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(3,iboun),lnodx(1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(lnodb_old(3,iboun),lnodb_old(2,iboun),lnodx(2)) ! Edge 3-2

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodb_old(1,iboun)
                   lnodb(2,kboun) = lnodb_old(3,iboun)
                   lnodb(3,kboun) = lnodx(1)                  
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodb_old(3,iboun)
                   lnodb(2,kboun) = lnodb_old(2,iboun)
                   lnodb(3,kboun) = lnodx(2)                                    
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                case ( BAR04 ) 
                   !
                   ! BAR04
                   !
                   call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(3,iboun),lnodx(1)) ! Edge 1-3
                   call mesh_multiplication_edge_node(lnodb_old(3,iboun),lnodb_old(4,iboun),lnodx(2)) ! Edge 3-4
                   call mesh_multiplication_edge_node(lnodb_old(4,iboun),lnodb_old(2,iboun),lnodx(3)) ! Edge 4-2

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodb_old(1,iboun)
                   lnodb(2,kboun) = lnodx(2)      
                   lnodb(3,kboun) = lnodx(1)                  
                   lnodb(4,kboun) = lnodb_old(3,iboun)               
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodx(2) 
                   lnodb(2,kboun) = lnodb_old(2,iboun)
                   lnodb(3,kboun) = lnodb_old(4,iboun)                               
                   lnodb(4,kboun) = lnodx(3)                                    
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                end select

             else

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(1,iboun)
                lnodb(2,kboun) = lnodb_old(2,iboun)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

             end if

          end do

       else

          do iboun = 1,nboun_old

             ielem = lelbo_old(iboun)
             pelty = abs(ltype_old(ielem))
             pblty = abs(ltypb_old(iboun))
             pnodb = element_type(pblty) % number_nodes

             select case ( pblty )

             case ( TRI03 ) 
                !
                ! TRI03
                !
                call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(2,iboun),lnodx(1)) ! Edge 1-2
                call mesh_multiplication_edge_node(lnodb_old(2,iboun),lnodb_old(3,iboun),lnodx(2)) ! Edge 2-3
                call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(3,iboun),lnodx(3)) ! Edge 1-3

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(1,iboun)
                lnodb(2,kboun) = lnodx(1)
                lnodb(3,kboun) = lnodx(3)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodx(1)
                lnodb(2,kboun) = lnodx(2)
                lnodb(3,kboun) = lnodx(3)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodx(1)
                lnodb(2,kboun) = lnodb_old(2,iboun)
                lnodb(3,kboun) = lnodx(2)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodx(2)
                lnodb(2,kboun) = lnodb_old(3,iboun)
                lnodb(3,kboun) = lnodx(3)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

             case ( TRI06 ) 
                !
                ! TRI06
                !
                call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(4,iboun),lnodx(1)) ! Edge 1-4
                call mesh_multiplication_edge_node(lnodb_old(4,iboun),lnodb_old(2,iboun),lnodx(2)) ! Edge 4-2
                call mesh_multiplication_edge_node(lnodb_old(2,iboun),lnodb_old(5,iboun),lnodx(3)) ! Edge 2-5
                call mesh_multiplication_edge_node(lnodb_old(5,iboun),lnodb_old(3,iboun),lnodx(4)) ! Edge 5-3
                call mesh_multiplication_edge_node(lnodb_old(3,iboun),lnodb_old(6,iboun),lnodx(5)) ! Edge 3-6
                call mesh_multiplication_edge_node(lnodb_old(6,iboun),lnodb_old(1,iboun),lnodx(6)) ! Edge 6-1

                call mesh_multiplication_face_node(lnodb_old(1:pnodb,iboun),kpoin,ifacg) 
                call mesh_multiplication_subface(1_ip,ifacg,pblty,lnodb_old(1:pnodb,iboun))        ! Faces: 6-4,4-5,6-5

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(1,iboun)
                lnodb(2,kboun) = lnodb_old(4,iboun) 
                lnodb(3,kboun) = lnodb_old(6,iboun) 
                lnodb(4,kboun) = lnodx(1)
                lnodb(5,kboun) = face_node(1,1)
                lnodb(6,kboun) = lnodx(6)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(4,iboun)
                lnodb(2,kboun) = lnodb_old(5,iboun) 
                lnodb(3,kboun) = lnodb_old(6,iboun) 
                lnodb(4,kboun) = face_node(2,1) 
                lnodb(5,kboun) = face_node(3,1)
                lnodb(6,kboun) = face_node(1,1)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(4,iboun)
                lnodb(2,kboun) = lnodb_old(2,iboun) 
                lnodb(3,kboun) = lnodb_old(5,iboun) 
                lnodb(4,kboun) = lnodx(2)
                lnodb(5,kboun) = lnodx(3)
                lnodb(6,kboun) = face_node(2,1)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)                            

                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(6,iboun)
                lnodb(2,kboun) = lnodb_old(5,iboun) 
                lnodb(3,kboun) = lnodb_old(3,iboun) 
                lnodb(4,kboun) = face_node(3,1) 
                lnodb(5,kboun) = lnodx(4)
                lnodb(6,kboun) = lnodx(5)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

             case ( QUA04 )
                !
                ! QUA04
                !
                ielem = lelbo_old(iboun)

                if( divide_boun(iboun) ) then
                   !
                   !   4           3
                   !   o-----*-----o
                   !   |     |     |
                   !   |     |     |
                   !   *-----*-----*  Boundary division Type 1
                   !   |     |     |
                   !   |     |     |
                   !   o-----*-----o
                   !   1           2
                   !
                   ipoi1 = lnodb_old(1,iboun)
                   ipoi2 = lnodb_old(2,iboun)
                   ipoi3 = lnodb_old(3,iboun)
                   ipoi4 = lnodb_old(4,iboun)
                   call mesh_multiplication_edge_node(ipoi1,ipoi2,lnodx(1)) ! Edge 1-2
                   call mesh_multiplication_edge_node(ipoi2,ipoi3,lnodx(2)) ! Edge 2-3
                   call mesh_multiplication_edge_node(ipoi3,ipoi4,lnodx(3)) ! Edge 3-4
                   call mesh_multiplication_edge_node(ipoi4,ipoi1,lnodx(4)) ! Edge 4-1
                   call mesh_multiplication_face_node([ipoi1,ipoi2,ipoi3,ipoi4],lnodx(5))

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodb_old(1,iboun)
                   lnodb(2,kboun) = lnodx(1)
                   lnodb(3,kboun) = lnodx(5)
                   lnodb(4,kboun) = lnodx(4)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodx(1)
                   lnodb(2,kboun) = lnodb_old(2,iboun)
                   lnodb(3,kboun) = lnodx(2)
                   lnodb(4,kboun) = lnodx(5)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodx(5)
                   lnodb(2,kboun) = lnodx(2)
                   lnodb(3,kboun) = lnodb_old(3,iboun)
                   lnodb(4,kboun) = lnodx(3)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                   kboun          = kboun + 1
                   lnodb(1,kboun) = lnodx(4)
                   lnodb(2,kboun) = lnodx(5)
                   lnodb(3,kboun) = lnodx(3)
                   lnodb(4,kboun) = lnodb_old(4,iboun)
                   ltypb(kboun)   = ltypb_old(iboun)
                   lboch(kboun)   = lboch_old(iboun)

                else
                   !                  ^
                   !   4           3  |
                   !   o-----*-----o  |
                   !   |     |     |  | Element normal direction
                   !   |     |     |  |
                   !   |     |     |    Boundary division Type 2
                   !   |     |     |
                   !   |     |     |
                   !   o-----*-----o
                   !   1           2
                   !
                   !
                   ! Correction boundary ordering required for kfl_elndi = 1
                   !
                   do ii = 1,4
                      !
                      ! Boundary points ordering
                      !
                      node1 = list_inodb_QUA04(ii,1)
                      node2 = list_inodb_QUA04(ii,2)
                      node3 = list_inodb_QUA04(ii,3)
                      node4 = list_inodb_QUA04(ii,4)

                      ipoi1 = lnodb_old(node1,iboun) ! BOT
                      ipoi2 = lnodb_old(node2,iboun) ! BOT
                      ipoi3 = lnodb_old(node3,iboun) ! TOP
                      ipoi4 = lnodb_old(node4,iboun) ! TOP

                      if(  any(lnods_old(1:4,ielem) == ipoi1) .and. any(lnods_old(1:4,ielem) == ipoi2) .and. &
                           any(lnods_old(5:8,ielem) == ipoi3) .and. any(lnods_old(5:8,ielem) == ipoi4) ) then

                         call mesh_multiplication_edge_node(ipoi1,ipoi2,lnodx(1)) ! Edge 1-2
                         call mesh_multiplication_edge_node(ipoi3,ipoi4,lnodx(2)) ! Edge 3-4

                         kboun          = kboun + 1
                         lnodb(1,kboun) = lnodb_old(node1,iboun)
                         lnodb(2,kboun) = lnodx(1)
                         lnodb(3,kboun) = lnodx(2)
                         lnodb(4,kboun) = lnodb_old(node4,iboun)
                         ltypb(kboun)   = ltypb_old(iboun)
                         lboch(kboun)   = lboch_old(iboun)

                         kboun          = kboun + 1
                         lnodb(1,kboun) = lnodx(1)
                         lnodb(2,kboun) = lnodb_old(node2,iboun)
                         lnodb(3,kboun) = lnodb_old(node3,iboun)
                         lnodb(4,kboun) = lnodx(2)
                         ltypb(kboun)   = ltypb_old(iboun)
                         lboch(kboun)   = lboch_old(iboun)

                      end if

                   end do

                end if

             case ( QUA09 ) 
                !
                ! QUA09
                !
                !   4           3
                !   o-----o-----o
                !   |     |     |
                !   |     |     |
                !   o-----o-----o 
                !   |     |     |
                !   |     |     |
                !   o-----o-----o
                !   1           2
                !
                call mesh_multiplication_edge_node(lnodb_old(1,iboun),lnodb_old(5,iboun),lnodx(1)) ! Edge 1-5
                call mesh_multiplication_edge_node(lnodb_old(5,iboun),lnodb_old(2,iboun),lnodx(2)) ! Edge 5-2
                call mesh_multiplication_edge_node(lnodb_old(2,iboun),lnodb_old(6,iboun),lnodx(3)) ! Edge 2-6
                call mesh_multiplication_edge_node(lnodb_old(6,iboun),lnodb_old(3,iboun),lnodx(4)) ! Edge 6-3
                call mesh_multiplication_edge_node(lnodb_old(3,iboun),lnodb_old(7,iboun),lnodx(5)) ! Edge 3-7
                call mesh_multiplication_edge_node(lnodb_old(7,iboun),lnodb_old(4,iboun),lnodx(6)) ! Edge 7-4
                call mesh_multiplication_edge_node(lnodb_old(4,iboun),lnodb_old(8,iboun),lnodx(7)) ! Edge 4-8
                call mesh_multiplication_edge_node(lnodb_old(8,iboun),lnodb_old(1,iboun),lnodx(8)) ! Edge 8-1
                call mesh_multiplication_face_node(lnodb_old(1:9,iboun),kpoin)
                ! 1
                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(1,iboun)
                lnodb(2,kboun) = lnodb_old(5,iboun)
                lnodb(3,kboun) = lnodb_old(9,iboun)
                lnodb(4,kboun) = lnodb_old(8,iboun)
                lnodb(5,kboun) = lnodx(1)
                lnodb(6,kboun) = kpoin + 1
                lnodb(7,kboun) = kpoin + 3
                lnodb(8,kboun) = lnodx(8)
                lnodb(9,kboun) = kpoin 
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)
                ! 2
                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(5,iboun)
                lnodb(2,kboun) = lnodb_old(2,iboun)
                lnodb(3,kboun) = lnodb_old(6,iboun)
                lnodb(4,kboun) = lnodb_old(9,iboun)
                lnodb(5,kboun) = lnodx(2)
                lnodb(6,kboun) = lnodx(3)
                lnodb(7,kboun) = kpoin + 4
                lnodb(8,kboun) = kpoin + 1
                lnodb(9,kboun) = kpoin + 2
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)
                ! 3
                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(8,iboun)
                lnodb(2,kboun) = lnodb_old(9,iboun)
                lnodb(3,kboun) = lnodb_old(7,iboun)
                lnodb(4,kboun) = lnodb_old(4,iboun)
                lnodb(5,kboun) = kpoin + 3
                lnodb(6,kboun) = kpoin + 6
                lnodb(7,kboun) = lnodx(6)
                lnodb(8,kboun) = lnodx(7)
                lnodb(9,kboun) = kpoin + 5
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)
                ! 4
                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(9,iboun)
                lnodb(2,kboun) = lnodb_old(6,iboun)
                lnodb(3,kboun) = lnodb_old(3,iboun)
                lnodb(4,kboun) = lnodb_old(7,iboun)
                lnodb(5,kboun) = kpoin + 4
                lnodb(6,kboun) = lnodx(4)
                lnodb(7,kboun) = lnodx(5)
                lnodb(8,kboun) = kpoin + 6
                lnodb(9,kboun) = kpoin + 7
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)               

             case ( POINT ) 
                !
                ! POINT
                !
                kboun          = kboun + 1
                lnodb(1,kboun) = lnodb_old(1,iboun)
                ltypb(kboun)   = ltypb_old(iboun)
                lboch(kboun)   = lboch_old(iboun)

             end select

          end do
       end if
       if( if_high_order ) then
          call memory_deallo(memor_dom,'FACE_DIVIDE','renpoi',face_divide)
       end if
       !
       ! Recompute LELBO: KBOUN => IBOUN => IELEM => Look for kelem
       !
       ! +-------+-------+      +---------------+
       ! |       |       |      |               |
       ! | kelem | kelem |      |               |
       ! |       |       |      |               |
       ! +-------+-------+  <=  |     ielem     |
       ! |       |       |      |               |
       ! | kelem | kelem |      |               |
       ! |       |       |      |               |
       ! +-------+-------+      +---------------+
       !   kboun   kboun              iboun
       !
       do kboun = 1,nboun
          iboun = perm_boun(kboun)
          ielem = lelbo_old(iboun)
          if( ielem > 0 ) then
             loop_kelem: do kelem = 1,nelem
                pnodb = nnode(abs(ltypb(kboun)))
                if( perm_elem(kelem) == ielem ) then
                   inodb = 0
                   pelty = abs(ltype(kelem))
                   do knode = 1,element_type(pelty) % number_nodes
                      kpoin = lnods(knode,kelem)
                      loop_knodb: do knodb = 1,pnodb
                         if( kpoin == lnodb(knodb,kboun) ) then
                            inodb = inodb + 1
                            if( inodb == pnodb ) then
                               lelbo(kboun) = kelem
                               exit loop_kelem
                            else
                               exit loop_knodb
                            end if
                         end if
                      end do loop_knodb
                   end do
                end if
             end do loop_kelem
          end if
       end do

       call memory_deallo(memor_dom,'PERM_ELEM','mesh_multiplication_divide_elements',perm_elem)
       call memory_deallo(memor_dom,'PERM_BOUN','mesh_multiplication_divide_elements',perm_boun)
       !
       ! Deallocate edge memory
       !
       call memory_deallo(memor_dom,'PEDGP','mesh_multiplication_divide_elements',pedgp)

       !--------------------------------------------------------------------
       !
       ! Node characteristic: hard to guess!!!
       !
       !--------------------------------------------------------------------

       !--------------------------------------------------------------------
       !
       ! LNINV_LOC
       !
       !--------------------------------------------------------------------

       if( kfl_mmpar == 2 ) then
          call mesh_multiplication_number_new_nodes(npoin_old)
       else
          call par_number_mesh(1_ip)
       end if

       !--------------------------------------------------------------------
       !
       ! Groups
       !
       !--------------------------------------------------------------------

       if( ngrou_dom > 0 .and. associated(lgrou_dom) ) then

          call memory_copy(memor_dom,'LGROU_DOM_OLD','mesh_multiplication_divide_elements',lgrou_dom,lgrou_dom_old,'DO_NOT_DEALLOCATE')
          call domain_memory_deallocate('LGROU_DOM')
          call domain_memory_allocate  ('LGROU_DOM')

          ! Old nodes
          do ipoin = 1,npoin_old
             lgrou_dom(ipoin) = lgrou_dom_old(ipoin)
          end do
          ! Edge nodes
          ipoin = npoin_old
          do iedgg = 1,nedgg
             if( divide_edge(iedgg) ) then
                ipoin = ipoin + 1
                node1 = ledgg(1,iedgg)
                node2 = ledgg(2,iedgg)
                lgrou_dom(ipoin) = min(lgrou_dom_old(node1),lgrou_dom_old(node2))
             end if
          end do
          ! Face nodes
          do ifacg = 1,nfacg
             if( divide_face(ifacg) ) then
                ipoin = ipoin + 1
                node1 = lfacg(1,ifacg)
                node2 = lfacg(2,ifacg)
                node3 = lfacg(3,ifacg)
                node4 = lfacg(4,ifacg)
                lgrou_dom(ipoin) = min(lgrou_dom_old(node1),lgrou_dom_old(node2),&
                     &                 lgrou_dom_old(node3),lgrou_dom_old(node4))
             end if
          end do
          ! Central nodes
          if( lnuty(QUA04) > 0 .or. lnuty(HEX08) > 0 ) then
             call memgen(1_ip,ngrou_dom,0_ip)
             do ielem = 1,nelem_old
                pelty = abs(ltype_old(ielem))
                if( ( pelty == QUA04 .or. pelty == HEX08 ) .and. lelch_old(ielem) /= ELINT .and. kfl_elndi == 0 ) then
                   pnode = nnode(pelty)
                   ipoin = ipoin + 1
                   do inode = 1,pnode
                      igrou = lgrou_dom_old(lnods_old(inode,ielem))
                      gisca(igrou) = gisca(igrou) + 1
                   end do
                   jgrou = 1
                   do igrou = 2,ngrou_dom
                      if( gisca(igrou) > gisca(jgrou) ) jgrou = igrou
                   end do
                   lgrou_dom(ipoin) = jgrou
                   do inode = 1,pnode
                      igrou = lgrou_dom_old(lnods_old(inode,ielem))
                      gisca(igrou) = 0
                   end do
                end if
             end do
             call memgen(3_ip,ngrou_dom,0_ip)
          end if

          call memory_deallo(memor_dom,'LGROU_DOM_OLD','mesh_multiplication_divide_elements',lgrou_dom_old)

       end if

       !--------------------------------------------------------------------
       !
       ! Sets: LESET, LBSET, LNSET
       !
       !--------------------------------------------------------------------

       if( neset > 0 ) then
          call memory_copy(memor_dom,'LESET_OLD','mesh_multiplication_divide_elements',leset,leset_old,'DO_NOT_DEALLOCATE')
       end if
       if( nbset > 0 ) then
          call memory_copy(memor_dom,'LBSET_OLD','mesh_multiplication_divide_elements',lbset,lbset_old,'DO_NOT_DEALLOCATE')
       end if
       if( nnset > 0 ) then
          call memory_copy(memor_dom,'LNSET_OLD','mesh_multiplication_divide_elements',lnset,lnset_old,'DO_NOT_DEALLOCATE')
       end if
       call domain_memory_reallocate('LESET')
       call domain_memory_reallocate('LBSET')
       call domain_memory_reallocate('LNSET')
       kelem = 0
       kboun = 0
       if( ndime == 2 ) then
          if( neset > 0 ) then
             do ielem = 1,nelem_old
                neldi = howdi(abs(ltype_old(ielem)))
                if( lelch_old(ielem) == ELINT .or. kfl_elndi == 1_ip ) neldi = howdi_ELINT
                do ii = 1,neldi
                   kelem = kelem + 1
                   leset(kelem) = leset_old(ielem)
                end do
             end do
          end if
          if( nbset > 0 ) then
             do iboun = 1,nboun_old
                nbodi = howdi(abs(ltypb_old(iboun)))
                if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
                do ii = 1,nbodi
                   kboun = kboun + 1
                   lbset(kboun) = lbset_old(iboun)
                end do
             end do
          end if
       else
          if( neset > 0 ) then
             do ielem = 1,nelem_old
                neldi = howdi(abs(ltype_old(ielem)))
                if( lelch_old(ielem) == ELINT .or. kfl_elndi == 1_ip ) neldi = howdi_ELINT
                do ii = 1,neldi
                   kelem = kelem + 1
                   leset(kelem) = leset_old(ielem)
                end do
             end do
          end if
          if( nbset > 0 ) then
             do iboun = 1,nboun_old
                nbodi = howdi(abs(ltypb_old(iboun)))
                if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
                do ii = 1,nbodi
                   kboun = kboun + 1
                   lbset(kboun) = lbset_old(iboun)
                end do
             end do
          end if
       end if
       if( nnset > 0 ) then
          do ipoin = 1,npoin_old
             lnset(ipoin) = lnset_old(ipoin)
          end do
       end if

       if( neset > 0 ) call memory_deallo(memor_dom,'LESET_OLD','mesh_multiplication_divide_elements',leset_old)
       if( nbset > 0 ) call memory_deallo(memor_dom,'LBSET_OLD','mesh_multiplication_divide_elements',lbset_old)
       if( nnset > 0 ) call memory_deallo(memor_dom,'LBSET_OLD','mesh_multiplication_divide_elements',lnset_old)

       !--------------------------------------------------------------------
       !
       ! Fields: XFIEL (node/element/boundary)
       !
       !--------------------------------------------------------------------

       if( nfiel > 0 ) then
          call memory_alloca(memor_dom,'XFIEL_OLD','mesh_multiplication_divide_elements',xfiel_old,nfiel)
          do ifiel = 1,nfiel
             if( kfl_field(1,ifiel) > 0 ) then

                if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then
                   nsteps = kfl_field(4,ifiel)
                else
                   nsteps = nsteps_fiel_ondemand
                end if

                if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                   call memory_alloca(memor_dom,'XFIEL_OLD % A','mesh_multiplication_divide_elements',xfiel_old(ifiel) % a,kfl_field(1,ifiel),npoin_old,nsteps,'DO_NOT_INITIALIZE')
                   do istep = 1,nsteps
                      do ipoin = 1,npoin_old
                         do idime = 1,kfl_field(1,ifiel)
                            xfiel_old(ifiel) % a(idime,ipoin,istep) = xfiel(ifiel) % a(idime,ipoin,istep)
                         end do
                      end do
                   end do

                else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                   call memory_alloca(memor_dom,'XFIEL_OLD % A','mesh_multiplication_divide_elements',xfiel_old(ifiel) % a,kfl_field(1,ifiel),nboun_old,nsteps,'DO_NOT_INITIALIZE')
                   do istep = 1,nsteps
                      do iboun = 1,nboun_old
                         do idime = 1,kfl_field(1,ifiel)
                            xfiel_old(ifiel) % a(idime,iboun,istep) = xfiel(ifiel) % a(idime,iboun,istep)
                         end do
                      end do
                   end do

                else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                   call memory_alloca(memor_dom,'XFIEL_OLD % A','mesh_multiplication_divide_elements',xfiel_old(ifiel) % a,kfl_field(1,ifiel),nelem_old,nsteps,'DO_NOT_INITIALIZE')
                   do istep = 1,nsteps
                      do ielem = 1,nelem_old
                         do idime = 1,kfl_field(1,ifiel)
                            xfiel_old(ifiel) % a(idime,ielem,istep) = xfiel(ifiel) % a(idime,ielem,istep)
                         end do
                      end do
                   end do
                end if
                call domain_memory_reallocate('XFIEL % A',NUMBER1=ifiel)
             end if
          end do
          do ifiel = 1,nfiel
             if( kfl_field(1,ifiel) > 0 ) then
                if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                   ! Old nodes
                   do istep = 1,nsteps
                      do ipoin = 1,npoin_old
                         do idime = 1,kfl_field(1,ifiel)
                            xfiel(ifiel) % a(idime,ipoin,istep) = xfiel_old(ifiel) % a(idime,ipoin,istep)
                         end do
                      end do
                   end do
                   ! Edge nodes
                   ipoin = npoin_old
                   !if( kfl_field(1,ifiel) == 0 ) then
                   do iedgg = 1,nedgg
                      if( divide_edge(iedgg) ) then
                         ipoin = ipoin + 1
                         node1 = ledgg(1,iedgg)
                         node2 = ledgg(2,iedgg)
                         do istep = 1,nsteps
                            do idime = 1,kfl_field(1,ifiel)
                               xfiel(ifiel) % a(idime,ipoin,istep) = 0.5_rp&
                                    *(xfiel_old(ifiel) % a(idime,node1,istep) + xfiel_old(ifiel) % a(idime,node2,istep) )
                            end do
                         end do
                      end if
                   end do
                   !else
                   !   do iedgg = 1,nedgg
                   !      ipoin = ipoin + 1
                   !      node1 = ledgg(1,iedgg)
                   !      node2 = ledgg(2,iedgg)
                   !      do idime = 1,kfl_field(1,ifiel)*nsteps
                   !         knode = 0
                   !         do inode = 1,2
                   !            xx(inode) = xfiel_old(ifiel) % a(idime,ledgg(inode,iedgg))
                   !            if( xx(inode) > 0.0_rp ) then
                   !               knode = knode + 1
                   !            end if
                   !         end do
                   !         if( knode == 2 ) &
                   !              xfiel(ifiel) % a(idime,ipoin) = 0.5_rp * ( xx(1) + xx(2) )
                   !      end do
                   !   end do
                   !end if
                   ! Face nodes
                   !if( kfl_field(1,ifiel) == 0 ) then
                   do ifacg = 1,nfacg
                      if( divide_face(ifacg) ) then
                         ipoin = ipoin + 1
                         node1 = lfacg(1,ifacg)
                         node2 = lfacg(2,ifacg)
                         node3 = lfacg(3,ifacg)
                         node4 = lfacg(4,ifacg)
                         do istep = 1,nsteps
                            do idime = 1,kfl_field(1,ifiel)
                               xfiel(ifiel) % a(idime,ipoin,istep) = 0.25_rp&
                                    &  * (   xfiel_old(ifiel) % a(idime,node1,istep) + xfiel_old(ifiel) % a(idime,node2,istep) &
                                    &      + xfiel_old(ifiel) % a(idime,node3,istep) + xfiel_old(ifiel) % a(idime,node4,istep) )
                            end do
                         end do
                      end if
                   end do
                   !else
                   !   do ifacg = 1,nfacg
                   !      ipoin = ipoin + 1
                   !      do idime = 1,kfl_field(1,ifiel)*nsteps
                   !         knode = 0
                   !         do inode = 1,4
                   !            xx(inode) = xfiel_old(ifiel) % a(idime,lfacg(inode,ifacg))
                   !            if( xx(inode) > 0.0_rp ) then
                   !               knode = knode + 1
                   !            end if
                   !         end do
                   !         if( knode == 4 ) &
                   !              xfiel(ifiel) % a(idime,ipoin) = 0.25_rp * ( xx(1) + xx(2) + xx(3) + xx(4) )
                   !      end do
                   !   end do
                   !end if
                   ! Central nodes
                   if( lnuty(QUA04) > 0 .or. lnuty(HEX08) > 0 ) then
                      do ielem = 1,nelem_old
                         pelty = abs(ltype_old(ielem))
                         if( ( pelty == QUA04 .or. pelty == HEX08 ) .and. lelch_old(ielem) /= ELINT .and. kfl_elndi == 0 ) then
                            ipoin = ipoin + 1
                            pnode = nnode(pelty)
                            do istep = 1,nsteps
                               do idime = 1,kfl_field(1,ifiel)
                                  xfiel(ifiel) % a(idime,ipoin,istep) = 0.0_rp
                               end do
                            end do
                            !if( kfl_field(1,ifiel) == 0 ) then
                            do istep = 1,nsteps
                               do inode = 1,pnode
                                  jpoin = lnods_old(inode,ielem)
                                  do idime = 1,kfl_field(1,ifiel)
                                     xfiel(ifiel) % a(idime,ipoin,istep) = xfiel(ifiel) % a(idime,ipoin,istep) &
                                          & + xfiel(ifiel) % a(idime,jpoin,istep)
                                  end do
                               end do
                            end do
                            dummr = 1.0_rp / real(pnode,rp)
                            do istep = 1,nsteps
                               do idime = 1,kfl_field(1,ifiel)
                                  xfiel(ifiel) % a(idime,ipoin,istep) = dummr * xfiel(ifiel) % a(idime,ipoin,istep)
                               end do
                            end do
                            !else
                            !   dummr = 1.0_rp / real(pnode)
                            !   do inode = 1,pnode
                            !      jpoin = lnods_old(inode,ielem)
                            !      do idime = 1,kfl_field(1,ifiel)*nsteps
                            !         xfiel(ifiel) % a(idime,ipoin) = xfiel(ifiel) % a(idime,ipoin) &
                            !              & + xfiel(ifiel) % a(idime,jpoin)
                            !      end do
                            !   end do
                            !   do idime = 1,kfl_field(1,ifiel)*nsteps
                            !
                            !      knode = 0
                            !      do inode = 1,pnode
                            !         xx(inode) = xfiel(ifiel) % a(idime,jpoin)
                            !         if( xx(inode) > 0.0_rp ) then
                            !            knode = knode + 1
                            !         end if
                            !      end do
                            !      if( knode == pnode ) then
                            !         xfiel(ifiel) % a(idime,ipoin) = dummr * xfiel(ifiel) % a(idime,ipoin)
                            !      else
                            !         xfiel(ifiel) % a(idime,ipoin) = 0.0_rp
                            !      end if
                            !   end do
                            !end if
                         end if
                      end do
                   end if

                else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                   kboun = 0
                   do iboun = 1,nboun_old
                      nbodi = howdi(abs(ltypb_old(iboun)))
                      if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
                      do ii = 1,nbodi
                         kboun = kboun + 1
                         do istep = 1,nsteps
                            do idime = 1,kfl_field(1,ifiel)
                               xfiel(ifiel) % a(idime,kboun,istep) = xfiel_old(ifiel) % a(idime,iboun,istep)
                            end do
                         end do
                      end do
                   end do

                else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                   kelem = 0
                   do ielem = 1,nelem_old
                      neldi = howdi(abs(ltype_old(ielem)))
                      if( lelch_old(ielem) == ELINT .or. kfl_elndi == 1_ip ) neldi = howdi_ELINT
                      do ii = 1,neldi
                         kelem = kelem + 1
                         do istep = 1,nsteps
                            do idime = 1,kfl_field(1,ifiel)
                               xfiel(ifiel) % a(idime,kelem,istep) = xfiel_old(ifiel) % a(idime,ielem,istep)
                            end do
                         end do
                      end do
                   end do
                else
                   call runend('SUBMSH: DIVISOR WITH ELEMENT FIELD NOT CODED')
                end if
                call memory_deallo(memor_dom,'XFIEL_OLD % A','mesh_multiplication_divide_elements',xfiel_old(ifiel) % a)
             end if
          end do
          call memory_deallo(memor_dom,'XFIEL_OLD','mesh_multiplication_divide_elements',xfiel_old)
       end if

       !-|CURVED| GEOMETRY FIELD MULTIPLICATION
       if(multiply_with_curvature.eq.1_ip) then
          beginDimCurvature = 1_ip
          numFieldCurvaturePerDimension = maxNumNodesGeometry!size(xfiel(curvatureField)%a,1,KIND=ip)/ndime
          beginDimCurvature = 1_ip
          !        do ielem = 1_ip,nelem
          !          if(geometry_order(ielem).gt.1_ip) then
          !            do idime = 1_ip,ndime
          !              endDimCurvature = beginDimCurvature+numFieldCurvaturePerDimension-1_ip
          !              xfiel(curvatureField)%a(beginDimCurvature:endDimCurvature,ielem) = subdivided_curved_geometry(:,idime,ielem)
          !              beginDimCurvature = endDimCurvature+1_ip
          !            end do
          !          end if
          !        end do

          do idime = 1_ip,ndime
             endDimCurvature = beginDimCurvature+numFieldCurvaturePerDimension-1_ip
             !          print*,'nelem',nelem,'  kelem',kelem, ' nelem_old',nelem_old
             !          print*,size(xfiel(curvatureField)%a(beginDimCurvature:endDimCurvature,:),1,KIND=ip),&
             !           size(xfiel(curvatureField)%a(beginDimCurvature:endDimCurvature,:),2,KIND=ip)
             !          print*,size(subdivided_curved_geometry(:,idime,:),1,KIND=ip),&
             !           size(subdivided_curved_geometry(:,idime,:),2,KIND=ip)
             xfiel(curvatureField)%a(beginDimCurvature:endDimCurvature,:,1) = subdivided_curved_geometry(:,idime,1:nelem)
             beginDimCurvature = endDimCurvature+1_ip
          end do

          !        do idime = 1,ndime
          !          do iaux = 1,numFieldCurvaturePerDimension
          !            xfiel(curvatureField)%a(iaux+((idime-1)*numFieldCurvaturePerDimension),:) = subdivided_curved_geometry(iaux,idime,:)
          !          end do
          !        end do
          deallocate(curved_geometry)
          deallocate(subdivided_curved_geometry)
       end if

       !--------------------------------------------------------------------
       !
       ! LNINV_LOC: put 0 temporarily on new nodes. LNINV_LOC will be
       ! later on calculated on par_submsh
       !
       !--------------------------------------------------------------------

       if( associated(lninv_loc) ) then
          call memory_resize(par_memor,'LNINV_LOC','submsh',lninv_loc,npoin)
       end if

       !--------------------------------------------------------------------
       !
       ! Boundary conditions
       !
       !--------------------------------------------------------------------

       !
       ! KFL_CODNO
       !
       if( kfl_icodn > 0 ) then
          allocate( kfl_codno_old(mcono,npoin_old) )
          do ipoin = 1,npoin_old
             do ii = 1,mcono
                kfl_codno_old(ii,ipoin) = kfl_codno(ii,ipoin)
             end do
          end do
          call domain_memory_reallocate('KFL_CODNO') ! KFL_CODNO
          ! Old nodes
          do ipoin = 1,npoin_old
             do ii = 1,mcono
                kfl_codno(ii,ipoin) = kfl_codno_old(ii,ipoin)
             end do
          end do
          mcod1 = mcodb + 1
          do ipoin = npoin_old+1,npoin
             do ii = 1,mcono
                kfl_codno(ii,ipoin) = mcod1
             end do
          end do
          ! Edge nodes
          ipoin = npoin_old
          do iedgg = 1,nedgg
             if( divide_edge(iedgg) ) then
                ipoin = ipoin + 1
                node1 = ledgg(1,iedgg)
                node2 = ledgg(2,iedgg)
                do ii = 1,1
                   code1 = kfl_codno_old(ii,node1)
                   code2 = kfl_codno_old(ii,node2)
                   if( code1 /= mcod1 .and. code2 /= mcod1 ) then
                      ! <GGU> ESTO SE HA COMENTADO EN LOS ULTIMOS COMMITS DEL TRUNK
                      !kfl_codno(ii,ipoin) = min(kfl_codno_old(ii,node1),kfl_codno_old(ii,node2))
                   end if
                end do
             end if
          end do
          ! Face nodes
          do ifacg = 1,nfacg
             if( divide_face(ifacg) ) then
                ipoin = ipoin + 1
                node1 = lfacg(1,ifacg)
                node2 = lfacg(2,ifacg)
                node3 = lfacg(3,ifacg)
                node4 = lfacg(4,ifacg)
                do ii = 1,1
                   code1 = kfl_codno_old(ii,node1)
                   code2 = kfl_codno_old(ii,node2)
                   code3 = kfl_codno_old(ii,node3)
                   code4 = kfl_codno_old(ii,node4)
                   if( code1 /= mcod1 .and. code2 /= mcod1 .and. code3 /= mcod1 .and. code4 /= mcod1 ) then
                      ! <GGU> ESTO SE HA COMENTADO EN LOS ULTIMOS COMMITS DEL TRUNK
                      !kfl_codno(ii,ipoin) = min(kfl_codno_old(ii,node1),kfl_codno_old(ii,node2),&
                      !     &                    kfl_codno_old(ii,node3),kfl_codno_old(ii,node4))
                   end if
                end do
             end if
          end do
          ! Central nodes
          if( lnuty(QUA04) > 0 .or. lnuty(HEX08) > 0 ) then
             do ielem = 1,nelem_old
                pelty = abs(ltype_old(ielem))
                if( ( pelty == QUA04 .or. pelty == HEX08 ) .and. lelch_old(ielem) /= ELINT .and. kfl_elndi == 0 ) then
                   ipoin = ipoin + 1
                   do ii = 1,1
                      kfl_codno(ii,ipoin) = mcodb + 1
                   end do
                end if
             end do
          end if

          deallocate( kfl_codno_old )

       end if
       !
       ! KFL_CODBO
       !
       if( kfl_icodb > 0 ) then
          allocate( kfl_codbo_old(nboun_old) )
          do iboun = 1,nboun_old
             kfl_codbo_old(iboun) = kfl_codbo(iboun)
          end do
          call domain_memory_reallocate('KFL_CODBO') ! KFL_CODBO
          kboun = 0
          if(      ndime == 1 ) then
             do iboun = 1,nboun_old
                kboun = kboun + 1
                kfl_codbo(kboun) = kfl_codbo_old(iboun)
             end do
          else if( ndime == 2 ) then
             do iboun = 1,nboun_old
                nbodi = 2_ip
                if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
                do ii = 1,nbodi
                   kboun = kboun + 1
                   kfl_codbo(kboun) = kfl_codbo_old(iboun)
                end do
             end do
          else
             do iboun = 1,nboun_old
                pblty = abs(ltypb_old(iboun))
                nbodi = howdi(pblty) !4_ip
                if( divide_boun(iboun) .eqv. .false. ) nbodi = howbo_ELINT
                do ii = 1,nbodi
                   kboun = kboun + 1
                   kfl_codbo(kboun) = kfl_codbo_old(iboun)
                end do
             end do
          end if
          deallocate( kfl_codbo_old )
       end if

       !--------------------------------------------------------------------
       !
       ! Interpolation
       !
       !--------------------------------------------------------------------

       if( 1 == 2 ) then

          linno => meshe(ndivi) % linno
          ! Old nodes
          do ipoin = 1,npoin_old
             linno(ipoin) % n = 1
             allocate( linno(ipoin) % l(linno(ipoin) % n) )
             linno(ipoin) % l(1) = ipoin
          end do
          ! Edge nodes
          ipoin = npoin_old
          do iedgg = 1,nedgg
             if( divide_edge(iedgg) ) then
                ipoin = ipoin + 1
                linno(ipoin) % n = 2
                allocate( linno(ipoin) % l(linno(ipoin) % n) )
                linno(ipoin) % l(1) = ledgg(1,iedgg)
                linno(ipoin) % l(2) = ledgg(2,iedgg)
             end if
          end do
          ! Face nodes
          do ifacg = 1,nfacg
             if( divide_face(ifacg) ) then
                ipoin = ipoin + 1
                linno(ipoin) % n = 4
                allocate( linno(ipoin) % l(linno(ipoin) % n) )
                linno(ipoin) % l(1) = lfacg(1,ifacg)
                linno(ipoin) % l(2) = lfacg(2,ifacg)
                linno(ipoin) % l(3) = lfacg(3,ifacg)
                linno(ipoin) % l(4) = lfacg(4,ifacg)
             end if
          end do
          ! Central nodes
          if( lnuty(QUA04) > 0 .or. lnuty(HEX08) > 0 ) then
             do ielem = 1,nelem_old
                pelty = abs(ltype_old(ielem))
                if( pelty == QUA04 .or. pelty == HEX08 .and. lelch_old(ielem) /= ELINT .and. kfl_elndi == 0 ) then
                   pnode = nnode(pelty)
                   linno(ipoin) % n = pnode
                   allocate( linno(ipoin) % l(linno(ipoin) % n) )
                   ipoin = ipoin + 1
                   do inode = 1,pnode
                      linno(ipoin) % l(inode) = lnods_old(inode,ielem)
                   end do
                end if
             end do
          end if
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_dom,'LMAST_OLD','mesh_multiplication_divide_elements',lmast_old)
       call memory_deallo(memor_dom,'LNOCH_OLD','mesh_multiplication_divide_elements',lnoch_old)
       call memory_deallo(memor_dom,'COORD_OLD','mesh_multiplication_divide_elements',coord_old)

       call memory_deallo(memor_dom,'LBOCH_OLD','mesh_multiplication_divide_elements',lboch_old)
       call memory_deallo(memor_dom,'LTYPB_OLD','mesh_multiplication_divide_elements',ltypb_old)
       call memory_deallo(memor_dom,'LELBO_OLD','mesh_multiplication_divide_elements',lelbo_old)
       call memory_deallo(memor_dom,'LNODB_OLD','mesh_multiplication_divide_elements',lnodb_old)

       call memory_deallo(memor_dom,'LMATE_OLD','mesh_multiplication_divide_elements',lmate_old)
       call memory_deallo(memor_dom,'LESUB_OLD','mesh_multiplication_divide_elements',lesub_old)
       call memory_deallo(memor_dom,'LNNOD_OLD','mesh_multiplication_divide_elements',lnnod_old)
       call memory_deallo(memor_dom,'LELCH_OLD','mesh_multiplication_divide_elements',lelch_old)
       call memory_deallo(memor_dom,'LTYPE_OLD','mesh_multiplication_divide_elements',ltype_old)
       call memory_deallo(memor_dom,'LNODS_OLD','mesh_multiplication_divide_elements',lnods_old)

       call memory_deallo(memor_dom,'DIVIDE_EDGE','mesh_multiplication_divide_elements',divide_edge)
       call memory_deallo(memor_dom,'DIVIDE_FACE','mesh_multiplication_divide_elements',divide_face)

    end if

    if( ISEQUEN ) then
       !
       ! LNINV_LOC: Node renumbering for sequential run
       !
       call memory_deallo(memor_dom,'LNINV_LOC','mesh_multiplication_divide_elements',lninv_loc)
       call memory_alloca(memor_dom,'LNINV_LOC','mesh_multiplication_divide_elements',lninv_loc,npoin,'IDENTITY')
       !
       ! LEINV_LOC: Element renumbering for sequential run
       !
       call memory_deallo(memor_dom,'LEINV_LOC','mesh_multiplication_divide_elements',leinv_loc)
       call memory_alloca(memor_dom,'LEINV_LOC','mesh_multiplication_divide_elements',leinv_loc,nelem,'IDENTITY')
       !
       ! LBINV_LOC: Boundary renumbering for sequential run
       !
       call memory_deallo(memor_dom,'LBINV_LOC','mesh_multiplication_divide_elements',lbinv_loc)
       call memory_alloca(memor_dom,'LBINV_LOC','mesh_multiplication_divide_elements',lbinv_loc,max(1_ip,nboun),'IDENTITY')
    end if

  end subroutine mesh_multiplication_divide_elements

  subroutine mesh_multiplication_number_new_nodes(npoin_old)

    use mod_parall,         only : PAR_COMM_MY_CODE_WM
    use mod_maths,          only : maths_Szudzik_pairing_function
    use mod_memory,         only : memory_copy
    use mod_communications, only : PAR_SUM

    integer(ip), intent(in) :: npoin_old
    integer(ip)             :: ipoin,node1,node2,iedgg,npoin_sum
    integer(ip), pointer    :: lninv_tmp(:)

    nullify(lninv_tmp)

    npoin_sum = npoin
    call PAR_SUM(npoin_sum,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)

    call memory_copy  (memor_dom,'LNINV_TMP','mesh_multiplication_number_new_nodes',lninv_loc,lninv_tmp)
    call memory_alloca(memor_dom,'LNINV_LOC','mesh_multiplication_number_new_nodes',lninv_loc,npoin)

    ! Old nodes
    do ipoin = 1,npoin_old
       lninv_loc(ipoin) = lninv_tmp(ipoin)
    end do
    ! Edge nodes
    ipoin = npoin_old
    do iedgg = 1,nedgg
       ipoin = ipoin + 1
       node1 = lninv_tmp(ledgg(1,iedgg))
       node2 = lninv_tmp(ledgg(2,iedgg))
       lninv_loc(ipoin) = npoin_sum + maths_Szudzik_pairing_function(node1,node2)
    end do
!!$    ! Face nodes
!!$    do ifacg = 1,nfacg
!!$       ipoin = ipoin + 1
!!$       node1 = lninv_tmp(lfacg(1,ifacg))
!!$       node2 = lninv_tmp(lfacg(2,ifacg))
!!$       node3 = lninv_tmp(lfacg(3,ifacg))
!!$       node4 = lninv_tmp(lfacg(4,ifacg))
!!$       ! order node1, node2, node3, node4
!!$       lninv_loc(ipoin) = min(lninv_loc_old(node1),lninv_loc_old(node2),&
!!$            &                 lninv_loc_old(node3),lninv_loc_old(node4))
!!$    end do
!!$    ! Central nodes
!!$    if( lnuty(QUA04) > 0 .or. lnuty(HEX08) > 0 ) then
!!$       call memgen(1_ip,ngrou_dom,0_ip)
!!$       do ielem = 1,nelem_old
!!$          pelty = abs(ltype_old(ielem))
!!$          if( pelty == QUA04 .or. pelty == HEX08 ) then
!!$             pnode = nnode(pelty)
!!$             ipoin = ipoin + 1
!!$             do inode = 1,pnode
!!$                igrou = lninv_loc_old(lnods_old(inode,ielem))
!!$                gisca(igrou) = gisca(igrou) + 1
!!$             end do
!!$             jgrou = 1
!!$             do igrou = 2,ngrou_dom
!!$                if( gisca(igrou) > gisca(jgrou) ) jgrou = igrou
!!$             end do
!!$             lninv_loc(ipoin) = jgrou
!!$             do inode = 1,pnode
!!$                igrou = lninv_loc_old(lnods_old(inode,ielem))
!!$                gisca(igrou) = 0
!!$             end do
!!$          end if
!!$       end do
!!$       call memgen(3_ip,ngrou_dom,0_ip)
!!$    end if
!!$
    call memory_deallo(memor_dom,'LNINV_TMP','mesh_multiplication_number_new_nodes',lninv_tmp)
    
  end subroutine mesh_multiplication_number_new_nodes

  subroutine mesh_multiplication_edge_node(node1,node2,nodex)
    !-----------------------------------------------------------------------
    !****f* domain/edgnod
    ! NAME
    !    edgnod
    ! DESCRIPTION
    !    Find the extra node NODEX on edge NODE1-NODE2
    ! OUTPUT
    ! USED BY
    !    Turnon
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: node1,node2
    integer(ip), intent(out) :: nodex
    integer(ip)              :: iedgg
    integer(ip)              :: iedgp

    if( node1 < node2 ) then
       iedgg = pedgp(node1)-1
       do iedgp = 1,pedgp(node1+1)-pedgp(node1)
          iedgg = iedgg + 1
          if( ledgg(2,iedgg) == node2 ) then
             nodex = ledgg(3,iedgg)
             return
          end if
       end do
    else
       iedgg = pedgp(node2)-1
       do iedgp = 1,pedgp(node2+1)-pedgp(node2)
          iedgg = iedgg + 1
          if( ledgg(2,iedgg) == node1 ) then
             nodex = ledgg(3,iedgg)
             return
          end if
       end do
    end if

    print*,'EDGE NOT FOUND= ',node1,node2
    call runend('EDGNOD: EDGE MIDDLE NODE NOT FOUND')

  end subroutine mesh_multiplication_edge_node
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-08-09
  !> @brief   Find the node corresponding to a subface
  !> @details Find the node corresponding to a subface
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_multiplication_subface(iface,ifacg,ielty,lnodb_loc)    

    integer(ip),           intent(in)  :: iface
    integer(ip),           intent(in)  :: ifacg
    integer(ip),           intent(in)  :: ielty
    integer(ip),           intent(in)  :: lnodb_loc(:)
    integer(ip)                        :: pnode,pnodf
    integer(ip)                        :: jj,kk
    integer(ip)                        :: lnodf(4)
    integer(ip)                        :: lnode(4)

    select case ( ielty )

    case ( QUA09 )
       
       do jj = 1,size(numQUA09)
          pnodf = numQUA09(jj)
          lnodf(1:pnodf) = lnodb_loc(subQUA09(1:pnodf,jj))
          call maths_heap_sort(2_ip,pnodf,lnodf)          
          sub1: do kk = 1,size(numQUA09)
             pnode = numQUA09(kk)
            if( pnode == pnodf ) then
               lnode(1:pnode) = face_divide(kk,ifacg) % l(1:pnode)
               call maths_heap_sort(2_ip,pnode,lnode)       
               if( all(lnode(1:pnodf)==lnodf(1:pnodf)) ) then
                  face_node(jj,iface) = lfacg(mnodm+1,ifacg)+kk-1
                  exit sub1
               end if
            end if
          end do sub1
       end do

    case ( TRI06 )

       do jj = 1,size(numTRI06)
          pnodf = numTRI06(jj)
          lnodf(1:pnodf) = lnodb_loc(subTRI06(1:pnodf,jj))
          call maths_heap_sort(2_ip,pnodf,lnodf)          
          sub2: do kk = 1,size(numTRI06)
             pnode = numTRI06(kk)
            if( pnode == pnodf ) then
               lnode(1:pnode) = face_divide(kk,ifacg) % l(1:pnode)
               call maths_heap_sort(2_ip,pnode,lnode)       
               if( all(lnode(1:pnodf)==lnodf(1:pnodf)) ) then
                  face_node(jj,iface) = lfacg(mnodm+1,ifacg)+kk-1
                  exit sub2
               end if
            end if
          end do sub2
       end do

    end select
    
  end subroutine mesh_multiplication_subface
  
  subroutine mesh_multiplication_face_node(lnode,nodex,facex)
    !-----------------------------------------------------------------------
    !****f* domain/facnod
    ! NAME
    !    facnod
    ! DESCRIPTION
    !    Find the extra node NODEX on face NODE1-NODE2
    ! OUTPUT
    ! USED BY
    !    Turnon
    !***
    !-----------------------------------------------------------------------

    integer(ip),           intent(in)  :: lnode(:)
    integer(ip),           intent(out) :: nodex
    integer(ip), optional, intent(out) :: facex
    integer(ip)                        :: ifacg,pnode
    integer(ip)                        :: i,j,k
    integer(ip)                        :: lnodf(size(lnode))
    !
    ! Order nodes
    !
    lnodf = lnode
    pnode = size(lnode)
    do i = 1,pnode-1
       do j = i+1,pnode
          if( lnodf(i) > lnodf(j) ) then
             k        = lnodf(i)
             lnodf(i) = lnodf(j)
             lnodf(j) = k
          end if
       end do
    end do
    !
    ! Find face
    !
    ifacg = 0
    nodex = 0
    do while( ifacg < nfacg )
       ifacg = ifacg + 1
       if( lnodf(1) == lfacg(1,ifacg) ) then
          if( lnodf(2) == lfacg(2,ifacg) ) then
             if( lnodf(3) == lfacg(3,ifacg) ) then
                if( lnodf(4) == lfacg(4,ifacg) ) then
                   if( pnode == 4 ) then
                      nodex = lfacg(mnodm+1,ifacg)
                      if( present(facex) ) facex = ifacg
                      return
                   else if( pnode == 9 ) then
                      if( lnodf(5) == lfacg(5,ifacg) ) then
                         if( lnodf(6) == lfacg(6,ifacg) ) then
                            if( lnodf(7) == lfacg(7,ifacg) ) then
                               if( lnodf(8) == lfacg(8,ifacg) ) then
                                  if( lnodf(9) == lfacg(9,ifacg) ) then
                                     nodex = lfacg(mnodm+1,ifacg)
                                     if( present(facex) ) facex = ifacg
                                     return
                                  end if
                               end if
                            end if
                         end if
                      end if
                   else if( pnode == 6 ) then
                      if( lnodf(5) == lfacg(5,ifacg) ) then
                         if( lnodf(6) == lfacg(6,ifacg) ) then
                            nodex = lfacg(mnodm+1,ifacg)
                            if( present(facex) ) facex = ifacg
                            return
                         end if
                      end if
                   else
                      stop 'face search not coded'
                   end if
                end if
             end if
          end if
       end if
    end do
    print*,'FACE NOT FOUND= ',lnode(:)
    !do ifacg =1, nfacg
    !   print*,'FACE= ',ifacg,' NODES= ',lfacg(1:4,ifacg)
    !end do
    call runend('FACNOD: COULD NOT FIND FACE MIDDLE NODE')

  end subroutine mesh_multiplication_face_node

  subroutine mesh_multiplication_output_info()
    !-----------------------------------------------------------------------
    !****f* domain/subdim
    ! NAME
    !    domain
    ! DESCRIPTION
    !    Create edge table
    ! OUTPUT
    !    NNEDG ... Number of edges
    !    LEDGG ... Edge table
    !    LEDGB ... Boundary edge table (when Parall is on)
    ! USED BY
    !    Turnon
    !***
    !-----------------------------------------------------------------------
    use def_kintyp
    use def_parame
    use def_domain
    use def_master
    use def_kermod
    use mod_memchk
    implicit none
    integer(ip) :: ipart,npoin_ok
    !
    ! Output geometry
    !
    npoin_ok = npoi3
    call PAR_SUM(npoin_ok)
    if( IMASTER ) then
       ioutp(1) = 0
       ioutp(3) = 0
       ioutp(2) = 0
       routp(1) = 0.0_rp
       routp(2) = 0.0_rp
       routp(3) = 0.0_rp
       ioutp(2) = npoin_ok
       routp(2) = real(npoin_ok,rp)
       do ipart = 2,npart+1
          ioutp(1) = ioutp(1) + nelem_tot(ipart)
          routp(1) = routp(1) + real(nelem_tot(ipart),rp)
          ioutp(3) = ioutp(3) + nboun_tot(ipart)
          routp(3) = routp(3) + real(nboun_tot(ipart),rp)
       end do
    else
       ioutp(1) = nelem
       ioutp(2) = npoin
       ioutp(3) = nboun
       routp(1) = real(nelem,rp)
       routp(2) = real(npoin,rp)
       routp(3) = real(nboun,rp)
    end if

    call livinf(97_ip,'MESH DIMENSIONS',0_ip)

  end subroutine mesh_multiplication_output_info

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-02-27
  !> @brief   Mesh dimension
  !> @details All gather mesh dimensions
  !>
  !-----------------------------------------------------------------------

  subroutine mesh_multiplication_dimensions()

    use mod_communications, only : PAR_GATHER

    integer(ip)          :: npoin_tmp

    npoin_tmp = npoi1+(npoi3-npoi2)+1
    call PAR_GATHER(npoin_tmp,npoin_tot)
    call PAR_GATHER(npoin    ,npoin_tot)
    call PAR_GATHER(nelem    ,nelem_tot)
    call PAR_GATHER(nboun    ,nboun_tot)

  end subroutine mesh_multiplication_dimensions

  subroutine mesh_multiplication_memory(itask)
    !-----------------------------------------------------------------------
    !****f* domain/submem
    ! NAME
    !    domain
    ! DESCRIPTION
    !    De/allocate memory for mesh division algorithm
    ! OUTPUT
    ! USED BY
    !    Turnon
    !***
    !-----------------------------------------------------------------------
    use def_kintyp
    use def_parame
    use def_domain
    use def_master
    use def_kermod
    use mod_memory
    use mod_graphs, only : graphs_elepoi_deallocate
    implicit none
    integer(ip), intent(in) :: itask
    integer(ip)             :: pelty
    
    if( INOTMASTER ) then

       select case ( itask )

       case ( 1_ip )
          !
          ! Define flags
          !
          if_high_order = .false.
          do pelty = 1,nelty
             if( lnuty(pelty) > 0 .and. element_type(pelty) % order > 1 ) if_high_order = .true.
          end do
          nullify(face_divide)
          !
          ! Allocate memory
          !
          if( .not. associated(lnlev) ) then
             call memory_alloca(memor_dom,'LNLEV','renpoi',lnlev,npoin)
          end if
          if( .not. associated(lelev) ) then
             call memory_alloca(memor_dom,'LELEV','renpoi',lelev,nelem)
          end if
          if( .not. associated(lblev) ) then
             call memory_alloca(memor_dom,'LBLEV','renpoi',lblev,nboun)
          end if

       case ( -1_ip )
          !
          ! LNLEV, LELEV, LBLEV: Deallocate memory
          !
          call memory_deallo(memor_dom,'LNLEV','renpoi',lnlev)
          call memory_deallo(memor_dom,'LELEV','renpoi',lelev)
          call memory_deallo(memor_dom,'LBLEV','renpoi',lblev)

       case ( 2_ip )
          !
          ! LEDGG: Deallocate memory
          !
          call memory_deallo(memor_dom,'LEDGG','renpoi',ledgg)
          call graphs_elepoi_deallocate(pelpo,lelpo,memor=memor_dom)
          !
          ! LFACG: Deallocate memory
          !
          if( nfacg > 0 ) then
             call memory_deallo(memor_dom,'LFACG','renpoi',lfacg)
          end if

       end select

    end if

  end subroutine mesh_multiplication_memory

  subroutine mesh_multiplication_node_codes()
    !-----------------------------------------------------------------------
    !****f* domain/subelm
    ! NAME
    !    domain
    ! DESCRIPTION
    !    Cancel codes for nodes which are not on the boundary to avoid
    !    false boundary condition when interpolating them
    ! OUTPUT
    ! USED BY
    !    submsh
    !***
    !-----------------------------------------------------------------------
    use def_kintyp
    use def_elmtyp
    use def_domain
    use def_master
    use def_kermod
    implicit none
    integer(ip) :: ipoin,ii

    if( ndivi > 0 .and. INOTMASTER ) then
       if( kfl_icodn > 0 ) then
          do ipoin = 1,npoin
             if( lpoty(ipoin) == 0 ) then
                do ii = 1,mcono
                   kfl_codno(ii,ipoin) = mcodb + 1
                end do
             end if
          end do
       end if
    end if

  end subroutine mesh_multiplication_node_codes

  subroutine mesh_multiplication_list_faces()
    !------------------------------------------------------------------------
    !****f* domain/lfaces
    ! NAME
    !    lfaces
    ! DESCRIPTION
    !    This routine computes the list of faces
    ! USED BY
    !    submsh
    !***
    !------------------------------------------------------------------------
    use def_parame
    use def_elmtyp
    use def_master
    use def_domain
    use mod_elmgeo, only : element_type
    use mod_memory, only : memory_alloca
   
    integer(ip) :: ielty,ielem,iface,inodb,ilist,pblty
    integer(ip) :: inode,jelem,jface,jelty,ipoin,porde
    integer(ip) :: pepoi,ielpo,pnodb,pelty
    logical(lg) :: equal_faces

    call messages_live('FACE TABLE')
    nfacg = 0
    nfaco = 0
    mnodm = max(4_ip,mnodb)
    
    if( INOTMASTER ) then

       if( lnuty(HEX08) /= 0 .or. lnuty(HEX27) /= 0 .or. lnuty(PYR05) /= 0 .or. lnuty(PEN06) /= 0 .or. lnuty(TET10) /= 0 ) then
          !
          ! Faces graph: Allocate memory for FACES
          !
          call memory_alloca(memor_dom,'FACEL','memgeo',facel,mnodb+1_ip,mface,nelem)
          !
          ! Construct and sort FACES
          !
          !$*OMP  PARALLEL DO SCHEDULE (STATIC)           &
          !$*OMP  DEFAULT (NONE)                          &
          !$*OMP  PRIVATE (ielem,ielty,iface,inodb,inode, &
          !$*OMP           pnodb)                         &
          !$*OMP  SHARED  (ltype,facel,lnods,nelem,nface, &
          !$*OMP           element_type,mnodb)
          !
          do ielem = 1,nelem
             ielty = abs(ltype(ielem))
             do iface = 1,element_type(ielty) % number_faces
                if(  element_type(ielty) % type_faces(iface) == QUA04 .or. &
                     element_type(ielty) % type_faces(iface) == QUA09 .or. &
                     element_type(ielty) % type_faces(iface) == TRI06 ) then
                   pnodb = element_type(ielty) % node_faces(iface)
                   do inodb = 1,pnodb
                      inode = element_type(ielty) % list_faces(inodb,iface)
                      facel(inodb,iface,ielem) = lnods(inode,ielem)
                   end do
                   facel(mnodb+1,iface,ielem) = 1
                   call sortin(pnodb,facel(1,iface,ielem))
                else
                   facel(mnodb+1,iface,ielem) = 0
                end if
             end do
          end do
          !
          ! Compute FACES
          !
          do ielem = 1,nelem                                            ! Compare the faces and
             ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
             do iface = 1,element_type(ielty) % number_faces
                if( facel(mnodb+1,iface,ielem) > 0 ) then
                   nfacg = nfacg + 1
                   ipoin = facel(1,iface,ielem)
                   ilist = 1
                   pepoi = pelpo(ipoin+1)-pelpo(ipoin)
                   ielpo = pelpo(ipoin)-1
                   do while( ilist <= pepoi )
                      ielpo = ielpo + 1
                      jelem = lelpo(ielpo)
                      if( jelem /= ielem ) then
                         jelty = abs(ltype(jelem))                              ! eliminate the repited faces
                         jface = 0
                         do while( jface /= element_type(jelty) % number_faces )
                            jface = jface + 1
                            if( facel(mnodb+1,jface,jelem) > 0 ) then
                               equal_faces = .true.
                               inodb = 0
                               do while( equal_faces .and. inodb /= element_type(jelty) % node_faces(jface) )
                                  inodb = inodb+1
                                  if( facel(inodb,iface,ielem) /= facel(inodb,jface,jelem) ) &
                                       equal_faces = .false.
                               end do
                               if( equal_faces ) then
                                  facel(mnodb+1,iface,ielem) =  jelem        ! Keep IELEM face
                                  facel(mnodb+1,jface,jelem) = -ielem        ! Elminate JELEM face
                                  facel(1,jface,jelem)       =  iface        ! Remember IFACE face
                                  jface = element_type(jelty) % number_faces ! Exit JFACE do
                                  ilist = pepoi                              ! Exit JELEM do
                               end if
                            end if
                         end do
                      end if
                      ilist = ilist + 1
                   end do
                end if
             end do
          end do
          !
          ! Allocate memory
          !
          call memory_alloca(memor_dom,'LFACG','memgeo',lfacg,mnodm+5,nfacg)

          nfacg = 0
          do ielem = 1,nelem                                            ! Compare the faces and
             ielty = abs(ltype(ielem))                                  ! eliminate the repeated faces
             do iface = 1,element_type(ielty) % number_faces
                if( facel(mnodb+1,iface,ielem) > 0 ) then
                   nfacg = nfacg + 1
                   pnodb = element_type(ielty) % node_faces(iface) 
                   do inode = 1,pnodb
                      lfacg(inode,nfacg) = facel(inode,iface,ielem)
                   end do
                   lfacg(mnodm+2,nfacg) = element_type(ielty) % type_faces(iface)
                   lfacg(mnodm+3,nfacg) = ielem
                   lfacg(mnodm+4,nfacg) = iface
                   facel(1,iface,ielem) = nfacg
                   pelty                = element_type(ielty) % type_faces(iface)
                   porde                = element_type(pelty) % order
                   nfaco(pelty)         = nfaco(pelty) + 1
                end if
             end do
          end do
          
          do ielem = 1,nelem
             ielty = abs(ltype(ielem))
             do iface = 1,element_type(ielty) % number_faces
                if( facel(mnodb+1,iface,ielem) < 0 ) then
                   jelem = -facel(mnodb+1,iface,ielem)
                   jface =  facel(1,iface,ielem)
                   facel(1,iface,ielem) = facel(1,jface,jelem)
                end if
             end do
          end do
          !
          ! Order faces
          !
          do iface = 1,nfacg
              pblty = lfacg(mnodm+2,iface)
              pnodb = element_type(pblty) % number_nodes
             call sortin(pnodb,lfacg(:,iface))
          end do

       end if

    end if

  end subroutine mesh_multiplication_list_faces

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet
  !> @date    October, 2019
  !> @brief   Divide interface elements / element normal division
  !> @details Divide edges, faces and boundaries for interface elements
  !>          or element normal division method.
  !>
  !>  2-d                  3-d
  !>
  !>                           8           7
  !>                           o-----*-----o
  !>                          /     /     /|
  !>                         *-----x-----* |       ^
  !>  4           3        5/     /    6/| |       |
  !>  o-----*-----o        o-----*-----o | |       |
  !>  |     |     |        |     |     | | |       | Element normal direction
  !>  |     |     |        |     |     | | o       | 2-d: 1-2-3-4
  !>  |     |     |        |     |     | |/3       | 3-d: 1-2-3-4-5-6-7-8
  !>  |     |     |        |     |     | *
  !>  |     |     |        |     |     |/
  !>  o-----*-----o        o-----*-----o
  !>  1           2        1           2
  !>
  !>
  !-----------------------------------------------------------------------

  subroutine mesh_multiplication_divide_ELINT(divide_edge,divide_face,divide_boun)

    use def_elmtyp, only : ELINT,ELFEM
    use mod_memory, only : memory_alloca

    implicit none

    logical(lg), pointer, intent(inout) :: divide_edge(:) !< Flag for edge division
    logical(lg), pointer, intent(inout) :: divide_face(:) !< Flag for face division
    logical(lg), pointer, intent(inout) :: divide_boun(:) !< Flag for boundary division

    integer(ip)                :: iedgg,ielem,ifacg,iboun
    integer(ip)                :: inode,knode
    integer(ip)                :: npoie,npoif        ! Number points edge and face
    integer(ip), pointer       :: nodee_ledgg(:)     ! List nodes edge global
    integer(ip), pointer       :: nodee_lelch(:)     ! List nodes edge element characteristic
    integer(ip), pointer       :: nodeb_lnodb(:)     ! List nodes edge boundary element
    integer(ip), pointer       :: nodef_lfacg(:)     ! List nodes face global
    integer(ip), pointer       :: nodef_lelch(:)     ! List nodes face element characteristic
    !
    ! Nullify pointers
    !
    nullify(nodee_ledgg)
    nullify(nodee_lelch)
    nullify(nodeb_lnodb)
    nullify(nodef_lfacg)
    nullify(nodef_lelch)
    !
    npoif = 4_ip  ! Number of points face
    npoie = 2_ip  ! Number of points edge
    !
    ! Allocate memory of nodes lists
    !
    call memory_alloca(memor_dom,'NODEE_LEDGG','mesh_multiplication_divide_elements',nodee_ledgg,npoie)
    call memory_alloca(memor_dom,'NODEE_LELCH','mesh_multiplication_divide_elements',nodee_lelch,npoie)
    call memory_alloca(memor_dom,'NODEF_LFACG','mesh_multiplication_divide_elements',nodef_lfacg,npoif)
    call memory_alloca(memor_dom,'NODEF_LELCH','mesh_multiplication_divide_elements',nodef_lelch,npoif)
    call memory_alloca(memor_dom,'NODEB_LNODB','mesh_multiplication_divide_elements',nodeb_lnodb,mnodb)
    !
    ! Nodes Edges
    !
    do iedgg = 1,nedgg
       nodee_ledgg(1:npoie) = ledgg(1:npoie,iedgg)
       do ielem = 1, nelem
          if ( lelch(ielem) == ELINT .or. (lelch(ielem) == ELFEM .and. kfl_elndi == 1) ) then
             if ( ndime == 2 ) then
                ! Edge 1-4
                nodee_lelch(1) = lnods(1,ielem)
                nodee_lelch(2) = lnods(4,ielem)
                call maths_heap_sort(2_ip,npoie,nodee_lelch(1:npoie))
                if ( all(nodee_ledgg == nodee_lelch) ) divide_edge(iedgg) = .false.
                ! Edge 2-3
                nodee_lelch(1) = lnods(2,ielem)
                nodee_lelch(2) = lnods(3,ielem)
                call maths_heap_sort(2_ip,npoie,nodee_lelch(1:npoie))
                if ( all(nodee_ledgg == nodee_lelch) ) divide_edge(iedgg) = .false.
             else
                do inode = 1,nnode(abs(ltype(ielem)))/2
                   knode  = inode + nnode(abs(ltype(ielem)))/2
                   ! Edge x - x
                   nodee_lelch(1) = lnods(inode,ielem)
                   nodee_lelch(2) = lnods(knode,ielem)
                   call maths_heap_sort(2_ip,npoie,nodee_lelch(1:npoie))
                   if ( all(nodee_ledgg == nodee_lelch) ) divide_edge(iedgg) = .false.
                end do
             end if
          end if
       end do

    end do
    !
    ! Nodes Faces
    !
    do ifacg = 1,nfacg
       nodef_lfacg(1:npoif) = lfacg(1:npoif,ifacg)
       do ielem = 1, nelem
          if ( lelch(ielem) == ELINT .or. (lelch(ielem) == ELFEM .and. kfl_elndi == 1) ) then
             ! Face 1
             nodef_lelch(1) = lnods(1,ielem)
             nodef_lelch(2) = lnods(2,ielem)
             nodef_lelch(3) = lnods(5,ielem)
             nodef_lelch(4) = lnods(6,ielem)
             call maths_heap_sort(2_ip,npoif,nodef_lelch(1:npoif))
             if ( all(nodef_lfacg == nodef_lelch) ) divide_face(ifacg) = .false.
             ! Face 2
             nodef_lelch(1) = lnods(2,ielem)
             nodef_lelch(2) = lnods(3,ielem)
             nodef_lelch(3) = lnods(6,ielem)
             nodef_lelch(4) = lnods(7,ielem)
             call maths_heap_sort(2_ip,npoif,nodef_lelch(1:npoif))
             if ( all(nodef_lfacg == nodef_lelch) ) divide_face(ifacg) = .false.
             ! Face 3
             nodef_lelch(1) = lnods(3,ielem)
             nodef_lelch(2) = lnods(4,ielem)
             nodef_lelch(3) = lnods(7,ielem)
             nodef_lelch(4) = lnods(8,ielem)
             call maths_heap_sort(2_ip,npoif,nodef_lelch(1:npoif))
             if ( all(nodef_lfacg == nodef_lelch) ) divide_face(ifacg) = .false.
             ! Face 4
             nodef_lelch(1) = lnods(1,ielem)
             nodef_lelch(2) = lnods(4,ielem)
             nodef_lelch(3) = lnods(5,ielem)
             nodef_lelch(4) = lnods(8,ielem)
             call maths_heap_sort(2_ip,npoif,nodef_lelch(1:npoif))
             if ( all(nodef_lfacg == nodef_lelch) ) divide_face(ifacg) = .false.
          end if
       end do

    end do
    !
    ! Nodes boundaries
    !
    do iboun = 1,nboun
       ielem = lelbo(iboun)
       nodeb_lnodb(1:mnodb) = lnodb(1:mnodb,iboun)
       call maths_heap_sort(2_ip,mnodb,nodeb_lnodb(1:mnodb))
       if( lelch(ielem) == ELINT .or. (lelch(ielem) == ELFEM .and. kfl_elndi == 1) ) then
          if( ndime == 2 ) then
             ! Edge 2-3
             nodee_lelch(1) = lnods(2,ielem)
             nodee_lelch(2) = lnods(3,ielem)
             call maths_heap_sort(2_ip,mnodb,nodee_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodee_lelch) ) divide_boun(iboun) = .false.
             ! Edge 1-4
             nodee_lelch(1) = lnods(1,ielem)
             nodee_lelch(2) = lnods(4,ielem)
             call maths_heap_sort(2_ip,mnodb,nodee_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodee_lelch) ) divide_boun(iboun) = .false.
          else
             ! Face 1
             nodef_lelch(1) = lnods(1,ielem)
             nodef_lelch(2) = lnods(2,ielem)
             nodef_lelch(3) = lnods(5,ielem)
             nodef_lelch(4) = lnods(6,ielem)
             call maths_heap_sort(2_ip,mnodb,nodef_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodef_lelch) ) divide_boun(iboun) = .false.
             ! Face 2
             nodef_lelch(1) = lnods(2,ielem)
             nodef_lelch(2) = lnods(3,ielem)
             nodef_lelch(3) = lnods(6,ielem)
             nodef_lelch(4) = lnods(7,ielem)
             call maths_heap_sort(2_ip,mnodb,nodef_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodef_lelch) ) divide_boun(iboun) = .false.
             ! Face 3
             nodef_lelch(1) = lnods(3,ielem)
             nodef_lelch(2) = lnods(4,ielem)
             nodef_lelch(3) = lnods(7,ielem)
             nodef_lelch(4) = lnods(8,ielem)
             call maths_heap_sort(2_ip,mnodb,nodef_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodef_lelch) ) divide_boun(iboun) = .false.
             ! Face 4
             nodef_lelch(1) = lnods(1,ielem)
             nodef_lelch(2) = lnods(4,ielem)
             nodef_lelch(3) = lnods(5,ielem)
             nodef_lelch(4) = lnods(8,ielem)
             call maths_heap_sort(2_ip,mnodb,nodef_lelch(1:mnodb))
             if ( all(nodeb_lnodb == nodef_lelch) ) divide_boun(iboun) = .false.
          end if
       end if

    end do
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_dom,'NODEE_LEDGG','mesh_multiplication_divide_elements',nodee_ledgg)
    call memory_deallo(memor_dom,'NODEE_LELCH','mesh_multiplication_divide_elements',nodee_lelch)
    call memory_deallo(memor_dom,'NODEB_LNODB','mesh_multiplication_divide_elements',nodeb_lnodb)
    call memory_deallo(memor_dom,'NODEF_LFACG','mesh_multiplication_divide_elements',nodef_lfacg)
    call memory_deallo(memor_dom,'NODEF_LELCH','mesh_multiplication_divide_elements',nodef_lelch)

  end subroutine mesh_multiplication_divide_ELINT

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-08-02
  !> @brief   Compute mid-point
  !> @details Compute mid-point using iso-parametric coordinates
  !> 
  !-----------------------------------------------------------------------
  
  subroutine coord_new_point(lnode,ielem,pelty,lnods_in,coord_in,x)
    
    integer(ip),            intent(in)  :: lnode(:)
    integer(ip),            intent(in)  :: ielem
    integer(ip),            intent(in)  :: pelty
    integer(ip),            intent(in)  :: lnods_in(:,:)
    real(rp),               intent(in)  :: coord_in(:,:)
    real(rp),               intent(out) :: x(:)
    real(rp)                            :: s(3)
    integer(ip)                         :: pnode,pdime,inode,ipoin,knode
    !
    ! Dimensions
    !
    pdime = element_type(pelty) % dimensions
    pnode = element_type(pelty) % number_nodes
    knode = size(lnode)
    !
    ! Average is-parametric coordinates
    !
    s = 0.0_rp
    do inode = 1,knode
       s(1:pdime) = s(1:pdime) + element_type(pelty) % coord(1:pdime,lnode(inode))
    end do
    s = s / real(knode,rp)

    call elmgeo_shapf_deriv_heslo(pdime,pnode,s,shapf)
    
    x = 0.0_rp
    do inode = 1,pnode
       ipoin      = lnods_in(inode,ielem)
       x(1:pdime) = x(1:pdime) + shapf(inode) * coord_in(1:pdime,ipoin)
    end do
    
  end subroutine coord_new_point

end module mod_mesh_multiplication
!> @}
