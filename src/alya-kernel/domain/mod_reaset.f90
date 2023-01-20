!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    mod_reaset.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read sets definitions
!> @details Read sets definitions.
!!          Sets are used to compute some global values, computed as
!!          volume integrals, boundary integrals or node-wise.
!!          It can be useful for example to compute the force exterted
!!          by the fluid on a body, by assigning a set to all the boundary elements
!!          that define the body. For each module with extension 'ext', the results
!!          are written in the following files:
!!          - Element sets: *-element.ext.set
!!          - Boundary sets: *-boundary.ext.set
!!          - Node sets: *-node.ext.set
!!          The output variables ares:
!!          \verbatim
!!          - NESET .......... Number of element sets (read)
!!          - NBSET .......... Number of boundary sets (read)
!!          - NNSET .......... Number of node sets (read)
!!          - LESET(NELEM) ... List of element sets (read)
!!          - LNSET(NBOUN) ... List of boundary sets (read)
!!          - LESEC(NESET) ... List of element set numbers
!!          - LBSEC(NBSET) ... List of boundary set numbers
!!          - LNSEC(NNSET) ... List of node set numbers
!!          \endverbatim
!> @}
!-----------------------------------------------------------------------

module mod_reaset
  
  use def_kintyp,        only : ip
  use def_master,        only : IPARALL
  use def_domain,        only : lesec
  use def_domain,        only : lbsec
  use def_domain,        only : neset
  use def_domain,        only : nbset
  use def_domain,        only : nnset_origi
  use def_domain,        only : neset_origi
  use def_domain,        only : nbset_origi
  use def_domain,        only : nperi
  use def_domain,        only : memor_dom
  use mod_messages,      only : messages_live
  use mod_domain,        only : domain_memory_allocate
  use mod_memory,        only : memory_alloca
  use mod_memory,        only : memory_deallo
  use mod_memory,        only : memory_size
  use mod_par_tools,     only : par_tools_merge_lists 
  use mod_maths_sort,    only : maths_heap_sort
  
  implicit none
  private

  public :: reaset_seq
  public :: reaset_set_connectivity

contains


  subroutine reaset_seq()
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    use mod_iofile
    use mod_ecoute, only :  ecoute
    implicit none
    integer(ip) :: ielem,ieset,iesta,iesto
    integer(ip) :: ibsta,ibsto,inset
    integer(ip) :: iboun,jbset,ipara,ipoin
    logical     :: set_processed


    set_processed = .FALSE.

    if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then
       !
       ! Initializations
       !
       neset       = 0                   ! There is no element set
       nbset       = 0                   ! There is no boundary set
       nnset       = 0                   ! There is no node set
       neset_origi = 0                   ! There is no element set
       nbset_origi = 0                   ! There is no boundary set
       nnset_origi = 0                   ! There is no node set
       !
       ! Read vectors
       !
       call ecoute('REASET')
       do while( words(1) /= 'SETS ' .and. words(1) /= 'SETSA' )
          call ecoute('REASET')
       end do
       if( exists('OFF  ') ) then
          do while( words(1) /= 'ENDSE' )
             call ecoute('reaset')
          end do
          return
       end if
       !.md<module>kernel
       !.md<input>case.dom.dat
       !.md<pos>3
       !.md<sec>
       !.md<0># Sets
       !.md<>
       !.md<code>
       !.md<0><b>SETS</b>
       !.md<field>SETS
       !.md<com>Sets are used to compute some global values, computed as
       !.md<com>volume integrals, boundary integrals or node-wise.
       !.md<com>It can be useful for example to compute the force exterted
       !.md<com>by the fluid on a body, by assigning a set to all the boundary elements
       !.md<com>that define the body. For each module with extension 'ext', the results
       !.md<com>are written in the following files:
       !.md<com>    -  Element sets: *-element.ext.set    
       !.md<com>    -  Boundary sets: *-boundary.ext.set  
       !.md<com>    -  Node sets: *-node.ext.set          
       !
       do while( words(1) /='ENDSE')
          call ecoute('reaset')

          if( words(1) =='ELEME' ) then
             !
             !.md<1>ELEMENTS: [ALL= int]
             !.md<2>int1 int2                                                                   $ Element, set number
             !.md<2>...
             !.md<1>END_ELEMENTS
             !.md<field>ELEMENTS
             !.md<com>This field defines the elements sets. It contains a list of elements (int1) and
             !.md<com>the set assigned to it (int2). The option ALL= int assigns automatically the
             !.md<com>set int to all the elements of the mesh.
             !
             set_processed = .TRUE.

             if( words(2) == 'ALL  ' ) then
                neset = 1
                call domain_memory_allocate('LESET')
                ieset = getint('ALL  ',0_ip,'#SET NUMBER')
                if( ieset == 0 ) call runend('REASET: WRONG SET NUMBER')
                do ielem = 1,nelem
                   leset(ielem)=ieset
                end do
                do while( words(1) /='ENDEL')
                   call ecoute('reaset')
                end do
             else
                call ecoute('reaset')
                if( words(1) /='ENDEL' ) then
                   neset = 1
                   call domain_memory_allocate('LESET')
                   do while( words(1) /='ENDEL')
                      if( words(1) =='FROM ' ) then
                         iesta=int(param(1))
                         iesto=int(param(2))
                         ieset=int(param(3))
                         do ielem=iesta,iesto
                            leset(ielem)=ieset
                         end do
                      else
                         ielem=int(param(1))
                         if(ielem<0.or.ielem>nelem) &
                              call runend('REASET: WRONG NUMBER OF ELEMENT '&
                              //adjustl(trim(intost(ielem))))
                         leset(ielem)=int(param(2))
                      end if
                      call ecoute('reaset')
                   end do
                end if
             end if

          else if( words(1) =='BOUND' ) then
             !
             !.md<1>BOUNDARIES: [ALL= int]
             !.md<2>int1 int2                                                                   $ Boundary, set number
             !.md<2>...
             !.md<1>END_BOUNDARIES
             !.md<field>BOUNDARIES
             !.md<com>This field defines the boundary sets. It contains a list of boundaries (int1) and
             !.md<com>the set assigned to it (int2). The option ALL= int assigns automatically the
             !.md<com>set int to all the boundaries of the mesh.
             !
             set_processed = .TRUE.

             if( words(2) == 'ALL  ' ) then
                nbset = 1
                jbset = getint('ALL  ',1_ip,'#SET NUMBER')
                call domain_memory_allocate('LBSET')                
                do iboun = 1,nboun
                   lbset(iboun) = jbset
                end do
                do while( words(1) /='ENDBO')
                   call ecoute('reaset')
                end do
             end if

             if( nbset == 0 ) call ecoute('reaset')
             if( words(1) /='ENDBO' ) then
                if( kfl_autbo == 1 ) then
                   call runend('BOUNDARY SET DEFINITION IMPOSSIBLE '&
                        //'WHEN USING AUTOMATIC BOUNDARY OPTION')
                end if
                nbset = 1
                call domain_memory_allocate('LBSET')                                
                do while( words(1) /='ENDBO')
                   if( words(1) =='FROM ' ) then
                      ibsta = int(param(1))
                      ibsto = int(param(2))
                      jbset = int(param(3))
                      do iboun = ibsta,ibsto
                         lbset(iboun) = jbset
                      end do
                   else
                      iboun = int(param(1))
                      if( iboun < 0 .or. iboun > nboun ) &
                           call runend('REASET: WRONG NUMBER OF BOUNDARY '&
                           //adjustl(trim(intost(iboun))))
                      lbset(iboun) = int(param(2))
                   end if
                   call ecoute('reaset')
                end do
             end if

          else if( words(1) =='NODES' ) then
             !
             !.md<1>NODES:
             !.md<2>int1 int2  $ Node_id Set_id
             !.md<2>...
             !.md<1>END_NODES
             !.md<field>NODES
             !.md<com>This field defines the node sets. Only one node per set
             !

             !-------------------------------------------------------------
             !
             ! LNSET: Nodes
             !
             !-------------------------------------------------------------

             call ecoute('reaset')

             if( words(1) /= 'ENDNO' ) then
                set_processed = .TRUE.
                nnset = 1
                call domain_memory_allocate('LNSET')                                
                if( words(2) == 'OLDFO' ) then
                   call messages_live('YOU ARE USING THE OLD FORMAT FOR DECLARING NODE SETS','WARNING')
                   do while( words(1) /='ENDNO')
                      inset = 0
                      do ipara = 1,nnpar
                         ipoin = int(param(ipara))
                         if( ipoin > 0 .and. ipoin <= npoin) then
                            inset        = inset + 1
                            lnset(ipoin) = inset
                         end if
                      end do
                      call ecoute('reaset')
                   end do
                else
                   do while( words(1) /='ENDNO')
                      ipoin        = int(param(1),ip)
                      inset        = int(param(2),ip)
                      lnset(ipoin) = inset
                      call ecoute('reaset')
                   end do
                end if
             end if
             
          end if
       end do
       !
       !.md<0>END_SETS
       !.md</code>
       !.md<>
       !

       if( .NOT. set_processed ) then
          call messages_live('Empty SETS...END_SETS in dom.dat. Likely missing ELEME, BOUND, or NODES.','WARNING')
       end if

       nnset_origi = nnset
       neset_origi = neset
       nbset_origi = nbset
       
    end if

  end subroutine reaset_seq

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-26
  !> @brief   Sets
  !> @details Master should allocate memory.
  !>          NESET, NBSET, NNSET must be the total number of sets
  !>          Compute:
  !>          LESEC
  !>          LBSEC
  !>          LNSEC
  !>
  !-----------------------------------------------------------------------

  subroutine reaset_set_connectivity()
            
    call reaset_element_set_connectivity()
    call reaset_boundary_set_connectivity()
    call reaset_node_set_connectivity()
    
    if( IPARALL ) then

       call par_tools_merge_lists(lesec,memor_dom,VARIABLE_NAME='LESEC')
       call par_tools_merge_lists(lbsec,memor_dom,VARIABLE_NAME='LBSEC')

       neset = memory_size(lesec)
       nbset = memory_size(lbsec)    
       !
       ! Very special case, node sets
       !
      call reaset_node_parall()

    end if
    
  end subroutine reaset_set_connectivity

  ! Count number of element sets and define the element set connectivity LESEC
  ! Local # body of ielem = LESEC(LESET(IELEM))

  subroutine reaset_element_set_connectivity

    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    implicit none

    integer(ip) :: ielem,jelem,keset,ieset
    
    if( neset_origi == 1 ) then
       neset = 0
       do ielem = 1,nelem
          if( leset(ielem) > 0 ) then
             neset = neset + 1
             keset = leset(ielem)
             do jelem = ielem,nelem
                if( leset(jelem) == keset ) leset(jelem) = -keset
             end do
          end if
       end do       
       call domain_memory_allocate('LESEC')                                
       ieset = 0
       do ielem = 1,nelem
          if( leset(ielem) < 0 ) then
             ieset = ieset + 1
             keset = leset(ielem)
             do jelem = ielem,nelem
                if( leset(jelem) == keset ) leset(jelem) = -keset
             end do
             lesec(ieset) = leset(ielem)
          end if
       end do
    end if

    if( associated(lesec) ) call maths_heap_sort(2_ip,neset,lesec)

  end subroutine reaset_element_set_connectivity

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-13
  !> @brief   Node set parallelization
  !> @details Node set parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine reaset_boundary_set_connectivity

    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    implicit none

    integer(ip) :: iboun,jboun,kbset,jbset

    if( nbset_origi == 1 ) then
       nbset = 0
       do iboun = 1,nboun
          if( lbset(iboun) > 0 ) then
             nbset = nbset + 1
             kbset = lbset(iboun)
             do jboun = iboun,nboun
                if( lbset(jboun) == kbset ) lbset(jboun) = -kbset
             end do
          end if
       end do
       call domain_memory_allocate('LBSEC')                                       
       jbset = 0
       do iboun = 1,nboun
          if( lbset(iboun) < 0 ) then
             jbset = jbset+1
             kbset = lbset(iboun)
             do jboun = iboun,nboun
                if( lbset(jboun) == kbset ) lbset(jboun) = -kbset
             end do
             lbsec(jbset) = lbset(iboun)
          end if
       end do       
       if( associated(lbsec) ) call maths_heap_sort(2_ip,nbset,lbsec)
    end if

  end subroutine reaset_boundary_set_connectivity

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Set connectivity
  !> @details Define the node set connectivity LNSEC
  !>          Node number= LBSEC(INSET)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine reaset_node_set_connectivity

    use def_kintyp,         only : ip
    use def_master,         only : ISEQUEN, intost
    use def_domain,         only : lnset
    use def_domain,         only : lnsec
    use def_domain,         only : nnset
    use def_domain,         only : nnset_origi
    use def_domain,         only : npoin
    use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
    use def_master,         only : THIS_NODE_IS_MINE
    use def_master,         only : npoi3
    implicit none

    integer(ip)          :: ipoin
    integer(ip), pointer :: take_min(:)
    integer(ip), pointer :: multi(:)
    integer(ip)          :: lnsec_size
    !
    ! Do we really have set... this is useful when repartitioning...
    !
    nullify(take_min,multi)
    if( nnset_origi == 1 ) then
       
       if( nperi == -1 ) then
          !
          ! Special treatment for periodic nodes. In fact, a node set could not
          ! be the own node of any subdomain.
          !  
          call memory_alloca(memor_dom,'TAKE_MIN','par_merge_parallel_list',take_min,npoin)
          call memory_alloca(memor_dom,'MULTI'   ,'par_merge_parallel_list',multi,npoin)

          do ipoin = 1,npoin
             if( lnset(ipoin) /= 0 ) then
                take_min(ipoin) = lnset(ipoin)
                multi(ipoin)    = 1
                lnset(ipoin)    = 0
             end if
          end do
          call PAR_INTERFACE_NODE_EXCHANGE(take_min,'SUM','IN MY CODE')
          call PAR_INTERFACE_NODE_EXCHANGE(multi   ,'SUM','IN MY CODE')

          do ipoin = 1,npoi3
             if( take_min(ipoin) /= 0 ) then
                lnset(ipoin) = take_min(ipoin) / multi(ipoin)
             end if
          end do
          
          call memory_deallo(memor_dom,'TAKE_MIN','par_merge_parallel_list',take_min)
          call memory_deallo(memor_dom,'MULTI'   ,'par_merge_parallel_list',multi)

       end if
       
       nnset = 0
       do ipoin = 1,npoin
          if( lnset(ipoin) > 0 .and. THIS_NODE_IS_MINE(ipoin) ) then
             nnset = nnset + 1
          end if
       end do
       
       call domain_memory_allocate('LNSEC')                                              
       nnset = 0
       lnsec_size = memory_size(lnsec)
       if( ISEQUEN ) then
          !
          ! In sequential, order the sets to have the same input as in parall
          !
          do ipoin = 1,npoin
             if( lnset(ipoin) > 0 ) then
                nnset = nnset + 1
                if ( lnset(ipoin) > lnsec_size ) call runend("ATTEMPING TO ACCESS LNSEC OF SIZE "//trim(intost(lnsec_size))//" AT INDEX "//trim(intost(lnset(ipoin)))//". CHECK IF NODE SETS ARE CONSECUTIVE NUMBERS")
                lnsec(lnset(ipoin)) = ipoin
             end if
          end do
       else
          !
          ! Parallel: order given my master
          !
          do ipoin = 1,npoin
             if( lnset(ipoin) > 0 .and. THIS_NODE_IS_MINE(ipoin) ) then
                nnset = nnset + 1
                if ( nnset > lnsec_size ) call runend("ATTEMPING TO ACCESS LNSEC OF SIZE "//trim(intost(lnsec_size))//" AT INDEX "//trim(intost(lnset(ipoin)))//". CHECK IF NODE SETS ARE CONSECUTIVE NUMBERS")
                lnsec(nnset) = ipoin
             end if
          end do
       end if
    end if
    
  end subroutine reaset_node_set_connectivity

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Setup node set communication
  !> @details Setup node set communication
  !> 
  !-----------------------------------------------------------------------

  subroutine reaset_node_parall()
    use def_master,         only : IMASTER
!    use def_master,         only : INOTMASTER
    use def_master,         only : INOTEMPTY
    use def_master,         only : lninv_loc
    use def_master,         only : THIS_NODE_IS_MINE
    use def_domain,         only : memor_dom
    use def_domain,         only : nnset,lnsec
    use def_domain,         only : npoin,lnset
    use mod_communications, only : PAR_GATHER
    use mod_communications, only : PAR_GATHERV
    use mod_communications, only : PAR_SUM
    use mod_communications, only : PAR_MIN
    use mod_communications, only : PAR_MAX
    use mod_maths,          only : maths_heap_sort
    use def_parall,         only : lnsec_par
    use def_parall,         only : nnset_par
    use mod_memory,         only : memory_copy
    use mod_memory,         only : memory_deallo
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_size
    use mod_parall,         only : PAR_CODE_SIZE

    integer(ip)          :: ii,inset,ipart,ipoin
    integer(ip), pointer :: lnsec_set_gat(:)
    integer(ip), pointer :: lnsec_nod_gat(:)
    integer(ip), pointer :: lnsec_nod(:)
    integer(ip), pointer :: lnsec_set(:)
    integer(4),  pointer :: nnset_par4(:)
    integer(4)           :: isize4
    integer(ip)          :: nnset_min
    integer(ip)          :: nnset_max
    integer(ip)          :: nnset_sum

    nullify(lnsec_set_gat)
    nullify(lnsec_nod_gat)
    nullify(lnsec_nod)
    nullify(lnsec_set)
    nullify(nnset_par4)
    
    nnset_sum = nnset
    call PAR_SUM(nnset_sum)
    if( nnset_sum > 0 ) then

       !if( INOTMASTER ) then
       !   nnset_min = minval(lnset,lnset/=0)
       !   nnset_max = maxval(lnset)
       !end if
       if( INOTEMPTY ) then
          nnset_min = minval(lnset,lnset/=0)
          nnset_max = maxval(lnset)
       else
          nnset_min =  huge(1_ip)
          nnset_max = -huge(1_ip)      
       end if
       call PAR_MIN(nnset_min)
       call PAR_MAX(nnset_max)
              
       if( (nnset_max-nnset_min+1 /= nnset_sum) .or. (nnset_min .ne. 1_ip) ) call runend('NODE SETS MUST BE CONSECUTIVE')

       if( IMASTER ) then
          nnset = nnset_sum
          call par_memset(1_ip)
          call domain_memory_allocate('LNSEC')                                                        
          call memory_alloca(memor_dom,'LNSEC_NOD_GAT','par_merge_parallel_list',lnsec_nod_gat,nnset_sum)
          call memory_alloca(memor_dom,'LNSEC_SET_GAT','par_merge_parallel_list',lnsec_set_gat,nnset_sum)
          nnset = 0
       else
          call memory_alloca(memor_dom,'LNSEC_NOD','par_merge_parallel_list',lnsec_nod,nnset)
          call memory_alloca(memor_dom,'LNSEC_SET','par_merge_parallel_list',lnsec_set,nnset)
          nnset = 0
          do ipoin = 1,npoin
             if( lnset(ipoin) > 0 .and. THIS_NODE_IS_MINE(ipoin) ) then
                nnset            = nnset + 1
                lnsec_set(nnset) = lnset(ipoin)     ! Global set number
                lnsec_nod(nnset) = lninv_loc(ipoin) ! Global node number
             end if
          end do
       end if
       !
       ! Gather
       ! NNSET_PAR: number of hosted sets by partitions
       ! LNSEC_PAR(index in communication order,1): partition
       ! LNSEC_PAR(index in communication order,2): Global node number
       ! LNSEC(IPOIN): set number
       !
       call PAR_GATHER(nnset,nnset_par)
       if( IMASTER ) then
          isize4 = int(memory_size(nnset_par),4)
          call memory_alloca(memor_dom,'NNSET_PAR4','par_merge_parallel_list',nnset_par4,isize4)
          nnset_par4 = int(nnset_par,4)
       end if
       call PAR_GATHERV(lnsec_set,lnsec_set_gat,nnset_par4)
       call PAR_GATHERV(lnsec_nod,lnsec_nod_gat,nnset_par4)

       if( IMASTER ) then
          call memory_deallo(memor_dom,'NNSET_PAR4','par_merge_parallel_list',nnset_par4)
          nnset = 0
          do ipart = 1,PAR_CODE_SIZE-1          
             do ii = 1,nnset_par(ipart)
                nnset              = nnset + 1
                inset              = lnsec_set_gat(nnset)
                ipoin              = lnsec_nod_gat(nnset)
                lnsec_par(nnset,1) = ipart                 ! Partition
                lnsec_par(nnset,2) = inset                 ! Global set number                   
                lnsec(inset)       = ipoin                 ! Global node number
             end do
          end do
          call memory_deallo(memor_dom,'LNSEC_NOD_GAT','par_merge_parallel_list',lnsec_nod_gat)
          call memory_deallo(memor_dom,'LNSEC_SET_GAT','par_merge_parallel_list',lnsec_set_gat)             
       else
          !
          ! Recover local numbering for LNSEC
          !
          call memory_deallo(memor_dom,'LNSEC_NOD','par_merge_parallel_list',lnsec_nod)
          call memory_deallo(memor_dom,'LNSEC_SET','par_merge_parallel_list',lnsec_set)
       end if

    end if

  end subroutine reaset_node_parall

end module mod_reaset
