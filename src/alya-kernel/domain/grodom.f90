!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Domain
!> @addtogroup Algebraic_Solver
!> @{
!> @file      grodom.f90
!> @author    houzeaux
!> @date      2018-11-06
!> @brief     Computing the groups for the defelated CG
!> @details   This routine computes the groups for the deflated CG.
!>            We have two possiblities:\n
!>            ITASK = 1 ... Sequential methods (done as a preprocess)
!>            ITASK = 2 ... Parallel methods
!> @} 
!-----------------------------------------------------------------------

subroutine grodom(itask)

  use def_kintyp,         only : ip,i1p,rp
  use def_master,         only : ISEQUEN,ISLAVE
  use def_master,         only : INOTMASTER,IMASTER
  use def_master,         only : READ_AND_RUN
  use def_master,         only : intost,kfl_paral
  use def_master,         only : npoin_par,npart
  use def_master,         only : lninv_loc
  use def_domain,         only : ndime,npoin
  use def_domain,         only : npoin_own
  use def_domain,         only : kfl_ngrou
  use def_domain,         only : ngrou_dom
  use def_domain,         only : ngrou_dom_target
  use def_domain,         only : coord
  use def_domain,         only : xfiel
  use def_domain,         only : lgrou_dom
  use def_domain,         only : r_dom
  use def_domain,         only : c_dom
  use def_domain,         only : memor_dom
  use def_domain,         only : ngrou_boxes_coarse
  use def_domain,         only : ngrou_boxes_fine
  use mod_parall,         only : PAR_COMM_MY_CODE_WM
  use def_parall,         only : num_repart_par
  use mod_memory,         only : memory_alloca  
  use mod_memory,         only : memory_deallo  
  use mod_memory,         only : memory_size  
  use mod_partitioning,   only : partitioning_frontal
  use mod_partitioning,   only : partitioning_frontal_materials
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_par_tools,      only : par_tools_gathered_graph
  use mod_partition_sfc,  only : partition_sfc
  use mod_messages,       only : messages_live
  use def_parall,         only : kfl_parseq_par
  use mod_parall,         only : PAR_PARALLEL_PARTITION 
  use mod_parall,         only : PAR_SEQUENTIAL_PARTITION 
  use mod_domain,         only : domain_memory_allocate
  use mod_domain,         only : domain_total_node_number
  use mod_std

  implicit none

  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ipart,iwarn
  integer(ip)             :: npoin_mesh,index_npoin
  integer(ip)             :: kgrou,igrou
  integer(ip), pointer    :: lninv_gat(:)
  integer(ip), pointer    :: lgrou_perm(:)
  integer(ip)             :: inul0(1)
  type(i1p),   pointer    :: ja_type(:)
  integer(ip)             :: boxes_coarse_sfc(3)
  integer(ip)             :: boxes_fine_sfc(3)
  !
  ! When shouldn't do anything
  !
  if( kfl_ngrou == 0 ) return
  !
  ! The number of groups calculated NGROU_DOM could be different from the target one NGROU_DOM_TARGET
  !
  nullify(ja_type)
  nullify(lninv_gat)
  nullify(lgrou_perm)
  
  if( itask == 1 .and. kfl_ngrou == -1 .and. IMASTER .and. .not. READ_AND_RUN() ) then

     !----------------------------------------------------------------------
     !
     ! Sequential Frontal approach, only called when Master has the whole
     ! mesh and if this is not a restart run
     !
     !----------------------------------------------------------------------
 
     call messages_live('COMPUTE GROUPS USING SEQUENTIAL FRONTAL METHOD')
     call domain_memory_allocate('LGROU_DOM')
     call partitioning_frontal(ngrou_dom,npoin,r_dom,c_dom,lgrou_dom)

  else if( itask == 2  .and. num_repart_par > 0 ) then

     !----------------------------------------------------------------------
     !
     ! Mesh has been repartitioned, no need to compute NGROU_DOM once more
     !
     !----------------------------------------------------------------------
     
     return
 
  else if( itask == 2 .and. num_repart_par == 0 ) then
     
     if( kfl_ngrou == -1  .and. kfl_parseq_par == PAR_SEQUENTIAL_PARTITION ) then

        !----------------------------------------------------------------------
        !
        ! Sequential Frontal approach
        !
        ! In parallel, groups were computed by Parall service.
        ! Just in case the master found another number of groups which does
        ! not corresponds to that read in readstr, NGROU_DOM
        ! is broadcast to slaves
        !
        !----------------------------------------------------------------------
        
        if( ISEQUEN ) then

           ngrou_dom = ngrou_dom_target
           call messages_live('COMPUTE GROUPS USING SEQUENTIAL FRONTAL METHOD')
           call domain_memory_allocate('LGROU_DOM')
           call partitioning_frontal(ngrou_dom,npoin,r_dom,c_dom,lgrou_dom)

           !call partitioning_frontal_materials(ngrou_dom,npoin,r_dom,c_dom,lgrou_dom)

       else

           call PAR_BROADCAST(ngrou_dom)
           return
           
        end if

     else if( kfl_ngrou == -2 .or. ( kfl_ngrou == -1 .and. kfl_parseq_par == PAR_PARALLEL_PARTITION ) ) then

        !----------------------------------------------------------------------
        !
        ! Parallel Frontal approach
        !
        ! NGROU_DOM is redefined just in case the frontal partitioning has
        ! found a different number of descent groups
        !
        !----------------------------------------------------------------------

        if( ISEQUEN ) then

           ngrou_dom = ngrou_dom_target
           call messages_live('COMPUTE GROUPS USING SEQUENTIAL FRONTAL METHOD')
           call domain_memory_allocate('LGROU_DOM')
           call partitioning_frontal(ngrou_dom,npoin,r_dom,c_dom,lgrou_dom)

        else
           !
           ! LGROU: allocate memory
           !
           call messages_live('COMPUTE GROUPS USING PARALLEL FRONTAL METHOD')
           if( INOTMASTER ) call domain_memory_allocate('LGROU_DOM')
           npoin_mesh = domain_total_node_number()
           !
           ! Construct groups: limit the number of groups to number of nodes
           !
           if( ngrou_dom >= npoin_mesh ) then

              do ipoin = 1,npoin
                 lgrou_dom(ipoin) = lninv_loc(ipoin)
              end do
              ngrou_dom = npoin_mesh

           else

              if( IMASTER ) then              
                 call memory_alloca(memor_dom,'LGROU_DOM','grodom',lgrou_dom,npoin_mesh)
              end if
              call par_tools_gathered_graph(ja_type,lninv_gat,MEMORY_COUNTER=memor_dom)

              if( IMASTER ) then
                 ngrou_dom = ngrou_dom_target
                 call partitioning_frontal(ngrou_dom,npoin_mesh,LPART=lgrou_dom,JA_TYPE=ja_type)
                 ngrou_dom = max(ngrou_dom,maxval(lgrou_dom(1:npoin)))
              end if
              call PAR_BROADCAST(ngrou_dom)
              if( IMASTER ) then
                 index_npoin = 0
                 do ipart = 1,npart
                    do ipoin = 1,npoin_par(ipart)
                       lninv_gat(index_npoin+ipoin) = lgrou_dom(lninv_gat(index_npoin+ipoin))
                    end do
                    call PAR_SEND_RECEIVE(npoin_par(ipart),0_ip,lninv_gat(index_npoin+1:),inul0,'IN MY CODE',ipart,'SYNCHRONOUS')
                    index_npoin = index_npoin + npoin_par(ipart)
                 end do
                 call memory_deallo(memor_dom,'LNINV_GAT','grodom',lninv_gat)        
                 call memory_deallo(memor_dom,'LGROU_DOM','grodom',lgrou_dom)        
                 call memory_deallo(memor_dom,'JA_TYPE'  ,'grodom',ja_type)
              else
                 if( npoin > 0 ) call PAR_SEND_RECEIVE(0_ip,npoin,inul0,lgrou_dom,'IN MY CODE',0_ip,'SYNCHRONOUS')
              end if

           end if
           
        end if

     else if( kfl_ngrou == -3 ) then

        !----------------------------------------------------------------------
        !
        ! 1 group per subdomaim. On non-own boundary nodes, take the group
        ! defined by my neighbors.
        !
        !----------------------------------------------------------------------

        call messages_live('COMPUTE GROUPS USING ONE GROUP PER PARTITION')
        ngrou_dom = max(1_ip,npart)
        if( INOTMASTER ) then
           call domain_memory_allocate('LGROU_DOM')
           if( ISLAVE ) then
              lgrou_dom(1:npoin_own) = kfl_paral
              call PAR_INTERFACE_NODE_EXCHANGE(lgrou_dom,'SUM')
           else
              lgrou_dom = 1
           end if
        end if

     else if( kfl_ngrou == -5 ) then

        !----------------------------------------------------------------------
        !
        ! 1 group per bunch of subdomains. On non-own boundary nodes, take the group
        ! defined by my neighbors.
        ! GROUP = floor((kfl_paral-1)/npart*ngrou_dom) + 1 
        !
        !----------------------------------------------------------------------

        call messages_live('COMPUTE GROUPS USING ONE GROUP PER PARTITION')

        ngrou_dom = min(ngrou_dom,npart)
        igrou     = floor(real(kfl_paral-1,rp)/real(npart,rp)*real(ngrou_dom,rp),kind=ip) + 1
        igrou     = max(1_ip,min(ngrou_dom,igrou))
        if( INOTMASTER ) then
           call domain_memory_allocate('LGROU_DOM')
           if( ISLAVE ) then
              lgrou_dom(1:npoin_own) = igrou
              call PAR_INTERFACE_NODE_EXCHANGE(lgrou_dom,'SUM')
           else
              lgrou_dom = 1
           end if
        end if

     else if( kfl_ngrou == -4 ) then

        !----------------------------------------------------------------------
        !
        ! Partitioning with SFC
        !
        !----------------------------------------------------------------------

        call messages_live('COMPUTE GROUPS USING A HILBERT SPACE FILLING CURVE IN PARALLEL')
        boxes_coarse_sfc = ngrou_boxes_coarse
        boxes_fine_sfc   = ngrou_boxes_fine

        if( INOTMASTER ) then
           call domain_memory_allocate('LGROU_DOM')
#ifdef __PGI
           if( ISEQUEN ) then
              call partition_sfc(lepar_=lgrou_dom,npart_sfc_=ngrou_dom,lenti_=coord,&
                   nenti_=npoin_own,sdim_=ndime,lweig_=std_int_1,boxes_coarse_=boxes_coarse_sfc,&
                   boxes_fine_=boxes_fine_sfc) 
           else
              call partition_sfc(lepar_=lgrou_dom,npart_sfc_=ngrou_dom,lenti_=coord,&
                   nenti_=npoin_own,sdim_=ndime,lweig_=std_int_1,boxes_coarse_=boxes_coarse_sfc,&
                   boxes_fine_=boxes_fine_sfc,lcorr_=std_real_1,PAR_COMM_=PAR_COMM_MY_CODE_WM) 
           end if
#else
           if( ISEQUEN ) then
              call partition_sfc(lgrou_dom,ngrou_dom,coord,npoin_own,ndime,  &
                   boxes_coarse_=boxes_coarse_sfc,boxes_fine_=boxes_fine_sfc)
           else
              call partition_sfc(lgrou_dom,ngrou_dom,coord,npoin_own,ndime,  &
                   boxes_coarse_=boxes_coarse_sfc,boxes_fine_=boxes_fine_sfc,&
                   PAR_COMM_=PAR_COMM_MY_CODE_WM)
           end if
#endif
           call PAR_INTERFACE_NODE_EXCHANGE(lgrou_dom,'SUM')
        end if

     else if( kfl_ngrou > 0 ) then

        !----------------------------------------------------------------------
        !
        ! Groups coming from a field
        !
        !----------------------------------------------------------------------

        call messages_live('COMPUTE GROUPS USING FIELD '//trim(intost(kfl_ngrou)))
        if( INOTMASTER ) then
           call domain_memory_allocate('LGROU_DOM')
           do ipoin = 1,npoin
              lgrou_dom(ipoin) = int(xfiel(kfl_ngrou) % a(1,ipoin,1))
           end do
           
        end if

     else

        return
        
     end if
     
  else

     return
     
  end if
  
  !----------------------------------------------------------------------------
  !
  ! Compact list of groups if necessary LGROU_DOM and recompute NGROU_DOM
  !
  !----------------------------------------------------------------------------
  !
  ! Check number of groups... could be unknown, e.g. if it is read from a field.
  !  
  if( itask /= 1 ) then
     ngrou_dom = 0
     if( npoin > 0 ) ngrou_dom = maxval(lgrou_dom)
     call PAR_MAX(ngrou_dom)
  end if
  !
  ! Mark non-empty groups
  !     
  call memory_alloca(memor_dom,'LGROU_PERM','grodom',lgrou_perm,ngrou_dom)
  if( npoin > 0 .and. associated(lgrou_dom) ) then
     do ipoin = 1,npoin
        igrou = lgrou_dom(ipoin)
        if( igrou > 0 ) lgrou_perm(igrou) = lgrou_perm(igrou) + 1
     end do
  end if
  
  if( itask /= 1 ) call PAR_SUM(lgrou_perm)
  !
  ! Remove empty groups and recompute number of groups
  !
  if( npoin > 0 .and. associated(lgrou_dom) ) then
     kgrou = 0
     do igrou = 1,ngrou_dom
        if( lgrou_perm(igrou) /= 0 ) then
           kgrou = kgrou + 1
           lgrou_perm(igrou) = kgrou
        end if
     end do
     ngrou_dom = kgrou
  end if
  
  if( itask /= 1 ) call PAR_MAX(ngrou_dom)
  !
  ! Reassign group if necessary
  ! When LGROU_DOM(IPOIN) = -1, node is not assigned any group
  ! and is understood to be prescribed
  !
  iwarn = 0
  if( npoin > 0 .and. associated(lgrou_dom) ) then
     do ipoin = 1,npoin
        igrou = lgrou_dom(ipoin)
        if( igrou > 0 ) then
           lgrou_dom(ipoin) = lgrou_perm(igrou)
        else if( igrou /= -1 ) then
           iwarn = 1
           lgrou_dom(ipoin) = ngrou_dom
        end if
     end do
  end if

  if( itask /=  1 ) call PAR_MAX(iwarn)
  call memory_deallo(memor_dom,'LGROU_PERM','grodom',lgrou_perm)
  !
  ! Warning message of the number of groups does not coincide
  ! with the one prescribed
  !
  if( ngrou_dom /= ngrou_dom_target ) &
       call messages_live('NUMBER OF GROUPS COMPUTED '//trim(intost(ngrou_dom))//&
       ' IS DIFFERENT FROM PRESCRIBED ONE '//trim(intost(ngrou_dom_target)),'WARNING')
  if( iwarn /= 0 ) &
       call messages_live('SOME NODES DID NOT HAVE ANY ASSIGNED GROUPS','WARNING')
  
end subroutine grodom
