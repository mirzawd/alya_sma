!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_clusters.f90
!> @author  houzeaux
!> @date    2020-05-28
!> @brief   Clusters
!> @details Create clusters from a grapg
!-----------------------------------------------------------------------

module mod_clusters

  use def_kintyp_basic,   only : ip,rp,lg
  use def_kintyp_mesh,    only : mesh_type
  use def_master,         only : INOTEMPTY
  use def_master,         only : IPARALL
  use def_master,         only : kfl_paral
  use def_master,         only : intost
  use def_domain
  !use def_domain,         only : memor_dom
  use mod_graphs,         only : graphs_element_element_graph
  use mod_graphs,         only : graphs_poipoi
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_resize
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_ALLGATHERV
  use mod_communications, only : PAR_GHOST_ELEMENT_EXCHANGE
  use mod_messages,       only : messages_live
  use mod_maths,          only : maths_heap_sort

  implicit none
  private

  integer(ip), parameter :: CLUSTERS_ON_ELEMENTS = 1
  integer(ip), parameter :: CLUSTERS_ON_NODES    = 2

  public :: test_clusters
  public :: clusters_volume
  public :: mask_from_field
  public :: clusters_basic
  public :: clusters

contains

  subroutine test_clusters()

!!$    use mod_postpr
!!$    use def_kermod
!!$    integer(ip), pointer :: lmask(:),legro(:)
!!$    integer(ip)          :: ielem,pelty,inode,pnode,ipoin,nclus,iclus
!!$    real(rp),    pointer :: regro(:)
!!$    integer(ip), pointer :: lclus(:)
!!$    
!!$    nullify(lmask)
!!$    nullify(legro)
!!$    nullify(regro)
!!$    nullify(lclus)
!!$    call memory_alloca(memor_dom,'LMASK','clusters_elements',lmask,nelem) 
!!$    call memory_alloca(memor_dom,'REGRO','clusters_elements',regro,nelem) 
!!$
!!$    do ielem = 1,nelem
!!$       lmask(ielem) = int(xfiel(1) % a(1,ielem,1),ip)      
!!$    end do   
!!$!    do ielem = 1,nelem
!!$!       lmask(ielem) = 0
!!$!       pelty = ltype(ielem)
!!$!       if (pelty > 0) then
!!$!          pnode = nnode(pelty)
!!$!          do inode = 1,pnode
!!$!             ipoin = lnods(inode,ielem)
!!$!             if( xfiel(1) % a(1,ipoin,1) > 0.0_rp ) lmask(ielem) = 1
!!$!          end do
!!$!       end if
!!$!    end do
!!$    if( nelem > 0 ) nclus = sum(lmask)
!!$    call PAR_SUM(nclus)
!!$    if(INOTSLAVE) print*,'COUNT= ',nclus
!!$    
!!$    !call memory_alloca(memor_dom,'LEGRO','clusters_elements',legro,nelem) 
!!$    !do ielem = 1,nelem
!!$    !   legro(ielem) = lmask(ielem)
!!$    !end do
!!$    call clusters(lmask,legro,nclus,meshe(ndivi),ON_ELEMENTS=.true.)
!!$    !if(nelem>0)print*,'caca=',maxval(legro)
!!$    
!!$    call memory_alloca(memor_dom,'REGRO','clusters_elements',lclus,nclus+1,lboun=0_ip)   
!!$    do ielem = 1,nelem
!!$       iclus = legro(ielem)
!!$       lclus(iclus) = lclus(iclus) + 1
!!$    end do
!!$    call PAR_SUM(lclus)
!!$    if( INOTSLAVE ) then
!!$       do iclus = 1,nclus
!!$          print*,'cluster: ',iclus,lclus(iclus)
!!$       end do
!!$    end if
!!$    
!!$    do ielem = 1,nelem
!!$       regro(ielem) = real(legro(ielem),rp)
!!$    end do
!!$    call postpr(regro,(/'CLUST','SCALA','NELEM'/),1_ip,1.0_rp)

  end subroutine test_clusters
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-17
  !> @brief   Return clusters
  !> @details Given a mask, return clusters by renumbering them
  !>
  !>          Unity test:
  !>
  !>          integer(ip), pointer :: lmask(:)
  !>          integer(ip), pointer :: legro(:)
  !>          integer(ip)          :: ielem,iboun
  !>          
  !>          nullify(legro,lmask)
  !>          allocate(lmask(nelem))
  !>          do ielem = 1,nelem
  !>             lmask(ielem) = 0
  !>          end do
  !>          do ielem = 1,nelem
  !>             lmask(ielem) = int(xfiel(1) % a(1,ielem,1),ip)
  !>          end do
  !>          call clusters(lmask,legro,meshe(ndivi))
  !>  
  !-----------------------------------------------------------------------

  subroutine clusters(lmask,legro,nclus,meshe,ia,ja,ON_ELEMENTS,ON_NODES,GRAPH_CRITERION)

    integer(ip),     pointer,           intent(inout) :: lmask(:)
    integer(ip),     pointer,           intent(inout) :: legro(:)
    integer(ip),                        intent(out)   :: nclus
    type(mesh_type),          optional, intent(inout) :: meshe     !< Mesh type
    integer(ip),     pointer, optional, intent(in)    :: ia(:)
    integer(ip),     pointer, optional, intent(in)    :: ja(:)
    logical(lg),              optional, intent(in)    :: ON_ELEMENTS
    logical(lg),              optional, intent(in)    :: ON_NODES
    character(len=*),         optional, intent(in)    :: GRAPH_CRITERION
    integer(ip)                                       :: ielem,kelem,istack,kclus
    integer(ip)                                       :: nenti_marked,nstack,iz,jelem,iclus
    integer(ip)                                       :: iclus_max,lclus,icont,nenti_2,iiter
    integer(ip)                                       :: nclus_tot,nclus_ini,nclus_end,nenti
    integer(ip)                                       :: jclus
    integer(ip),     pointer                          :: ia_loc(:)
    integer(ip),     pointer                          :: ja_loc(:)
    integer(ip),     pointer                          :: global_num_loc(:)
    integer(ip),     pointer                          :: lstack(:)
    integer(ip),     pointer                          :: nclus_gat(:)
    integer(ip),     pointer                          :: permr(:)
    integer(ip),     pointer                          :: permr_gat(:)
    integer(ip),     pointer                          :: legro_send(:)
    integer(ip),     pointer                          :: legro_recv(:)
    integer(ip),     pointer                          :: permr_ptr(:)
    integer(ip)                                       :: wherein
    character(50)                                     :: w_GRAPH_CRITERION

    nullify(ia_loc)
    nullify(ja_loc)
    nullify(lstack)
    nullify(nclus_gat)
    nullify(permr)
    nullify(permr_gat)
    nullify(legro_send)
    nullify(legro_recv)
    wherein = CLUSTERS_ON_ELEMENTS
    iiter   = 0
    !
    ! Options
    !
    if( present(GRAPH_CRITERION) ) then
       w_GRAPH_CRITERION = trim(GRAPH_CRITERION)
    else
       w_GRAPH_CRITERION = 'INCLUDING HALO'
    end if

    !-----------------------------------------------------------------
    !
    ! Compute graphs if needed
    !
    !-----------------------------------------------------------------

    if( present(ia) .and. present(ja) ) then

       ia_loc  => ia
       ja_loc  => ja
       nenti_2 =  memory_size(ia)-1
       call runend('MOD_CLUSTERS: NUMBER OF ENTITIES SHOULD BE GIVEN')

    else if( present(meshe) ) then

       if( present(ON_ELEMENTS) ) then
          if( ON_ELEMENTS ) wherein = CLUSTERS_ON_ELEMENTS
       end if
       if( present(ON_NODES) ) then
          if( ON_NODES ) wherein = CLUSTERS_ON_NODES
       end if

       if( wherein == CLUSTERS_ON_ELEMENTS ) then
          nenti_2        =  meshe % nelem_2
          nenti          =  meshe % nelem
          global_num_loc => meshe % leinv_loc
          if( INOTEMPTY ) & 
               call graphs_element_element_graph(meshe,'SHARING FACES',trim(w_GRAPH_CRITERION),'INCLUDING DIAGONAL',ia_loc,ja_loc,memor=memor_dom)
       else
          call runend('CLUSTERS: NOT CODED')
          nenti_2        =  meshe % npoin_2
          nenti          =  meshe % npoin
          global_num_loc => meshe % lninv_loc
          if( INOTEMPTY ) & 
               call graphs_poipoi(meshe,ia_loc,ja_loc,IA_NAME='IA_LOC',JA_NAME='JA_LOC',INCLUDE_HALOS=.true.,memor=memor_dom)
       end if
    else
       call runend('CLUSTERS: SOME ARGUMENTS ARE MISSING')
    end if

    if( wherein == CLUSTERS_ON_ELEMENTS ) then
       call messages_live('IDENTIFY CLUSTERS ON ELEMENTS')
    else
       call messages_live('IDENTIFY CLUSTERS ON NODES')
    end if

    !-----------------------------------------------------------------
    !
    ! Mark initial clusters
    !
    !-----------------------------------------------------------------

    if( .not. associated(legro) ) call memory_alloca(memor_dom,'LEGRO','clusters_elements',legro,nenti_2)
    if( memory_size(lmask) < nenti_2 ) then
       call memory_resize(memor_dom,'LMASK','clusters_elements',lmask,nenti_2)
    end if

    call PAR_GHOST_ELEMENT_EXCHANGE(lmask,'SUBSTITUTE','IN MY CODE')

    do ielem = 1,nenti_2
       if( lmask(ielem) /= 0 ) legro(ielem) = 1
    end do

    if( nenti_2 > 0 ) then
       nenti_marked = count(legro(1:nenti_2)/=0)
    else
       nenti_marked = 0           
    end if

    !-----------------------------------------------------------------
    !
    ! Mark elements recursively
    !
    !-----------------------------------------------------------------

    call memory_alloca(memor_dom,'LSTACK'    ,'clusters_elements',lstack,nenti_2)
    nclus = 0
    kelem = 0

    do while( kelem /= nenti_marked ) 

       nclus = nclus + 1
       ielem = 1
       do while( lmask(ielem) == 0 )
          ielem = ielem + 1
          if( ielem > nenti_2 ) goto 10
       end do

       nstack       = 1
       lstack(1)    = ielem
       legro(ielem) = nclus
       lmask(ielem) = 0
       istack       = 0        
       kelem        = kelem + 1

       do while( istack /= nstack )

          istack = istack + 1   
          ielem  = lstack(istack)

          do iz = ia_loc(ielem),ia_loc(ielem+1)-1
             jelem = ja_loc(iz)
             if( lmask(jelem) /= 0 ) then
                legro(jelem)   = nclus
                nstack         = nstack + 1
                lstack(nstack) = jelem
                kelem          = kelem + 1
                lmask(jelem)   = 0
             end if
          end do
       end do

    end do
    call memory_deallo(memor_dom,'LSTACK','clusters_elements',lstack)
    nclus_ini = nclus
    call PAR_SUM(nclus_ini)
    nclus_end = nclus

    if( nclus_ini /= 0 ) then

       if( IPARALL ) then 
          !
          ! Renumber groups in a lexicographical way: PERMR
          !
          call memory_alloca(memor_dom,'NCLUS_GAT' ,'clusters_elements',nclus_gat,PAR_WORLD_SIZE,LBOUN=0_ip)
          call memory_alloca(memor_dom,'PERMR'     ,'clusters_elements',permr,nclus)

          call PAR_ALLGATHER(nclus,nclus_gat)
          kclus = sum(nclus_gat(1:kfl_paral-1))
          do iclus = 1,nclus
             permr(iclus) = kclus+iclus
          end do
          do ielem = 1,nenti_2
             iclus = legro(ielem)
             if( iclus /= 0 ) legro(ielem) = permr(iclus)
          end do

          nclus_tot = sum(nclus_gat)      
          call memory_alloca(memor_dom,'PERMR_GAT','clusters_elements',permr_gat,nclus_tot+1,lboun=0_ip)

          permr_ptr => permr_gat(1:)
          call PAR_ALLGATHERV(permr,permr_ptr,nclus_gat)

          call memory_deallo(memor_dom,'NCLUS_GAT','clusters_elements',nclus_gat)
          call memory_deallo(memor_dom,'PERMR'    ,'clusters_elements',permr)

          call memory_alloca(memor_dom,'LEGRO_SEND','clusters_elements',legro_send,nenti_2)
          call memory_alloca(memor_dom,'LEGRO_RECV','clusters_elements',legro_recv,nenti_2)

          !-----------------------------------------------------------------
          !
          ! Parallelization, take always maximum group
          !
          !-----------------------------------------------------------------

          icont = 1
          do while( icont /= 0 )

             icont = 0
             iiter = iiter + 1

             !if( iiter < 10 ) then
             !   call postpr_right_now('LEGR'//trim(intost(iiter)),'SCALA','NELEM',legro)
             !else
             !   call postpr_right_now('LEG'//trim(intost(iiter)) ,'SCALA','NELEM',legro)
             !end if

             do ielem = 1,nenti_2
                legro_send(ielem) = legro(ielem)
                legro_recv(ielem) = legro(ielem)
             end do

             call PAR_GHOST_ELEMENT_EXCHANGE(legro_recv,'MAX','IN MY CODE') ! I receive in my halo

             do ielem = nenti+1,nenti_2             
                kclus            = legro(ielem)
                lclus            = legro_recv(ielem)
                iclus            = permr_gat(kclus)
                jclus            = permr_gat(lclus)            
                iclus_max        = max(iclus,jclus,kclus,lclus)
                permr_gat(iclus) = iclus_max
                permr_gat(kclus) = iclus_max
                permr_gat(jclus) = iclus_max
                permr_gat(lclus) = iclus_max
             end do

             do ielem = 1,nenti_2
                iclus = legro(ielem)
                if( iclus /= permr_gat(iclus) ) then
                   icont = 1
                   legro(ielem) = permr_gat(iclus)
                end if
             end do

             call PAR_MAX(icont,'IN MY CODE')

          end do

       else
          nclus_tot = nclus
          call memory_alloca(memor_dom,'PERMR_GAT','clusters_elements',permr_gat,nclus+1,lboun=0_ip)

       end if

       !-----------------------------------------------------------------
       !
       ! Renumber clusters
       !
       !-----------------------------------------------------------------
       !
       ! Renumber in the order of appearance
       !
       permr_gat = 0
       do ielem = 1,nenti
          iclus = legro(ielem)
          permr_gat(iclus) = iclus
       end do
       call PAR_MAX(permr_gat)
       kclus = 0
       do iclus = 1,nclus_tot
          if( permr_gat(iclus) > 0 ) then
             kclus = kclus + 1
             permr_gat(iclus) = kclus
          end if
       end do

       do ielem = 1,nenti_2
          iclus = legro(ielem)
          legro(ielem) = permr_gat(iclus)
       end do

       nclus_end = count(permr_gat/=0)
       nclus     = nclus_end
       !
       ! Order according to global numbering
       !
       call memory_alloca(memor_dom,'PERMR','clusters_elements',permr,nclus+1,lboun=0_ip)
       do ielem = 1,nenti
          iclus = legro(ielem)
          if( iclus > 0 ) permr(iclus) = max(permr(iclus),global_num_loc(ielem))
       end do
       call PAR_MAX(permr)
       do iclus = 1,nclus
          permr_gat(iclus) = iclus
       end do
       call maths_heap_sort(2_ip,nclus,permr(1:),ivo1=permr_gat(1:))
       do iclus = 1,nclus
          permr(permr_gat(iclus)) = iclus 
       end do
       do ielem = 1,nenti_2
          iclus        = legro(ielem)
          legro(ielem) = permr(iclus) 
       end do
       call memory_deallo(memor_dom,'PERMR','clusters_elements',permr)

    end if
    !
    ! Messages, summary
    !
    if( IPARALL ) then
       call messages_live('   -> CLUSTERS COMPUTED WITH '//trim(intost(iiter))//' ITERATIONS')
       call messages_live('   -> INITIAL # OF CLUSTERS '//trim(intost(nclus_ini)))
       call messages_live('   -> FINAL   # OF CLUSTERS '//trim(intost(nclus_end)))
    else
       call messages_live('   -> # OF CLUSTERS '//trim(intost(nclus_ini)))
    end if
    !
    ! DEallocate
    !
    if( .not. ( present(ia) .and. present(ja) ) ) then       
       call memory_deallo(memor_dom,'IA_LOC','clusters_elements',ia_loc)
       call memory_deallo(memor_dom,'JA_LOC','clusters_elements',ja_loc)
    end if
    call memory_deallo(memor_dom,'PERMR_GAT' ,'clusters_elements',permr_gat)
    call memory_deallo(memor_dom,'LEGRO_SEND','clusters_elements',legro_send)
    call memory_deallo(memor_dom,'LEGRO_RECV','clusters_elements',legro_recv)


10  continue
    !call postpr_right_now('LEGRO','SCALA','NELEM',legro)

  end subroutine clusters

  
  !-----------------------------------------------------------------------
  !> 
  !> @author  Margarida Moragues
  !> @date    2020-03-20
  !> @brief   Computes clusters volume
  !> @details Given a list of clusters it computes its volume
  !>
  !-----------------------------------------------------------------------

  subroutine clusters_volume(legro,nclus,volume)

    use def_domain
    use def_master
    use mod_func
    use mod_integrals

    implicit none
    integer(ip),           pointer,       intent(in)    :: legro(:)
    integer(ip),                          intent(in)    :: nclus
    real(rp),                             intent(inout) :: volume(:)
    integer(ip)                                         :: ipoin, iclus
    type(func_ptr)                                      :: array_func(1,1)
    type(field_arrays)                                  :: array_fields(1)
    real(rp)                                            :: integrals(1,nclus)
    
    call func_initialization(array_func)

    allocate(array_fields(1)%a(1,npoin))
    do ipoin = 1,npoin
       array_fields(1) % a(1,ipoin) = 1.0_rp
    end do

    call integrals_volume(array_fields,array_func,legro,integrals)
    
    deallocate(array_fields(1)%a)

    do iclus = 1, nclus
       volume(iclus) = integrals(1,iclus)
    end do
    
  end subroutine clusters_volume

  !-----------------------------------------------------------------------
  !> 
  !> @author  Margarida Moragues
  !> @date    2020-03-20
  !> @brief   Computes mask from field
  !> @details Given a 1D field and a threshold, it computes a mask (field of 0 and 1).
  !>
  !-----------------------------------------------------------------------

  subroutine mask_from_field(lmask,field,field_threshold,ON_ELEMENTS,ON_NODES)
     implicit none
     integer(ip),              pointer,  intent(inout) :: lmask(:)
     real(rp),                 pointer,  intent(in)    :: field(:)
     real(rp),                           intent(in)    :: field_threshold
     logical(lg),              optional, intent(in)    :: ON_ELEMENTS
     logical(lg),              optional, intent(in)    :: ON_NODES
     integer(ip)                                       :: wherein
     integer(ip)                                       :: ipoin,ielem,igaus,inode,pelty,pnode,pgaus
     real(rp)                                          :: lmask_aux,gpcon(mgaus)
 
     wherein = CLUSTERS_ON_ELEMENTS
     if( present(ON_ELEMENTS) ) then
        if( ON_ELEMENTS ) wherein = CLUSTERS_ON_ELEMENTS
     else if( present(ON_NODES) ) then
        if( ON_NODES ) wherein = CLUSTERS_ON_NODES
     else 
        call runend('MOD_CLUSTERS: SPECIFY IF THE MASK SHOULD BE BUILD ON ELEMENTS OR ON NODES')
     end if
 
     if( wherein == CLUSTERS_ON_ELEMENTS ) then
        do ielem=1,nelem
           lmask_aux = 0.0_rp
           pelty = ltype(ielem)
           if (pelty > 0) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              do igaus = 1, pgaus
                 gpcon(igaus) = 0.0_rp
                 do inode = 1, pnode
                    ipoin = lnods(inode,ielem)
                    gpcon(igaus) = gpcon(igaus) + elmar(pelty) % shape(inode,igaus) * field(ipoin)
                 end do
                 lmask_aux = lmask_aux + gpcon(igaus)
              end do
              lmask_aux = lmask_aux / real(pgaus,rp)
              if( lmask_aux > field_threshold ) then
                 lmask(ielem) = 1
              else
                 lmask(ielem) = 0
              end if
           end if
        end do
     else
        do ipoin=1,npoin
           if( field(ipoin) > field_threshold ) then
              lmask(ipoin) = 1
           else
              lmask(ipoin) = 0
           end if
        end do
     end if
  end subroutine mask_from_field
 
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-28
  !> @brief   Define clusters
  !> @details Define a basic cluster from a graph. This is a sequential
  !>          task
  !> 
  !-----------------------------------------------------------------------

  subroutine clusters_basic(nn,ia,ja,nclus,lclus,lmask)

    integer(ip),                    intent(in)    :: nn
    integer(ip),           pointer, intent(in)    :: ia(:)
    integer(ip),           pointer, intent(in)    :: ja(:)
    integer(ip),                    intent(out)   :: nclus
    integer(ip),           pointer, intent(inout) :: lclus(:)
    logical(lg), optional, pointer, intent(inout) :: lmask(:)
    integer(ip),           pointer                :: lstack(:)
    integer(ip)                                   :: ii,jj,kk
    integer(ip)                                   :: nn_marked
    integer(ip)                                   :: istack,nstack,iz

    if( nn <= 0 ) return
    
    nullify(lstack)

    if( .not. associated(lclus) ) then
       call memory_alloca(memor_dom,'LCLUS','clusters_elements',lclus,nn,'DO_NOT_INITIALIZE')
    end if  
    
    if( present(lmask) ) then
       do ii = 1,nn
          if( lmask(ii) ) then
             lclus(ii) = 1
          else
             lclus(ii) = 0
          end if
       end do
    else
       do ii = 1,nn
          lclus(ii) = 1
       end do
    end if
    nn_marked = count(lclus>0)
    ! 
    ! Mark elements recursively
    !
    call memory_alloca(memor_dom,'LSTACK','clusters_elements',lstack,nn)
    
    nclus = 0
    kk    = 0

    do while( kk /= nn_marked ) 

       nclus = nclus + 1
       ii = 1
       do while( lclus(ii) <= 0 )
          ii = ii + 1
          if( ii > nn ) goto 10
       end do

       nstack    = 1
       lstack(1) = ii
       lclus(ii) = -nclus
       istack    = 0        
       kk        = kk + 1

       do while( istack /= nstack )

          istack = istack + 1   
          ii     = lstack(istack)

          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             if( lclus(jj) > 0 ) then
                lclus(jj)      = -nclus
                nstack         =  nstack + 1
                lstack(nstack) =  jj
                kk             =  kk + 1
             end if
          end do
       end do

    end do

    do ii = 1,nn
       lclus(ii) = abs(lclus(ii))
    end do
    
10  continue
    
    call memory_deallo(memor_dom,'LSTACK','clusters_elements',lstack)

  end subroutine clusters_basic

end module mod_clusters
!> @}
