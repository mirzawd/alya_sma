!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_parall_openmp.f90
!> @author  houzeaux
!> @date    2018-12-11
!> @brief   OpenMP module
!> @details All you need to run OpenMP and OmpSs, of course!
!>          
!-----------------------------------------------------------------------
module mod_parall_openmp

  use def_kintyp,     only : ip,rp,lg
  use def_kintyp,     only : ompss_domain
  use def_master,     only : INOTMASTER
  use def_master,     only : intost
  use def_domain,     only : mesh_type
  use def_domain,     only : ngaus,memor_dom
  use mod_memory,     only : memory_alloca
  use mod_memory,     only : memory_deallo
  use mod_graphs,     only : graphs_eleele
  use mod_graphs,     only : graphs_eleele_deallocate
  use mod_graphs,     only : graphs_domlis
  use mod_graphs,     only : graphs_coloring_greedy
  use mod_graphs,     only : graphs_coloring
  use mod_graphs,     only : graphs_coloring_deallocate
  use mod_graphs,     only : graphs_deallocate
  use mod_maths,      only : maths_heap_sort

  use mod_parall,     only : par_memor
  use mod_parall,     only : par_omp_num_blocks             ! Size of blocks to compute chunks n/par_omp_num_blocks
  use mod_parall,     only : par_omp_granularity            ! Granularity to compute par_omp_num_blocks=par_omp_num_threads*par_omp_granularity
  use mod_parall,     only : par_omp_nelem_chunk            ! Element loop chunk size
  use mod_parall,     only : par_omp_npoin_chunk            ! Node loop chunk size
  use mod_parall,     only : par_omp_nboun_chunk            ! Boundary loop chunk size
  use mod_parall,     only : par_omp_num_threads            ! Number of openmp threads
  use mod_parall,     only : par_omp_coloring_alg           ! Coloring algorithm
  use mod_parall,     only : par_omp_partition_alg          ! Partitioning algorithm for OmpSs (use same nomenclature for mesh partitioning)
  use mod_parall,     only : par_omp_num_colors             ! Element: Number of colors
  use mod_parall,     only : par_omp_list_colors            ! Element: Element colors
  use mod_parall,     only : par_omp_ia_colors              ! Element: Linked list IA for colors
  use mod_parall,     only : par_omp_ja_colors              ! Element: Linked list JA for colors
  use mod_parall,     only : par_omp_nboun_num_colors       ! Boundary: Number of colors
  use mod_parall,     only : par_omp_nboun_list_colors      ! Boundary: Element colors
  use mod_parall,     only : par_omp_nboun_ia_colors        ! Boundary: Linked list IA for colors
  use mod_parall,     only : par_omp_nboun_ja_colors        ! Boundary: Linked list JA for colors
  use mod_parall,     only : PAR_HYBRID_OFF         
  use mod_parall,     only : PAR_OPENMP_COLORING    
  use mod_parall,     only : PAR_OPENMP_NO_COLORING 
  use mod_parall,     only : PAR_OMPSS              
  use mod_parall,     only : par_hybrid
  use mod_messages,   only : livinf
  use mod_alya2metis, only : alya2metis_METIS_PartGraph
  use mod_std
  implicit none

  private

  public :: parall_openmp_partition_and_adjacency_ompss
  public :: parall_openmp_chunk_sizes
  public :: parall_openmp_coloring
  public :: parall_openmp_adjacency_ompss_unity_test
  
contains

  !-----------------------------------------------------------------------
  !
  !> @author  Antoni Artigues
  !> @brief   Partition a graph and creates ompss graph adjacency
  !> @details Partition the element graph using and creates ompss graph adjacency
  !
  !-----------------------------------------------------------------------

  subroutine parall_openmp_partition_and_adjacency_ompss(chunk_size,meshe,ompss_domains,ON_BOUNDARIES)

    integer(ip),          intent(in)  :: chunk_size    !< Number of blocks 
    type(mesh_type),      intent(in)  :: meshe    !< Mesh type
    TYPE(ompss_domain), pointer, intent(inout)   :: ompss_domains(:) !< Inverse permutation new = invpe(old) and adjacency

    integer(ip)                       :: nelem          
    integer(ip)                       :: npoin          
    integer(ip)                       :: mnode          
    integer(ip), pointer              :: lnods(:,:) 
    integer(ip), pointer              :: lnnod(:)       
    integer(ip), pointer              :: ltype(:)       
    logical(lg), optional             :: ON_BOUNDARIES

    integer(ip)                       :: ielem
    integer(ip)                       :: nedge       
    integer(ip)                       :: medge      
    integer(ip)                       :: kfl_weigh 
    integer(ip)                       :: ipart
    integer(ip)                       :: mepoi
    integer(ip)                       :: wedge(1)   
    integer(ip), pointer              :: pelel(:)    
    integer(ip), pointer              :: lelel(:)    
    integer(ip), pointer              :: pelpo(:)    
    integer(ip), pointer              :: lelpo(:)    
    integer(ip), pointer              :: wvert(:)
    integer(ip), pointer              :: lepar(:)
    integer(ip), pointer              :: lnpar(:)
    integer(ip), pointer              :: npoin_dom(:)
    integer(ip), pointer              :: lneigh_dom(:)

    integer(ip), pointer              :: nelem_part(:)
    integer(ip), pointer              :: pelem_part(:)


    !--ADJACENCY VARS
    integer(ip), pointer              :: neighDomOmpss(:)
    integer(ip)                       :: inter
    integer(ip)                       :: front
    integer(ip)                       :: inode,dom2,ii,jj,kk
    integer(ip)                       :: ndomi,domin,dom1
    integer(ip), pointer              :: domli(:) => null()
    integer(ip)                       :: npart
    !--END ADJACENCY VARS
    logical(lg)                       :: if_on_boundaries

    if( present(ON_BOUNDARIES) ) then
       if_on_boundaries = ON_BOUNDARIES
    else
       if_on_boundaries = .false.
    end if
    
    if( if_on_boundaries ) then
       nelem  = meshe % nboun
       npoin  = meshe % npoin
       mnode  = meshe % mnodb
       lnods => meshe % lnodb
       lnnod => meshe % lnnob
       ltype => meshe % ltypb
    else
       nelem  = meshe % nelem
       npoin  = meshe % npoin
       mnode  = meshe % mnode
       lnods => meshe % lnods
       lnnod => meshe % lnnod
       ltype => meshe % ltype
    end if

    if( nelem      <= 0     ) return
    if( chunk_size <= 0     ) call runend('graphs_partition_and_adjacency_ompss: PRESCRIBE AN ELEMENT CHUNK SIZE')
    if( chunk_size  > nelem ) call runend('graphs_partition_and_adjacency_ompss: ELEMENT CHUNK IS TOO HIGH')

    npart = max(nelem / chunk_size,1_ip)

    nullify(pelpo)
    nullify(lelpo)
    nullify(pelel)
    nullify(lelel)
    nullify(wvert)
    nullify(lepar)
    nullify(lnpar)
    nullify(npoin_dom)
    nullify(lneigh_dom)
    nullify(nelem_part)
    nullify(pelem_part)
    !
    ! Weights
    !
    allocate(wvert(nelem))
    allocate(lepar(nelem))
    do ielem = 1,nelem
       wvert(ielem) = ngaus(abs(ltype(ielem)))
    end do
    kfl_weigh = 2
    wedge     = 0
    !
    ! Compute element graph
    !
#ifdef __PGI
    call graphs_eleele(&
         nelem,npoin,mnode,mepoi,lnods,lnnod,&
         pelpo,lelpo,nedge,medge,pelel,lelel)
#else
    call graphs_eleele(&
         nelem,npoin,mnode,mepoi,lnods,lnnod,&
         pelpo,lelpo,nedge,medge,pelel,lelel,memor=memor_dom)
#endif
    !
    ! Partitioning
    !
    call alya2metis_METIS_PartGraph(npart,nelem,pelel,lelel,wvert,lepar)
    !
    ! Compute adjacency graph following ompss type
    !
    allocate( ompss_domains( npart ) )
    allocate( nelem_part(npart) )
    allocate( pelem_part(npart) )
    nelem_part = 0 
    pelem_part = 0

    do ielem = 1,nelem
       ipart = lepar(ielem)
       nelem_part(ipart) = nelem_part(ipart) + 1
    end do
    do ipart = 1,npart
       nullify  ( ompss_domains( ipart ) % elements ) 
       allocate ( ompss_domains( ipart ) % elements (nelem_part(ipart)))
       ompss_domains( ipart ) % elemIdx = 1
    end do

    do ielem = 1,nelem
       ipart             = lepar(ielem)
       ompss_domains( ipart ) % elements (ompss_domains( ipart ) % elemIdx) = ielem
       ompss_domains( ipart ) % elemIdx = ompss_domains( ipart ) % elemIdx + 1       
    end do
    !
    ! Create a half-matrix with the neighbourhood domains (including main diagonal)
    ! 
    allocate(neighDomOmpss((npart*(npart+1))/2) )
    allocate(npoin_dom(npart))
    do dom1 = 1,npart
       do dom2 = 1,dom1
          neighDomOmpss( (dom1*(dom1-1))/2 + dom2) = 0
       end do
       npoin_dom(dom1) = 0
    end do


    call memory_alloca(par_memor,'DOMLI','par_domgra' , domli , mepoi )
    allocate(lnpar(npoin))
    allocate(lneigh_dom(npart))


    front = 0
    inter = 0

    do inode = 1,npoin
       call graphs_domlis( pelpo, lelpo, inode, lepar, ndomi, domli )
       if( ndomi == 1 ) then
          !
          ! This is an internal node
          !
          inter            = inter + 1
          domin            = domli(1)
          lnpar(inode) = domin
          npoin_dom(domin) = npoin_dom(domin) + 1      
       else
          !
          ! This is a boundary node
          ! neighDomOmpss(i,j)
          !      if i.ne.j   boundary nodes shared by domains i and j
          !      if i.eq.j   boundary nodes of domain i
          !
          front = front + 1
          do ii = 1,ndomi
             dom1             = domli(ii)
             lnpar(inode) = 0
             npoin_dom(dom1)  = npoin_dom(dom1) + 1
             do jj = 1,ii
                dom2 = domli(jj)
                if( dom1 > dom2 ) then
                   kk = (dom1*(dom1-1))/2 + dom2
                else
                   kk = (dom2*(dom2-1))/2 + dom1
                end if
                neighDomOmpss(kk) = neighDomOmpss(kk) + 1          
             end do
          end do
       end if
    end do
    call memory_deallo(par_memor,'DOMLI','par_domgra' , domli )
    !
    ! Count the number of neighbours per domain 
    !
    do dom1 = 1, npart
       lneigh_dom(dom1) = 0
    enddo
    do dom1 = 1, npart
       kk = (dom1*(dom1-1))/2
       do dom2 = 1, dom1-1
          if( neighDomOmpss(kk + dom2) /= 0 ) then
             lneigh_dom(dom1) = lneigh_dom(dom1) + 1
             lneigh_dom(dom2) = lneigh_dom(dom2) + 1
          end if
       end do
    end do
    !
    ! Build the domain interconnection graph     
    !   
    do dom1 = 1,npart
       lneigh_dom(dom1) = lneigh_dom(dom1) + 1
       nullify  ( ompss_domains( dom1 ) % neighbours )
       allocate ( ompss_domains( dom1 ) % neighbours(lneigh_dom(dom1)) )
       ompss_domains( dom1 ) % neighbours(1) = dom1
       ompss_domains( dom1 ) % neighIdx = 2
    end do
    !
    ! Create and order adjacencies 
    !  
    do dom1= 1, npart
       kk = (dom1*(dom1-1))/2
       do dom2= 1, dom1-1
          if( neighDomOmpss(kk + dom2) /= 0 ) then            
             ompss_domains( dom1 ) % neighbours (ompss_domains( dom1 ) % neighIdx) = dom2
             ompss_domains( dom2 ) % neighbours (ompss_domains( dom2 ) % neighIdx) = dom1
             ompss_domains( dom1 ) % neighIdx = ompss_domains( dom1 ) % neighIdx + 1
             ompss_domains( dom2 ) % neighIdx = ompss_domains( dom2 ) % neighIdx + 1
          end if
       end do
    end do
    do dom1 = 1,npart
       call maths_heap_sort(2_ip,lneigh_dom(dom1),ompss_domains( dom1 ) % neighbours)
    end do
    !
    ! END Compute adjacency graph following ompss type
    !
    !
    !
    ! Deallocate memory
    !
    call memory_deallo(par_memor,'pelpo        ','graphs_partition_and_adjacency_ompss',pelpo)
    call memory_deallo(par_memor,'lelpo        ','graphs_partition_and_adjacency_ompss',lelpo)
    call memory_deallo(par_memor,'pelel        ','graphs_partition_and_adjacency_ompss',pelel)
    call memory_deallo(par_memor,'lelel        ','graphs_partition_and_adjacency_ompss',lelel)
    call memory_deallo(par_memor,'wvert        ','graphs_partition_and_adjacency_ompss',wvert)
    call memory_deallo(par_memor,'lepar        ','graphs_partition_and_adjacency_ompss',lepar)
    call memory_deallo(par_memor,'nelem_part   ','graphs_partition_and_adjacency_ompss',nelem_part)
    call memory_deallo(par_memor,'pelem_part   ','graphs_partition_and_adjacency_ompss',pelem_part)
    call memory_deallo(par_memor,'lnpar        ','graphs_partition_and_adjacency_ompss',lnpar)
    call memory_deallo(par_memor,'npoin_dom    ','graphs_partition_and_adjacency_ompss',npoin_dom)
    call memory_deallo(par_memor,'lneigh_dom   ','graphs_partition_and_adjacency_ompss',lneigh_dom)
    call memory_deallo(par_memor,'neighDomOmpss','graphs_partition_and_adjacency_ompss',neighDomOmpss)

  end subroutine parall_openmp_partition_and_adjacency_ompss

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   XSR to COO format
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  subroutine parall_openmp_output_ompss_domain(ompss_domains,meshe)

    TYPE(ompss_domain), intent(in)  :: ompss_domains(:)
    type(mesh_type),    intent(in)  :: meshe    !< Mesh type

    integer(ip), DIMENSION(:),  ALLOCATABLE   :: ia
    integer(ip), DIMENSION(:),  ALLOCATABLE   :: ja
    integer(ip), DIMENSION(:),  ALLOCATABLE   :: nodesMark
    real(rp),    DIMENSION(:,:),ALLOCATABLE   :: coordAux
    integer(ip)                               :: NUM_DOMAINS
    integer(ip)                               :: NUM_NEIGHBOURS, nodeSize, numNodesMarket
    integer(ip)                               :: idx1, idx2, I, J, ielem, pelty, pnode, inode, ipoin
    real(rp)                                  :: x,y,z

    write(*,*) 'LLEGOO1\n'
    !--IA AND JA CALCULATION
    NUM_DOMAINS = SIZE(ompss_domains)
    allocate(ia(NUM_DOMAINS+1))

    NUM_NEIGHBOURS = 0
    DO I = 1, NUM_DOMAINS
       NUM_NEIGHBOURS = NUM_NEIGHBOURS + SIZE(ompss_domains(I) % neighbours)      
    END DO
    allocate(ja(NUM_NEIGHBOURS))
    NUM_NEIGHBOURS = 0
    idx1 = 1
    idx2 = 1
    ia(1) = 1
    DO I = 1, NUM_DOMAINS
       NUM_NEIGHBOURS = SIZE(ompss_domains(I) % neighbours)      
       idx2 = ia(I)
       ia(I+1) = idx1 + NUM_NEIGHBOURS
       idx1 = idx1 + NUM_NEIGHBOURS
       DO J = 1, NUM_NEIGHBOURS
          ja(idx2) = ompss_domains(I) % neighbours(J)
          idx2 = idx2 + 1 
       END DO
    END DO



    !--COORDS CALCULATION
    allocate(coordAux(meshe % ndime, NUM_DOMAINS))

    nodeSize = SIZE(meshe % coord,2)
    allocate(nodesMark(nodeSize))
    DO I = 1, NUM_DOMAINS
       !--marcar los nodos del subdom
       nodesMark = 0
       NUM_NEIGHBOURS = SIZE(ompss_domains(I) % neighbours)
       DO J = 1, NUM_NEIGHBOURS
          ielem = ompss_domains(I) % neighbours(J)
          pelty = meshe % ltype(ielem)
          if( pelty > 0 ) then
             pnode = meshe % lnnod(ielem)
             do inode = 1,pnode
                ipoin = meshe % lnods(inode,ielem)
                nodesMark(ipoin) = 1
             end do
          end if
       END DO

       !--calcular media de los nodos marcados
       x = 0
       y = 0
       z = 0
       numNodesMarket = 0
       DO inode = 1, nodeSize
          IF (nodesMark(inode) == 1) THEN
             numNodesMarket = numNodesMarket + 1
             if (meshe % ndime == 1) then
                x = x + meshe % coord(1,inode)
             end if
             if (meshe % ndime == 2) then
                x = x + meshe % coord(1,inode)
                y = y + meshe % coord(2,inode)       
             end if
             if (meshe % ndime == 3) then
                x = x + meshe % coord(1,inode)
                y = y + meshe % coord(2,inode)
                z = z + meshe % coord(3,inode)
             end if
          END IF
       END DO
       x = x / real(numNodesMarket,rp)
       y = y / real(numNodesMarket,rp)
       z = z / real(numNodesMarket,rp)
       if (meshe % ndime == 1) then
          coordAux(1,I) = x
       end if
       if (meshe % ndime == 2) then
          coordAux(1,I) = x
          coordAux(2,I) = y      
       end if
       if (meshe % ndime == 3) then
          coordAux(1,I) = x
          coordAux(2,I) = y
          coordAux(3,I) = z  
       end if

    END DO

    !write(*,*) 'LLEGOO4\n'
    !call graphs_output(NUM_DOMAINS, ia, ja, coordAux, 2672_ip, 'grafo_deps_dominio10')

  end subroutine parall_openmp_output_ompss_domain

  !------------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    27/10/2015
  !> @brief   Openm chunk sizes
  !> @details If they are not prescribed by the user, compute the chunk
  !>          size for elements, boundaries and nodes. For elements,
  !>          the folloiwng formula is used:
  !>
  !>          OMP PARALLEL DO SCHEDULE (DYNAMIC,n/par_omp_num_blocks) &
  !>          OMP PRIVATE  (i) &
  !>          OMP SHARED   (n)
  !>          do i = 1,n
  !>          ...
  !>          end do
  !>                                       n
  !>          Chunk size is:   -------------------------
  !>                           num_threads x granularity  (=num_blocks)
  !------------------------------------------------------------------------

  subroutine parall_openmp_chunk_sizes(meshe)
    type(mesh_type), intent(in) :: meshe !< Mesh
    integer(ip)                 :: kelem
    integer(ip)                 :: icolo

    if( par_omp_num_threads <= 0 ) return

    if( par_omp_granularity < 1 ) then
       call runend('PAR_OMP_NUMBER_BLOCKS: WRONG GRANULARITY')
    else

       par_omp_num_blocks = max(1_ip,par_omp_num_threads) * par_omp_granularity

       if( par_omp_nelem_chunk == 0 ) then
          if( par_hybrid == PAR_OPENMP_COLORING ) then
             !
             ! With coloring: average number of elements
             !
             if( INOTMASTER ) then
                kelem = 0
                do icolo = 1,par_omp_num_colors
                   kelem = kelem + par_omp_ia_colors(icolo+1)-par_omp_ia_colors(icolo)
                end do
                kelem = kelem / par_omp_num_colors
                par_omp_nelem_chunk = kelem / par_omp_num_blocks
             else
                par_omp_nelem_chunk = meshe % nelem / par_omp_num_blocks
             end if
          else
             !
             ! No coloring
             !
             par_omp_nelem_chunk = meshe % nelem / par_omp_num_blocks
          end if
       end if

       if( par_omp_nboun_chunk == 0 ) then
          if( par_hybrid == PAR_OPENMP_COLORING ) then
             !
             ! With coloring: average number of elements
             !
             if( INOTMASTER ) then
                kelem = 0
                do icolo = 1,par_omp_nboun_num_colors
                   kelem = kelem + par_omp_nboun_ia_colors(icolo+1)-par_omp_nboun_ia_colors(icolo)
                end do
                kelem = kelem / par_omp_nboun_num_colors
                par_omp_nboun_chunk = kelem / par_omp_num_blocks
             else
                par_omp_nboun_chunk = meshe % nboun / par_omp_num_blocks
             end if
          else
             !
             ! No coloring
             !
             par_omp_nboun_chunk = meshe % nboun / par_omp_num_blocks
          end if
       end if

       !if( par_omp_nboun_chunk == 0 ) par_omp_nboun_chunk = meshe % nboun / par_omp_num_blocks
       if( par_omp_npoin_chunk == 0 ) par_omp_npoin_chunk = meshe % npoin / par_omp_num_blocks
       !
       ! Impose minimum and maximum values for chunks
       !
       par_omp_nelem_chunk = min( par_omp_nelem_chunk , meshe % nelem / par_omp_num_threads )
       par_omp_nboun_chunk = min( par_omp_nboun_chunk , meshe % nboun / par_omp_num_threads )
       par_omp_npoin_chunk = min( par_omp_npoin_chunk , meshe % npoin / par_omp_num_threads )

       par_omp_nelem_chunk = max( par_omp_nelem_chunk , 1_ip )
       par_omp_nboun_chunk = max( par_omp_nboun_chunk , 1_ip )
       par_omp_npoin_chunk = max( par_omp_npoin_chunk , 1_ip )

    end if

  end subroutine parall_openmp_chunk_sizes

  !------------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    20/05/2014
  !> @brief   Coloring
  !> @details Color the elements and create a liked list for the colored
  !>          elements
  !------------------------------------------------------------------------

  subroutine parall_openmp_coloring(meshe,mepoi_opt,pelpo_opt,lelpo_opt,ON_BOUNDARIES)

#ifdef __PGI
#define MEMPGI )
#else
#define MEMPGI ,memor=memor_dom)
#endif

    type(mesh_type),                     intent(in)    :: meshe         !< Mesh
    integer(ip),     optional,           intent(inout) :: mepoi_opt     !< Max nb of nodes per element
    integer(ip),     optional, pointer,  intent(inout) :: pelpo_opt(:)  !< Element-node linked list
    integer(ip),     optional, pointer,  intent(inout) :: lelpo_opt(:)  !< Element-node linked list
    logical(lg),     optional, intent(in)              :: ON_BOUNDARIES !< If coloring should be carried out on boundaries
    integer(ip)                                        :: nedge
    integer(ip)                                        :: medge
    integer(ip)                                        :: mepoi
    integer(ip)                                        :: ielem
    integer(ip)                                        :: nenti
    integer(ip),     pointer                           :: lelel(:)
    integer(ip),     pointer                           :: pelel(:)
    integer(ip),     pointer                           :: lelpo(:)
    integer(ip),     pointer                           :: pelpo(:)
    logical(lg)                                        :: if_on_boundaries
    !
    ! Nullify pointers
    !
    nullify(lelel,pelel)
    nullify(lelpo,pelpo)

    if( present(ON_BOUNDARIES) ) then
       if_on_boundaries = ON_BOUNDARIES
    else
       if_on_boundaries = .false.       
    end if
    if( if_on_boundaries ) then
       nenti = meshe % nboun
    else
       nenti = meshe % nelem
    end if
    if( present(mepoi_opt) ) mepoi = mepoi_opt

    if( INOTMASTER ) then
       if( par_omp_num_threads /= 0 .and. nenti /= 0 ) then

          if( if_on_boundaries ) then
             call livinf(0_ip,'COLOR BOUNDARIES FOR OPENMP',0_ip)
          else
             call livinf(0_ip,'COLOR ELEMENTS FOR OPENMP',0_ip)
          end if
          !
          ! Element-element graph
          !
          if( present(lelpo_opt) ) lelpo => lelpo_opt
          if( present(pelpo_opt) ) pelpo => pelpo_opt

          if( if_on_boundaries ) then
             call graphs_eleele(&
                  meshe % nboun,meshe % npoin,meshe % mnodb,mepoi,&
                  meshe % lnodb,meshe % lnnob,pelpo,lelpo,nedge,  &
                  medge,pelel,lelel MEMPGI
          else
             call graphs_eleele(&
                  meshe % nelem,meshe % npoin,meshe % mnode,mepoi,&
                  meshe % lnods,meshe % lnnod,pelpo,lelpo,nedge,  &
                  medge,pelel,lelel MEMPGI
          end if
          if( present(mepoi_opt) ) mepoi_opt = mepoi
          !
          ! Color elements
          !
          if( par_omp_coloring_alg == 0 ) then
             if( if_on_boundaries ) then
                call graphs_coloring_greedy(&
                     meshe % nboun,pelel,lelel,par_omp_nboun_list_colors,&
                     par_omp_nboun_num_colors,par_omp_nboun_ia_colors,&
                     par_omp_nboun_ja_colors,&
                     IA_NAME='PAR_OMP_NBOUN_IA_COLORS',JA_NAME='PAR_OMP_NBOUN_JA_COLORS',&
                     memor=memor_dom)
             else
                call graphs_coloring_greedy(&
                     meshe % nelem,pelel,lelel,par_omp_list_colors,&
                     par_omp_num_colors,par_omp_ia_colors,&
                     par_omp_ja_colors,&
                     IA_NAME='PAR_OMP_IA_COLORS',JA_NAME='PAR_OMP_JA_COLORS',&
                     memor=memor_dom)
             end if
          else
             if( if_on_boundaries ) then
                call graphs_coloring(&
                     meshe % nboun,pelel,lelel,par_omp_nboun_list_colors,&
                     par_omp_nboun_num_colors,par_omp_nboun_ia_colors,&
                     par_omp_nboun_ja_colors,&
                     IA_NAME='PAR_OMP_NBOUN_IA_COLORS',JA_NAME='PAR_OMP_NBOUN_JA_COLORS',&
                     memor=memor_dom)
             else
                call graphs_coloring(&
                     meshe % nelem,pelel,lelel,par_omp_list_colors,&
                     par_omp_num_colors,par_omp_ia_colors,&
                     par_omp_ja_colors,&
                     IA_NAME='PAR_OMP_IA_COLORS',JA_NAME='PAR_OMP_JA_COLORS',&
                     memor=memor_dom)
             end if
          end if
          !
          ! Deallocate graphs
          !
          if( if_on_boundaries ) then
             call graphs_coloring_deallocate(par_omp_nboun_list_colors,memor=memor_dom)
          else
             call graphs_coloring_deallocate(par_omp_list_colors,memor=memor_dom)
          end if

          call graphs_eleele_deallocate(PELEL=pelel,LELEL=lelel,memor=memor_dom)
          
          if( .not. present(pelpo_opt) ) call graphs_eleele_deallocate(PELPO=pelpo,memor=memor_dom)
          if( .not. present(lelpo_opt) ) call graphs_eleele_deallocate(LELPO=lelpo,memor=memor_dom)
          !
          ! Output
          !
          if( if_on_boundaries ) then
             call livinf(0_ip,'NUMBER OF COLORS FOUND (BOUNDARIES)= '//trim(intost(par_omp_nboun_num_colors)),0_ip)
          else
             call livinf(0_ip,'NUMBER OF COLORS FOUND (ELEMENTS)= '//trim(intost(par_omp_num_colors)),0_ip)
          end if

       else

          if( if_on_boundaries ) then
             par_omp_nboun_num_colors = 1
             call memory_alloca(par_memor,'PAR_OMP_NBOUN_IA_COLORS','par_omp_coloring',par_omp_nboun_ia_colors,par_omp_nboun_num_colors+1_ip)
             call memory_alloca(par_memor,'PAR_OMP_NBOUN_JA_COLORS','par_omp_coloring',par_omp_nboun_ja_colors,meshe % nboun)
             par_omp_nboun_ia_colors(1) = 1
             par_omp_nboun_ia_colors(2) = meshe % nboun + 1
             do ielem = 1,meshe % nboun
                par_omp_nboun_ja_colors(ielem) = ielem
             end do
          else
             par_omp_num_colors = 1
             call memory_alloca(par_memor,'PAR_OMP_IA_COLORS','par_omp_coloring',par_omp_ia_colors,par_omp_num_colors+1_ip)
             call memory_alloca(par_memor,'PAR_OMP_JA_COLORS','par_omp_coloring',par_omp_ja_colors,meshe % nelem)
             par_omp_ia_colors(1) = 1
             par_omp_ia_colors(2) = meshe % nelem + 1
             do ielem = 1,meshe % nelem
                par_omp_ja_colors(ielem) = ielem
             end do
          end if

       end if

    end if

  end subroutine parall_openmp_coloring

!!$    subroutine parall_openmp_test_coloring(nvert,ia,ja,lcolo,ncolo,ia_color,ja_color)
!!$    integer(ip), intent(in)                       :: nvert
!!$    integer(ip), intent(in),    pointer           :: ia(:)
!!$    integer(ip), intent(in),    pointer           :: ja(:)
!!$    integer(ip), intent(inout), pointer           :: lcolo(:)
!!$    integer(ip), intent(out),            optional :: ncolo
!!$    integer(ip), intent(inout), pointer, optional :: ia_color(:)
!!$    integer(ip), intent(inout), pointer, optional :: ja_color(:)
!!$    integer(ip)                                   :: ivert,jvert,kvert,mvert
!!$    integer(ip)                                   :: icolo,jcolo,izver,jzver
!!$    integer(ip)                                   :: istack,ihuge,nvert0
!!$    integer(ip)                                   :: nstack,ivert_last,color
!!$    integer(ip)                                   :: mz,hh
!!$    integer(ip), pointer                          :: lcolo_loc(:)
!!$    
!!$    if( .not. associated(lcolo) ) call memory_alloca(memor_dom,'LCOLO','graphs_coloring',lcolo,nvert)
!!$
!!$    hh    = huge(1_ip)
!!$    lcolo = hh
!!$    mz = 0
!!$    do ivert = 1,nvert
!!$       mz = max(mz,ia(ivert+1)-ia(ivert))
!!$    end do
!!$    allocate(lcolo_loc(mz))
!!$    !
!!$    ! Color linked list
!!$    !
!!$    if( present(ncolo) .and. present(ia_color) .and. present(ja_color) ) then
!!$       call memory_alloca(memor_dom,'IA_COLOR','graphs_coloring',ia_color,ncolo+1)
!!$       call memory_alloca(memor_dom,'JA_COLOR','graphs_coloring',ja_color,nvert)
!!$       do ivert = 1,nvert
!!$          lcolo_loc = hh
!!$          do izver = ia(ivert),ia(ivert+1)-1
!!$             icolo = minval(lcolo(ja(ia(ivert):ia(ivert+1)-1)
!!$             if( lcolo(ivert) == icolo ) then
!!$                kvert           = jvert + ia_color(icolo) 
!!$                jvert           = jvert + 1
!!$                ja_color(kvert) = ivert
!!$             end if
!!$          end do
!!$          ia_color(icolo+1) = ia_color(icolo) + jvert
!!$       end do
!!$    end if

  subroutine parall_openmp_adjacency_ompss_unity_test(ompss_domains,list_domains,num_pack,npoin_in,lnnod_in,lnods_in)

    use mod_parall, only : typ_list_elements_par
    use def_master, only : IMASTER
    use mod_memory, only : memory_size
    
    type(ompss_domain),          pointer, intent(in) :: ompss_domains(:) !< Inverse permutation new = invpe(old) and adjacency
    type(typ_list_elements_par), pointer, intent(in) :: list_domains(:)
    integer(ip),                 pointer, intent(in) :: num_pack(:)
    integer(ip),                          intent(in) :: npoin_in
    integer(ip),                 pointer, intent(in) :: lnnod_in(:)
    integer(ip),                 pointer, intent(in) :: lnods_in(:,:)
    integer(ip)                                      :: isubd,ii,ineig,ipoin
    integer(ip)                                      :: ielem,inode
    integer(ip)                                      :: ipack,jsubd,num_neigh
    integer(ip)                                      :: num_subd
    logical(lg),                 pointer             :: lmark(:)
    logical(lg),                 pointer             :: smark(:)

    nullify(lmark,smark)

    if( .not. associated(ompss_domains) ) then
       return
    else
       num_subd = size(ompss_domains,KIND=ip)
    end if

    if( IMASTER )      return
    if( num_subd < 1 ) return

    allocate(lmark(npoin_in))
    allocate(smark(num_subd))

    do isubd = 1,num_subd
       !
       ! Mark nodes of this subdomain
       !
       lmark = .false.
       do ipack = 1,num_pack(isubd)
          do ii = 1,memory_size(list_domains(isubd) % packs(ipack) % l)
             ielem = list_domains(isubd) % packs(ipack) % l(ii)
             if( ielem > 0 ) then
                do inode = 1,lnnod_in(ielem)
                   ipoin = lnods_in(inode,ielem)
                   lmark(ipoin) = .true.
                end do
             end if
          end do
       end do
       !
       ! Mark independent subdomains
       !
       smark = .true.
       num_neigh = size(ompss_domains(isubd) % neighbours,KIND=ip)
       do ii = 1,num_neigh
          ineig = ompss_domains(isubd) % neighbours(ii)
          smark(ineig) = .false.
       end do
       !
       ! Loop over independent subdomains
       !
       do jsubd = 1,num_subd
          if( smark(isubd) .and. isubd /= jsubd ) then

             do ipack = 1,num_pack(jsubd)
                do ii = 1,memory_size(list_domains(jsubd) % packs(ipack) % l)
                   ielem = list_domains(jsubd) % packs(ipack) % l(ii)
                   if( ielem > 0 ) then
                      do inode = 1,lnnod_in(ielem)
                         ipoin = lnods_in(inode,ielem)
                         if( lmark(ipoin) ) call runend('NODE ALREADY MARKED!')
                      end do
                   end if
                end do
             end do

          end if
       end do

    end do

    deallocate(lmark,smark)

  end subroutine parall_openmp_adjacency_ompss_unity_test

end module mod_parall_openmp
!> @}
