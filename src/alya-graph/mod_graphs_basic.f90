!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Graphs_Toolbox
!> Toolbox for graphs operations
!> @{
!> @name    ToolBox for graphs and renumbering
!> @file    mod_graphs.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for graphs and renumbering.
!> @details ToolBox for graphs and renumbering. Uses METIS_NodeND,
!>          (Node dissection) for renumbering
!
!-----------------------------------------------------------------------

module mod_graphs_basic

  use def_kintyp_basic,       only : ip,rp,lg,i1p,i2p
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_memory,             only : memory_resize
  use mod_memory,             only : memory_size
  use mod_maths,              only : maths_heap_sort
  use mod_alya2metis,         only : alya2metis_METIS_NodeND
  use mod_alya2metis,         only : alya2metis_initialization
  use mod_strings,            only : integer_to_string
  use mod_optional_argument,  only : optional_argument
  use mod_maths,              only : maths_heap_sort

  implicit none
  private 
  
  type TeliTree
     integer(ip)          :: nn                !< Number of nodes
     integer(ip), pointer :: descendiente(:)   !< Descendiente del nodo i  
     integer(ip), pointer :: nant(:)           !< N. antecesores del nodo i  
     integer(ip), pointer :: Tnant(:)          !< N total de antecesores de i 
     type(i1p),   pointer :: antecesores(:)    !< Antecesores del nodo i         
  end type TeliTree

  public :: graphs_permut_metis_postordering
  public :: graphs_permut_metis_postordering_deallocate
  public :: graphs_postorder
  public :: graphs_rengra
  public :: graphs_comper
  public :: graphs_iniper
  public :: graphs_permut
  public :: graphs_subgra
  public :: graphs_number_to_linked_list
  
contains

  !------------------------------------------------------------------------
  !
  !> @brief   Remove diagonal from a graph
  !> @details Remove the diagonal from a graph.
  !>          It can be useful for example to feed METIS which does not
  !>          require the diagonal to renumber or to partition a graph 
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_offdia(itask,nn,nnz,ia,ja,iaren,jaren,memor)

    integer(ip),          intent(in)            :: itask
    integer(ip),          intent(in)            :: nn
    integer(ip),          intent(in)            :: nnz
    integer(ip), pointer, intent(inout)         :: ia(:)
    integer(ip), pointer, intent(inout)         :: ja(:)
    integer(ip), pointer, intent(inout), optional :: iaren(:)
    integer(ip), pointer, intent(inout), optional :: jaren(:)
    integer(8),           intent(inout)           :: memor(2)
    integer(ip), pointer                        :: iare2(:)
    integer(ip), pointer                        :: jare2(:)
    integer(ip)                                 :: ii,iiz,jjz

    select case ( itask )

    case (  1_ip ) 
       !
       ! Allocate and take off diagonal
       !
       if( nnz <= nn ) return

       if( present(iaren) ) then
          call memory_alloca(memor,'IAREN','graphs_offdia',iaren,nn+1)
          call memory_alloca(memor,'JAREN','graphs_offdia',jaren,max(1_ip,nnz-nn))
          iare2 => iaren
          jare2 => jaren
       else
          call memory_alloca(memor,'IAREN','graphs_offdia',iare2,nn+1)
          call memory_alloca(memor,'JAREN','graphs_offdia',jare2,max(1_ip,nnz-nn))          
       end if

       iare2(1) = 1
       do ii = 2,nn+1
          iare2(ii) = ia(ii)-ii+1
       end do

       jjz = 0
       do ii = 1,nn
          do iiz = ia(ii),ia(ii+1)-1
             if( ja(iiz) /= ii ) then
                jjz = jjz + 1
                jare2(jjz) = ja(iiz)
             end if
          end do
       end do

       if( .not. present(iaren) ) then
          call memory_deallo(memor,'IA','graphs_offdia',ia)
          call memory_deallo(memor,'JA','graphs_offdia',ja)  
          call memory_alloca(memor,'IA','graphs_offdia',ia,nn+1)
          call memory_alloca(memor,'JA','graphs_offdia',ja,max(1_ip,nnz-nn))
          do ii = 1,nn+1
             ia(ii) = iare2(ii)
          end do
          do iiz = 1,max(1_ip,nnz-nn)
             ja(iiz) = jare2(iiz)
          end do
          call memory_deallo(memor,'IARE2','graphs_offdia',iare2)
          call memory_deallo(memor,'JARE2','graphs_offdia',jare2)          
       end if

    case ( -1_ip )  
       !
       ! Deallocate graph
       !
       if( present(iaren) ) then
          call memory_deallo(memor,'IAREN','graphs_offdia',iaren)
          call memory_deallo(memor,'JAREN','graphs_offdia',jaren)
       end if

    end select

  end subroutine graphs_offdia

  !-----------------------------------------------------------------------
  !
  !> @brief   Permutes a graph
  !> @details Permutes a graph containing its diagonal \n
  !!          1. Take off diag \n
  !!          2. Compute permutation arrays permr and invpr \n
  !!          3. Renumber the graph using nested dissection (METIS) \n
  !!          4. Order the graph \n
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_permut(nn,nnz,ia,ja,permr,invpr,PERMR_NAME,INVPR_NAME,memor)

    integer(ip),               intent(in)             :: nn         !< Number of rows
    integer(ip),               intent(in)             :: nnz        !< Size of graph
    integer(ip),      pointer, intent(inout)          :: permr(:)   !< Permutation:  NEW = INVPR(OLD)
    integer(ip),      pointer, intent(inout)          :: invpr(:)   !< Inverse perm: OLD = PERMR(NEW)
    integer(ip),      pointer, intent(inout)          :: ia(:)      !< Graph pointer
    integer(ip),      pointer, intent(inout)          :: ja(:)      !< Graph list
    character(len=*),          intent(in),   optional :: PERMR_NAME
    character(len=*),          intent(in),   optional :: INVPR_NAME
    integer(8),                intent(inout)          :: memor(2)
    integer(ip),      pointer                         :: iadiag(:)  !< Graph pointer
    integer(ip),      pointer                         :: jadiag(:)  !< Graph list
    integer(ip)                                       :: ii

    nullify(iadiag)
    nullify(jadiag)
    !
    ! Permutation arrays
    !
    if( .not. associated(permr) ) then
       if( present(PERMR_NAME) ) then
          call memory_alloca(memor,trim(PERMR_NAME),'graphs_permut',permr,nn)
       else
          call memory_alloca(memor,'PERMR','graphs_permut',permr,nn)
       end if
    end if
    if( .not. associated(invpr) ) then
       if( present(INVPR_NAME) ) then
          call memory_alloca(memor,trim(INVPR_NAME),'graphs_permut',invpr,nn)
       else
          call memory_alloca(memor,'INVPR','graphs_permut',invpr,nn)
       end if
    end if
    !
    ! Initialize permutation: useful if Parall is off
    !    
    do ii = 1,nn
       permr(ii) = ii
       invpr(ii) = ii
    end do
    !
    ! If matrix is diagonal do not do anuthing
    !
    if( nnz == nn ) return
    !
    ! IAREN, JAREN: Compute graph without diagonal from IA, JA
    !
    call graphs_offdia(1_ip,nn,nnz,ia,ja,iadiag,jadiag,memor=memor)
    !
    ! Renumber graph using METIS (if available via Parall service)
    ! using Nested Bisection
    !
    if( nnz == nn ) return
    !
    ! In METIS call, permutation and inverse permutation are inverted!!!!!!
    ! This is why the call to METIS_NodeNd is upside down
    !
    call alya2metis_initialization()
    call alya2metis_METIS_NodeND(nn,iadiag,jadiag,permr,invpr)
    !
    ! Deallocate memory for non-diagonal graph
    !
    call graphs_offdia(-1_ip,nn,nnz,ia,ja,iadiag,jadiag,memor=memor)
    !
    ! Renumber and sort graph
    !
    call graphs_rengra(nn,nnz,permr,invpr,ia,ja,memor=memor)

  end subroutine graphs_permut

  !------------------------------------------------------------------------
  !
  !> @brief   Renumber a graph
  !> @details This subroutine renumber a graph (IA, JA)  using permutation 
  !!          arrays (PERMR, INVPR) and order the entries of JA
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_rengra(nn,nnz,permr,invpr,ia,ja,memor)

    integer(ip),          intent(in)    :: nn
    integer(ip),          intent(in)    :: nnz
    integer(ip), pointer, intent(in)    :: permr(:)
    integer(ip), pointer, intent(in)    :: invpr(:)
    integer(ip), pointer, intent(inout) :: ia(:)
    integer(ip), pointer, intent(inout) :: ja(:)
    integer(8),           intent(inout) :: memor(2)
    integer(ip), pointer                :: iacpy_rengra(:)
    integer(ip), pointer                :: jacpy_rengra(:)
    integer(ip)                         :: ii,jj,kk,iiz,jjz,ll,nnz2
    !
    ! Initialization
    !
    nullify(iacpy_rengra)
    nullify(jacpy_rengra)
    if( nnz == 0 ) then
       nnz2 = size(ja,KIND=ip)
    else
       nnz2 = nnz
    end if
    !
    ! Copy graph in temporary graphs IACPY_RENGRA, JACPY_RENGRA
    !
    call memory_alloca(memor,'IACPY_RENGRA','graphs_rengra',iacpy_rengra,nn+1_ip)
    call memory_alloca(memor,'JACPY_RENGRA','graphs_rengra',jacpy_rengra,nnz2)
    !
    ! Copy graph to IACPY_RENGRA, JACPY_RENGRA 
    !
    do ii = 1,nn+1
       iacpy_rengra(ii) = ia(ii)
    end do
    do iiz = 1,nnz2
       jacpy_rengra(iiz) = ja(iiz)
    end do
    !
    ! IA(II) is new number of connections of II
    !
    do ii = 1,nn
       jj     = invpr(ii)  ! old = invpr(new)
       ia(ii) = iacpy_rengra(jj+1) - iacpy_rengra(jj)
    end do
    !
    ! Renumber
    !
    kk    = ia(1)
    ia(1) = 1 
    do ii = 2,nn+1
       ll     = ia(ii)
       ia(ii) = ia(ii-1) + kk
       kk     = ll
    end do

    do ii = 1,nn
       jj = invpr(ii)
       kk = ia(ii)
       do jjz = iacpy_rengra(jj),iacpy_rengra(jj+1)-1
          ja(kk) = permr(jacpy_rengra(jjz))
          kk     = kk + 1
       end do
    end do
    !
    ! Order graph
    !
    do ii = 1,nn
       ll = ia(ii+1)-ia(ii)
       if( ll > 0 ) call heapsorti1(2_ip,ll,ja(ia(ii)))
    end do
    !
    ! Deallocate temporary graphs IACPY_RENGRA, JACPY_RENGRA
    !
    call memory_deallo(memor,'JACPY_RENGRA','graphs_rengra',jacpy_rengra)
    call memory_deallo(memor,'IACPY_RENGRA','graphs_rengra',iacpy_rengra)

  end subroutine graphs_rengra

  subroutine graphs_iniper(nn,permr,invpr,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_iniper
    ! NAME
    !    graphs_iniper
    ! DESCRIPTION
    !    This subroutine allocates and initializes permutation arrays
    ! INPUT
    !    NN ................ Number of rows
    ! OUTPUT
    !    PERMR(NN) ......... Permutation:  Identity
    !    INVPR(NN) ......... Inverse perm: Identity
    ! USED BY
    !***
    !------------------------------------------------------------------------

    integer(ip),          intent(in)    :: nn
    integer(ip), pointer, intent(inout)   :: permr(:)
    integer(ip), pointer, intent(inout)   :: invpr(:)
    integer(8),           intent(inout)   :: memor(2)
    integer(ip)                         :: ii
    !
    ! Permutation arrays
    !
    call memory_alloca(memor,'PERMR','graphs_iniper',permr,nn)
    call memory_alloca(memor,'INVPR','graphs_iniper',invpr,nn)
    !
    ! Initialize permutation: useful if Parall is off
    !    
    do ii = 1,nn
       permr(ii) = ii
       invpr(ii) = ii
    end do

  end subroutine graphs_iniper

  subroutine graphs_postorder(nn,ia,ja,permr,invpr,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_iniper
    ! NAME
    !    graphs_iniper
    ! DESCRIPTION
    !    This subroutine allocates and initializes permutation arrays
    ! INPUT
    !    NN ................ Number of rows
    ! OUTPUT
    !    PERMR(NN) ......... Permutation:  Identity
    !    INVPR(NN) ......... Inverse perm: Identity
    ! USED BY
    !***
    !------------------------------------------------------------------------
    
    integer(ip),          intent(in)    :: nn
    integer(ip), pointer, intent(in)    :: ia(:)
    integer(ip), pointer, intent(in)    :: ja(:)
    integer(ip), pointer, intent(inout) :: permr(:)
    integer(ip), pointer, intent(inout) :: invpr(:)
    integer(8),           intent(inout) :: memor(2)
    type(TeliTree)                      :: tree
    integer(ip)                         :: node,last
    integer(ip), pointer                :: iL(:),jL(:)
    integer(ip), pointer                :: iU(:),jU(:)
    !
    ! Symbolic factorization
    !
    nullify(iL,jL,iU,jU)
    call graphs_symbolical_LU(nn,iA,jA,iL,jL,iU,jU,MEMORY_COUNTER=memor)
    !
    ! Build Elimination tree
    !
    call graphs_BuildEliminationTree(nn,iL,jL,tree,memor=memor)

    node = nn    ! Nodo en la cima del arbol
    last = 1     ! Inicio de la nueva numeracion
    !
    ! Mientras queden sub-arboles independientes 
    !
    do while( node >= 1 )

       call graphs_generatePostOrder1(tree,node,invpr,last,memor=memor)
       !
       ! Look for top of next sub-tree 
       !
       node = node - 1

       loop_node: do while ( node >= 1 )
          if( tree % descendiente(node) >= 1 ) then
             node = node - 1
          else
             exit loop_node
          end if
       end do loop_node

    end do

    if ( last /= nn + 1 ) call runend('WRONG TREE')

    do node = 1,nn
       permr(invpr(node)) = node
    end do
    !
    ! Deallocate memory
    !
    call graphs_symbolical_LU_deallocate(iL,jL,iU,jU,MEMORY_COUNTER=memor)
    call graphs_deallocatetree(tree,memor=memor)

  end subroutine graphs_postorder

  subroutine graphs_BuildEliminationTree(nn,iL,jL,tree,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_iniper
    ! NAME
    !    graphs_iniper
    ! DESCRIPTION
    !    Computes the elimination tree of the factorisation.
    ! INPUT
    !   NN ...... Number of rows
    !   IL,JL ... L-matrix graph
    ! OUTPUT
    !   TREE .... Elimination tree
    ! USED BY
    !***
    !------------------------------------------------------------------------
    
    integer(ip),          intent(in)    :: nn
    integer(ip), pointer, intent(in)    :: iL(:)
    integer(ip), pointer, intent(in)    :: jL(:)
    type(TeliTree),       intent(inout) :: tree
    integer(8),           intent(inout) :: memor(2)
    integer(ip)                         :: ii,jj,kk
    integer(ip), pointer                :: descendiente(:)
    integer(ip), pointer                :: nant(:)
    integer(ip), pointer                :: Tnant(:)
    type(i1p),   pointer                :: antecesores(:)

    tree % nn = nn
    nullify(tree % descendiente)
    nullify(tree % nant)
    nullify(tree % Tnant)
    nullify(tree % antecesores)

    call memory_alloca(memor,'TREE % DESCENDIENTE','graphs_BuildEliminationTree',tree % descendiente,nn)
    call memory_alloca(memor,'TREE % NANT'        ,'graphs_BuildEliminationTree',tree % nant        ,nn)
    call memory_alloca(memor,'TREE % TNANT'       ,'graphs_BuildEliminationTree',tree % Tnant       ,nn)
    call memory_alloca(memor,'TREE % ANTECESORES' ,'graphs_BuildEliminationTree',tree % antecesores ,nn)
    
    descendiente => tree % descendiente
    nant         => tree % nant
    Tnant        => tree % Tnant
    antecesores  => tree % antecesores

    do ii = 1,nn
       do kk = iL(ii),iL(ii+1)-1
          jj = jL(kk)
          if( descendiente(jj) <= 0 ) descendiente(jj) = ii
       end do
    end do
    !
    ! Count the amount of predecessors
    ! and the accumulated load for each node
    !
    do ii = 1,nn
       nant(ii)  = 0
       Tnant(ii) = 1
    end do

    do ii = 1,nn
       kk = descendiente(ii)
       if( kk > 0 ) then
          nant(kk)  = nant(kk)  + 1
          Tnant(kk) = Tnant(kk) + Tnant(ii)
       end if
    end do
    !
    ! Allocate space for the predecessor list
    !
    do ii = 1,nn
       kk = nant(ii)
       if( kk > 0 ) then
          call memory_alloca(memor,'TREE % ANTECESORES % L','graphs_BuildEliminationTree', antecesores(ii) % l,kk)
          !allocate(antecesores(ii) % l(kk))
       else
          nullify( antecesores(ii) % l )
       end if
    end do
    !
    ! Fill the predecessor list for each node
    !
    do ii = 1,nn
       nant(ii) = 0
    end do

    do ii = 1,nn
       kk = descendiente(ii)
       if( kk > 0 ) then
          nant(kk) = nant(kk) + 1
          antecesores(kk) % l(nant(kk)) = ii
       end if
    end do

  end subroutine graphs_BuildEliminationTree

  subroutine graphs_generatePostOrder1(tree,node,invpr,last,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_generatePostOrder1
    ! NAME
    !    graphs_generatePostOrder1
    ! DESCRIPTION
    !  Generates a postorder with maximum locality. It returns the number
    !   of ordered nodes starting form "node".
    ! INPUT
    !   tree        Elimination tree
    !   node        Starting node
    !   invp        Inverse permutation array
    !   last        Number of ordered nodes
    ! OUTPUT
    ! USED BY
    !***
    !------------------------------------------------------------------------

    type(TeliTree),       intent(in)    :: tree
    integer(ip),          intent(in)    :: node
    integer(ip), pointer, intent(inout) :: invpr(:)
    integer(ip),          intent(inout) :: last
    integer(8)                          :: memor(2)
    integer(ip)                         :: ii,kk,active_node,n_act
    integer(ip)                         :: nn
    integer(ip), pointer                :: nant(:)
    type(i1p),   pointer                :: antecesores(:)
    integer(ip), pointer                :: iwa1(:)
    logical(lg), pointer                :: mask(:)

    nant         => tree % nant
    antecesores  => tree % antecesores
    nn           =  tree % nn
    nullify(iwa1)
    nullify(mask)

    call memory_alloca(memor,'iwa1','graphs_generatePostOrder1',iwa1,nn)
    call memory_alloca(memor,'mask','graphs_generatePostOrder1',mask,nn)       

    do ii = 1,nn
       mask(ii) = .false.
    end do

    active_node = 1
    iwa1(active_node) = node

    do while( active_node > 0 )

       n_act = iwa1(active_node)

       if ( .not. mask(n_act) ) then

          mask(n_act) = .true.
          do ii = nant(n_act),1,-1
             active_node = active_node + 1
             kk = antecesores(n_act) % l(ii)
             iwa1(active_node) = kk
          end do

       else

          active_node = active_node - 1
          invpr(last) = n_act
          last = last + 1

       end if

    end do

    call memory_deallo(memor,'MASK','graphs_generatePostOrder1',mask)
    call memory_deallo(memor,'IWA1','graphs_generatePostOrder1',iwa1)

  end subroutine graphs_generatePostOrder1

  subroutine graphs_deallocatetree(tree,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_deallocatetree
    ! NAME
    !    graphs_deallocatetree
    ! DESCRIPTION
    !    Deallocates the tree
    ! INPUT
    !    TREE ................ tree
    ! OUTPUT
    ! USED BY
    !***
    !------------------------------------------------------------------------
   
    type(TeliTree)  :: tree
    integer(8)      :: memor(2)

    call memory_deallo(memor,'TREE % DESCENDIENTE','graphs_deallocatetree',tree % descendiente)
    call memory_deallo(memor,'TREE % NANT'        ,'graphs_deallocatetree',tree % nant)
    call memory_deallo(memor,'TREE % TNANT'       ,'graphs_deallocatetree',tree % Tnant)
    call memory_deallo(memor,'TREE % ANTECESORES' ,'graphs_deallocatetree',tree % antecesores)
    
  end subroutine graphs_deallocatetree

  subroutine graphs_comper(nn,permr_1,invpr_1,permr_2,invpr_2,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_comper
    ! NAME
    !    graphs_comper
    ! DESCRIPTION
    !    This subroutine composes permutations
    !                       II
    !                       ||
    !                       \/
    !    II=INVPR_1(KK)     KK     KK=PERMR_1(II)
    !                       ||
    !                       \/
    !    KK=INVPR_2(JJ)     JJ     JJ=PERMR_2(KK)
    !
    !    JJ = PERMR_2(PERMR_1(II))
    !    II = INVPR_1(INVPR_2(jj))
    !
    ! INPUT
    !    PERMR_1, INVPR_1 ............ First permutation
    !    PERMR_2, INVPR_2 ............ Second permutation
    ! OUTPUT
    !    PERMR_1, INVPR_1 ............ Final permutation
    ! USED BY
    !***
    !------------------------------------------------------------------------

    integer(ip),          intent(in)    :: nn
    integer(ip), pointer, intent(inout) :: permr_1(:)
    integer(ip), pointer, intent(inout) :: invpr_1(:)
    integer(ip), pointer, intent(in)    :: permr_2(:)
    integer(ip), pointer, intent(in)    :: invpr_2(:)
    integer(ip), pointer                :: permr_cpy(:)
    integer(ip), pointer                :: invpr_cpy(:)
    integer(8),           intent(inout) :: memor(2)
    integer(ip)                         :: ii,jj,kk
    !
    ! Copy permutation
    !
    nullify(permr_cpy,invpr_cpy)
    call graphs_iniper(nn,permr_cpy,invpr_cpy,memor=memor)
    do ii = 1,nn
       permr_cpy(ii) = permr_1(ii)
       invpr_cpy(ii) = invpr_1(ii)
    end do
    !
    ! II = OLD, JJ = NEW, KK = INTERMEDIATE
    !
    do ii = 1,nn
       kk = permr_cpy(ii)
       jj = permr_2(kk)
       permr_1(ii) = jj 
    end do
    do jj = 1,nn
       kk = invpr_2(jj)
       ii = invpr_cpy(kk)
       invpr_1(jj) = ii 
    end do
    !
    ! Deallocate
    !
    call graphs_deaper(permr_cpy,invpr_cpy,memor=memor)

  end subroutine graphs_comper

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-17
  !> @brief   Symbolical factorization
  !> @details Symbolical factorization using a CSR format graph
  !> 
  !-----------------------------------------------------------------------
  
  subroutine graphs_symbolical_LU(nbrows,iA,jA,iL,jL,iU,jU,rfillin,MEMORY_COUNTER)

    integer(ip),                       intent(in)    :: nbrows
    integer(ip),           pointer,    intent(in)    :: iA(:)
    integer(ip),           pointer,    intent(in)    :: jA(:)
    integer(ip),           pointer,    intent(inout) :: iL(:),jL(:)
    integer(ip),           pointer,    intent(inout) :: iU(:),jU(:)
    real(rp),    optional,             intent(out)   :: rfillin
    integer(8),  optional,             intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                      :: i,j,k,col,nxt,totalL,totalU
    integer(ip)                                      :: rowIA_ini,rowIA_diag,rowIA_fin  
    integer(ip)                                      :: minColL,maxColL,minColU,maxColU
    integer(ip), pointer                             :: firstL(:),firstU(:),seen(:),nextU(:),iwa(:)   
    integer(ip), pointer                             :: nzL(:),nzU(:)       !Local arrays for L/U matrices 
    type(i1p),   pointer                             :: ptL(:),ptU(:)
    integer(8)                                       :: memor(2)

    memor = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    if( associated(iL) .or. associated(jL) .or. associated(iU) .or. associated(jU) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
       !
       ! Alloc local work arrays
       !
       nullify(firstL)
       nullify(firstU)
       nullify(seen)
       nullify(nextU)
       nullify(iwa)
       nullify(nzL)
       nullify(nzU)
       nullify(ptL)
       nullify(ptU)

       call memory_alloca(memor,'FIRSTL','graphs_symbolical_LU',firstL,nbrows)
       call memory_alloca(memor,'FIRSTU','graphs_symbolical_LU',firstU,nbrows)
       call memory_alloca(memor,'SEEN'  ,'graphs_symbolical_LU',seen,nbrows)
       call memory_alloca(memor,'NEXTU' ,'graphs_symbolical_LU',nextU,nbrows)
       call memory_alloca(memor,'IWA'   ,'graphs_symbolical_LU',iwa,nbrows)
       call memory_alloca(memor,'NZL'   ,'graphs_symbolical_LU',nzL,nbrows)
       call memory_alloca(memor,'NZU'   ,'graphs_symbolical_LU',nzU,nbrows)
       call memory_alloca(memor,'PTL'   ,'graphs_symbolical_LU',ptL,nbrows)
       call memory_alloca(memor,'PTU'   ,'graphs_symbolical_LU',ptU,nbrows)
       !
       ! Initialize the fill-in links and local L/U arrays
       !
       do i= 1, nbrows
          firstL(i) = -1
          firstU(i) = -1
          seen(i)   = -1
          nzL(i)    =  0
          nzU(i)    =  0
          nullify( ptL(i) % l , ptU(i) % l )   
       end do
       !
       ! Initialize total number of non-zero elements in L/U
       !
       totalL = 0
       totalU = 0
       !
       ! Main loop in rows
       !
       do i= 1, nbrows
          minColL = nbrows
          maxColL = -1
          minColU = nbrows
          maxColU = -1
          !
          !For all the elements in the i-th row of the matrix A
          !
          rowIA_ini  = iA(i)
          rowIA_fin  = iA(i+1)-1
          rowIA_diag = iA(i+1)
          !
          ! Find a diagonal element in the i-th row of the matrix A
          !
          loop1: do k = rowIA_ini, rowIA_fin
             if( jA(k) == i ) then
                rowIA_diag = k
                exit loop1
             end if
          end do loop1
          !
          ! For all the elements in the i-th row of the matrix A before diagonal
          !
          do j = rowIA_ini,rowIA_diag-1
             col = jA(j)
             if( seen(col) /= i ) then
                seen(col) = i
                minColL   = min(minColL,col)
                maxColL   = max(maxColL,col)
                !
                ! Compute Reachable Set from L
                !
                col = firstL(col)
                loop2:do while (col /= -1)
                   if( seen(col) /= i ) then
                      minColL   = min(minColL,col)
                      maxColL   = max(maxColL,col)
                      seen(col) = i
                      col = firstL(col)
                   else
                      exit loop2
                   end if
                end do loop2
             end if
          end do
          !
          ! For all the elements in the i-th row of the matrix A after diagonal
          !
          do j = rowIA_diag+1, rowIA_fin
             col       = jA(j)
             seen(col) = i
             minColU   = min(minColU,col)
             maxColU   = max(maxColU,col)
          end do
          !
          ! Compute Reachable Set from U
          ! 
          k = firstU(i)
          do while (k /= -1)
             !
             ! For all the elements in the k-th row of the matrix U without the diagonal element
             !
             do j= 2, nzU(k)
                col = ptU(k) % l(j)
                if( col > i ) then
                   minColU = min(minColU,col)
                   maxColU = max(maxColU,col)
                   seen(col) = i
                end if
             end do
             k = nextU(k)
          end do
          !
          ! For all the non-zero elements of the matrix L. L matrix is stored without the diagonal element
          !
          nxt = 0
          do j= minColL, maxColL
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
                if( firstL(j) == -1) firstL(j) = i
             end if
          end do
          !
          ! Allocate space for L values and copy iwa to ptL 
          !
          totalL = totalL + nxt
          nzL(i) = nxt
          call memory_alloca(memor,'PTL % L','graphs_symbolical_LU',ptL(i) % l,nxt,'DO_NOT_INITIALIZE')
          do k= 1, nxt
             ptL(i) % l(k) = iwa(k)
          end do
          !
          ! For all the non-zero elements of the matrix U. The diagonal element must be included
          !
          iwa(1) = i       !The diagonal element 
          nxt    = 1
          do j= minColU, maxColU
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
             end if
          end do
          !
          ! Put the proper link for future fill-in generation 
          !
          if( nxt > 2 ) then
             col         = iwa(2)
             nextU(i)    = firstU(col)
             firstU(col) = i
          end if
          !
          ! Allocate space for U values and copy iwa to ptU 
          !
          totalU = totalU + nxt
          nzU(i) = nxt 
          call memory_alloca(memor,'PTU % L','graphs_symbolical_LU',ptU(i) % l,nxt,'DO_NOT_INITIALIZE')
          do k= 1, nxt
             ptU(i) % l(k) = iwa(k)
          end do
       end do
       !
       ! Put L and U in CSR format 
       !
       call memory_alloca(memor,'IL','graphs_symbolical_LU',iL,nbrows+1_ip)
       call memory_alloca(memor,'IU','graphs_symbolical_LU',iU,nbrows+1_ip)
       call memory_alloca(memor,'JL','graphs_symbolical_LU',jL,totalL)
       call memory_alloca(memor,'JU','graphs_symbolical_LU',jU,totalU)
       !
       ! CSR format for the L matrix
       !
       iL(1) = 1
       do i= 1, nbrows
          iL(i+1) = iL(i) + nzL(i)
          nxt = iL(i) - 1
          do k= 1, nzL(i)
             jL(nxt+k) = ptL(i) % l(k)
          end do
       end do
       !
       ! CSR format for the U matrix
       !
       iU(1) = 1
       do i= 1, nbrows
          iU(i+1) = iU(i) + nzU(i)

          nxt = iU(i) - 1
          do k= 1, nzU(i)
             jU(nxt+k) =  ptU(i) % l(k)
          end do
       end do

       call memory_deallo(memor,'PTU'   ,'graphs_symbolical_LU',ptU)
       call memory_deallo(memor,'PTL'   ,'graphs_symbolical_LU',ptL)
       call memory_deallo(memor,'NZU'   ,'graphs_symbolical_LU',nzU)
       call memory_deallo(memor,'NZL'   ,'graphs_symbolical_LU',nzL)
       call memory_deallo(memor,'IWA'   ,'graphs_symbolical_LU',iwa)
       call memory_deallo(memor,'NEXTU' ,'graphs_symbolical_LU',nextU)
       call memory_deallo(memor,'SEEN'  ,'graphs_symbolical_LU',seen)
       call memory_deallo(memor,'FIRSTU','graphs_symbolical_LU',firstU)
       call memory_deallo(memor,'FIRSTL','graphs_symbolical_LU',firstL)

       if( present(rfillin) ) then
          rfillin = (real(iL(nbrows+1)-1,rp)+real(iU(nbrows+1)-1,rp))/real(iA(nbrows+1)-1,rp)
       end if

       if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

    end if

  end subroutine graphs_symbolical_LU

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-17
  !> @brief   Deallocate symbolical factorization
  !> @details Deallocate symbolical factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine graphs_symbolical_LU_deallocate(iL,jL,iU,jU,MEMORY_COUNTER)
 
    integer(ip),           pointer, intent(inout) :: iL(:),jL(:)
    integer(ip),           pointer, intent(inout) :: iU(:),jU(:)
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                    :: memor(2)

    memor = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    call memory_deallo(memor,'iL','graphs_symbolical_LU',iL)
    call memory_deallo(memor,'jL','graphs_symbolical_LU',jL)
    call memory_deallo(memor,'iU','graphs_symbolical_LU',iU)
    call memory_deallo(memor,'jU','graphs_symbolical_LU',jU)

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

  end subroutine graphs_symbolical_LU_deallocate

  subroutine graphs_deaper(permr,invpr,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_deaper
    ! NAME
    !    graphs_deaper
    ! DESCRIPTION
    !    Deallocate permutation arrays
    ! USED BY
    !    par_partit
    !***
    !------------------------------------------------------------------------

    integer(ip), pointer, intent(inout) :: permr(:)
    integer(ip), pointer, intent(inout) :: invpr(:)
    integer(8),           intent(inout)          :: memor(2)

    call memory_deallo(memor,'INVPR','graphs_permut',invpr)   
    call memory_deallo(memor,'PERMR','graphs_permut',permr)  

  end subroutine graphs_deaper

    !-----------------------------------------------------------------------
  !
  !> @brief   Compute a graph
  !> @details Compute the list of faces
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------
  
  subroutine graphs_permut_metis_postordering_deallocate(permr,invpr,memor_opt,PERMR_NAME,INVPR_NAME)
    
    integer(ip), pointer, intent(inout), optional :: permr(:)   !< Permutation:  OLD = PERMR(NEW)
    integer(ip), pointer, intent(inout), optional :: invpr(:)   !< Inverse perm: NEW = INVPR(OLD)    
    integer(8),           intent(inout)           :: memor_opt(2)
    character(len=*),     intent(in),    optional :: PERMR_NAME
    character(len=*),     intent(in),    optional :: INVPR_NAME
    character(len=:),     allocatable             :: my_permr_name
    character(len=:),     allocatable             :: my_invpr_name

    my_permr_name = optional_argument('PERMR_POST',PERMR_NAME)
    my_invpr_name = optional_argument('INVPR_POST',INVPR_NAME)
    
    if( present(permr) ) then
       call memory_deallo(memor_opt,my_permr_name,'graphs_permut_metis_postordering',permr)       
    end if
    if( present(invpr) ) then
       call memory_deallo(memor_opt,my_invpr_name,'graphs_permut_metis_postordering',invpr)
    end if

    if( allocated(my_permr_name) ) deallocate(my_permr_name)
    if( allocated(my_invpr_name) ) deallocate(my_invpr_name)
    
  end subroutine graphs_permut_metis_postordering_deallocate

  subroutine graphs_permut_metis_postordering(nn,nnz,ia,ja,permr,invpr,PERMR_NAME,INVPR_NAME,memor)
    integer(ip),               intent(in)             :: nn         !< Number of rows
    integer(ip),               intent(in)             :: nnz        !< Size of graph
    integer(ip),      pointer, intent(inout)          :: permr(:)   !< Permutation:  OLD = PERMR(NEW)
    integer(ip),      pointer, intent(inout)          :: invpr(:)   !< Inverse perm: NEW = INVPR(OLD)
    integer(ip),      pointer, intent(inout)          :: ia(:)      !< Graph pointer
    integer(ip),      pointer, intent(inout)          :: ja(:)      !< Graph list
    character(len=*),          intent(in),   optional :: PERMR_NAME
    character(len=*),          intent(in),   optional :: INVPR_NAME
    integer(8),                intent(inout)          :: memor(2)
    integer(ip),      pointer                         :: permr_2(:)        
    integer(ip),      pointer                         :: invpr_2(:)        

    if( nnz /= size(ja,KIND=ip) ) call runend('GRAPHS_PERMUT_METIS_POSTORDERING: WRONG GRAPH SIZE')
    !
    ! Permute using METIS
    !
    call graphs_permut(nn,nnz,ia,ja,permr,invpr,PERMR_NAME,INVPR_NAME,memor=memor)
    !
    ! Postordering 
    !
    ! - Copy second permutation using postordering: PERMR_2, INVPR_2
    ! - Renumber the graph IA, JA
    !
    nullify(permr_2,invpr_2)
    call graphs_iniper(nn,permr_2,invpr_2,memor=memor)
    call graphs_postorder(nn,ia,ja,permr_2,invpr_2,memor=memor)
    call graphs_rengra(nn,nnz,permr_2,invpr_2,ia,ja,memor=memor)
    !
    ! Compose permutation PERMR, INVPR and PERMR_2, INVPR_2
    !
    call graphs_comper(nn,permr,invpr,permr_2,invpr_2,memor=memor)
    call graphs_deaper(permr_2,invpr_2,memor=memor)

  end subroutine graphs_permut_metis_postordering

  !------------------------------------------------------------------------
  !
  !> @brief   Compute a subgraph from a graph
  !> @details This subroutine renumber a graph (IA, JA)  using permutation 
  !>          arrays (PERMR, INVPR) and order the entries of JA
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_subgra(nn,nnz,permr,invpr,ia,ja,memor)

    implicit none  
    integer(ip),          intent(in)    :: nn
    integer(ip),          intent(inout) :: nnz
    integer(ip), pointer, intent(in)    :: permr(:)
    integer(ip), pointer, intent(in)    :: invpr(:)
    integer(ip), pointer, intent(inout) :: ia(:)
    integer(ip), pointer, intent(inout) :: ja(:)
    integer(8),           intent(inout) :: memor(2)
    integer(ip), pointer                :: iacpy_subgra(:) 
    integer(ip), pointer                :: jacpy_subgra(:) 
    integer(ip)                         :: ii,jj,kk,iiz,jjz,ll,pp,mm,nnz1
    !
    ! Initialization
    !
    nullify(iacpy_subgra)
    nullify(jacpy_subgra)
    !
    ! Copy graph in temporary graphs IACPY_SUBGRA, JACPY_SUBGRA
    !
    call memory_alloca(memor,'IACPY_SUBGRA','graphs_subgra',iacpy_subgra,nn+1_ip)
    call memory_alloca(memor,'JACPY_SUBGRA','graphs_subgra',jacpy_subgra,nnz)
    !
    ! Copy graph to IACPY_SUBGRA, JACPY_SUBGRA 
    !
    do ii = 1,nn+1
       iacpy_subgra(ii) = ia(ii)
       ia(ii)           = 0
    end do
    do iiz = 1,nnz
       jacpy_subgra(iiz) = ja(iiz)
    end do
    !
    ! IA(II) is new number of connections of II
    ! II, JJ: OLD
    ! KK, LL: NEW => NEW = INVPR(OLD) 
    !                    = 0 if vertex does not exist in new graph
    !                OLD = PERMR(NEW)
    !
    mm = 0
    do ii = 1,nn
       kk = invpr(ii)
       if( kk /= 0 ) then
          mm = mm + 1
          do jjz = iacpy_subgra(ii),iacpy_subgra(ii+1)-1
             jj = jacpy_subgra(jjz)
             ll = invpr(jj)
             if( ll /= 0 ) ia(kk) = ia(kk) + 1
          end do
       end if
    end do
    !
    ! IA is now total number of edges on new subgraph for each node
    !
    call graphs_number_to_linked_list(mm,ia)
    !
    ! Compute subgraph
    !
    nnz  = 0
    nnz1 = 0
    do kk = 1,mm
       ii = permr(kk)
       pp = 0
       do jjz = iacpy_subgra(ii),iacpy_subgra(ii+1)-1
          jj = jacpy_subgra(jjz)
          ll = invpr(jj)
          if( ll /= 0 ) then    
             ja( ia(kk) + pp ) = ll
             pp = pp + 1
             nnz1 = nnz1 + 1
          end if
       end do
    end do
    !
    ! Order graph
    !
    do ii = 1,mm
       ll = ia(ii+1)-ia(ii)
       if( ll > 0 ) call maths_heap_sort(2_ip,ll,ja(ia(ii):))
    end do
    nnz = ia(mm+1)-1
    !
    ! Deallocate temporary graphs IACPY_SUBGRA, JACPY_SUBGRA
    !
    call memory_deallo(memor,'JACPY_SUBGRA','graphs_subgra',jacpy_subgra)
    call memory_deallo(memor,'IACPY_SUBGRA','graphs_subgra',iacpy_subgra)
    
    call memory_resize(memor,'IA','graphs_subgra',ia,mm+1_ip)             
    call memory_resize(memor,'JA','graphs_subgra',ja,nnz)             

  end subroutine graphs_subgra

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/01/1016
  !> @brief   Transform a sum to a linked list
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  pure subroutine graphs_number_to_linked_list(nn,ia)
    integer(ip), intent(in)    :: nn
    integer(ip), intent(inout) :: ia(nn+1)
    integer(ip)                :: ii,kk,ll

    kk    = ia(1)
    ia(1) = 1 
    do ii = 2,nn+1
       ll     = ia(ii)
       ia(ii) = ia(ii-1) + kk
       kk     = ll
    end do

  end subroutine graphs_number_to_linked_list
  
end module mod_graphs_basic
!> @}
