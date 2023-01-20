!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_ras.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   RAS preconditioner.
!> @details RAS. If you want to use RAS as a solver, you have to 
!>          select richardson solver and RAS as a preconditioner.
!>          By default the preconditioner takes one domain per MPI
!>          subdomains. But inside subdomain, a multiplicative Schwarz
!>          can be used as well with smaller subdomains by giving as
!>          input LGROU, the list of nodes for each group. There, any overlap
!>          is possible. An addition array PGROU enables to select the
!>          responsible nodes of each subdomains, just like in the usual
!>          RAS.
!>         
!-----------------------------------------------------------------------

module def_ras

  use def_kintyp,             only : ip,rp,i1p
  use def_mat,                only : mat
  use def_mat_csr,            only : mat_csr
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_memory,             only : memory_copy
  use mod_memory,             only : memory_size
  use def_direct_solvers,     only : direct_factorization
  use def_iterative_solvers,  only : ddm
  use mod_alya_direct_solver, only : alya_direct_solver
  use mod_graphs_basic,       only : graphs_permut_metis_postordering
  use mod_graphs_basic,       only : graphs_subgra
  implicit none
  private

  type                         :: dom_typ
     type(alya_direct_solver)  :: dir
     integer(ip),  pointer     :: ia(:)       ! Graph
     integer(ip),  pointer     :: ja(:)       ! Graph
     integer(ip),  pointer     :: perm(:)     ! Permutation
     integer(ip),  pointer     :: invp(:)     ! Inverse permutation
     integer(ip),  pointer     :: perm_ras(:) ! Permutation wrt to original problem
     integer(ip),  pointer     :: invp_ras(:) ! Inverse permutation wrt to original problem
     type(mat_csr)             :: a           ! Matrix
   contains
     procedure,           pass :: init   => init_dom
     procedure,           pass :: deallo => deallo_dom
  end type dom_typ
  type, extends(ddm)           :: ras
     type(dom_typ), pointer    :: dom(:)
     integer(ip)               :: ndom
     integer(ip),   pointer    :: pgrou(:)
   contains
     procedure,           pass :: init     
     procedure,           pass :: alloca   
     procedure,           pass :: deallo   
     procedure,           pass :: solve    
     procedure,           pass :: set    
     procedure,           pass :: setup    
     procedure,           pass :: symbolical     
     procedure,           pass :: numerical     
     procedure,           pass :: reordering     
     procedure,           pass :: cleaning     
     procedure,           pass :: partial_cleaning     
     procedure,           pass :: subsystems        ! Construct subsystems fro original one
  end type ras

  character(3), parameter :: vacal = 'ras'

  public :: ras
    
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init_dom(self)

    class(dom_typ), intent(inout) :: self
    
    call self % dir % init()
    call self % a   % init()
    nullify(self % ia)
    nullify(self % ja)
    nullify(self % perm)
    nullify(self % invp)
    nullify(self % perm_ras)
    nullify(self % invp_ras)
        
  end subroutine init_dom
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(ras), intent(inout) :: self

    nullify(self % dom)
    nullify(self % pgrou)
    self % ndom = 1

  end subroutine init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Partial cleaning
  !> @details Partial cleaning
  !> 
  !-----------------------------------------------------------------------

  subroutine partial_cleaning(self)
    
    class(ras), intent(inout) :: self
    integer(ip)               :: idom

    do idom = 1,self % ndom
       call self % dom(idom) % dir % deallo_numerical()
    end do
    
  end subroutine partial_cleaning
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Partial cleaning
  !> @details Partial cleaning
  !> 
  !-----------------------------------------------------------------------

  subroutine cleaning(self)
    
    class(ras), intent(inout) :: self
    integer(ip)               :: idom

    do idom = 1,self % ndom
       call self % dom(idom) % dir % deallo_numerical ()
       call self % dom(idom) % dir % deallo_symbolical()
       call self % dom(idom) % deallo(self % memor)
    end do
    if( associated(self % dom) ) deallocate(self % dom)
    
  end subroutine cleaning
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    
    class(ras), intent(inout) :: self

    call self % cleaning()
    call memory_deallo(self % memor,'SELF % PGROU',vacal,self % pgrou)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_dom(self,memor)
    
    class(dom_typ), intent(inout) :: self
    integer(8),     intent(inout) :: memor(2)
    
    call memory_deallo(memor,'SELF % IA'      ,vacal,self % ia)
    call memory_deallo(memor,'SELF % JA'      ,vacal,self % ja)
    call memory_deallo(memor,'SELF % PERM'    ,vacal,self % perm)
    call memory_deallo(memor,'SELF % INVP'    ,vacal,self % invp)
    call memory_deallo(memor,'SELF % PERM_RAS',vacal,self % perm_ras)
    call memory_deallo(memor,'SELF % INVP_RAS',vacal,self % invp_ras)

    call self % a % deallo(memor)
    
  end subroutine deallo_dom

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Symbolical factorization
  !> @details Symbolical factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine symbolical(self,a)
    
    class(ras), intent(inout) :: self
    class(*),   intent(in)    :: a    !< Matrix A
    integer(ip)               :: idom

    if( self % info(1_ip) ) write(*,*) 'RAS: symbolical factorization'
    do idom = 1,self % ndom
       call self % dom(idom) % dir % factor_symbolical(       &
            self % dom(idom) % a      ,                       &
            self % dom(idom) % ia     , self % dom(idom) % ja )
    end do
    
  end subroutine symbolical
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Numerical factorization
  !> @details Numerical factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine numerical(self,a)
    
    class(ras), intent(inout) :: self
    class(*),   intent(in)    :: a    !< Matrix A
    integer(ip)               :: idom

    if( self % info(1_ip) ) write(*,*) 'RAS: numerical factorization'
    select type ( a )
    class is ( mat_csr )
       do idom = 1,self % ndom
          call self % dom(idom) % dir % factor_numerical(         &
               self % dom(idom) % a    ,                          &
               self % dom(idom) % invp , self % dom(idom) % perm, &
               self % dom(idom) % ia   , self % dom(idom) % ja    )
       end do
    end select
    
  end subroutine numerical
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)
    
    class(ras), intent(inout) :: self
    integer(ip)               :: idom
    
    select case ( self % input % kfl_block_ras )
    case ( 0_ip ) ; self % ndom  = self % input % num_subd_ras
    case default  ; self % ndom  = self % input % num_subd_ras
    end select

    allocate(self % dom(self % ndom))
    do idom = 1,self % ndom
       call self % dom(idom) % init()
    end do

  end subroutine alloca
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Reordering
  !> @details Reordering
  !> 
  !-----------------------------------------------------------------------

  subroutine reordering(self,a)
    
    class(ras), intent(inout) :: self          !< Solver
    class(mat), intent(in)    :: a             !< Matrix A
    integer(ip)               :: n,nz,i,idom

    if( self % info(1_ip) ) write(*,*) 'RAS: reordering'

    select type ( a )
    class is ( mat_csr )

       do idom = 1,self % ndom          
          n  = self % dom(idom) % a % nrows 
          nz = self % dom(idom) % a % nz

          call memory_alloca(self % memor,'SELF % DOM % PERM',vacal,self % dom(idom) % perm,n     )
          call memory_alloca(self % memor,'SELF % DOM % INVP',vacal,self % dom(idom) % invp,n     )
          call memory_alloca(self % memor,'SELF % DOM % IA'  ,vacal,self % dom(idom) % ia  ,n+1_ip)
          call memory_alloca(self % memor,'SELF % DOM % JA'  ,vacal,self % dom(idom) % ja  ,nz    )
          do i = 1,n+1
             self % dom(idom) % ia(i) = self % dom(idom) % a % ia(i)
          end do
          do i = 1,nz
             self % dom(idom) % ja(i) = self % dom(idom) % a % ja(i)
          end do
          call graphs_permut_metis_postordering(&
               n                            , nz                        , &
               self % dom(idom) % ia        , self % dom(idom) % ja     , &
               self % dom(idom) % perm      , self % dom(idom) % invp   , &
               memor=self % memor)

       end do

    end select

  end subroutine reordering
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Set
  !> @details Set preconditioner
  !>          If lgrou is provided, we are doing schwarz inside subdomains:
  !>          lgrou is overlapping: multiplicative overlapping schwarz
  !>          lgrou is disjoint:    non overlapping multiplicative schwarz
  !>          If pgrou is present, then only one group is responsible
  !>          for one node, just like the restricted additive schwarz
  !>
  !-----------------------------------------------------------------------

  subroutine set(self,a,lgrou,pgrou)

    class(ras),                      intent(inout) :: self            !< Solver
    class(*),                        intent(inout) :: a               !< Matrix A
    type(i1p),    pointer, optional, intent(inout) :: lgrou(:)
    integer(ip),  pointer, optional, intent(inout) :: pgrou(:)

    if( present(pgrou) ) then
       call memory_copy(self % memor,'PGROU',vacal,pgrou,self % pgrou,'DO_NOT_DEALLOCATE')
    end if
    
    select type ( a )
    class is ( mat_csr )
       call self % subsystems  (a,lgrou)
       call self % reordering  (a)
       call self % symbolical  (a)
    end select
    
  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Setup
  !> @details Setup and reorder to minimize fillin
  !> 
  !-----------------------------------------------------------------------

  subroutine subsystems(self,a,lgrou)

    class(ras),                        intent(inout) :: self          !< Solver
    class(mat_csr),                    intent(inout) :: a             !< Matrix A
    type(i1p),      pointer, optional, intent(inout) :: lgrou(:)
    integer(ip)                                      :: idofn
    integer(ip)                                      :: nras,idom
    integer(ip)                                      :: iras,igrou
    integer(ip)                                      :: n,i,k
    !
    ! Allocate domains
    !
    if( present(lgrou) ) then
       self % input % num_subd_ras = max(1_ip,memory_size(lgrou))
    else
       self % input % num_subd_ras = 1
    end if

    call self % alloca()
    
    select case ( self % input % kfl_block_ras )
    case ( 0_ip ) ; nras  = 1          ; idofn = a % ndof1       
    case default  ; nras  = a % ndof1  ; idofn = 1
    end select
    !
    ! Define each subsystem for each domain
    !
    idom = 0
    n    = a % nrows
    if( self % input % num_subd_ras > 1 ) then
       !
       ! Extended RAS
       !
       if( self % info(1_ip) ) write(*,*) 'RAS: mini-subdomain graphs'
       do igrou = 1,self % input % num_subd_ras
          do iras = 1,nras
             idom = idom + 1          
             call memory_alloca(self % memor,'SELF % DOM % PERM_RAS',vacal,self % dom(idom) % perm_ras,n)
             call memory_alloca(self % memor,'SELF % DOM % INVP_RAS',vacal,self % dom(idom) % invp_ras,n)

             do k = 1,memory_size(lgrou(igrou) % l)
                i                              = lgrou(igrou) % l(k)
                self % dom(idom) % perm_ras(k) = i
                self % dom(idom) % invp_ras(i) = k
             end do

             self % dom(idom) % a % nrows = memory_size(lgrou(igrou) % l)
             self % dom(idom) % a % ncols = memory_size(lgrou(igrou) % l)
             self % dom(idom) % a % ndof1 = idofn
             self % dom(idom) % a % ndof2 = idofn
             self % dom(idom) % a % nz    = a % nz
             
             call memory_copy(self % memor,'SELF % DOM % IA',vacal,a % ia,self % dom(idom) % a % ia,'DO_NOT_DEALLOCATE')
             call memory_copy(self % memor,'SELF % DOM % JA',vacal,a % ja,self % dom(idom) % a % ja,'DO_NOT_DEALLOCATE')

             call graphs_subgra(&
                  n,&
                  self % dom(idom) % a % nz  ,self % dom(idom) % perm_ras,&
                  self % dom(idom) % invp_ras,self % dom(idom) % a % ia  ,&
                  self % dom(idom) % a % ja  ,memor=self % memor)

          end do
       end do
       
    else
       !
       ! Full RAS
       !      
       do iras = 1,nras
          
          idom                         = idom + 1
          self % dom(idom) % a % nrows = a % nrows
          self % dom(idom) % a % ncols = a % nrows
          self % dom(idom) % a % ndof1 = idofn
          self % dom(idom) % a % ndof2 = idofn
          self % dom(idom) % a % nz    = a % nz
          
          call memory_copy  (self % memor,'SELF % DOM % IA',vacal,a % ia,self % dom(idom) % a % ia,'DO_NOT_DEALLOCATE')
          call memory_copy  (self % memor,'SELF % DOM % JA',vacal,a % ja,self % dom(idom) % a % ja,'DO_NOT_DEALLOCATE')
          
       end do
       
    end if

  end subroutine subsystems
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Setup
  !> @details Setup and reorder to minimize fillin
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(ras),            intent(inout) :: self          !< Solver
    class(mat),            intent(inout) :: a             !< Matrix A
    integer(ip)                          :: idom

    select type ( a )
    class is ( mat_csr ) 
       !
       ! Copy matrices
       !
       if( self % input % kfl_block_ras == 0 ) then
          !
          ! NDOFN degrees of freedom
          !
          do idom = 1,self % ndom         
             call self % dom(idom) % a % alloca_matrix(MEMORY_COUNTER=self % memor)
             call self % dom(idom) % a % copy(a,self % dom(idom) % invp_ras)
          end do
       else
          !
          ! Degrees of freedom treated separately
          !          
          stop 2
       end if
    end select
    !
    ! Numerical factorization of block matrices
    !
    call self % numerical(a)

  end subroutine setup 
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Solve
  !> @details Factorize, solve and partial clean
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)

    class(ras),               intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A
    integer(ip)                             :: i,idom,nn
    integer(ip)                             :: k,l,n,ndof
    real(rp),        pointer                :: b_ras(:)
    real(rp),        pointer                :: x_ras(:)

    nullify(x_ras,b_ras)
    
    if( self % input % kfl_block_ras == 0 ) then
       if( self % input % num_subd_ras == 1 ) then
          call self % dom(1) % dir % solve(b,x,self % dom(1) % a,self % dom(1) % invp)
       else
          n    = self % n
          nn   = self % nn
          ndof = self % ndof
          call memory_alloca(self % memor,'B_RAS',vacal,b_ras,nn)
          call memory_alloca(self % memor,'X_RAS',vacal,x_ras,nn)
          do idom = 1,self % input % num_subd_ras
             !
             ! Global to local system
             !
             do i = 1,self % dom(idom) % a % nrows
                k = self % dom(idom) % perm_ras(i)
                do l = 1,ndof
                   b_ras((i-1)*ndof+l) = b((k-1)*ndof+l)
                end do
             end do
             call self % dom(1) % dir % solve(b_ras,x_ras,self % dom(idom) % a,self % dom(idom) % invp)
             !
             ! Local to global system
             !
             if( associated(self % pgrou) ) then
                do i = 1,self % dom(idom) % a % nrows
                   k = self % dom(idom) % perm_ras(i)
                   if( self % pgrou(k) == idom ) then
                      do l = 1,ndof
                         x((k-1)*ndof+l) = x_ras((i-1)*ndof+l)
                      end do
                   end if                
                end do
             else                
                do i = 1,self % dom(idom) % a % nrows
                   k = self % dom(idom) % perm_ras(i)
                   do l = 1,ndof
                      x((k-1)*ndof+l) = x_ras((i-1)*ndof+l)
                   end do
                end do
             end if
          end do
          call memory_deallo(self % memor,'B_RAS',vacal,b_ras)
          call memory_deallo(self % memor,'X_RAS',vacal,x_ras)
       end if
    end if
    
  end subroutine solve

end module def_ras
!> @}
