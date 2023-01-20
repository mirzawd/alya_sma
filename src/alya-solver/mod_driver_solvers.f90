!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_all_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   All solvers
!> @details All solvers
!>         
!-----------------------------------------------------------------------

module mod_driver_solvers

  use def_mat,             only : mat
  use def_preconditioners, only : preconditioner_setup
  use mod_strings,         only : integer_to_string
  use def_all_solvers
  use def_solver
  implicit none
  private

  public :: solver_init
  public :: solver_alloca
  public :: solver_dim
  public :: solver_setup
  public :: solver_preprocess
  public :: solver_solve
  public :: solver_postprocess
  public :: solver_deallo

  public :: solver_all

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solve a complete system
  !> @details Solve a complete system
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_all(self,b,x,a)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    real(rp),         pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                intent(inout) :: a             !< Matrix A
    
    call solver_preprocess (self,b,x,a)
    call solver_solve      (self,b,x,a)
    call solver_postprocess(self,b,x,a)

  end subroutine solver_all
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solve a complete system
  !> @details Solve a complete system
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_solve(self,b,x,a)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    real(rp),         pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                intent(inout) :: a             !< Matrix A
    
    call self % solve(b,x,a)

  end subroutine solver_solve
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Set solver
  !> @details Allocate solver
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_alloca(solve,name)

    class(solver), pointer, intent(inout) :: solve
    integer(ip),            intent(in)    :: name

    if( associated(solve) ) then

       call runend('SOLVER_ALLOCA: SOLVER ALREADY ALLOCATED')

    else

       select case ( name )

       case ( SOL_SOLVER_CG                  ) ; allocate( cg                :: solve )
       case ( SOL_SOLVER_DEFLATED_CG         ) ; allocate( dcg               :: solve ) 
       case ( SOL_SOLVER_BICGSTAB            ) ; allocate( bicgstab          :: solve )
       case ( SOL_SOLVER_GMRES               ) ; allocate( gmres             :: solve )
       case ( SOL_SOLVER_MINRES              ) ; allocate( minres            :: solve )
       case ( SOL_SOLVER_JACOBI              ) ; allocate( jacobi            :: solve )
       case ( SOL_SOLVER_RICHARDSON          ) ; allocate( richardson        :: solve )
       case ( SOL_SOLVER_SSOR                ) ; allocate( ssor              :: solve )
       case ( SOL_SOLVER_SOR                 ) ; allocate( sor               :: solve )
       case ( SOL_SOLVER_GAUSS_ELIMINATION   ) ; allocate( gauss_elimination :: solve )
       case ( SOL_SOLVER_APPROX_INVERSE      ) ; allocate( approx_inv        :: solve )
       case ( SOL_SOLVER_LINELET             ) ; allocate( linelet           :: solve )
       case ( SOL_SOLVER_DIRECT              ) ; allocate( lu_factorization  :: solve )
       case ( SOL_DIRECT_SOLVER_SKYLINE_ALYA ) ; allocate( cholesky          :: solve )
       case ( SOL_SOLVER_RAS                 ) ; allocate( ras               :: solve )
       case default ; call runend('SOLVER_ALLOCA: UNKNOWN SOLVER '//integer_to_string(name))

       end select

    end if

  end subroutine solver_alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solver process
  !> @details Solver process
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solver_natural(self,b)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    integer(ip)                              :: i,idof,k

    if( self % input % kfl_bvnat /= 0 ) then
       do i = 1,self % nn
          do idof = 1,self % ndof
             k = (i-1)*self % ndof + idof
             if( self % arrays % fixno(k) <= 0 ) then
                b(k) = b(k) + self % arrays % bvnat(k)
             end if
          end do
       end do       
    end if
    
  end subroutine solver_natural
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solver process
  !> @details Solver process
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solver_dirichlet(self,b,x,a)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    real(rp),         pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                intent(inout) :: a             !< Matrix A
    integer(ip)                              :: i,idof,k
    
    if( associated(self % arrays % fixno) ) then
       
       select case ( self % input % kfl_iffix )
          
       case ( 1_ip )
          !
          ! Fixity and value prescribed 
          !
          if( associated(self % arrays % bvess) ) then
             call a % dirichlet(self % arrays % fixno,self % arrays % bvess,b)
             do i = 1,self % nn
                do idof = 1,self % ndof
                   k = (i-1)*self % ndof + idof
                   if( self % arrays % fixno(k) > 0 ) then
                      x(k) = self % arrays % bvess(k)
                   end if
                end do
             end do
          else
             call a % dirichlet(self % arrays % fixno)
             do i = 1,self % nn
                do idof = 1,self % ndof
                   k = (i-1)*self % ndof + idof
                   if( self % arrays % fixno(k) > 0 ) then
                      b(k) = 0.0_rp
                      x(k) = 0.0_rp
                   end if
                end do
             end do
          end if
          
       case ( 2_ip )
          !
          ! Fixity and 0 value
          !
          call a % dirichlet(self % arrays % fixno)
          do i = 1,self % nn
             do idof = 1,self % ndof
                k = (i-1)*self % ndof + idof
                if( self % arrays % fixno(k) > 0 ) then
                   b(k) = 0.0_rp
                   x(k) = 0.0_rp
                end if
             end do
          end do
          
       end select
    end if

  end subroutine solver_dirichlet
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solver process
  !> @details Solver process
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solver_preprocess(self,b,x,a)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    real(rp),         pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                intent(inout) :: a             !< Matrix A
    !
    ! Exchange RHS
    !
    if( self % input % kfl_exchange == 1 ) then
       select type ( self )
       class is ( iterative_solver ) ; call self % parallel_exchange(b)
       end select
    end if
    !
    ! Impose Dirichlet conditions if required
    !
    call solver_dirichlet(self,b,x,a)
    !
    ! Impose Natural conditions if required
    !
    call solver_natural(self,b)
    !
    ! Setup solver (numerical factorization for direct solvers)
    !
    call solver_setup(self,a)
    !
    ! Preconditioner setup
    !
    select type ( self )
    class is ( iterative_solver ) ; call preconditioner_setup(self,a)
    end select
    
  end subroutine solver_preprocess
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Solver postprocess
  !> @details Solver postprocess: deallocate solver here?
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solver_postprocess(self,b,x,a)
    
    class(solver),             intent(inout) :: self          !< Solver
    real(rp),         pointer, intent(inout) :: b(:)          !< RHS b
    real(rp),         pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                intent(inout) :: a             !< Matrix A
    !
    ! Output convergence
    !
    continue

  end subroutine solver_postprocess

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Set up all solvers
  !> @details Set up all solvers
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_dim(self,a)

    class(solver),              intent(inout) :: self          !< Solver
    class(mat),                 intent(inout) :: a             !< Matrix A

    call self % dim(a)
    
  end subroutine solver_dim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Set up all solvers
  !> @details Set up all solvers
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_setup(self,a)

    class(solver),              intent(inout) :: self          !< Solver
    class(mat),                 intent(inout) :: a             !< Matrix A
    
    self % n    = a % nrows
    self % ndof = a % ndof1
    self % nn   = self % n * self % ndof
    
    if( associated(self % comm) ) then
       self % nint = self % comm % npoi1
       self % nown = self % comm % npoi3
    else
       self % nint = a % nrows
       self % nown = a % nrows
    end if

    call self % dim  (a)
    call self % setup(a)
    
  end subroutine solver_setup

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Initialize common variables of the class
  !> @details Initialization 
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_init(self)
    
    class(solver), intent(inout) :: self

    nullify(self % comm)
    self % memor             = 0_8
    self % iwrite            = .false.
    self % n                 = 0
    self % nown              = 0
    self % nint              = 0
    self % nn                = 0
    self % ndof              = 0
    self % comm_is_allocated = .true.

    call self % input  % init()
    call self % output % init()
    call self % arrays % init()
    call self % init()

    select type ( self )
    class is ( iterative_solver ) ; nullify(self % preco)
    end select
    
  end subroutine solver_init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Deallocate solver
  !> @details Deallocate solver
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_deallo(self)
    
    class(solver), intent(inout) :: self
    integer(ip)                  :: i
    
    if( associated(self % comm) ) then
       select case ( self % comm_is_allocated )
       case ( .true.  ) ; call self % comm % deallo()
       case ( .false. ) ; nullify(self % comm)
       end select       
    end if
    
    call self % deallo()

    select type ( self )
    class is ( iterative_solver ) 
       if( associated(self % preco) ) then
          do i = 1,size(self % preco)
             call self % preco(i) % deallo()
          end do
          deallocate(self % preco)
          nullify(self % preco)
       end if
    end select
    
  end subroutine solver_deallo

end module mod_driver_solvers
!> @}
