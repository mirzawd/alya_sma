!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_iterative_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Krylov solvers
!> @details Krylov solvers
!>         
!-----------------------------------------------------------------------

module def_iterative_solvers
  
  use def_kintyp_basic,                  only : ip,rp,lg
  use def_mat,                           only : mat
  use def_mat_dia,                       only : mat_dia
  use def_mat_csr,                       only : mat_csr
  use def_solvers,                       only : solver
  use def_solver,                        only : SOL_RIGHT_PRECONDITIONING
  use def_solver,                        only : SOL_LEFT_PRECONDITIONING
  use def_solver,                        only : SOL_LEFT_RIGHT_PRECONDITIONING
  use def_solver,                        only : SOL_SOLVER_DIAGONAL
  use mod_memory_basic,                  only : memory_alloca
  use mod_memory_basic,                  only : memory_deallo
  use mod_memory_basic,                  only : memory_size
  use mod_maths_basic,                   only : maths_inverse
  use mod_optional_argument,             only : optional_argument
  implicit none

  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  character(21), parameter :: vacal = 'def_iterative_solvers'

  private

  !----------------------------------------------------------------------
  !
  ! Iterative solver abstract class
  !
  !----------------------------------------------------------------------

  type :: preconditioner
     class(solver),             pointer                    :: p                       ! Should be allocated
     class(mat),                pointer                    :: a                       ! Preconditioner
     logical(lg)                                           :: a_is_allocated          ! If matrix is allocated or points
   contains
     procedure,                             pass           :: init   => init_preco
     procedure,                             pass           :: deallo => deallo_preco
  end type preconditioner
  type, extends(solver), abstract :: iterative_solver
     type(preconditioner),      pointer                    :: preco(:)
   contains
     procedure,                             pass           :: init_all                ! Common initializations
     procedure,                             pass           :: deallo_all              ! Common deallocate
     procedure,                             pass           :: left_preconditioning    ! Preconditioning L z = r
     procedure,                             pass           :: preconditioning         ! Preconditioning L z = A R^-1 r
     procedure,                             pass           :: right_preconditioning   ! Preconditioning   r = R z
     procedure,                             pass           :: start                   ! Starting operations
     procedure,                             pass           :: end                     ! Ending operations
     procedure,                             pass           :: tolerance               ! Compute tolerance
     procedure,                             pass           :: outcvg                  ! Output convergence
     procedure,                             pass           :: init_iterations         ! Initialize iterations
  end type iterative_solver
  
  !-----------------------------------------------------------------------
  !
  ! Different methods: Krylov, DDM, stationary
  !
  !-----------------------------------------------------------------------

  type, abstract, extends(iterative_solver) :: krylov
  end type krylov
  
  type, abstract, extends(krylov)           :: krylov_symmetric
  end type krylov_symmetric
  
  type, abstract, extends(krylov)           :: krylov_unsymmetric
  end type krylov_unsymmetric
  
  type, abstract, extends(iterative_solver) :: stationary
   contains
     procedure (stationary_mv), pass, deferred :: mv          ! Matrix vector product used for right preconditioning
  end type stationary

  type, abstract, extends(iterative_solver) :: ddm
  end type ddm

  abstract interface          
     subroutine stationary_mv(self,x,y,a)
       import                                    :: rp
       import                                    :: stationary
       import                                    :: mat
       class(stationary),          intent(in)    :: self          !< Solver
       real(rp),          pointer, intent(in)    :: x(:)          !< RHS b
       real(rp),          pointer, intent(inout) :: y(:)          !< Solve x
       class(mat),                 intent(in)    :: a             !< Matrix A
     end subroutine stationary_mv
  end interface
  
  
  public :: iterative_solver
  public :: preconditioner
  public :: krylov_symmetric
  public :: krylov_unsymmetric
  public :: stationary
  public :: ddm

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Initialize common variables of the class
  !> @details Initialization 
  !> 
  !-----------------------------------------------------------------------

  subroutine init_all(self)
    
    class(iterative_solver), intent(inout) :: self

    !call self % init_solver()
    !nullify(self % preco  )

  end subroutine init_all

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Initialize common variables of the class
  !> @details Initialization 
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_all(self)
    
    class(iterative_solver), intent(inout) :: self
    
    !call self % deallo_solver()
       
    !if( associated(self % preco) ) then
    !   do i = 1,size(self % preco)
    !      call self % preco(i) % deallo()
    !   end do
    !   deallocate(self % preco)
    !   nullify(self % preco)
    !end if

  end subroutine deallo_all

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Preconditioning
  !> @details Precondition L z = r
  !>          If z is not present, result is copied in r
  !> 
  !-----------------------------------------------------------------------

  subroutine left_preconditioning(self,r,z)

    class(iterative_solver),                    intent(inout) :: self          !< Solver
    real(rp),                          pointer, intent(in)    :: r(:)          !< Input r
    real(rp),                optional, pointer, intent(inout) :: z(:)          !< Output z
    real(rp),                          pointer                :: w(:)
    integer(ip)                                               :: i

    if( associated(self % preco) ) then

       nullify(w)
       if( present(z) ) then
          w => z
       else if( self % nn > 0 ) then
          allocate(w(self % nn))
       end if

       select case ( self % input % kfl_preco )

       case ( SOL_SOLVER_DIAGONAL )

          select type ( v => self % preco(1) % a )
          class is ( mat_dia ) ; call v % mv(r,w)
          class default        ; call self % preco(1) % p % solve(r,w,v)
          end select

       case default

          call self % preco(1) % p % solve(r,w,self % preco(1) % a)

       end select

       if( .not. present(z) .and. associated(w) ) then
          do i = 1,self % nn
             r(i) = w(i)
          end do
          deallocate(w)
       end if
       
    else

       if( present(z) ) then
          do i = 1,self % nn
             z(i) = r(i)
          end do
       end if

    end if

  end subroutine left_preconditioning

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Preconditioning
  !> @details Solver system: L z = A R^-1 r
  !>                           w = A R^-1 r
  !>          If A is not present, solve simply L z = r
  !> 
  !-----------------------------------------------------------------------

  subroutine preconditioning(self,r,z,w,a)

    class(iterative_solver),          intent(inout) :: self          !< Solver
    real(rp),                pointer, intent(in)    :: r(:)          !< RHS b
    real(rp),                pointer, intent(inout) :: z(:)          !< Solve x
    real(rp),      optional, pointer, intent(inout) :: w(:)          !< Working
    class(mat),    optional,          intent(in)    :: a             !< Matrix A
    integer(ip)                                     :: i
    
    if( present(a) ) then
       
       if( present(w) ) then
          
          select case ( self % input % kfl_leftr )

          case ( SOL_RIGHT_PRECONDITIONING ) 

             call self % left_preconditioning(r,w)  ! R z = r => z = R^-1 r
             call self % parallel_mv(a,w,z)         ! z   = A R^-1 r
             
          case ( SOL_LEFT_RIGHT_PRECONDITIONING )
             
             call self % left_preconditioning(r,z)  ! R z = r => z = R^-1 r
             call self % parallel_mv(a,z,w)         ! w   = A z   = A R^-1 r
             call self % left_preconditioning(w,z)  ! z   = L-1 w = L-1 A R^-1 r
          
          case ( SOL_LEFT_PRECONDITIONING )
             
             call self % parallel_mv(a,r,w)         ! w = A r
             call self % left_preconditioning(w,z)  ! z = L-1 w = L-1 A r
             
          end select

       else

          call runend('PRECONDITIONING: PROVIDE WORKING ARRAY W')
          
       end if
       
    else
      
       select case ( self % input % kfl_leftr )
       case ( SOL_LEFT_PRECONDITIONING , &
            & SOL_LEFT_RIGHT_PRECONDITIONING )
          call self % left_preconditioning(r,z)     ! L z = r
       case default
          do i = 1,self % nn
             z(i) = r(i)                            ! z = r
          end do
       end select
       
    end if
    
  end subroutine preconditioning

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Preconditioning
  !> @details Precondition r = R z
  !> 
  !-----------------------------------------------------------------------

  subroutine right_preconditioning(self,r,z)

    class(iterative_solver),                    intent(inout) :: self          !< Solver
    real(rp),                          pointer, intent(in)    :: r(:)          !< RHS b
    real(rp),                optional, pointer, intent(inout) :: z(:)          !< Solve x
    real(rp),                          pointer                :: w(:)          !< Solve x
    integer(ip)                                               :: i
    class(mat),                        pointer                :: a
    
    if(  associated(r) .and. &
         self % input % kfl_leftr == SOL_RIGHT_PRECONDITIONING      .or. &
         self % input % kfl_leftr == SOL_LEFT_RIGHT_PRECONDITIONING ) then

       nullify(w)
       if( present(z) ) then
          w => z
       else if( self % nn > 0 ) then
          allocate(w(self % nn))
       end if
       a => self % preco(1) % a
       
       select case ( self % input % kfl_preco )
          
       case ( SOL_SOLVER_DIAGONAL )
          
          select type ( v => self % preco(1) % a )
          class is ( mat_dia ) ; call v % solve(r,w) ! Diagonal preconditioner is inverse diag
          class is ( mat     ) ; call v % mv   (r,w) ! Just a MV product
          end select

       case default

          select type ( v => self % preco(1) % p )
          class is ( stationary ) ; call v % mv(r,w,a)
          class default           ; call runend('NOT CODED')
          end select          
          
       end select

       do i = 1,self % nn
          r(i) = w(i)
       end do

       if( .not. present(z) .and. associated(w) ) deallocate(w)

    end if
    
  end subroutine right_preconditioning

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Ending operation of solver
  !> @details Ending operation of solver
  !> 
  !-----------------------------------------------------------------------

  subroutine end(self,b,x,a,r)
    
    class(iterative_solver),                    intent(inout) :: self          !< Solver
    real(rp),                          pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),                          pointer, intent(inout) :: x(:)          !< Solve x
    real(rp),                optional, pointer, intent(inout) :: r(:)          !< Solve x
    class(mat),                                 intent(in)    :: a             !< Matrix A
    real(rp),                          pointer                :: y(:)
    real(rp)                                                  :: rnorm
    integer(ip)                                               :: nn

    !----------------------------------------------------------------------
    !
    ! Right preconditioning: Solve R x = x'
    !
    !----------------------------------------------------------------------

    select case ( self % input % kfl_leftr )
    case ( SOL_RIGHT_PRECONDITIONING, SOL_LEFT_RIGHT_PRECONDITIONING )
       call self % left_preconditioning(x)
    end select
    
    !----------------------------------------------------------------------
    !
    ! Final residual r = b - Ax
    !
    !----------------------------------------------------------------------

    if( present(r) ) then 
       rnorm = self % parallel_L2norm(r)
    else
       nullify(y)
       nn = a % ndof1 * a % nrows 
       call memory_alloca(self % memor,'Y','end',y,nn)
       call self % parallel_residual(a,x,b,y) 
       rnorm = self % parallel_L2norm(y)
       call memory_deallo(self % memor,'Y','end',y)
    end if

    self % output % resid_final = rnorm * self % output % bnorm_inv
    
  end subroutine end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Starting operation of solver
  !> @details Starting operation of solver
  !>
  !>
  !>          Outputs:
  !>
  !>          r ............................ residual r = b-Ax
  !>          z ............................ preconditioned residual z = L^-1 r
  !>          self % output % bnorm ........ RHS norm ||b||
  !>          self % output % bnorp ........ preconditioned RHS norm ||b||'
  !>                                         (b,L^-1 b )^1/2 for symemtric solvers (CG, etc.)
  !>                                         ||L^-1 b||      for unsymmetric solvers
  !>          self % output % bnorm_inv .... 1/bnorm
  !>          self % output % bnorp_inv .... 1/bnorp
  !>          self % output % resid_init ... normalized residual norm || r || / || b ||
  !>             
  !> 
  !-----------------------------------------------------------------------

  subroutine start(self,b,x,a,r,z,ierro)

    class(iterative_solver),          intent(inout) :: self          !< Solver
    real(rp),                pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),                pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),                       intent(in)    :: a             !< Matrix A
    real(rp),                pointer, intent(inout) :: r(:)          !< Residual
    real(rp),     optional,  pointer, intent(inout) :: z(:)          !< Preconditioned residual
    integer(ip),                      intent(out)   :: ierro         !< Error status
    integer(ip)                                     :: n,nn,ndof,i
    real(rp),                pointer                :: w(:)
    real(rp)                                        :: rnorm

    !----------------------------------------------------------------------
    !
    ! Allocate and initialize
    !    
    !----------------------------------------------------------------------
    
    nullify(w)
      
    call cputim(self % output % time_init)
    self % output % time_final  = self % output % time_init
    self % output % time_old    = self % output % time_init
    
    ierro                       = 0
    n                           = a % nrows
    ndof                        = a % ndof1
    nn                          = n * ndof    
    self % output % iters       = 0
    self % output % resip_init  = 1.0_rp
    self % output % resip_old   = 1.0_rp
    self % output % resip_final = 1.0_rp
    self % output % xorth       = 0.0_rp

    !----------------------------------------------------------------------
    !
    ! w = L^-1 b ... Preconditioned RHS
    !    
    !----------------------------------------------------------------------

    call memory_alloca(self % memor,'W','solve',w,nn)
    call self % preconditioning(b,w)    
    
    !----------------------------------------------------------------------
    !
    ! bnorm = || b || 
    ! bnorp = (b , L^-1 b )^1/2 for symemtric solvers (CG, etc.)
    !       = || L^-1 b ||      for unsymmetric solvers
    !
    !----------------------------------------------------------------------

    select type ( self )
    class is ( krylov_symmetric ) ; self % output % bnorp = sqrt(self % parallel_dot(b,w))
    class default                 ; self % output % bnorp = sqrt(self % parallel_dot(w,w))
    end select 
    self % output % bnorm     = self % parallel_L2norm(b)
    self % output % bnorp_inv = maths_inverse(self % output % bnorp,0.0_rp)
    self % output % bnorm_inv = maths_inverse(self % output % bnorm,0.0_rp)

    !----------------------------------------------------------------------
    !
    ! bnorm = 0 ... trivial solution x = 0
    !
    !----------------------------------------------------------------------
        
    if( self % output % bnorm <= self % input % bnorm_min .and. self % input % solco >= 0.0_rp ) then
       do i = 1,nn
          x(i) = 0.0_rp
        end do
        self % output % resip_init  = 0.0_rp
        self % output % resip_final = 0.0_rp
        self % output % resip_old   = 0.0_rp
        ierro                       = -1
        if( self % input % kfl_cvgso == 1 ) call self % outcvg(b,x,a,r)
        goto 10
     end if
     
     !----------------------------------------------------------------------
     !
     ! r             = b - Ax .............. Residual
     ! resid_initial = || r || / || b || ... Normalized residual norm
     !
     !----------------------------------------------------------------------

     call self % parallel_residual(a,x,b,r) 
     rnorm = self % parallel_L2norm(r)
     self % output % resid_init = rnorm * self % output % bnorm_inv
     call self % outcvg(b,x,a,r)

     !----------------------------------------------------------------------
     !
     ! Right preconditioning: initial x0 = R x
     !
     !----------------------------------------------------------------------

     call self % right_preconditioning(x)
     
     !----------------------------------------------------------------------
     !
     ! z = L^{-1} r ... Preconditioned residual
     !
     !----------------------------------------------------------------------

     if( present(z) ) call self % preconditioning(r,z)  

10   continue
     
     call memory_deallo(self % memor,'W','solve',w)

  end subroutine start

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-04
  !> @brief   Tolerance
  !> @details Compute solver tolerance according to the strategy
  !>          toler = ||b'|| * eps ... Stopping criterion
  !> 
  !-----------------------------------------------------------------------

  function tolerance(self) result(toler)

    class(iterative_solver), intent(in) :: self          !< Solver
    real(rp)                            :: eps
    real(rp)                            :: toler
    
    select case ( self % input % kfl_adres )
    case ( 1_ip ) ; eps = max( self % output % resid_init * self % input % adres , self % input % solmi )
    case default  ; eps = self % input % solco
    end select

    toler = eps * self % output % bnorp

  end function tolerance
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Initialize preconditioner
  !> @details Initialization of the preconditioners
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_preco(preco)
    
    class(preconditioner), intent(inout) :: preco

    nullify(preco % p)
    nullify(preco % a)
    preco % a_is_allocated = .false.
    
  end subroutine init_preco

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Deallocate preconditioner
  !> @details Deallocate preconditioner
  !> 
  !-----------------------------------------------------------------------
  
  recursive subroutine deallo_preco(preco)
    
    class(preconditioner), intent(inout) :: preco
    
    if( preco % a_is_allocated .and. associated(preco % a) ) then
       call preco % a % deallo()
       if( associated(preco % a) ) deallocate(preco % a)
    end if
    
    if( associated(preco % p) ) then
       !select type ( v => preco % p )
       !class is ( iterative_solver )
       !   if( associated(v % preco) ) then
       !      if( associated(v % preco(1) % p) ) then
       !         call v % preco(1) % p % deallo()
       !      end if
       !   end if
       !end select
       call preco % p % deallo()
       deallocate(preco % p)
    end if
    
  end subroutine deallo_preco
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Output convergence
  !> @details Output convergence
  !>          If required, compute RNORM = || b - Ax || / || b || 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine outcvg(self,b,x,a,r)

    class(iterative_solver),          intent(inout) :: self
    real(rp),                pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),                pointer, intent(in)    :: x(:)          !< Solve x
    class(mat),                       intent(in)    :: a             !< Matrix A
    real(rp),    optional,   pointer, intent(inout) :: r(:)          !< Residual
    real(rp),                pointer                :: y(:)
    real(rp)                                        :: tdiff,telap
    real(rp)                                        :: rnorm
    
    if( self % input % kfl_cvgso /= 0 ) then

       if( self % input % kfl_exres /= 0 ) then
          
          if ( self % input % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
             !
             ! Right preconditioning
             !
             rnorm = self % output % resip_final
          else
             !
             ! Left preconditioning
             !
             if( present(r) ) then  
                rnorm = self % parallel_L2norm(r)
             else
                nullify(y)
                call memory_alloca(self % memor,'Y',vacal,y,a % nrows * a % ndof1)
                call self % parallel_residual(a,x,b,y)
                rnorm = self % parallel_L2norm(y)             
                call memory_deallo(self % memor,'Y',vacal,y)
             end if
             rnorm = rnorm * self % output % bnorm_inv
          end if

       else
          rnorm = 0.0_rp
       end if

       if( self % iwrite ) then
          self % output % time_old = self % output % time_final
          call cputim(self % output % time_final)       
          tdiff = self % output % time_final - self % output % time_old
          telap = self % output % time_final - self % output % time_init
          write(self % input % lun_cvgso,1) & 
               self % output % iters,       & ! 1. Iteration number
               self % output % resip_final, & ! 2. Preconditioned residual norm
               rnorm,                       & ! 3. Non-preconditioned residual norm
               tdiff,                       & ! 4. Time since last iteration
               telap,                       & ! 5. Time since beginning
               self % output % xorth          ! 6. Orthogonality
       end if
    end if
    
1   format(i7,20(1x,e16.8E3))

  end subroutine outcvg


  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Initital residuals and tolerance
  !> @details Parallel sum
  !>
  !----------------------------------------------------------------------
  
  subroutine init_iterations(self,resid) 

    class(iterative_solver), intent(inout) :: self
    real(rp),                intent(in)    :: resid
    
    self % output % resip_init  = resid * self % output % bnorp_inv
    self % output % resip_old   = self % output % resip_init 
    self % output % resip_final = self % output % resip_init 
    self % output % toler       = self % tolerance()

  end subroutine init_iterations
  
end module def_iterative_solvers
!> @}

