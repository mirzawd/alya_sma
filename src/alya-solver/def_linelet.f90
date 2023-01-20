!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_linelet.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Linelet solver
!> @details Classical loop over linelets
!>          do iline = 1,self % nline
!>            do ipoin = self % lline(iline),self % lline(iline+1)-1
!>               jpoin = self % lrenup(ipoin)
!>               ...
!>            end do
!>          end do
!>         
!-----------------------------------------------------------------------

module def_linelet

  use def_kintyp,            only : ip,rp,lg,i1p
  use def_mat,               only : mat
  use def_mat_den,           only : mat_den
  use def_mat_csr,           only : mat_csr
  use def_mat_dia,           only : mat_dia
  use def_mat_tri,           only : mat_tri
  use def_direct_solvers,    only : direct_solver
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_maths_basic,       only : maths_inverse
  use mod_strings
  implicit none
  private
  
  type, extends(direct_solver) :: linelet
     integer(ip)             :: nline
     integer(ip),  pointer   :: lline(:)
     logical(lg),  pointer   :: inlin(:)
     integer(ip),  pointer   :: perm(:)
     type(mat_tri)           :: tri
     type(mat_dia)           :: dia
     real(rp),     pointer   :: xx(:,:)
     real(rp),     pointer   :: bb(:,:)
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
     procedure,         pass :: set
  end type linelet

  character(7), parameter :: vacal = 'linelet'

  public :: linelet
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(linelet), intent(inout) :: self
        
    call self % tri % init()
    call self % dia % init()
    
    self % nline = 0_ip
    nullify(self % lline)  
    nullify(self % inlin)  
    nullify(self % perm)
    nullify(self % xx)
    nullify(self % bb)

  end subroutine init  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    
    class(linelet), intent(inout) :: self

    call self % tri % deallo()
    call self % dia % deallo()
    call memory_deallo(self % memor,'XX'   ,vacal,self % xx   )
    call memory_deallo(self % memor,'BB'   ,vacal,self % bb   )
    call memory_deallo(self % memor,'INLIN',vacal,self % inlin)
       
  end subroutine deallo
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)

    class(linelet), intent(inout) :: self          !< Solver
    integer(ip)                   :: ndof,nl

    if( self % n > 0 ) then
       ndof = self % ndof
       if( self % nline > 0 ) then
          nl = self % lline(self % nline+1)-1
       else
          nl = 0
       end if
       call memory_alloca(self % memor,'XX'   ,vacal,self % xx   ,ndof,nl)
       call memory_alloca(self % memor,'BB'   ,vacal,self % bb   ,ndof,nl)       
       call memory_alloca(self % memor,'INLIN',vacal,self % inlin,self % n)       
       self % tri % nrows = nl
       self % tri % ndof1 = self % ndof
       call self % tri % alloca()
    end if
    
  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Set
  !> @details Set preconditioner
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,nline,lline,perm)

    class(linelet),        intent(inout) :: self          !< Solver
    integer(ip),           intent(in)    :: nline
    integer(ip),  pointer, intent(in)    :: lline(:)
    integer(ip),  pointer, intent(in)    :: perm(:)

    self % nline  =  nline
    self % lline  => lline
    self % perm   => perm

  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Set up solver
  !> @details Set up solver. Fill in matrix
  !>
  !>          ia   i   ib
  !>          o----o---o
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(linelet),                     intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A
    integer(ip)                                       :: il,i,k,r
    integer(ip)                                       :: ia,ib,ndof,idof
    real(rp),                  pointer                :: val(:,:)
    !
    ! Setup and allocate
    !
    call self % alloca()
    !
    ! Diagonal
    !
    call a    % diag             (self % dia)
    call self % parallel_exchange(self % dia)
    call self % dia % inverse    ()    
    !
    ! Dimensions
    !
    ndof = self % ndof

    allocate(val(ndof,ndof))
    r = 0

    do il = 1,self % nline
       !
       ! Beginning of linelet
       !
       i                                    = self % perm(self % lline(il))
       ib                                   = self % perm(self % lline(il)+1)
       r                                    = r + 1
       val                                  = a % get_val(i,i)       
       self % tri % va(1:ndof,1:ndof,2,r)   = val(1:ndof,1:ndof)
       val                                  = a % get_val(i,ib)       
       self % tri % va(1:ndof,1:ndof,3,r)   = val(1:ndof,1:ndof)
       !
       ! Inside a linelet
       !
       ia = i       
       do k = self % lline(il)+1,self % lline(il+1)-2
          i                                    = self % perm(k)
          ib                                   = self % perm(k+1)
          r                                    = r + 1
          val                                  = a % get_val(i,ia)
          self % tri % va(1:ndof,1:ndof,1,r)   = val(1:ndof,1:ndof)
          val                                  = a % get_val(i,i)
          self % tri % va(1:ndof,1:ndof,2,r)   = val(1:ndof,1:ndof)
          val                                  = a % get_val(i,ib)
          self % tri % va(1:ndof,1:ndof,3,r)   = val(1:ndof,1:ndof)             
          ia                                   = i
       end do
       !
       ! End of linelet
       !
       if( self % lline(il+1) - self % lline(il) > 1 ) then
          i                                    = self % perm(self % lline(il+1)-1)
          ia                                   = self % perm(self % lline(il+1)-2)
          r                                    = r + 1
          val                                  = a % get_val(i,i)       
          self % tri % va(1:ndof,1:ndof,2,r)   = val(1:ndof,1:ndof)
          val                                  = a % get_val(i,ia)       
          self % tri % va(1:ndof,1:ndof,1,r)   = val(1:ndof,1:ndof)
       end if
    end do
    !
    ! Prescribe diagonal on interface nodes 
    !
    do il = 1,self % nline
       do k = self % lline(il),self % lline(il+1)-1
          i               = self % perm(k)
          self % inlin(i) = .true.
          if( i > self % nint ) then
             self % tri % va(1:ndof,1:ndof,1,i) = 0.0_rp
             self % tri % va(1:ndof,1:ndof,2,i) = 0.0_rp 
             self % tri % va(1:ndof,1:ndof,3,i) = 0.0_rp
             do idof = 1,ndof
                self % tri % va(idof,idof,2,i) = maths_inverse(self % dia % va(idof,i),0.0_rp)
             end do 
          end if
       end do
    end do

  end subroutine setup
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Solve
  !> @details Solve x = A'^-1 b ... unsymmetric case
  !>                x = L L^t b ... symmetric case
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)
    
    class(linelet),           intent(inout) :: self       !< Solver
    real(rp),        pointer, intent(in)    :: b(:)       !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)       !< Solve x
    class(mat),               intent(in)    :: a          !< Matrix A
    integer(ip)                             :: i,k,idof
    integer(ip)                             :: ndof

    if( self % n > 0 ) then
       !
       ! Interface nodes
       !
       ndof = self % ndof
       do i = 1,self % n
          if( .not. self % inlin(i) ) then
             do idof = 1,self % ndof
                x((i-1)*self % ndof+idof) =  b((i-1)*self % ndof+idof) * self % dia % vA(idof,i)
             end do
          end if
       end do
       !
       ! Interior nodes
       !
       if( self % nline > 0 ) then
          do k = 1,self % lline(self % nline+1)-1
             i = self % perm(k)
             do idof = 1,ndof
                self % bb(idof,k) = b((i-1)*ndof+idof)
             end do
          end do
          call self % tri % solve(self % bb,self % xx)
          do k = 1,self % lline(self % nline+1)-1
             i = self % perm(k)
             do idof = 1,ndof
                x((i-1)*ndof+idof) = self % xx(idof,k)
             end do
          end do
       end if
    end if
    
  end subroutine solve
  
end module def_linelet
!> @}

