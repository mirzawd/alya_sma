!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_blk.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   BLK Matrix
!> @details BLK Matrix class
!-----------------------------------------------------------------------

module def_mat_blk

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat,               only : BLK_FORMAT
  use def_mat_csr,           only : mat_csr
  use def_mat_coo,           only : mat_coo
  use def_mat_ell,           only : mat_ell
  use def_mat_den,           only : mat_den
  use def_mat_dia,           only : mat_dia
  use mod_optional_argument, only : optional_argument
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end

  implicit none

  private 
  
  type, extends(mat) :: mat_blk
     class(mat),   pointer :: vA(:,:)
   contains
     procedure,               pass    :: init 
     procedure,               pass    :: alloca 
     procedure,               pass    :: deallo
     procedure,               pass    :: assign
     procedure,               pass    :: diag
     procedure,               pass    :: assign_blk
     procedure,               pass    :: norm           ! Matrix norm
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: mv_row         ! Single row MV product
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
   end type mat_blk

  character(11), parameter :: vacal = 'def_mat_blk'
  
  public :: mat_blk

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   Initialization 
  !> @details Initialization 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init(self)

    class(mat_blk), intent(inout) :: self

    call self % init_mat()
    nullify(self % vA)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
    
  subroutine alloca(self,param,MEMORY_COUNTER)
    
    class(mat_blk),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nrows_loc
    integer(ip)                              :: ncols_loc

    nrows_loc = optional_argument(self % nrows,param,1_ip)
    ncols_loc = optional_argument(self % ncols,param,2_ip)
    
    if( nrows_loc * ncols_loc > 0 ) then
       self % nrows = nrows_loc  
       self % ncols = ncols_loc  
    end if
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
    
  subroutine deallo(self,MEMORY_COUNTER)
    
    class(mat_blk),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Assign an array matrix to mat class vA
  !> @details Assign a rank-1 array to mat class vA. Nothing to do in this 
  !>          particular case of a block matrix
  !> 
  !-----------------------------------------------------------------------
    
  subroutine assign(self,a)!,i,j)
    class(mat_blk),                      intent(inout) :: self
    real(rp),      pointer,  contiguous, intent(in)    :: a(:)
    !integer(ip),    optional,            intent(in)    :: i
    !integer(ip),    optional,            intent(in)    :: j

    !select type ( a )
    !class is ( mat_csr )
    !   !self % vA(i:i,j:j) => a 
    !class is ( mat_coo )
    !   !self % vA(i:i,j:j) => a
    !class is ( mat_den )
    !   !self % vA(i:i,j:j) => a
    !end select
    
  end subroutine assign

  subroutine assign_blk(self,a,i,j)
    class(mat_blk),                      intent(inout) :: self
    class(*),       target,  contiguous, intent(in)    :: a(:)
    integer(ip),                         intent(in)    :: i
    integer(ip),                         intent(in)    :: j

    select type ( a )
    class is ( mat_csr ) ; self % vA(i:i,j:j) => a
    class is ( mat_coo ) ; self % vA(i:i,j:j) => a
    class is ( mat_ell ) ; self % vA(i:i,j:j) => a
    class is ( mat_den ) ; self % vA(i:i,j:j) => a
    class is ( mat_dia ) ; self % vA(i:i,j:j) => a
    end select
    
  end subroutine assign_blk

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal
  !> @details Compute the diagonal matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diag(self,dia,diagonal,row,val)
    class(mat_blk),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)

!!$    if( present(diagonal) ) then
!!$       do ii = 1,self % nrows
!!$          zd = self % ia(ii) - 1
!!$          jj    = 0
!!$          do while( jj /= ii )
!!$             zd = zd + 1
!!$             jj    = self % ja(zd)
!!$          end do
!!$          do idofn = 1,self % ndof2
!!$             diagonal(idofn,ii) = self % vA(idofn,idofn,zd)
!!$          end do
!!$       end do
!!$    else if( present(dia) ) then
!!$       select type ( dia )
!!$       class is ( mat_blk )
!!$          do ii = 1,self % nrows
!!$             zd = self % ia(ii) - 1
!!$             jj    = 0
!!$             do while( jj /= ii )
!!$                zd = zd + 1
!!$                jj    = self % ja(zd)
!!$             end do
!!$             do idofn = 1,self % ndof2
!!$                dia % vA(idofn,ii) = self % vA(idofn,idofn,zd)
!!$             end do
!!$          end do
!!$       end select
!!$    end if
    
  end subroutine diag

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   y = b - Ax
  !> @details Compute residual y = b - Ax
  !> 
  !-----------------------------------------------------------------------
  
  subroutine residual_1(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    
  end subroutine residual_1
  
  subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:,:)                 !< RHS
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    
  end subroutine residual_2
      
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   Parallel mv
  !> @details Parallel mv
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mv_11(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
 
    if( associated(xx) .and. associated(yy) ) then
       !ndofn = self % ndofn
       !call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)       
    end if
    
  end subroutine mv_11
  
  subroutine mv_12(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
 
    if( associated(xx) .and. associated(yy) ) then
       !ndofn = self % ndofn
       !call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)       
    end if
    
  end subroutine mv_12
  
  subroutine mv_22(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    
    if( associated(xx) .and. associated(yy) ) then
       !ndofn = self % ndofn
       !call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)       
    end if
    
  end subroutine mv_22

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Matrix norm
  !> @details Norm of a matrix
  !>
  !----------------------------------------------------
      
  function norm(self,wnorm) result(anorm)

    class(mat_blk), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm

  end function norm

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Mv on a single row
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_blk),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof

  end function mv_row

    !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   (L+D) product
  !> @details Lower matrix vector product
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine mv_lower(self,xx,yy,n1,n2,INITIALIZATION)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,nn1,nn2
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if

  end subroutine mv_lower
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   (U+D) product
  !> @details Upper matrix vector product
  !> 
  !-----------------------------------------------------------------------

  pure subroutine mv_upper(self,xx,yy,n1,n2,INITIALIZATION)
    
    class(mat_blk), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,nn1,nn2
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
    
  end subroutine mv_upper

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Get value
  !> @details Get value at position i,j
  !> 
  !-----------------------------------------------------------------------

  pure function get_val(self,i,j) result(a)
     
    class(mat_blk), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    
  end function get_val

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine dirichlet(self,fixno,bvess,rhs,n1,n2)
    
    class(mat_blk),                    intent(inout) :: self
    integer(ip),    optional, pointer, intent(in)    :: fixno(:)
    real(rp),       optional, pointer, intent(in)    :: bvess(:)
    real(rp),       optional, pointer, intent(inout) :: rhs(:)
    integer(ip),    optional,          intent(in)    :: n1            
    integer(ip),    optional,          intent(in)    :: n2
    integer(ip)                                      :: nn1,nn2
    integer(ip)                                      :: ndofc,ndofr
    real(rp),       allocatable                      :: diag(:)
    logical(lg)                                      :: if_fixno
    
    nn1      = optional_argument(1_ip,n1)
    nn2      = optional_argument(self % nrows,n2)
    ndofr    = self % ndof2
    ndofc    = self % ndof1
    if_fixno = .true.
    
    allocate(diag(ndofr))
    deallocate(diag)

  end subroutine dirichlet
  
end module def_mat_blk
!> @}

