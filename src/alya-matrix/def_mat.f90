!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   Matrix
!> @details Matrix type. All formats are block formats
!-----------------------------------------------------------------------

module def_mat

  use def_kintyp_basic, only : ip,rp,lg
  implicit none

  integer(ip),  parameter :: CSR_FORMAT   = 1
  integer(ip),  parameter :: COO_FORMAT   = 2
  integer(ip),  parameter :: ELL_FORMAT   = 3
  integer(ip),  parameter :: DEN_FORMAT   = 4
  integer(ip),  parameter :: BLK_FORMAT   = 5
  integer(ip),  parameter :: DIA_FORMAT   = 6
  integer(ip),  parameter :: TRI_FORMAT   = 7
  integer(ip),  parameter :: SKY_FORMAT   = 8
  integer(ip),  parameter :: BND_FORMAT   = 9
  integer(ip),  parameter :: OMP_CHUNK    = 1000
  integer(ip),  parameter :: OMP_STATIC   = 2
  integer(ip),  parameter :: OMP_GUIDED   = 3
  integer(ip),  parameter :: OMP_DYNAMIC  = 4

  private

  type, abstract :: mat     
     integer(ip)                           :: nrows
     integer(ip)                           :: ncols
     integer(ip)                           :: ndof1
     integer(ip)                           :: ndof2
     integer(ip)                           :: kfl_format
   contains
     procedure(init),       pass, deferred :: init                  ! Initialize
     procedure(alloca),     pass, deferred :: alloca                ! Allocate
     procedure(deallo),     pass, deferred :: deallo                ! Deallocate
     procedure(assign),     pass, deferred :: assign                ! Initialize
     procedure(diag),       pass, deferred :: diag                  ! Compute diagonal
     procedure(norm),       pass, deferred :: norm                  ! Matrix norm
     procedure(dirichlet),  pass, deferred :: dirichlet             ! Impose a Dirichlet condition
     procedure(get_val),    pass, deferred :: get_val               ! Get a value (i,j)
     procedure(mv_row),     pass, deferred :: mv_row                ! MV product on on single row
     procedure(mv_lower),   pass, deferred :: mv_lower              ! (L+D) product
     procedure(mv_upper),   pass, deferred :: mv_upper              ! (U+D) product
     procedure(mv_11),      pass, deferred :: mv_11                 ! MV product on (:)
     procedure(mv_12),      pass, deferred :: mv_12                 ! MV product on (:)
     procedure(mv_22),      pass, deferred :: mv_22                 ! MV product on (:,:)
     generic                               :: mv =>        &        ! Sparse matrix-vector product    
          &                                   mv_11,       &
          &                                   mv_12,       &
          &                                   mv_22
     procedure(residual_1), pass, deferred :: residual_1            ! Residual p:)
     procedure(residual_2), pass, deferred :: residual_2            ! Residual p:,:)
     generic                               :: residual =>  &        ! Residual
          &                                   residual_1,  &
          &                                   residual_2
     procedure,   pass                     :: init_mat
  end type mat

  abstract interface

     subroutine init(self)
       import                    :: mat
       class(mat), intent(inout) :: self
     end subroutine init

     subroutine alloca(self,param,MEMORY_COUNTER)
       import                                           :: mat
       import                                           :: ip  
       class(mat),                        intent(inout) :: self
       integer(ip),    optional,          intent(in)    :: param(:)
       integer(8),     optional,          intent(inout) :: MEMORY_COUNTER(2)
     end subroutine alloca

     subroutine deallo(self,MEMORY_COUNTER)
       import                                           :: mat
       import                                           :: ip  
       class(mat),                        intent(inout) :: self
       integer(8),     optional,          intent(inout) :: MEMORY_COUNTER(2)
     end subroutine deallo

     subroutine assign(self,a)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       class(mat),                        intent(inout) :: self
       real(rp),     pointer, contiguous, intent(in)    :: a(:) 
     end subroutine assign

     subroutine diag(self,dia,diagonal,row,val)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       class(mat),                        intent(in)    :: self
       class(*),       optional,          intent(inout) :: dia
       real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
       integer(ip),    optional,          intent(in)    :: row
       real(rp),       optional,          intent(out)   :: val(:)
     end subroutine diag

     subroutine mv_11(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
       real(rp),                 pointer, intent(inout) :: yy(:)                   !< Output vector
       integer(ip),    optional,          intent(in)    :: n1                      !< Starting node
       integer(ip),    optional,          intent(in)    :: n2                      !< Final node
       logical(lg),    optional,          intent(in)    :: INITIALIZATION          !< If array should be initialized
       logical(lg),    optional,          intent(in)    :: OPENMP                  !< If OpenMP should be used
       integer(ip),    optional,          intent(in)    :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
       integer(ip),    optional,          intent(in)    :: OPENMP_SCHEDULE         !< OpenMP scheduling
       logical(lg),    optional,          intent(in)    :: TRANSPOSE               !< Transpose product
    end subroutine mv_11

     subroutine mv_12(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
       real(rp),                 pointer, intent(inout) :: yy(:,:)                 !< Output vector
       integer(ip),    optional,          intent(in)    :: n1                      !< Starting node
       integer(ip),    optional,          intent(in)    :: n2                      !< Final node
       logical(lg),    optional,          intent(in)    :: INITIALIZATION          !< If array should be initialized
       logical(lg),    optional,          intent(in)    :: OPENMP                  !< If OpenMP should be used
       integer(ip),    optional,          intent(in)    :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
       integer(ip),    optional,          intent(in)    :: OPENMP_SCHEDULE         !< OpenMP scheduling
       logical(lg),    optional,          intent(in)    :: TRANSPOSE               !< Transpose product
    end subroutine mv_12

     subroutine mv_22(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:,:)                 !< Input vector
       real(rp),                 pointer, intent(inout) :: yy(:,:)                 !< Output vector
       integer(ip),    optional,          intent(in)    :: n1                      !< Starting node
       integer(ip),    optional,          intent(in)    :: n2                      !< Final node
       logical(lg),    optional,          intent(in)    :: INITIALIZATION          !< If array should be initialized
       logical(lg),    optional,          intent(in)    :: OPENMP                  !< If OpenMP should be used
       integer(ip),    optional,          intent(in)    :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
       integer(ip),    optional,          intent(in)    :: OPENMP_SCHEDULE         !< OpenMP scheduling
       logical(lg),    optional,          intent(in)    :: TRANSPOSE               !< Transpose product
     end subroutine mv_22

     real(rp) pure function mv_row(self,xx,n,ndof) 
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
       integer(ip),                       intent(in)    :: n                       !< Node
       integer(ip),                       intent(in)    :: ndof                    !< dof
     end function mv_row

     subroutine residual_1(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
       real(rp),                 pointer, intent(inout) :: yy(:)                   !< Output vector
       real(rp),                 pointer, intent(in)    :: bb(:)                   !< RHS
       integer(ip),    optional,          intent(in)    :: n1                      !< Starting node
       integer(ip),    optional,          intent(in)    :: n2                      !< Final node
       logical(lg),    optional,          intent(in)    :: INITIALIZATION          !< If array should be initialized
       logical(lg),    optional,          intent(in)    :: OPENMP                  !< If OpenMP should be used
       integer(ip),    optional,          intent(in)    :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
       integer(ip),    optional,          intent(in)    :: OPENMP_SCHEDULE         !< OpenMP scheduling
     end subroutine residual_1

     subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),                        intent(in)    :: self
       real(rp),                 pointer, intent(in)    :: xx(:,:)                 !< Input vector
       real(rp),                 pointer, intent(inout) :: yy(:,:)                 !< Output vector
       real(rp),                 pointer, intent(in)    :: bb(:,:)                 !< RHS
       integer(ip),    optional,          intent(in)    :: n1                      !< Starting node
       integer(ip),    optional,          intent(in)    :: n2                      !< Final node
       logical(lg),    optional,          intent(in)    :: INITIALIZATION          !< If array should be initialized
       logical(lg),    optional,          intent(in)    :: OPENMP                  !< If OpenMP should be used
       integer(ip),    optional,          intent(in)    :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
       integer(ip),    optional,          intent(in)    :: OPENMP_SCHEDULE         !< OpenMP scheduling
     end subroutine residual_2

     function norm(self,wnorm) result(anorm)
       import                                           :: mat
       import                                           :: rp            
       class(mat),                        intent(in)    :: self
       character(1),                      intent(in)    :: wnorm
       real(rp)                                         :: anorm       
     end function norm

     pure subroutine mv_lower(self,xx,yy,n1,n2,INITIALIZATION)       
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),     intent(in)                       :: self
       real(rp),       intent(in),    pointer           :: xx(:)                   !< Input vector
       real(rp),       intent(inout), pointer           :: yy(:)                   !< Output vector
       integer(ip),    intent(in),            optional  :: n1                      !< Starting node
       integer(ip),    intent(in),            optional  :: n2                      !< Final node
       logical(lg),    intent(in),            optional  :: INITIALIZATION          !< If array should be initialized
     end subroutine mv_lower
     
     pure subroutine mv_upper(self,xx,yy,n1,n2,INITIALIZATION)       
       import                                           :: mat
       import                                           :: ip     
       import                                           :: rp     
       import                                           :: lg     
       class(mat),     intent(in)                       :: self
       real(rp),       intent(in),    pointer           :: xx(:)                   !< Input vector
       real(rp),       intent(inout), pointer           :: yy(:)                   !< Output vector
       integer(ip),    intent(in),            optional  :: n1                      !< Starting node
       integer(ip),    intent(in),            optional  :: n2                      !< Final node
       logical(lg),    intent(in),            optional  :: INITIALIZATION          !< If array should be initialized
     end subroutine mv_upper

     pure function get_val(self,i,j) result(a)  
       import                                           :: mat
       import                                           :: ip          
       import                                           :: rp          
       class(mat),     intent(in)                       :: self
       integer(ip),    intent(in)                       :: i
       integer(ip),    intent(in)                       :: j
       real(rp)                                         :: a(self % ndof2,self % ndof1)
     end function get_val     

     pure subroutine dirichlet(self,fixno,bvess,rhs,n1,n2)
       import                                           :: mat
       import                                           :: ip             
       import                                           :: rp             
       class(mat),                        intent(inout) :: self
       integer(ip),    optional, pointer, intent(in)    :: fixno(:)
       real(rp),       optional, pointer, intent(in)    :: bvess(:)
       real(rp),       optional, pointer, intent(inout) :: rhs(:)
       integer(ip),    optional,          intent(in)    :: n1            
       integer(ip),    optional,          intent(in)    :: n2
     end subroutine dirichlet
     
  end interface

  public :: mat

  public :: COO_FORMAT  
  public :: CSR_FORMAT  
  public :: ELL_FORMAT  
  public :: DEN_FORMAT
  public :: BLK_FORMAT
  public :: DIA_FORMAT
  public :: TRI_FORMAT
  public :: SKY_FORMAT
  public :: BND_FORMAT

  public :: OMP_CHUNK  
  public :: OMP_STATIC 
  public :: OMP_GUIDED 
  public :: OMP_DYNAMIC 

contains

  subroutine init_mat(self)

    class(mat), intent(inout) :: self

    self % kfl_format = -1
    self % nrows      =  0
    self % ncols      =  0 

  end subroutine init_mat

end module def_mat
!> @}
