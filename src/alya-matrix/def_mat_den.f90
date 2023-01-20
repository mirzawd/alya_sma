!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_demat.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   Sparse matrix
!> @details Sparse matrix type
!-----------------------------------------------------------------------

module def_mat_den

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_mat,               only : DEN_FORMAT
  use def_mat,               only : OMP_CHUNK  
  use def_mat,               only : OMP_STATIC 
  use def_mat,               only : OMP_GUIDED 
  use def_mat,               only : OMP_DYNAMIC
  use mod_optional_argument, only : optional_argument
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end
  use mod_iofile_basic,      only : iofile_available_unit 
  use mod_strings,           only : upper_case
  use mod_strings,           only : add_extension
  implicit none
  private
  
  !----------------------------------------------------------------------
  !
  ! Sparse matrix 
  !
  !----------------------------------------------------------------------
  
  type, extends(mat) :: mat_den
     real(rp),    contiguous, pointer :: vA(:,:,:,:)
  contains
     procedure,               pass    :: init           ! Initialization           
     procedure,               pass    :: alloca         ! Allocate 
     procedure,               pass    :: deallo         ! Deallocate 
     procedure,               pass    :: assign         ! Assign 
     procedure,               pass    :: diag           ! Diagonal 
     procedure,               pass    :: norm           ! Norm of a matrix 
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: solve          ! Solve a system
     procedure,               pass    :: output         ! Output 
     procedure,               pass    :: csr2den        ! CSR format to dense format 
     procedure,               pass    :: symmetry       ! Check symmetry
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: mv_row         ! Single row MV product
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
  end type mat_den

  character(11), parameter :: vacal = 'def_mat_den'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_den
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Assign an array matrix to mat class vA
  !> @details Assign a rank-1 array to mat class vA
  !> 
  !-----------------------------------------------------------------------
    
  subroutine assign(self,a)
    class(mat_den),                       intent(inout) :: self
    real(rp),        pointer, contiguous, intent(in)    :: a(:)
       
    self % vA(1:self % ndof1,1:self % ndof2,1:self % nrows,1:self % ncols) => a
    
  end subroutine assign

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------
    
  subroutine init(self)
    class(mat_den), intent(inout) :: self

    self % kfl_format = DEN_FORMAT
    self % nrows      = 0
    self % ncols      = 0
    self % ndof1      = 1
    self % ndof2      = 1
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
    
    class(mat_den),            intent(inout)  :: self
    integer(ip),    optional,  intent(in)     :: param(:)
    integer(8),     optional,  intent(inout)  :: MEMORY_COUNTER(2)
    integer(8)                                :: memor_loc(2)
    integer(ip)                               :: nrows_loc,ncols_loc
    integer(ip)                               :: ndof1_loc,ndof2_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    ndof1_loc = optional_argument(self % ndof1,param,1_ip)
    ndof2_loc = optional_argument(self % ndof2,param,2_ip)
    nrows_loc = optional_argument(self % nrows,param,3_ip)
    ncols_loc = optional_argument(self % ncols,param,4_ip)
    
    if( nrows_loc * ncols_loc > 0 ) then
       call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,ndof2_loc,nrows_loc,ncols_loc)        
       self % ndof1 = ndof1_loc   
       self % ndof2 = ndof2_loc   
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
    
    class(mat_den),          intent(inout)  :: self
    integer(8),   optional,  intent(inout)  :: MEMORY_COUNTER(2)
    integer(8)                              :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % VA',vacal,self % vA)
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal
  !> @details Compute the diagonal matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diag(self,dia,diagonal,row,val)
    class(mat_den),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,idofn

    if( present(diagonal) ) then
       
       do ii = 1,self % nrows
          do idofn = 1,self % ndof2
             diagonal(idofn,ii) = self % vA(idofn,idofn,ii,ii)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia )
          do ii = 1,self % nrows
             do idofn = 1,self % ndof2
                dia % vA(idofn,ii) = self % vA(idofn,idofn,ii,ii)
             end do
          end do
       end select
       
    else if( present(row) .and. present(val) ) then
       
       ii = row
       do ii = 1,self % nrows
          do idofn = 1,self % ndof2
             val(idofn) = self % vA(idofn,idofn,ii,ii)
          end do
       end do
       
    end if
    
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
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    integer(ip)                                     :: ndofr,ii,ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofc = self % ndof1
       ndofr = self % ndof2
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       do ii = 1,self % nrows * ndofr
          yy(ii) = bb(ii) - yy(ii)
       end do
    end if
    
  end subroutine residual_1
  
  subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:,:)                 !< RHS
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    integer(ip)                                     :: ndofr,ii,ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofc = self % ndof1
       ndofr = self % ndof2
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)       
       do ii = 1,self % nrows
          yy(1:ndofr,ii) = bb(1:ndofr,ii) - yy(1:ndofr,ii)
       end do
    end if
    
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
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofr
    integer(ip)                                     :: ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof2
       ndofc = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if
    
  end subroutine mv_11
  
  subroutine mv_12(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofr
    integer(ip)                                     :: ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof2
       ndofc = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if
    
  end subroutine mv_12
  
  subroutine mv_22(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofr
    integer(ip)                                     :: ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof2
       ndofc = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if
    
  end subroutine mv_22
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   Parallel mv
  !> @details Parallel mv
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)

    class(mat_den), intent(in)                    :: self
    integer(ip),    intent(in)                    :: ndofr                   !< Number of dofs
    integer(ip),    intent(in)                    :: ndofc                   !< Number of dofs
    real(rp),       intent(in)                    :: xx(ndofc,*)             !< Input vector
    real(rp),       intent(inout)                 :: yy(ndofr,*)             !< Output vector
    integer(ip),    intent(in),          optional :: n1                      !< Starting node
    integer(ip),    intent(in),          optional :: n2                      !< Final node
    logical(lg),    intent(in),          optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),          optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),          optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),          optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),          optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                   :: ii,jj
    integer(ip)                                   :: idofc,idofr
    logical(lg)                                   :: use_openmp
    logical(lg)                                   :: do_initialize
    integer(ip)                                   :: my_schedule
    integer(ip)                                   :: my_chunk,nn1,nn2
    !
    ! Optional arguments
    !
    use_openmp    = optional_argument(.false.,OPENMP)
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    my_chunk      = optional_argument(OMP_CHUNK,OPENMP_CHUNK)
    my_schedule   = optional_argument(OMP_STATIC,OPENMP_SCHEDULE)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    !
    ! Initialize solution
    !
    if( do_initialize ) then
       do ii = nn1,nn2
          yy(:,ii) = 0.0_rp
       end do
    end if

    if( optional_argument(.false.,TRANSPOSE) ) then

       do jj = 1,self % ncols
          do ii = nn1,nn2
             do idofc = 1,ndofc
                do idofr = 1,ndofr
                   yy(idofc,jj) = yy(idofc,jj) + self % vA(idofr,idofc,jj,ii) * xx(idofr,ii)
                end do
             end do
          end do
       end do

    else

       if( use_openmp ) then
          !$OMP PARALLEL   DO                                      &
          !$OMP SCHEDULE ( STATIC )                                &
          !$OMP DEFAULT  ( NONE )                                  &
          !$OMP SHARED   ( self, nn1, nn2, xx, yy, ndofc, ndofr  ) &
          !$OMP PRIVATE  ( ii, jj, idofc )
          do jj = 1,self % ncols
             do ii = nn1,nn2
                do idofc = 1,ndofc
                   yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(idofc,1:ndofr,ii,jj) * xx(idofc,jj)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          do jj = 1,self % ncols
             do ii = nn1,nn2
                do idofc = 1,ndofc
                   yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(idofc,1:ndofr,ii,jj) * xx(idofc,jj)
                end do
             end do
          end do
       end if

    end if

  end subroutine mv_go

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Matrix norm
  !> @details Norm of a matrix
  !>
  !----------------------------------------------------
      
  function norm(self,wnorm) result(anorm)

    class(mat_den), intent(in)  :: self
    character(1),   intent(in)  :: wnorm
    real(rp)                    :: anorm
    real(rp)                    :: dummr
    integer(ip)                 :: ii,jj
    integer(ip)                 :: idof1,idof2
    
    anorm = 0.0_rp
    if( wnorm == 'i' .or. wnorm == 'I' ) then
       
       do ii = 1,self % ncols
          do idof2 = 1,self % ndof2
             dummr = 0.0_rp
             do jj = 1,self % nrows
                do idof1 = 1,self % ndof1
                   dummr = dummr + abs(self % vA(idof1,idof2,ii,jj))
                end do
             end do
             anorm = max(anorm,dummr)
          end do
       end do
       
    else if( wnorm == '1' ) then
       
       do jj = 1,self % ncols
          do idof1 = 1,self % ndof1
             dummr = 0.0_rp
             do ii = 1,self % nrows
                do idof2 = 1,self % ndof2
                   dummr = dummr + abs(self % vA(idof1,idof2,ii,jj))
                end do
             end do
             anorm = max(anorm,dummr)
          end do
       end do

    else if( wnorm == 'f' .or. wnorm == 'F' ) then

       do jj = 1,self % ncols
          do ii = 1,self % nrows
             do idof2 = 1,self % ndof2
                do idof1 = 1,self % ndof1
                   anorm = anorm + self % vA(idof1,idof2,ii,jj)*self % vA(idof1,idof2,ii,jj)
                end do
             end do
          end do
       end do
       anorm = sqrt(anorm)
       
    else
       call runend('MOD_MATHS: THIS NORM IS NOT CODED')
    end if
    
  end function norm

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-17
  !> @brief   Solve Ax=b
  !> @details Solve linear equations by gaussian
  !>          elimination with partial pivoting. Taken from:
  !>          https://www.osti.gov/servlets/purl/6636181
  !>          https://en.wikipedia.org/wiki/Gaussian_elimination#:~:text=In%20mathematics%2C%20Gaussian%20elimination%2C%20also,the%20corresponding%20matrix%20of%20coefficients.
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solve(self,rhs,x,istat,PIVOTING,TRANSPOSE)

    class(mat_den),            intent(in)     :: self
    real(rp),      target,     intent(inout)  :: rhs(:)
    real(rp),                  intent(inout)  :: x(:)
    integer(ip),   optional,   intent(out)    :: istat
    logical(lg),   optional                   :: PIVOTING
    logical(lg),   optional                   :: TRANSPOSE
    real(rp),      contiguous, pointer        :: a(:,:)
    real(rp),      contiguous, pointer        :: b(:)
    integer(ip)                               :: n,m
    integer(ip)                               :: i,k,j,h,pivot
    real(rp)                                  :: f
    real(rp),      allocatable                :: row(:)

    n = self % nrows
    m = self % ncols

    if( present(istat) ) istat=0

    select case ( n )

    case ( 1_ip )

       if( self % vA(1,1,1,1) < epsil ) then
          if( present(istat) ) istat=1
       else
          x(1) = rhs(1) / self % vA(1,1,1,1)
       end if

    case default

       allocate(a,source = self % vA(1,1,:,:))
       allocate(b,source = rhs(1:n))

       if( optional_argument(.true.,PIVOTING) ) then

          !-----------------------------------------------------------------
          !
          ! With partial pivoting
          !
          !-----------------------------------------------------------------

          allocate(row(m))
          
          do k = 1,n-1
             pivot = maxloc(abs(a(k:n,k)),dim=1)+k-1
             if( pivot /= k ) then
                row(:)     = a(k,:)
                a(k,:)     = a(pivot,:)
                a(pivot,:) = row(:)
                f          = b(k)
                b(k)       = b(pivot)
                b(pivot)   = f
             end if
             !
             ! Do for all rows below pivot
             !
             do i = k+1,n
                f       = a(i,k)  / a(k,k)
                a(i,k:) = a(i,k:) - f * a(k,k:)
                b(i)    = b(i)    - f * b(k)
             end do
          end do
          !
          ! Back substitution
          !
          x(n) = b(n) / a(n,n)
          do i = n-1,1,-1
             x(i) = ( b(i) - dot_product(a(i,i+1:n), x(i+1:n)) ) / a(i,i)
          end do

          deallocate(row)

       else

          !-----------------------------------------------------------------
          !
          ! Without pivoting
          !
          !-----------------------------------------------------------------
          !
          ! Elimination
          !
          do k = 1,n-1
             do i = k+1,n
                f       = a(i,k)  / a(k,k)
                a(i,k:) = a(i,k:) - f * a(k,k:)
                b(i)    = b(i)    - f * b(k)
             end do
          end do
          !
          ! Back substitution
          !
          x(n) = b(n) / a(n,n)
          do i = n-1,1,-1
             x(i) = ( b(i) - dot_product(a(i,i+1:n), x(i+1:n)) ) / a(i,i)
          end do

       end if

       deallocate(a,b)

    end select

  end subroutine solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Mv on a single row
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_den),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: ll,jj

    yy = 0.0_rp
    do jj = 1,self % ncols
       do ll = 1,self % ndof2
          yy = yy + self % vA(ll,ndof,n,jj) * xx((jj-1)*self % ndof2+ll)
       end do
    end do
    
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
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if

    do jj = 1,self % ncols
       do jd = 1,self % ndof2
          jt = (jj-1) * self % ndof2+jd 
          do ii = nn1,nn2
             if( jj <= ii ) then
                do id = 1,self % ndof1
                   if( jd <= id ) then
                      it     = (ii-1) * self % ndof1 + id
                      yy(it) = yy(it) + self % vA(id,jd,ii,jj) * xx(jt)
                   end if
                end do
             end if
          end do
       end do
    end do

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
    
    class(mat_den), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    do jj = 1,self % ncols
       do jd = 1,self % ndof2
          jt = (jj-1) * self % ndof2+jd
          do ii = nn1,nn2
             if( jj >= ii ) then
                do id = 1,self % ndof1
                   if( jd >= id ) then
                      it     = (ii-1) * self % ndof1 + id
                      yy(it) = yy(it) + self % vA(id,jd,ii,jj) * xx(jt)
                   end if
                end do
             end if
          end do
       end do
    end do
    
  end subroutine mv_upper

    !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-08
  !> @brief   Output
  !> @details Output matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine output(self,FMT,FILENAME,PERM)

    class(mat_den),                   intent(in) :: self
    character(*),                     intent(in) :: FMT
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: ii,jj
    integer(ip)                                  :: idof1,idof2,ndof1,ndof2
    integer(ip)                                  :: ndofn,nrows,ncols
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc
    !
    ! Dimensions
    !
    ndof1 = self % ndof1
    ndof2 = self % ndof2
    ndofn = self % ndof1 * self % ndof2
    nrows = self % nrows
    ncols = self % ncols
    !
    ! Unit and files
    !
    unit4        = iofile_available_unit(90_ip) ! TODO: possible int8 to int4 conversion

    select case ( upper_case(FMT) )

    case ( 'MATRIX MARKET' , 'MTX' , 'MM' )


    case ( 'PAJEK NET' )

       !-----------------------------------------------------------------
       !
       ! Matrix market
       !
       !-----------------------------------------------------------------

    case ( 'DENSE' )

       !-----------------------------------------------------------------
       !
       ! Dense
       !
       !-----------------------------------------------------------------

       filename_loc = optional_argument('matrix-dense.txt',FILENAME)
       call add_extension(filename_loc,'txt')
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       do ii = 1,nrows
          do idof2 = 1,ndof2
             do jj = 1,ncols
                do idof1 = 1,ndof1
                   write(unit4,'(1x,e12.6,a)',advance='no') self % vA(idof1,idof2,ii,jj),' ,'
                end do
             end do
             write(unit4,*)
          end do
       end do
       close(unit4)

    case ( 'GID' )

       !-----------------------------------------------------------------
       !
       ! GID
       !
       !-----------------------------------------------------------------

    end select

  end subroutine output

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-08
  !> @brief   CSR to dense
  !> @details CSR to dense
  !>           
  !-----------------------------------------------------------------------

  subroutine csr2den(self,csr) 

    class(mat_den),  intent(inout) :: self
    type(mat_csr),   intent(in)    :: csr
    integer(ip)                    :: ndof1,ndof2,nrows,ncols
    integer(ip)                    :: ii,idof2,jj,kk,ll,idof1
    !
    ! Dimensions
    !
    ndof1        = csr % ndof1
    ndof2        = csr % ndof2
    nrows        = csr % nrows
    ncols        = csr % get_ncols()
    self % nrows = nrows
    self % ncols = ncols
    self % ndof1 = ndof1
    self % ndof2 = ndof2

    call self % alloca()
    
    do ii = 1,nrows
       do idof2 = 1,ndof2
          do jj = 1,ncols
             kk = 0
             ll = csr % ia(ii)-1
             do while( ll < csr % ia(ii+1)-1 .and. kk == 0 )
                ll = ll + 1
                if( csr % ja(ll) == jj ) kk = ll
             end do
             do idof1 = 1,ndof1
                if( kk /= 0 ) &
                     self % Va(idof1,idof2,ii,jj) = csr % vA(idof1,idof2,kk)
             end do
          end do
       end do
    end do
       
  end subroutine csr2den

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-08
  !> @brief   Check symmetry
  !> @details Check symmetry. Returns a negative value if symmetric
  !>          TOLERANCE > 0: A unsymmetric if abs(A_ij-A_ji) > TOLERANCE
  !>          TOLERANCE < 0: A unsymmetric if abs(A_ij-A_ji) > TOLERANCE * Amax
  !>                         with Amax = max_ij( abs(A_ij) )
  !>           
  !-----------------------------------------------------------------------

  function symmetry(self,TOLERANCE) result(sym)

    class(mat_den),           intent(in) :: self
    real(rp),       optional, intent(in) :: TOLERANCE
    integer(ip)                          :: ii,jj
    integer(ip)                          :: idof1,idof2
    real(rp)                             :: toler,vamax,sym

    sym = -1.0_rp

    if( associated(self % vA) ) then

       if( self % ndof1 /= self % ndof2 ) then
          sym = 1.0_rp
       else         
          toler = optional_argument(1.0e-12_rp,TOLERANCE)
          vamax = max(epsil,maxval(abs(self % vA)))
          if( toler < 0.0_rp ) toler = abs(toler) * vamax
          do ii = 1,self % nrows
             do jj = ii+1,self % nrows
                do idof2 = 1,self % ndof2
                   do idof1 = idof2+1,self % ndof1
                      if( abs(self % vA(idof1,idof2,ii,jj)-self % vA(idof2,idof1,jj,ii)) > toler ) then
                         sym = max(sym,abs(self % vA(idof1,idof2,ii,jj)-self % vA(idof2,idof1,jj,ii)) / vamax)
                      end if
                   end do
                end do
             end do
          end do
       end if

    end if

  end function symmetry

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Get value
  !> @details Get value at position i,j
  !> 
  !-----------------------------------------------------------------------

  pure function get_val(self,i,j) result(a)
     
    class(mat_den), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)

    if( i <= self % nrows .and. j <= self % ncols ) then
       a(:,:) = self % vA(:,:,i,j)
    else
       a = 0.0_rp
    end if
    
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
    
    class(mat_den),                    intent(inout) :: self
    integer(ip),    optional, pointer, intent(in)    :: fixno(:)
    real(rp),       optional, pointer, intent(in)    :: bvess(:)
    real(rp),       optional, pointer, intent(inout) :: rhs(:)
    integer(ip),    optional,          intent(in)    :: n1            
    integer(ip),    optional,          intent(in)    :: n2
    integer(ip)                                      :: ii,idofn,jdofn,nn1,nn2
    integer(ip)                                      :: jj
    integer(ip)                                      :: ndofc,ndofr
    logical(lg)                                      :: if_fixno
    
    nn1      = optional_argument(1_ip,n1)
    nn2      = optional_argument(self % nrows,n2)
    ndofr    = self % ndof2
    ndofc    = self % ndof1
    if_fixno = .true.
    
    do ii = nn1,nn2
       do idofn = 1,ndofr
          if( present(fixno) ) then
             if( fixno((ii-1)*ndofr+idofn)  > 0 ) then
                if_fixno = .true.
             else
                if_fixno = .false.
             end if
          end if
          if( if_fixno ) then
             do jj = 1,self % nrows
                do jdofn = 1,ndofc
                   if( jj /= ii .and. idofn /= jdofn ) then
                      self % vA(jdofn,idofn,ii,jj) = 0.0_rp
                   end if
                   self % vA(idofn,jdofn,jj,ii) = 0.0_rp  
                end do
             end do
          end if
       end do
    end do
    
  end subroutine dirichlet
  
end module def_mat_den
!> @}
