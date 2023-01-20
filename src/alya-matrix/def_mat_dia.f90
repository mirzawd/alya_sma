!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_dia.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   CSR Matrix
!> @details CSR Matrix class
!-----------------------------------------------------------------------

module def_mat_dia

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat,               only : DIA_FORMAT
  use def_mat,               only : OMP_CHUNK  
  use def_mat,               only : OMP_STATIC 
  use def_mat,               only : OMP_GUIDED 
  use def_mat,               only : OMP_DYNAMIC 
  use mod_optional_argument, only : optional_argument
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end

  implicit none

  private

  type, extends(mat) :: mat_dia
     real(rp),    contiguous, pointer :: vA(:,:)
   contains
     procedure,               pass    :: init           ! Initialize the class
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix 
     procedure,               pass    :: diag           ! Diagonal
     procedure,               pass    :: norm           ! Norm of a matrix 
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: solve          ! Solve a system
     procedure,               pass    :: solve_inv      ! Solve inverse system
     procedure,               pass    :: inverse        ! Inverse diagonal
     procedure,               pass    :: scale          ! Scale
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: mv_row         ! Multiply one row
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
  end type mat_dia

  character(11), parameter :: vacal = 'def_mat_dia'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)

  public :: mat_dia

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

    class(mat_dia), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = DIA_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    nullify(self % vA)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate: parameters:
  !>          PARAM(1:2) = NDOF1, NROWS
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self,param,MEMORY_COUNTER)

    class(mat_dia),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: ndof1_loc,nrows_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    ndof1_loc = optional_argument(self % ndof1,param,1_ip)
    nrows_loc = optional_argument(self % nrows,param,2_ip)

    call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,nrows_loc)        
    self % ndof1 = ndof1_loc
    self % ndof2 = ndof1_loc
    self % nrows = nrows_loc       

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

    class(mat_dia),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    call memory_deallo(memor_loc,'SELF % VA',vacal,self % vA)

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Assign an array matrix to mat class vA
  !> @details Assign a rank-1 array to mat class vA
  !> 
  !-----------------------------------------------------------------------

  subroutine assign(self,a)
    class(mat_dia),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)

    self % vA(1:self % ndof1,1:self % nrows) => a

  end subroutine assign

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   y = b - Ax
  !> @details Compute residual y = b - Ax
  !> 
  !-----------------------------------------------------------------------

  subroutine residual_1(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)

    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    integer(ip)                                     :: ndofr,ii

    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       do ii = 1,self % nrows * ndofr
          yy(ii) = bb(ii) - yy(ii)
       end do
    end if

  end subroutine residual_1

  subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)

    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    real(rp),       intent(in),    pointer          :: bb(:,:)                 !< RHS
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    integer(ip)                                     :: ndofr,ii

    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)       
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

    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofn

    if( associated(xx) .and. associated(yy) ) then
       ndofn = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if

  end subroutine mv_11

  subroutine mv_12(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)

    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofn

    if( associated(xx) .and. associated(yy) ) then
       ndofn = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if

  end subroutine mv_12

  subroutine mv_22(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)

    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:,:)                 !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:,:)                 !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),            optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),            optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),            optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),            optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                     :: ndofn

    if( associated(xx) .and. associated(yy) ) then
       ndofn = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
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

  pure subroutine mv_go(self,xx,yy,n1,n2,ndofn,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)

    class(mat_dia), intent(in)                    :: self
    integer(ip),    intent(in)                    :: ndofn                   !< Number of dofs
    real(rp),       intent(in)                    :: xx(ndofn,*)             !< Input vector
    real(rp),       intent(inout)                 :: yy(ndofn,*)             !< Output vector
    integer(ip),    intent(in),          optional :: n1                      !< Starting node
    integer(ip),    intent(in),          optional :: n2                      !< Final node
    logical(lg),    intent(in),          optional :: INITIALIZATION          !< If array should be initialized
    logical(lg),    intent(in),          optional :: OPENMP                  !< If OpenMP should be used
    integer(ip),    intent(in),          optional :: OPENMP_CHUNK            !< Chunks for dynamic scheduling of OpenMP
    integer(ip),    intent(in),          optional :: OPENMP_SCHEDULE         !< OpenMP scheduling
    logical(lg),    intent(in),          optional :: TRANSPOSE               !< Transpose product
    integer(ip)                                   :: ii,jj
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

    if( do_initialize ) then
       do ii = nn1,nn2
          do jj = 1,self % ndof1
             yy(jj,ii) = self % vA(jj,ii) * xx(jj,ii)
          end do
       end do
    else
       do ii = nn1,nn2
          do jj = 1,self % ndof1
             yy(jj,ii) = yy(jj,ii) + self % vA(jj,ii) * xx(jj,ii)
          end do
       end do
    end if

  end subroutine mv_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal
  !> @details Compute the diagonal matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diag(self,dia,diagonal,row,val)

    class(mat_dia),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,idofn

    if( present(diagonal) ) then

       do ii = 1,self % nrows
          do idofn = 1,self % ndof1
             diagonal(idofn,ii) = self % vA(idofn,ii)
          end do
       end do

    else if( present(dia) ) then

       select type ( dia )
       class is ( mat_dia)
          do ii = 1,self % nrows
             do idofn = 1,self % ndof1
                dia % vA(idofn,ii) = self % vA(idofn,ii)
             end do
          end do
       end select

    else if( present(row) .and. present(val) ) then

       ii = row
       do idofn = 1,self % ndof1
          val(idofn) = self % vA(idofn,ii)
       end do

    end if

  end subroutine diag

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Matrix norm
  !> @details Norm of a matrix
  !>
  !----------------------------------------------------

  function norm(self,wnorm) result(anorm)

    class(mat_dia), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    integer(ip)                   :: ii,idof1

    anorm = 0.0_rp

    if( self % nrows > 0 ) then

       if( wnorm == 'i' .or. wnorm == 'I' ) then

          do ii = 1,self % nrows
             do idof1 = 1,self % ndof1
                anorm = max(anorm,abs(self % vA(idof1,ii)))
             end do
          end do

       else if( wnorm == '1' ) then

          do ii = 1,self % nrows
             do idof1 = 1,self % ndof1
                anorm = max(anorm,abs(self % vA(idof1,ii)))
             end do
          end do

       else if( wnorm == 'f' .or. wnorm == 'F' ) then

          do ii = 1,self % nrows
             do idof1 = 1,self % ndof1
                anorm = anorm + self % vA(idof1,ii)*self % vA(idof1,ii)
             end do
          end do
          anorm = sqrt(anorm)

       else
          call runend('MOD_MATHS: THIS NORM IS NOT CODED')
       end if

    end if

  end function norm

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-17
  !> @brief   Solve Ax=b
  !> @details Solve diagonal system
  !> 
  !-----------------------------------------------------------------------

  pure subroutine solve(self,rhs,x,istat)

    class(mat_dia),            intent(in)     :: self
    real(rp),      target,     intent(inout)  :: rhs(:)
    real(rp),                  intent(inout)  :: x(:)
    integer(ip),   optional,   intent(out)    :: istat
    integer(ip)                               :: i,j,id

    if( present(istat) ) then
       istat = 0
       do i = 1,self % nrows
          do id = 1,self % ndof1
             if( abs(self % vA(id,i)) == 0.0_rp ) then
                istat = 1
                return
             end if
          end do
       end do       
    end if
    
    select case ( self % ndof1 )

    case ( 1_ip )

       do i = 1,self % nrows
          x(i) = rhs(i) / self % vA(1,i)
       end do

    case default

       j = 0
       do i = 1,self % nrows
          do id = 1,self % ndof1
             j    = j + 1
             x(j) = rhs(j) / self % vA(id,i)
          end do
       end do

    end select

  end subroutine solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-17
  !> @brief   Solve x=Ab
  !> @details Solve inverse diagonal system
  !> 
  !-----------------------------------------------------------------------

  pure subroutine solve_inv(self,rhs,x,istat)

    class(mat_dia),            intent(in)     :: self
    real(rp),      target,     intent(inout)  :: rhs(:)
    real(rp),                  intent(inout)  :: x(:)
    integer(ip),   optional,   intent(out)    :: istat
    integer(ip)                               :: i,j,id

    select case ( self % ndof1 )

    case ( 1_ip )

       do i = 1,self % nrows
          x(i) = rhs(i) * self % vA(1,i)
       end do

    case default

       j = 0
       do i = 1,self % nrows
          do id = 1,self % ndof1
             j    = j + 1
             x(j) = rhs(j) * self % vA(id,i)
          end do
       end do

    end select

  end subroutine solve_inv

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Inverse diagonal matrix 
  !> @details Inverse diagonal matrix 
  !>
  !----------------------------------------------------

  subroutine inverse(self) 

    class(mat_dia), intent(inout) :: self
    integer(ip)                   :: ii,idof1

    do ii = 1,self % nrows
       do idof1 = 1,self % ndof1
          if( self % vA(idof1,ii) /= 0.0_rp ) then
             self % vA(idof1,ii) = 1.0_rp / self % vA(idof1,ii)
          else
             self % vA(idof1,ii) = 0.0_rp
          end if
       end do
    end do

  end subroutine inverse

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Scale
  !> @details Scale a vector
  !> 
  !-----------------------------------------------------------------------

  pure subroutine scale(self,xx,row,dof)

    class(mat_dia),                    intent(in)    :: self
    real(rp),       optional,          intent(out)   :: xx(:)
    integer(ip),    optional,          intent(in)    :: row
    integer(ip),    optional,          intent(in)    :: dof
    integer(ip)                                      :: ii,kk,idof1

    if( present(row) ) then
       if( present(dof) ) then
          kk     = (row-1) * self % ndof1 + dof
          xx(kk) = xx(kk)  * self % vA(dof,row)
       else
          do idof1 = 1,self % ndof1
             kk     = (row-1) * self % ndof1 + idof1
             xx(kk) = xx(kk)  * self % vA(idof1,row)
          end do
       end if
    else
       kk = 0
       do ii = 1,self % nrows
          do idof1 = 1,self % ndof1
             kk     = (ii-1) * self % ndof1+idof1
             xx(kk) = xx(kk) * self % vA(idof1,ii)
          end do
       end do
    end if
    
  end subroutine scale

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Scale
  !> @details Scale a vector
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_dia),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: ii
    
    ii = (n-1)  * self % ndof1 + ndof
    yy = xx(ii) * self % vA(ndof,n)
    
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
    
    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized

    call mv_go(self,xx,yy,n1,n2,self % ndof1,INITIALIZATION)
    
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
    
    class(mat_dia), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    
    call mv_go(self,xx,yy,n1,n2,self % ndof1,INITIALIZATION)
    
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
     
    class(mat_dia), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    integer(ip)                :: idof

    a = 0.0_rp
    if( i == j ) then
       do idof = 1,self % ndof1
          a(idof,idof) = self % vA(idof,i)
       end do
    end if
    
  end function get_val

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Nothing to do!
  !> 
  !-----------------------------------------------------------------------

  pure subroutine dirichlet(self,fixno,bvess,rhs,n1,n2)
    
    class(mat_dia),                    intent(inout) :: self
    integer(ip),    optional, pointer, intent(in)    :: fixno(:)
    real(rp),       optional, pointer, intent(in)    :: bvess(:)
    real(rp),       optional, pointer, intent(inout) :: rhs(:)
    integer(ip),    optional,          intent(in)    :: n1            
    integer(ip),    optional,          intent(in)    :: n2
    
  end subroutine dirichlet

end module def_mat_dia
!> @}
