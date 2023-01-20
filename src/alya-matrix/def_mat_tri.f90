!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_tri.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   TRI Matrix
!> @details TRI Matrix class
!>          SELF % NDOF1
!>          SELF % NDOF2
!>          SELF % NROWS
!>          SELF % VA(SELF % NDOF1,SELF % NDOF2,3,SELF % NROWS)
!-----------------------------------------------------------------------

module def_mat_tri

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat,               only : TRI_FORMAT
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
 
  type, extends(mat) :: mat_tri
     real(rp),    contiguous, pointer :: vA(:,:,:,:)
   contains
     procedure,               pass    :: init           ! Initialize the class
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix
     procedure,               pass    :: diag           ! Compute diagonal
     procedure,               pass    :: solve          ! Solve system
     procedure,               pass    :: norm           ! Norm of a matrix 
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
  end type mat_tri

  character(11), parameter :: vacal = 'def_mat_tri'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_tri
  
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

    class(mat_tri), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = TRI_FORMAT
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
    
    class(mat_tri),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: ndof1_loc
    integer(ip)                              :: ndof2_loc
    integer(ip)                              :: nrows_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nrows_loc = optional_argument(self % nrows,param,1_ip)
    ndof1_loc = optional_argument(self % ndof1,param,2_ip)
    ndof2_loc = optional_argument(self % ndof2,param,3_ip)
    
    call memory_alloca(memor_loc,'SPMAT % VA',vacal,self % vA,ndof1_loc,ndof2_loc,3_ip,nrows_loc)        
    self % nrows = nrows_loc       
    self % ndof1 = ndof1_loc
    self % ndof2 = ndof2_loc
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
    
    class(mat_tri),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SPMAT % VA',vacal,self % vA)
    
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
    class(mat_tri),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)
  
    self % vA(1:self % ndof1,1:self % ndof2,1:3,1:self % nrows) => a
    
  end subroutine assign
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine diag(self,dia,diagonal,row,val)
    
    class(mat_tri),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,nrows,ndofn,idofn

    nrows = self % nrows
    ndofn = self % ndof2
    
    if( present(diagonal) ) then
       
       do ii = 1,nrows
          do idofn = 1,ndofn
             diagonal(idofn,ii) = self % vA(idofn,idofn,2,ii)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia )
          
          if( .not. associated(dia % vA) ) call dia % alloca((/ndofn,nrows/)) 
          do ii = 1,nrows
             do idofn = 1,ndofn
                dia % vA(idofn,ii) = self % vA(idofn,idofn,2,ii)
             end do
          end do
          
       end select
       
    else if( present(row) .and. present(val) ) then

       ii = row
       do idofn = 1,ndofn
          val(idofn) = self % vA(idofn,idofn,2,ii)
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
    
    class(mat_tri), intent(in)                      :: self
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
    integer(ip)                                     :: ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof2
       ndofc = self % ndof1
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       do ii = 1,self % nrows * ndofr
          yy(ii) = bb(ii) - yy(ii)
       end do
    end if
    
  end subroutine residual_1
  
  subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_tri), intent(in)                      :: self
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
    integer(ip)                                     :: ndofc
    
    if( associated(xx) .and. associated(yy) ) then
       ndofr = self % ndof2
       ndofc = self % ndof1
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
  !> O
  !-----------------------------------------------------------------------
  
  subroutine mv_11(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_tri), intent(in)                      :: self
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
       ndofc = self % ndof1
       ndofr = self % ndof2
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if
    
  end subroutine mv_11
  
  subroutine mv_12(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_tri), intent(in)                      :: self
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
       ndofc = self % ndof1
       ndofr = self % ndof2
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)       
    end if
    
  end subroutine mv_12
  
  subroutine mv_22(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_tri), intent(in)                      :: self
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
       ndofc = self % ndof1
       ndofr = self % ndof2
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

    class(mat_tri), intent(in)                    :: self
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
    integer(ip)                                   :: ii,idofc
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
    if( nn1 == 1 ) then
       do idofc = 1,ndofc
          yy(1:ndofr,1) = yy(1:ndofr,1) &
               + self % vA(idofc,1:ndofr,2,1) * xx(idofc,1  ) &
               + self % vA(idofc,1:ndofr,3,1) * xx(idofc,2  )
       end do
    else
       nn1 = nn1 - 1
    end if
    if( nn2 == self % nrows ) then
       do idofc = 1,ndofc
          yy(1:ndofr,nn2) = yy(1:ndofr,nn2) &
               + self % vA(idofc,1:ndofr,1,nn2) * xx(idofc,nn2-1) &
               + self % vA(idofc,1:ndofr,2,nn2) * xx(idofc,nn2  )
       end do
    else
       nn2 = nn2 - 1
    end if

    if( use_openmp ) then
       !$OMP PARALLEL   DO                                 &
       !$OMP SCHEDULE ( STATIC )                           &
       !$OMP DEFAULT  ( NONE )                             &
       !$OMP SHARED   ( self, nn1, nn2, xx, yy, idofc     ) &
       !$OMP PRIVATE  ( ii, ndofr, ndofc )
       do ii = nn1+1,nn2-1
          do idofc = 1,ndofc
             yy(1:ndofr,ii) = yy(1:ndofr,ii) &
                  + self % vA(idofc,1:ndofr,1,ii) * xx(idofc,ii-1) &
                  + self % vA(idofc,1:ndofr,2,ii) * xx(idofc,ii  ) &
                  + self % vA(idofc,1:ndofr,3,ii) * xx(idofc,ii+1) 
          end do
       end do
       !$OMP END PARALLEL DO
    else
       do ii = nn1+1,nn2-1
          do idofc = 1,ndofc
             yy(1:ndofr,ii) = yy(1:ndofr,ii) &
                  + self % vA(idofc,1:ndofr,1,ii) * xx(idofc,ii-1) &
                  + self % vA(idofc,1:ndofr,2,ii) * xx(idofc,ii  ) &
                  + self % vA(idofc,1:ndofr,3,ii) * xx(idofc,ii+1) 
          end do
       end do
    end if

  end subroutine mv_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Solve A x= d
  !> @details Solve tridiagonal system
  !>          https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  !>          For i = 2 To N
  !>             W = A(i) / B(i - 1)
  !>             B(i) = B(i) - W * C(i - 1)
  !>             D(i) = D(i) - W * D(i - 1)
  !>         Next i
  !>         X(N) = D(N) / B(N)
  !>         For i = N - 1 To 1 Step -1
  !>             X(i) = (D(i) - C(i) * X(i + 1)) / B(i)
  !>         Next i
  !> 
  !-----------------------------------------------------------------------
    
  subroutine solve(self,dd,xx)
    class(mat_tri),          intent(inout)   :: self
    real(rp),       pointer, intent(in)      :: dd(:,:)
    real(rp),       pointer, intent(inout)   :: xx(:,:)
    integer(ip)                              :: ii,idofn,nn,ndofn
    real(rp),       allocatable              :: ww(:)
    real(rp),       allocatable              :: d2(:,:)
    real(rp),       allocatable              :: bb(:,:)

    if( self % nrows > 0 ) then
       nn    = self % nrows
       ndofn = self % ndof2
       allocate(d2(ndofn,nn))
       allocate(bb(ndofn,nn))
       allocate(ww(ndofn))
       do ii = 1,nn
          do idofn = 1,ndofn
             bb(idofn,ii) = self % vA(idofn,idofn,2,ii)
             d2(idofn,ii) = dd(idofn,ii)
          end do
       end do       
       do ii = 2,nn
          do idofn = 1,ndofn
             ww(idofn)    = self % vA(idofn,idofn,1,ii) / bb(idofn,ii-1)
             bb(idofn,ii) = bb(idofn,ii) - ww(idofn) * self % vA(idofn,idofn,3,ii-1)
             d2(idofn,ii) = d2(idofn,ii) - ww(idofn) * d2(idofn,ii-1)
          end do
       end do
       do idofn = 1,ndofn
          xx(idofn,nn) = d2(idofn,nn) / bb(idofn,nn)
       end do
       do ii = nn-1,1,-1
          do idofn = 1,ndofn
             xx(idofn,ii) = (d2(idofn,ii) - self % vA(idofn,idofn,3,ii) * xx(idofn,ii+1)) / bb(idofn,ii)
          end do
       end do
       deallocate(d2,bb)
    end if

  end subroutine solve

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Matrix norm
  !> @details Norm of a matrix
  !>
  !----------------------------------------------------
      
  function norm(self,wnorm) result(anorm)

    class(mat_tri), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    real(rp)                      :: dummr
    integer(ip)                   :: ii
    integer(ip)                   :: idof1,idof2

    anorm = 0.0_rp
    if( self % nrows > 0 ) then

       if( wnorm == 'i' .or. wnorm == 'I' ) then

          do ii = 1,self % nrows
             do idof2 = 1,self % ndof2
                dummr = 0.0_rp
                do idof1 = 1,self % ndof1
                   dummr = dummr + abs(self % vA(idof1,idof2,1,ii))
                   dummr = dummr + abs(self % vA(idof1,idof2,2,ii))
                   dummr = dummr + abs(self % vA(idof1,idof2,3,ii))
                end do
                anorm = max(anorm,dummr)
             end do
          end do

       else if( wnorm == '1' ) then
          
          do idof2 = 1,self % ndof2
             dummr = 0.0_rp
             do idof1 = 1,self % ndof1
                dummr = dummr + abs(self % vA(idof1,idof2,1,2))
                dummr = dummr + abs(self % vA(idof1,idof2,2,1))
             end do
             anorm = max(anorm,dummr)          
             dummr = 0.0_rp
             do idof1 = 1,self % ndof1
                dummr = dummr + abs(self % vA(idof1,idof2,2,self % nrows))
                dummr = dummr + abs(self % vA(idof1,idof2,3,self % nrows-1))
             end do
             anorm = max(anorm,dummr)
          end do
          
          do ii = 2,self % nrows-1
             do idof2 = 1,self % ndof2
                dummr = 0.0_rp
                do idof1 = 1,self % ndof1
                   dummr = dummr + abs(self % vA(idof1,idof2,1,ii+1))
                   dummr = dummr + abs(self % vA(idof1,idof2,2,ii))
                   dummr = dummr + abs(self % vA(idof1,idof2,3,ii-1))
                end do
                anorm = max(anorm,dummr)
             end do
          end do

       else if( wnorm == 'f' .or. wnorm == 'F' ) then

          do ii = 1,self % nrows
             do idof2 = 1,self % ndof2
                do idof1 = 1,self % ndof1
                   anorm = anorm + self % vA(idof1,idof2,1,ii)*self % vA(idof1,idof2,1,ii)
                   anorm = anorm + self % vA(idof1,idof2,2,ii)*self % vA(idof1,idof2,2,ii)
                   anorm = anorm + self % vA(idof1,idof2,3,ii)*self % vA(idof1,idof2,3,ii)
                end do
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
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Mv on a single row
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_tri),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: ll

    if( n == 1 ) then
       do ll = 1,self % ndof2
          yy =   self % vA(ll,ndof,2,1) * xx((1-1)*self % ndof1+ll) &
               + self % vA(ll,ndof,3,1) * xx((2-1)*self % ndof1+ll)
       end do
    else if( n == self % nrows ) then
       do ll = 1,self % ndof2
          yy =   self % vA(ll,ndof,1,n) * xx((n-2)*self % ndof1+ll) &
               + self % vA(ll,ndof,2,n) * xx((n-1)*self % ndof1+ll)
       end do
    else
       do ll = 1,self % ndof2
          yy =   self % vA(ll,ndof,1,n) * xx((n-2)*self % ndof1+ll) &
               + self % vA(ll,ndof,2,n) * xx((n-1)*self % ndof1+ll) &
               + self % vA(ll,ndof,3,n) * xx((n-0)*self % ndof1+ll) 
       end do    
    end if

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
    
    class(mat_tri), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,jt1
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    if( nn1 == 1 ) then
       ii = nn1
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = 1,id
             jt  = (ii-1) * self % ndof1 + jd
             jt1 = (ii-2) * self % ndof1 + jd
             yy(it) = yy(it) &
                  + self % vA(jd,id,2,ii) * xx(jt)
          end do
       end do       
    else
       nn1 = nn1 - 1
    end if
    
    if( nn2 == self % nrows ) then
       ii = nn2
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = 1,id
             jt  = (ii-1) * self % ndof1 + jd
             jt1 = (ii-2) * self % ndof1 + jd
             yy(it) = yy(it) &
                  + self % vA(jd,id,1,ii) * xx(jt1) &
                  + self % vA(jd,id,2,ii) * xx(jt)
          end do
       end do       
    else
       nn2 = nn2 - 1
    end if
    
    do ii = nn1+1,nn2-1
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = 1,id
             jt  = (ii-1) * self % ndof1 + jd
             jt1 = (ii-2) * self % ndof1 + jd
             yy(it) = yy(it) &
               + self % vA(jd,id,1,ii) * xx(jt1) &
               + self % vA(jd,id,2,ii) * xx(jt)
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
    
    class(mat_tri), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize
   
    integer(ip)                                     :: ii,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,jt1
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    if( nn1 == 1 ) then
       ii = nn1
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = id,self % ndof2
             jt  = (ii-1) * self % ndof1 + jd
             jt1 = (ii  ) * self % ndof1 + jd
             yy(it) = yy(it) &
                  + self % vA(jd,id,2,ii) * xx(jt)  &
                  + self % vA(jd,id,3,ii) * xx(jt1)
          end do
       end do       
    else
       nn1 = nn1 - 1
    end if
    
    if( nn2 == self % nrows ) then
       ii = nn2
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = id,self % ndof2
             jt  = (ii-1) * self % ndof1 + jd
             yy(it) = yy(it) &
                  + self % vA(jd,id,2,ii) * xx(jt)
          end do
       end do       
    else
       nn2 = nn2 - 1
    end if
    
    do ii = nn1+1,nn2-1
       do id = 1,self % ndof1
          it = (ii-1) * self % ndof1 + id
          do jd = id,self % ndof2
             jt  = (ii-1) * self % ndof1 + jd
             jt1 = (ii  ) * self % ndof1 + jd
             yy(it) = yy(it) &
               + self % vA(jd,id,2,ii) * xx(jt)  &
               + self % vA(jd,id,3,ii) * xx(jt1)
          end do
       end do
    end do
    
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
     
    class(mat_tri), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)

    if( i <= self % nrows ) then
       if(       j == i-1 ) then
          a(:,:) = self % vA(:,:,1,i)
       else if ( j == i   ) then
          a(:,:) = self % vA(:,:,2,i)
       else if ( j == i+1 ) then
          a(:,:) = self % vA(:,:,3,i)
       else
          a(:,:) = 0.0_rp
       end if
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
    
    class(mat_tri),                    intent(inout) :: self
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

end module def_mat_tri
!> @}
