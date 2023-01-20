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
!> @brief   BND Matrix
!> @details BND Matrix class for band matrices
!-----------------------------------------------------------------------

module def_mat_bnd

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_den,           only : mat_den
  use def_mat,               only : BND_FORMAT
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
 
  type, extends(mat) :: mat_bnd
     logical(lg)                      :: symme          ! Symmetry
     integer(ip)                      :: bandw          ! Bandwidth
     integer(ip)                      :: id             ! Diagonal place
     real(rp),    contiguous, pointer :: vA(:,:)        ! Values
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
     procedure,               pass    :: factorize      ! Compute L or L and U
     procedure,               pass    :: output         ! Output matrix
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
  end type mat_bnd

  character(11), parameter :: vacal = 'def_mat_bnd'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_bnd
  
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

    class(mat_bnd), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = BND_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    self % bandw      = 0
    self % symme      = .false.
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
    
    class(mat_bnd),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: ndof1_loc
    integer(ip)                              :: ndof2_loc
    integer(ip)                              :: nrows_loc
    integer(ip)                              :: bandw_loc
    integer(ip)                              :: symmi_loc
    logical(lg)                              :: symme_loc
    
    symmi_loc = log2int(self % symme)
    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nrows_loc = optional_argument(self % nrows,param,1_ip)
    ndof1_loc = optional_argument(self % ndof1,param,2_ip)
    ndof2_loc = optional_argument(self % ndof2,param,3_ip)
    bandw_loc = optional_argument(self % bandw,param,4_ip)
    symmi_loc = optional_argument(symmi_loc   ,param,5_ip)
        
    self % nrows = nrows_loc       
    self % ndof1 = ndof1_loc
    self % ndof2 = ndof2_loc
    self % bandw = bandw_loc
    self % symme = int2log(symmi_loc)
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    self % id    = (self % bandw-1)/2+1
    
    if( self % symme ) then
       call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,self % id,nrows_loc*ndof1_loc*ndof2_loc)
    else
       call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,bandw_loc,nrows_loc*ndof1_loc*ndof2_loc)
    end if

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
    
    class(mat_bnd),            intent(inout) :: self
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
    class(mat_bnd),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)
  
    self % vA(1:self % bandw,1:self % ndof1 * self % ndof2 * self % nrows) => a
    
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
    
    class(mat_bnd),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,nrows,ndofn,idofn,id,bandw

    nrows = self % nrows
    ndofn = self % ndof2
    bandw = self % bandw
    id    = (bandw-1)/2+1
    
    if( present(diagonal) ) then
       
       do ii = 1,nrows
          do idofn = 1,ndofn
             diagonal(idofn,ii) = self % vA(id,ii)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia )
          
          if( .not. associated(dia % vA) ) call dia % alloca((/ndofn,nrows/)) 
          do ii = 1,nrows
             do idofn = 1,ndofn
                dia % vA(idofn,ii) = self % vA(id,ii)
             end do
          end do
          
       end select
       
    else if( present(row) .and. present(val) ) then

       ii = row
       do idofn = 1,ndofn
          val(idofn) = self % vA(id,ii)
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
    
    class(mat_bnd), intent(in)                      :: self
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
    
    class(mat_bnd), intent(in)                      :: self
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
    
    class(mat_bnd), intent(in)                      :: self
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
    
    class(mat_bnd), intent(in)                      :: self
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
    
    class(mat_bnd), intent(in)                      :: self
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

    class(mat_bnd), intent(in)                    :: self
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
    integer(ip)                                   :: ii,idofc,jj,kk,id
    integer(ip)                                   :: bandw,jmin,jmax
    logical(lg)                                   :: use_openmp
    logical(lg)                                   :: do_initialize
    integer(ip)                                   :: my_schedule
    integer(ip)                                   :: my_chunk,nn1,nn2,nn
    !
    ! Optional arguments
    !
    use_openmp    = optional_argument(.false.,OPENMP)
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    my_chunk      = optional_argument(OMP_CHUNK,OPENMP_CHUNK)
    my_schedule   = optional_argument(OMP_STATIC,OPENMP_SCHEDULE)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    nn            = self % nrows
    id            = (self % bandw-1)/2+1
    if( self % symme ) then
       bandw = id
    else
       bandw = self % bandw
    end if
    !
    ! Initialize solution
    !
    if( do_initialize ) then
       do ii = nn1,nn2
          yy(:,ii) = 0.0_rp
       end do
    end if

    if( ndofr == 1 .and. ndofc == 1 ) then
       do ii = nn1,nn2
          jmin = max( 1_ip, 1+id-ii)
          jmax = min(bandw,nn+id-ii)
          do jj = jmin,jmax
             kk = ii-id+jj
             yy(1,ii) = yy(1,ii) + self % vA(jj,ii) * xx(1,kk)
          end do
       end do
       if( self % symme ) then
          do ii = max(1_ip,nn1-id+1),min(nn,nn2+id-1)
             jmin = max( 1_ip, 1+id-ii)
             jmax = min(bandw,nn+id-ii)
             do jj = jmin,jmax
                kk = ii-id+jj
                if( ii /= kk .and. kk >= nn1 .and. kk <= nn2 ) &
                     yy(1,kk) = yy(1,kk) + self % vA(jj,ii) * xx(1,ii)
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
      
  subroutine factorize(self,l)

    class(mat_bnd), intent(in)    :: self
    class(mat),     intent(inout) :: l
    integer(ip)                   :: i,j,k,n,id,bandw
    integer(ip)                   :: ki,ji,jk,ik,ij
    integer(ip)                   :: maxi,maxk,mink,mini
    real(rp)                      :: s


    n         = self % nrows
    bandw     = self % bandw
    id        = (bandw-1)/2+1
    if( self % symme ) bandw = id
    
    call l % init()
    l % nrows = n
    l % ncols = n

    select type ( l )
    class is ( mat_bnd )
       l % symme = self % symme
       l % bandw = self % bandw
    end select

    call l % alloca()

    select type ( l )
    class is ( mat_den )

       !-----------------------------------------------------------------
       !
       ! L is a dense matrix
       !
       !-----------------------------------------------------------------

       do k = 1,n
          s = dot_product(l % va(1,1,k,1:k-1),l % va(1,1,k,1:k-1))
          l % va(1,1,k,k) = sqrt(self % va(id,k)-s)       
          mink = max(k+1,k+1-id+1)
          maxk = min(n,k+id-1)
          do j = mink,maxk 
             jk = id +k-j
             l % va(1,1,j,k) = self % va(jk,j)
             mini = max(1_ip,k-id+1,j-id+1)
             maxi = min(k-1,j+id-1,k+id-1)
             l % va(1,1,j,k) = ( l % va(1,1,j,k) &
                  - dot_product(l % va(1,1,j,mini:maxi),l % va(1,1,k,mini:maxi)) ) / l % va(1,1,k,k)
          end do
       end do
       
    class is ( mat_bnd )
       
       !-----------------------------------------------------------------
       !
       ! L is a band matrix
       !
       ! do k = 1,n
       !    Lkk = ( Akk -sum_i=1^k-1 Lki^2 )^1/2
       !    do j = k+1,n
       !       Ljk = ( Ajk - sum_i*1^k-1 LjiLki ) / Lkk
       !    end do
       ! end do
       !
       !-----------------------------------------------------------------

       do k = 1,n
          s            = dot_product(l % va(1:id-1,k),l % va(1:id-1,k))
          l % va(id,k) = sqrt(self % va(id,k)-s)
          mink = max(k+1,k+1-id+1)
          maxk = min(n,k+id-1)
          do j = mink,maxk            
             jk   = id + k-j
             mini = max(1_ip,k-id+1,j-id+1)
             maxi = min(k-1,j+id-1,k+id-1)
             s    = 0.0_rp
             do i = mini,maxi
                ji = id + i-j
                ki = id + i-k
                s  = s + l % va(ji,j) * l % va(ki,k)
             end do
             l % va(jk,j) = ( self % va(jk,j) -s ) / l % va(id,k)          
          end do
       end do
       
    end select

  end subroutine factorize

  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Solve
  !> @details Solve system
  !>
  !----------------------------------------------------
      
  subroutine solve(self,l,b,x)

    class(mat_bnd),          intent(in)    :: self
    class(mat),              intent(in)    :: l
    real(rp),       pointer, intent(in)    :: b(:,:)
    real(rp),       pointer, intent(inout) :: x(:,:)
    integer(ip)                            :: i,j,k,n,bandw
    integer(ip)                            :: id,kt,it,ij
    integer(ip)                            :: jmin,jmax
    real(rp),       pointer                :: y(:,:)
    real(rp)                               :: s

    n     = self % nrows
    bandw = self % bandw
    id    = (bandw-1)/2+1

    allocate(y(1,n))
    !
    ! L y = b
    ! L^t x = y
    !
    select type ( l )
    class is ( mat_den )
       
       !-----------------------------------------------------------------
       !
       ! L is a dense matrix
       !
       !-----------------------------------------------------------------
       
       y(1,1) = b(1,1) / l % va(1,1,1,1)
       do i = 2,n
          s      = dot_product(l % va(1,1,i,1:i-1),y(1,1:i-1))
          y(1,i) = (b(1,i)-s)/l % va(1,1,i,i)
       end do
       x(1,n) = y(1,n) / l % va(1,1,n,n)
       do i = n-1,1,-1
          s      = dot_product(l % va(1,1,i+1:n,i),x(1,i+1:n))
          x(1,i) = (y(1,i)-s)/l % va(1,1,i,i)
       end do
       
    class is ( mat_bnd )
       
       !-----------------------------------------------------------------
       !
       ! L is a band matrix
       !
       !-----------------------------------------------------------------
    
       y(1,1) = b(1,1) / l % va(id,1)
       do i = 2,n
          s = 0.0_rp
          jmin = max(1_ip,1-i+id)
          jmax = min(id-1,n-i+id)
          do ij = jmin,jmax
             j = i - id + ij
             s = s + l % va(ij,i) * y(1,j)
          end do
          y(1,i) = (b(1,i)-s)/l % va(id,i)
       end do
       x(1,n) = y(1,n) / l % va(id,n)
       do i = n-1,1,-1
          s = 0.0_rp
          do j = max(i+1,id+i-bandw),min(id+i-1,n)
             k = id + i-j
             s = s + l % va(k,j) * x(1,j)
          end do
          x(1,i) = (y(1,i)-s)/l % va(id,i)
       end do
       
    end select

    deallocate(y)

  end subroutine solve

  pure function inband(ij,bandw)
    integer(ip), intent(in) :: ij
    integer(ip), intent(in) :: bandw
    logical(lg)             :: inband

    if( ij >= 1 .and. ij <= bandw ) then
       inband = .true.
    else
       inband = .false.
    end if
    
  end function inband
  
  !----------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2019
  !> @brief   Matrix norm
  !> @details Norm of a matrix
  !>
  !----------------------------------------------------
      
  function norm(self,wnorm) result(anorm)

    class(mat_bnd), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    real(rp)                      :: dummr
    integer(ip)                   :: ii
    integer(ip)                   :: idof1,idof2

!!$    anorm = 0.0_rp
!!$    if( self % nrows > 0 ) then
!!$
!!$       if( wnorm == 'i' .or. wnorm == 'I' ) then
!!$
!!$          do ii = 1,self % nrows
!!$             do idof2 = 1,self % ndof2
!!$                dummr = 0.0_rp
!!$                do idof1 = 1,self % ndof1
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,1,ii))
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,2,ii))
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,3,ii))
!!$                end do
!!$                anorm = max(anorm,dummr)
!!$             end do
!!$          end do
!!$
!!$       else if( wnorm == '1' ) then
!!$          
!!$          do idof2 = 1,self % ndof2
!!$             dummr = 0.0_rp
!!$             do idof1 = 1,self % ndof1
!!$                dummr = dummr + abs(self % vA(idof1,idof2,1,2))
!!$                dummr = dummr + abs(self % vA(idof1,idof2,2,1))
!!$             end do
!!$             anorm = max(anorm,dummr)          
!!$             dummr = 0.0_rp
!!$             do idof1 = 1,self % ndof1
!!$                dummr = dummr + abs(self % vA(idof1,idof2,2,self % nrows))
!!$                dummr = dummr + abs(self % vA(idof1,idof2,3,self % nrows-1))
!!$             end do
!!$             anorm = max(anorm,dummr)
!!$          end do
!!$          
!!$          do ii = 2,self % nrows-1
!!$             do idof2 = 1,self % ndof2
!!$                dummr = 0.0_rp
!!$                do idof1 = 1,self % ndof1
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,1,ii+1))
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,2,ii))
!!$                   dummr = dummr + abs(self % vA(idof1,idof2,3,ii-1))
!!$                end do
!!$                anorm = max(anorm,dummr)
!!$             end do
!!$          end do
!!$
!!$       else if( wnorm == 'f' .or. wnorm == 'F' ) then
!!$
!!$          do ii = 1,self % nrows
!!$             do idof2 = 1,self % ndof2
!!$                do idof1 = 1,self % ndof1
!!$                   anorm = anorm + self % vA(idof1,idof2,1,ii)*self % vA(idof1,idof2,1,ii)
!!$                   anorm = anorm + self % vA(idof1,idof2,2,ii)*self % vA(idof1,idof2,2,ii)
!!$                   anorm = anorm + self % vA(idof1,idof2,3,ii)*self % vA(idof1,idof2,3,ii)
!!$                end do
!!$             end do
!!$          end do
!!$          anorm = sqrt(anorm)
!!$
!!$       else
!!$          call runend('MOD_MATHS: THIS NORM IS NOT CODED')
!!$       end if
!!$
!!$    end if

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
    
    class(mat_bnd),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: ll

!!$    if( n == 1 ) then
!!$       do ll = 1,self % ndof2
!!$          yy =   self % vA(ll,ndof,2,1) * xx((1-1)*self % ndof1+ll) &
!!$               + self % vA(ll,ndof,3,1) * xx((2-1)*self % ndof1+ll)
!!$       end do
!!$    else if( n == self % nrows ) then
!!$       do ll = 1,self % ndof2
!!$          yy =   self % vA(ll,ndof,1,n) * xx((n-2)*self % ndof1+ll) &
!!$               + self % vA(ll,ndof,2,n) * xx((n-1)*self % ndof1+ll)
!!$       end do
!!$    else
!!$       do ll = 1,self % ndof2
!!$          yy =   self % vA(ll,ndof,1,n) * xx((n-2)*self % ndof1+ll) &
!!$               + self % vA(ll,ndof,2,n) * xx((n-1)*self % ndof1+ll) &
!!$               + self % vA(ll,ndof,3,n) * xx((n-0)*self % ndof1+ll) 
!!$       end do    
!!$    end if

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
    
    class(mat_bnd), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,jt1
    
!!$    do_initialize = optional_argument(.true. ,INITIALIZATION)
!!$    nn1           = optional_argument(1_ip,n1)
!!$    nn2           = optional_argument(self % nrows,n2)
!!$    
!!$    if( do_initialize ) then
!!$       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
!!$          yy(ii) = 0.0_rp
!!$       end do
!!$    end if
!!$ 
!!$    if( nn1 == 1 ) then
!!$       ii = nn1
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = 1,id
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             jt1 = (ii-2) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$                  + self % vA(jd,id,2,ii) * xx(jt)
!!$          end do
!!$       end do       
!!$    else
!!$       nn1 = nn1 - 1
!!$    end if
!!$    
!!$    if( nn2 == self % nrows ) then
!!$       ii = nn2
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = 1,id
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             jt1 = (ii-2) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$                  + self % vA(jd,id,1,ii) * xx(jt1) &
!!$                  + self % vA(jd,id,2,ii) * xx(jt)
!!$          end do
!!$       end do       
!!$    else
!!$       nn2 = nn2 - 1
!!$    end if
!!$    
!!$    do ii = nn1+1,nn2-1
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = 1,id
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             jt1 = (ii-2) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$               + self % vA(jd,id,1,ii) * xx(jt1) &
!!$               + self % vA(jd,id,2,ii) * xx(jt)
!!$          end do
!!$       end do
!!$    end do

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
    
    class(mat_bnd), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize
   
    integer(ip)                                     :: ii,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,jt1
    
!!$    do_initialize = optional_argument(.true. ,INITIALIZATION)
!!$    nn1           = optional_argument(1_ip,n1)
!!$    nn2           = optional_argument(self % nrows,n2)
!!$    
!!$    if( do_initialize ) then
!!$       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
!!$          yy(ii) = 0.0_rp
!!$       end do
!!$    end if
!!$ 
!!$    if( nn1 == 1 ) then
!!$       ii = nn1
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = id,self % ndof2
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             jt1 = (ii  ) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$                  + self % vA(jd,id,2,ii) * xx(jt)  &
!!$                  + self % vA(jd,id,3,ii) * xx(jt1)
!!$          end do
!!$       end do       
!!$    else
!!$       nn1 = nn1 - 1
!!$    end if
!!$    
!!$    if( nn2 == self % nrows ) then
!!$       ii = nn2
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = id,self % ndof2
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$                  + self % vA(jd,id,2,ii) * xx(jt)
!!$          end do
!!$       end do       
!!$    else
!!$       nn2 = nn2 - 1
!!$    end if
!!$    
!!$    do ii = nn1+1,nn2-1
!!$       do id = 1,self % ndof1
!!$          it = (ii-1) * self % ndof1 + id
!!$          do jd = id,self % ndof2
!!$             jt  = (ii-1) * self % ndof1 + jd
!!$             jt1 = (ii  ) * self % ndof1 + jd
!!$             yy(it) = yy(it) &
!!$               + self % vA(jd,id,2,ii) * xx(jt)  &
!!$               + self % vA(jd,id,3,ii) * xx(jt1)
!!$          end do
!!$       end do
!!$    end do
    
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
     
    class(mat_bnd), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)

!!$    if( i <= self % nrows ) then
!!$       if(       j == i-1 ) then
!!$          a(:,:) = self % vA(:,:,1,i)
!!$       else if ( j == i   ) then
!!$          a(:,:) = self % vA(:,:,2,i)
!!$       else if ( j == i+1 ) then
!!$          a(:,:) = self % vA(:,:,3,i)
!!$       else
!!$          a(:,:) = 0.0_rp
!!$       end if
!!$    else
!!$       a = 0.0_rp
!!$    end if
    
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
    
    class(mat_bnd),                    intent(inout) :: self
    integer(ip),    optional, pointer, intent(in)    :: fixno(:)
    real(rp),       optional, pointer, intent(in)    :: bvess(:)
    real(rp),       optional, pointer, intent(inout) :: rhs(:)
    integer(ip),    optional,          intent(in)    :: n1            
    integer(ip),    optional,          intent(in)    :: n2
    integer(ip)                                      :: nn1,nn2
    integer(ip)                                      :: ndofc,ndofr
    real(rp),       allocatable                      :: diag(:)
    logical(lg)                                      :: if_fixno
    
!!$    nn1      = optional_argument(1_ip,n1)
!!$    nn2      = optional_argument(self % nrows,n2)
!!$    ndofr    = self % ndof2
!!$    ndofc    = self % ndof1
!!$    if_fixno = .true.
!!$    
!!$    allocate(diag(ndofr))
!!$    deallocate(diag)

  end subroutine dirichlet

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-10-18
  !> @brief   Real to logical
  !> @details Convert a real to a logical
  !> 
  !-----------------------------------------------------------------------

  function int2log(int)
    
    integer(ip), intent(in) :: int
    logical(lg)             :: int2log

    if( int <= 0 ) then
       int2log = .false.
    else
       int2log = .true.
    end if

  end function int2log
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-10-18
  !> @brief   Real to logical
  !> @details Convert a real to a logical
  !> 
  !-----------------------------------------------------------------------

  function log2int(logi)
    
    logical(lg), intent(in) :: logi
    integer(ip)             :: log2int

    if( logi ) then
       log2int = 1_ip
    else
       log2int = 0_ip
    end if

  end function log2int

    !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-08
  !> @brief   Output
  !> @details Output matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine output(self,FMT,FILENAME,PERM)

    class(mat_bnd),                   intent(in) :: self
    character(*),                     intent(in) :: FMT
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: ii,jj,kk,bandw
    integer(ip)                                  :: ndofn,nrows,ncols
    integer(ip)                                  :: ndof1,ndof2,id,ll
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc
    real(rp),      allocatable                   :: vA(:)
    !
    ! Dimensions
    !
    ndof1 = self % ndof1
    ndof2 = self % ndof2
    ndofn = self % ndof1 * self % ndof2
    bandw = self % bandw
    nrows = self % nrows
    id    = self % id
    if( nrows <= 0 ) return
    if( self % symme ) bandw = id
    !
    ! Unit and files
    !
    unit4 = iofile_available_unit()

    allocate(vA(nrows))
    select case ( upper_case(FMT) )

    case ( 'DENSE' )

       !-----------------------------------------------------------------
       !
       ! Dense
       !
       !-----------------------------------------------------------------

       filename_loc = optional_argument('matrix-band.txt',FILENAME)
       call add_extension(filename_loc,'txt')
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       do ii = 1,nrows
          vA = 0.0_rp
          do kk = 1,bandw    
             jj = kk-id+ii
             if( jj >= 1 .and. jj <= nrows ) &
                  vA(jj) = self % vA(kk,ii)
          end do
          if( self % symme ) then
             do jj = ii+1,min(nrows,ii+bandw)
                ll = ii-jj+id
                if( ll >= 1 .and. ll<= bandw ) &
                     vA(jj) = self % vA(ll,jj)
             end do
          end if
          do jj = 1,nrows-1
             write(unit4,'(1x,e12.6,a)',advance='no') vA(jj),' ,'
          end do
          write(unit4,'(1x,e12.6,a)',advance='no') vA(nrows)             
          write(unit4,*)
          
       end do
       close(unit4)

    end select
    deallocate(vA)

  end subroutine output

end module def_mat_bnd
!> @}
