!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_csr.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   SKY Matrix
!> @details SKY Matrix class
!>          Required to allocate:
!>          SELF % NZ
!>          SELF % NDOF1
!>          SELF % NDOF2
!>          SELF % NROWS
!>          SELF % IA(SELF % NROWS+1)
!>          SELF % JA(SELF % NZ)
!>          SELF % VA(SELF % NDOF1,SELF % NDOF2,SELF % NZ)
!-----------------------------------------------------------------------

module def_mat_sky

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_mat,               only : SKY_FORMAT
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

  type, extends(mat) :: mat_sky
     logical(lg)                      :: symme
     integer(ip)                      :: nskyl
     integer(ip),             pointer :: iskyl(:)
     integer(ip),             pointer :: idiag(:)     
     real(rp),    contiguous, pointer :: vA(:)
   contains
     procedure,               pass    :: init           ! Initialize the class
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix           
     procedure,               pass    :: set_nrows      ! Set rows      
     procedure,               pass    :: get_nrows      ! Get rows      
     procedure,               pass    :: set_ncols      ! Set columns      
     procedure,               pass    :: get_ncols      ! Get columns
     procedure,               pass    :: bandwidth      ! Compute bandwidth and profile
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: alloca_matrix  ! Allocate matrix only
     procedure,               pass    :: deallo_matrix  ! Deallocate matrix only
     procedure,               pass    :: output         ! Output matrix
     procedure,               pass    :: diag           ! Compute diagonal
     procedure,               pass    :: diagz          ! Compute diagonal position
     procedure,               pass    :: norm           ! Matrix norm
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: symmetry       ! Determine the symmetry
     procedure,               pass    :: csr2sky        ! Convert CSR to skyline format
     procedure,               pass    :: mv_row         ! Single row MV product
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
  end type mat_sky

  character(11), parameter :: vacal = 'def_mat_sky'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_sky
  
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

    class(mat_sky), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = SKY_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    self % nskyl      = 0
    self % symme      = .true.
    nullify(self % iskyl)
    nullify(self % idiag)
    nullify(self % vA)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   CSR to skyline
  !> @details CSR to skyline
  !> 
  !-----------------------------------------------------------------------
    
  subroutine csr2sky(self,csr,ia,ja,SYMMETRIC)

    class(mat_sky),                    intent(inout) :: self
    class(mat_csr), optional,          intent(in)    :: csr
    integer(ip),    optional, pointer, intent(in)    :: ia(:)
    integer(ip),    optional, pointer, intent(in)    :: ja(:)
    logical(lg),    optional,          intent(in)    :: SYMMETRIC
    integer(ip),              pointer                :: ia_loc(:)
    integer(ip),              pointer                :: ja_loc(:)
    integer(ip)                                      :: n,ndof,i,iz,j,k,l
    integer(ip)                                      :: idof,jdof,kskyl
    integer(ip)                                      :: id,jd,it,jt

    if( present(csr) ) then
       n            = csr % nrows
       ndof         = csr % ndof1
       self % nrows = csr % nrows
       self % ndof1 = csr % ndof1
       self % ndof2 = csr % ndof2
    else
       n            = self % nrows
       ndof         = self % ndof1
    end if
    self % symme = optional_argument(self % symme,SYMMETRIC)

    call self % alloca()

    if( present(csr) ) then
       ia_loc => csr % ia
       ja_loc => csr % ja
    else
       ia_loc => ia
       ja_loc => ja
    end if
    !
    ! Skyline format
    !
    do i = 1,n*ndof+1
       self % iskyl(i) = n*ndof
    end do

    do i = 1,n
       do iz = ia_loc(i),ia_loc(i+1)-1
          j = ja_loc(iz)  
          if( i >= j ) then
             do idof = 1,ndof 
                k = (i-1)*ndof+idof+1
                do jdof = 1,ndof
                   l = (j-1)*ndof+jdof
                   if( l < self % iskyl(k) ) self % iskyl(k) = l
                end do
             end do
          end if
       end do
    end do

    self % nskyl    = 1
    self % iskyl(1) = 1

    if( self % symme ) then
       ! 
       ! For the symmetric case, do not need idiag
       !
       do k = 1,n*ndof
          kskyl             = k - self % iskyl(k+1) + 1
          self % nskyl      = self % nskyl + kskyl
          self % iskyl(k+1) = self % nskyl
       end do

    else
       !
       ! For the nonsymmetric case, set idiag 
       !
       do k = 1,n*ndof
          kskyl             = k - self % iskyl(k+1)
          self % idiag(k)   = self % nskyl + kskyl
          kskyl             = 2 * kskyl + 1  
          self % nskyl      = self % nskyl + kskyl
          self % iskyl(k+1) = self % nskyl
       end do

    end if
    self % nskyl = self % nskyl - 1
    !
    ! Copy CSR matrix
    !
    if( present(csr) ) then

       if( .not. associated(self % va) ) call self % alloca_matrix()

       if( ndof == 1 ) then

          if( self % symme ) then
             do i = 1,csr % nrows
                do iz = csr % ia(i),csr % ia(i+1)-1
                   j = csr % ja(iz)
                   if( j <= i ) then
                      kskyl            = self % iskyl(i+1) - 1 - (i-j)
                      self % va(kskyl) = self % va(kskyl) + csr % va(1,1,iz)
                   end if
                end do
             end do
          else
             do i = 1,csr % nrows
                do iz = csr % ia(i),csr % ia(i+1)-1
                   j = csr % ja(iz)
                   if( i < j ) then
                      kskyl            = self % iskyl(j+1) + (i-j)
                      self % va(kskyl) = self % va(kskyl)  + csr % va(1,1,iz)
                   else
                      kskyl            = self % idiag(i)  + (j-i)
                      self % va(kskyl) = self % va(kskyl) + csr % va(1,1,iz)
                   end if
                end do
             end do
          end if

       else

          if( self % symme ) then
             do i = 1,csr % nrows
                do iz = csr % ia(i),csr % ia(i+1)-1
                   j = csr % ja(iz)
                   if( j <= i ) then
                      kskyl            = self % iskyl(i+1) - 1 - (i-j)
                      self % va(kskyl) = self % va(kskyl) + csr % va(1,1,iz)
                   end if
                end do
             end do
          else
             do i = 1,csr % nrows
                do iz = csr % ia(i),csr % ia(i+1)-1
                   j = csr % ja(iz)
                   do id = 1,ndof
                      it = (i-1)*ndof+id
                      do jd = 1,ndof
                         jt = (j-1)*ndof+jd                         
                         if( it < jt ) then
                            kskyl            = self % iskyl(jt+1) + (it-jt)
                            self % va(kskyl) = self % va(kskyl)  + csr % va(1,1,iz)
                         else
                            kskyl            = self % idiag(it)  + (jt-it)
                            self % va(kskyl) = self % va(kskyl) + csr % va(1,1,iz)
                         end if
                      end do
                   end do
                end do
             end do
          end if
          
       end if
    end if

  end subroutine csr2sky

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
    
  subroutine alloca(self,param,MEMORY_COUNTER)
    
    class(mat_sky),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nz_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    
    call memory_alloca(memor_loc,'SELF % ISKYL',vacal,self % iskyl,self % ndof1*self % nrows+1_ip)
    call memory_alloca(memor_loc,'SELF % IDIAG',vacal,self % idiag,self % ndof1*self % nrows)
    call memory_alloca(memor_loc,'SELF % VA'   ,vacal,self % vA   ,self % nskyl)

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine alloca

  subroutine alloca_matrix(self,param,MEMORY_COUNTER)
    
    class(mat_sky),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nskyl

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nskyl     = optional_argument(self % nskyl   ,param,1_ip)
    
    call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,nskyl)
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine alloca_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
    
  subroutine deallo(self,MEMORY_COUNTER)
    
    class(mat_sky),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % ISKYL',vacal,self % iskyl)
    call memory_deallo(memor_loc,'SELF % IDIAG',vacal,self % idiag)
    call memory_deallo(memor_loc,'SELF % VA'   ,vacal,self % vA)
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo

  subroutine deallo_matrix(self,MEMORY_COUNTER)
    
    class(mat_sky),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % VA',vacal,self % vA)
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Assign an array matrix to mat class vA
  !> @details Assign a rank-1 array to mat class vA
  !> 
  !-----------------------------------------------------------------------
    
  subroutine assign(self,a)
    class(mat_sky),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)
  
    !self % vA(1:self % ndof1,1:self % ndof2,1:self % nz) => a
    
  end subroutine assign

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Set number of rows
  !> @details Set number of rows
  !> 
  !-----------------------------------------------------------------------
    
  subroutine set_nrows(self)
    class(mat_sky), intent(inout) :: self

  end subroutine set_nrows
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Set number of columns
  !> @details Set number of columns
  !> 
  !-----------------------------------------------------------------------
    
  subroutine set_ncols(self)
    class(mat_sky), intent(inout) :: self

 
  end subroutine set_ncols
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Get number of rows
  !> @details Get number of rows
  !> 
  !-----------------------------------------------------------------------
    
  pure function get_nrows(self) result(nrows)
    class(mat_sky), intent(in) :: self
    integer(ip)                :: nrows

    nrows = self % nrows

  end function get_nrows

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Get number of columns
  !> @details Get number of columns
  !> 
  !-----------------------------------------------------------------------
    
  pure function get_ncols(self) result(ncols)
    class(mat_sky), intent(in) :: self
    integer(ip)                :: ncols

    ncols = self % nrows
 
  end function get_ncols

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the bandwidth and profile of the matrix
  !> @details Compute the bandwidth and profile of the matrix
  !>
  !-----------------------------------------------------------------------

  subroutine bandwidth(self,bandw,profi)

    class(mat_sky),       intent(inout) :: self
    integer(ip),          intent(out)   :: bandw
    real(rp),             intent(out)   :: profi
    integer(ip)                         :: ii,jj,band,iz

    bandw = 0
    do ii = 1,self % nrows
       iz      = self % iskyl(ii+1)-self % iskyl(ii)
       jj      = ii - iz + 1
       band    = abs(ii-jj)
       bandw   = max(bandw,band)
       profi   = profi + real(band,rp)
    end do

  end subroutine bandwidth

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   y = b - Ax
  !> @details Compute residual y = b - Ax
  !> 
  !-----------------------------------------------------------------------
  
  subroutine residual_1(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_sky), intent(in)                      :: self
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
    
    class(mat_sky), intent(in)                      :: self
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
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mv_11(self,xx,yy,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE,TRANSPOSE)
    
    class(mat_sky), intent(in)                      :: self
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
    
    class(mat_sky), intent(in)                      :: self
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
    
    class(mat_sky), intent(in)                      :: self
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

    class(mat_sky), intent(in)                    :: self
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
    integer(ip)                                   :: ii,jj,iz,ki,kj
    integer(ip)                                   :: kk,ll,kskyl
    real(rp)                                      :: raux1,raux
    real(rp)                                      :: raux2,raux3
    logical(lg)                                   :: use_openmp
    logical(lg)                                   :: do_initialize
    integer(ip)                                   :: do_schedule
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

    if( self % ndof1 == 1 ) then
       if( self % symme ) then

          if( nn1 == 1 .and. nn2 == self % nrows ) then

             do ii = 1,self % nrows
                iz = self % iskyl(ii+1)-self % iskyl(ii)
                jj = ii - iz
                do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
                   jj       = jj + 1
                   yy(1,ii) = yy(1,ii)  + self % vA(kskyl) * xx(1,jj)
                   if( jj /= ii ) &
                        yy(1,jj) = yy(1,jj) + self % vA(kskyl) * xx(1,ii)
                end do
             end do

          else

             do ii = 1,self % nrows
                iz = self % iskyl(ii+1)-self % iskyl(ii)
                jj = ii - iz
                do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
                   jj       = jj + 1
                   if( ii >= nn1 .and. ii <= nn2 ) then
                      yy(1,ii) = yy(1,ii)  + self % vA(kskyl) * xx(1,jj)
                   end if
                   if( jj >= nn1 .and. jj /= ii .and. jj <= nn2 ) then
                      yy(1,jj) = yy(1,jj) + self % vA(kskyl) * xx(1,ii)
                   end if
                end do
             end do

          end if
       else
          if( nn1 == 1 .and. nn2 == self % nrows ) then
             do ii = 1,self % nrows
                jj = ii-(self % idiag(ii)-self % iskyl(ii))
                do kk = 1,self % idiag(ii)-self % iskyl(ii)+1
                   kskyl    = self % idiag(ii)  + (jj-ii)
                   yy(1,ii) = yy(1,ii) + self % va(kskyl) * xx(1,jj)
                   jj       = jj + 1
                end do
                jj = ii-(self % iskyl(ii+1)-1-self % idiag(ii))          
                do kskyl    = self % idiag(ii)+1,self % iskyl(ii+1)-1
                   yy(1,jj) = yy(1,jj) + self % va(kskyl) * xx(1,ii)
                   jj       = jj + 1
                end do
             end do
          else
             do ii = 1,self % nrows
                jj = ii-(self % idiag(ii)-self % iskyl(ii))
                do kk = 1,self % idiag(ii)-self % iskyl(ii)+1
                   kskyl    = self % idiag(ii)  + (jj-ii)
                   if( ii >= nn1 .and. ii <= nn2 ) &
                        yy(1,ii) = yy(1,ii) + self % va(kskyl) * xx(1,jj)
                   jj       = jj + 1
                end do
                jj = ii-(self % iskyl(ii+1)-1-self % idiag(ii))          
                do kskyl    = self % idiag(ii)+1,self % iskyl(ii+1)-1
                   if( jj >= nn1 .and. jj <= nn2 ) &
                        yy(1,jj) = yy(1,jj) + self % va(kskyl) * xx(1,ii)
                   jj       = jj + 1
                end do
             end do
          end if
       end if
    else
       stop 77
    end if

  end subroutine mv_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine dirichlet(self,fixno,bvess,rhs,n1,n2)
    
    class(mat_sky),                    intent(inout) :: self
    integer(ip),    optional, pointer, intent(in)    :: fixno(:)
    real(rp),       optional, pointer, intent(in)    :: bvess(:)
    real(rp),       optional, pointer, intent(inout) :: rhs(:)
    integer(ip),    optional,          intent(in)    :: n1            
    integer(ip),    optional,          intent(in)    :: n2
    integer(ip)                                      :: ii,idofn,jdofn,nn1,nn2
    integer(ip)                                      :: izdod,iz,jj,kk,jz
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

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal position
  !> @details Diagonal position in matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diagz(self,diagonal)
    
    class(mat_sky),          intent(in)    :: self
    integer(ip),    pointer, intent(inout) :: diagonal(:)
    integer(ip)                            :: ii,jj,iz

  end subroutine diagz
     
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine diag(self,dia,diagonal,row,val)
    
    class(mat_sky),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,jj,iz,idofn
    integer(ip)                                      :: nrows,ndofn,it

    nrows = self % nrows
    ndofn = self % ndof2
    
    if( present(diagonal) ) then

       if( self % symme ) then
          do ii = 1,nrows
             iz = self % iskyl(ii+1)-1
             do idofn = 1,ndofn
                it = (iz-1)*ndofn*ndofn+(idofn-1)*ndofn+idofn
                diagonal(idofn,ii) = self % vA(it)
             end do
          end do
       else
          do ii = 1,nrows
             iz = self % idiag(ii)
             do idofn = 1,ndofn
                it = (iz-1)*ndofn*ndofn+(idofn-1)*ndofn+idofn
                diagonal(idofn,ii) = self % vA(it)
             end do
          end do
       end if
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia )
          
          if( .not. associated(dia % vA) ) call dia % alloca((/ndofn,nrows/)) 
          if( self % symme ) then
             do ii = 1,nrows
                iz = self % iskyl(ii+1)-1
                do idofn = 1,ndofn
                   it = (iz-1)*ndofn*ndofn+(idofn-1)*ndofn+idofn
                   dia % vA(idofn,ii) = self % vA(it)
                end do
             end do
          else
             do ii = 1,nrows
                iz = self % idiag(ii)
                do idofn = 1,ndofn
                   it = (iz-1)*ndofn*ndofn+(idofn-1)*ndofn+idofn
                   dia % vA(idofn,ii) = self % vA(it)
                end do
             end do
          end if
          
       end select
       
    else if( present(row) .and. present(val) ) then
       
       ii = row
       if( self % symme ) then
          iz = self % iskyl(ii+1)-1
       else
          iz = self % idiag(ii)          
       end if
       do idofn = 1,ndofn
          it = (iz-1)*ndofn*ndofn+(idofn-1)*ndofn+idofn
          val(idofn) = self % vA(it)
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

    class(mat_sky), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    real(rp)                      :: dummr
    integer(ip)                   :: ii,jj,iz,ncols
    integer(ip)                   :: idof1,idof2
    real(rp),       pointer       :: aa(:,:)


  end function norm

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-08
  !> @brief   Output
  !> @details Output matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine output(self,FMT,FILENAME,PERM)

    class(mat_sky),                   intent(in) :: self
    character(*),                     intent(in) :: FMT
    character(*),  optional,          intent(in) :: FILENAME 
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: ndof1,ndof2
    integer(ip)                                  :: ndofn,nrows,ncols,jt
    integer(ip)                                  :: ii,idof2,jj,kk,ll
    integer(ip)                                  :: col,kskyl,idof1,irow
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc
    real(rp),      allocatable                   :: vA(:)
    !
    ! Dimensions
    !
    ndof1 = self % ndof1
    ndof2 = self % ndof2
    ndofn = self % ndof1 * self % ndof2
    nrows = self % nrows
    ncols = self % get_ncols()
    !
    ! Unit and files
    !
    unit4 = iofile_available_unit(90_ip) 

    select case ( upper_case(FMT) )

    case ( 'DENSE' )

       !-----------------------------------------------------------------
       !
       ! Dense
       !
       !-----------------------------------------------------------------

       filename_loc = optional_argument('matrix-dense.txt',FILENAME)
       call add_extension(filename_loc,'txt')
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       
       allocate(vA(nrows))
       
       if( self % symme ) then

          if( ndof1 /= 1 .or. ndof2 /= 1 ) stop 1
          do ii = 1,nrows
             vA  = 0.0_rp
             ll  = self % iskyl(ii+1)-self % iskyl(ii)               
             ! column number of the first non zero in row i
             col = ii-ll
             do kk = 1,ll
                col     = col + 1
                vA(col) = self % vA(self % iskyl(ii)+kk-1)
             end do
             do jj = 1,nrows
                if( ii /= jj ) then
                   ll = self % iskyl(jj+1)-self % iskyl(jj)
                   jt = jj - ll 
                   loop_kk: do kk = 1,ll
                      jt = jt + 1
                      if( jt == ii ) then
                         vA(jj) = self % vA(self % iskyl(jj)+kk-1)
                         exit loop_kk
                      end if
                   end do loop_kk
                end if
             end do
             do jj = 1,nrows-1
                write(unit4,'(1x,e12.6,a)',advance='no') vA(jj),' ,'
             end do
             write(unit4,'(1x,e12.6,a)',advance='no') vA(jj)             
             write(unit4,*)
          end do

       else
          
          do ii = 1,nrows
             vA = 0.0_rp
             !
             ! Lower part
             !
             ! First column of row ii
             jj = ii-(self % idiag(ii)-self % iskyl(ii))
             do kk = 1,self % idiag(ii)-self % iskyl(ii)+1
                kskyl  = self % idiag(ii)  + (jj-ii)
                vA(jj) = self % va(kskyl)
                jj     = jj + 1
             end do
             ! Traverse down the columns to find i column
             do jj = ii+1,nrows
                ! First row of column jj
                irow = jj-(self % iskyl(jj+1)-1-self % idiag(jj))
                loop_kk2: do kskyl = self % idiag(jj)+1,self % iskyl(jj+1)-1
                   if( ii == irow ) then
                      vA(jj) = self % va(kskyl)
                      exit loop_kk2
                   end if
                   irow = irow + 1
                end do loop_kk2
             end do
             do jj = 1,nrows-1
                write(unit4,'(1x,e12.6,a)',advance='no') vA(jj),' ,'
             end do
             write(unit4,'(1x,e12.6,a)',advance='no') vA(jj)
             write(unit4,*)
          end do

          
       end if
       
       deallocate(vA)
       close(unit4)

    end select

  end subroutine output

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

    class(mat_sky),           intent(in) :: self
    real(rp),       optional, intent(in) :: TOLERANCE
    real(rp)                             :: sym

    sym = -1.0_rp
    if( self % symme ) then
    else
    end if
    
  end function symmetry

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Mv on a single row
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_sky),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof

    yy = 0.0_rp

    if( self % ndof1 == 1 ) then
       if( self % symme ) then
          
       end if
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
    
    class(mat_sky), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: iz,kskyl
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
    if( self % symme ) then
        do ii = nn1,nn2
          iz = self % iskyl(ii+1)-self % iskyl(ii)
          jj = ii - iz
          do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
             jj       = jj + 1
             yy(ii) = yy(ii)  + self % vA(kskyl) * xx(jj)
          end do
       end do      
    else
       do ii = nn1,nn2
          iz = self % iskyl(ii+1)-self % iskyl(ii)
          jj = ii - iz
          do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
             jj       = jj + 1
             if( jj <= ii ) &
                  yy(ii) = yy(ii)  + self % vA(kskyl) * xx(jj)
          end do
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
    
    class(mat_sky), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: iz,kskyl
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
    if( self % symme ) then
       do ii = nn1,nn2
          iz = self % iskyl(ii+1)-self % iskyl(ii)
          jj = ii - iz
          do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
             jj     = jj + 1
             yy(ii) = yy(ii)  + self % vA(kskyl) * xx(jj)
          end do
       end do      
    else
       do ii = nn1,nn2
          iz = self % iskyl(ii+1)-self % iskyl(ii)
          jj = ii - iz
          do kskyl = self % iskyl(ii),self % iskyl(ii+1)-1
             jj       = jj + 1
             if( jj >= ii ) &
                  yy(ii) = yy(ii)  + self % vA(kskyl) * xx(jj)
          end do
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
     
    class(mat_sky), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    
  end function get_val

end module def_mat_sky
!> @}
