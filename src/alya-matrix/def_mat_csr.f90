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
!> @brief   CSR Matrix
!> @details CSR Matrix class
!>          Required to allocate:
!>          SELF % NZ
!>          SELF % NDOF1
!>          SELF % NDOF2
!>          SELF % NROWS
!>          SELF % IA(SELF % NROWS+1)
!>          SELF % JA(SELF % NZ)
!>          SELF % VA(SELF % NDOF1,SELF % NDOF2,SELF % NZ)
!-----------------------------------------------------------------------

module def_mat_csr

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat,               only : CSR_FORMAT
  use def_mat,               only : OMP_CHUNK  
  use def_mat,               only : OMP_STATIC 
  use def_mat,               only : OMP_GUIDED 
  use def_mat,               only : OMP_DYNAMIC 
  use mod_optional_argument, only : optional_argument
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end
  use mod_iofile_basic,      only : iofile_available_unit
  use mod_strings,           only : upper_case
  use mod_strings,           only : add_extension
  implicit none

  private

  type, extends(mat) :: mat_csr
     integer(ip)                      :: nz
     integer(ip),             pointer :: iA(:)
     integer(ip),             pointer :: jA(:)     
     real(rp),    contiguous, pointer :: vA(:,:,:)
   contains
     procedure,               pass    :: init           ! Initialize the class
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix           
     procedure,               pass    :: set_nrows      ! Set rows      
     procedure,               pass    :: get_nrows      ! Get rows      
     procedure,               pass    :: set_ncols      ! Set columns      
     procedure,               pass    :: get_ncols      ! Get columns
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: bandwidth      ! Compute bandwidth and profile
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: alloca_matrix  ! Allocate matrix only
     procedure,               pass    :: deallo_matrix  ! Deallocate matrix only
     procedure,               pass    :: output         ! Output matrix
     procedure,               pass    :: diag           ! Compute diagonal
     procedure,               pass    :: diagz          ! Compute diagonal position
     procedure,               pass    :: scale          ! Scale a matrix
     procedure,               pass    :: norm           ! Matrix norm
     procedure,               pass    :: symmetry       ! Determine the symmetry
     procedure,               pass    :: mv_row         ! Single row MV product
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: copy           ! Copy a matrix
     procedure,               pass    :: edge           ! Find edge number corresponding to position (i,j)
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
     procedure,               pass    :: size_to_linked_list ! Transform ia
     procedure,               pass    :: equal               ! Check if two matrices are equal
  end type mat_csr

  character(11), parameter :: vacal = 'def_mat_csr'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_csr
  
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

    class(mat_csr), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = CSR_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    self % nz         = 0
    nullify(self % iA)
    nullify(self % jA)
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
    
    class(mat_csr),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nz_loc,ndof1_loc,ndof2_loc
    integer(ip)                              :: nrows_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nz_loc    = optional_argument(self % nz   ,param,1_ip)
    ndof1_loc = optional_argument(self % ndof1,param,2_ip)
    ndof2_loc = optional_argument(self % ndof2,param,3_ip)
    nrows_loc = optional_argument(self % nrows,param,4_ip)
    
    if( nz_loc > 0 ) then
       call memory_alloca(memor_loc,'SELF % IA',vacal,self % iA,nrows_loc+1_ip)
       call memory_alloca(memor_loc,'SELF % JA',vacal,self % jA,nz_loc)
       call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,ndof2_loc,nz_loc)        
       self % nz    = nz_loc       
       self % ndof1 = ndof1_loc
       self % ndof2 = ndof2_loc       
       self % nrows = nrows_loc       
    end if
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine alloca

  subroutine alloca_matrix(self,param,MEMORY_COUNTER)
    
    class(mat_csr),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nz_loc,ndof1_loc,ndof2_loc
    integer(ip)                              :: nrows_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nz_loc    = optional_argument(self % nz   ,param,1_ip)
    ndof1_loc = optional_argument(self % ndof1,param,2_ip)
    ndof2_loc = optional_argument(self % ndof2,param,3_ip)
    nrows_loc = optional_argument(self % nrows,param,4_ip)
    
    if( nz_loc > 0 ) then
       call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,ndof2_loc,nz_loc) 
    end if
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
    
    class(mat_csr),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % IA',vacal,self % iA)
    call memory_deallo(memor_loc,'SELF % JA',vacal,self % jA)
    call memory_deallo(memor_loc,'SELF % VA',vacal,self % vA)
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo

  subroutine deallo_matrix(self,MEMORY_COUNTER)
    
    class(mat_csr),            intent(inout) :: self
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
    class(mat_csr),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)
  
    self % vA(1:self % ndof1,1:self % ndof2,1:self % nz) => a
    
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
    class(mat_csr), intent(inout) :: self

    if( associated(self % iA) ) then
       self % nrows = size(self % IA)-1
    else
       self % nrows = 0
    end if
       
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
    class(mat_csr), intent(inout) :: self

    if( associated(self % jA) ) then
       self % ncols = maxval(self % jA)
    else
       self % ncols = 0
    end if
       
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
    class(mat_csr), intent(in) :: self
    integer(ip)                :: nrows

    if( self % nrows > 0 ) then
       nrows = self % nrows
    else
       if( associated(self % iA) ) then
          nrows = int(size(self % IA),ip)-1_ip
       else
          nrows = 0
       end if
    end if
    
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
    class(mat_csr), intent(in) :: self
    integer(ip)                :: ncols

    if( self % ncols > 0 ) then
       ncols = self % ncols
    else
       if( associated(self % jA) ) then
          ncols = maxval(self % jA)
       else
          ncols = 0
       end if
    end if
    
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

    class(mat_csr),       intent(inout) :: self
    integer(ip),          intent(out)   :: bandw
    real(rp),             intent(out)   :: profi
    integer(ip)                         :: ii,jj,band,iz,bandloc,ipmax,jpmax  

    bandw =  0_ip
    profi =  0.0_rp

    if( associated(self % ia) ) then

       do ii = 1,size(self % iA)-1

          bandloc = 0_ip
          !
          ! Loop on neighbors
          !
          do iz = self % iA(ii),self % iA(ii+1)-1
             jj = self % jA(iz)
             if( ii /= jj ) then
                band = abs(jj-ii)
                !
                ! Test bandwidth
                !
                if( band > bandw ) then
                   bandw = band
                   ipmax = ii
                   jpmax = jj
                endif
                !
                ! Test profile
                !
                if( jj < ii ) bandloc = max(bandloc,band)
             end if
          end do
          !
          ! Accumulate profile
          !
          profi = profi + real(bandloc,rp)

       end do

    end if

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
    
    class(mat_csr), intent(in)                      :: self
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
    
    class(mat_csr), intent(in)                      :: self
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
    
    class(mat_csr), intent(in)                      :: self
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
    
    class(mat_csr), intent(in)                      :: self
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
    
    class(mat_csr), intent(in)                      :: self
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

    class(mat_csr), intent(in)                    :: self
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
    integer(ip)                                   :: col,kk,ll
    real(rp)                                      :: raux1,raux
    real(rp)                                      :: raux2,raux3
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
       
       do ii = nn1,nn2
          do jj  = self % iA(ii),self % iA(ii+1)-1
             col = self % jA(jj)
             do ll = 1,ndofc
                do kk = 1,ndofr
                   yy(ll,col) = yy(ll,col) + self % vA(kk,ll,jj) * xx(kk,ii)
                end do
             end do
          end do
       end do
       
    else
       if( ndofr == 1 .and. ndofc == 1 ) then
          !
          ! NDOF=1
          !
          if( use_openmp ) then
             if( my_schedule == OMP_STATIC ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( STATIC )                           &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = nn1,nn2
                   do jj = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( GUIDED )                           &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = nn1,nn2
                   do jj = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                 &
                !$OMP SCHEDULE ( DYNAMIC, my_chunk )                &
                !$OMP DEFAULT  ( NONE )                             &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj )
                do ii = nn1,nn2
                   do jj = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * xx(1,col)
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = nn1,nn2
                do jj = self % iA(ii),self % iA(ii+1)-1
                   col      = self % jA(jj)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * xx(1,col)
                end do
             end do
          end if

       else if( ndofr == 2 .and. ndofc == 2 ) then
          !
          ! NDOF=2
          !
          if( use_openmp ) then
             if( my_schedule == OMP_STATIC ) then
                !$OMP PARALLEL  DO                                    &
                !$OMP SCHEDULE ( STATIC )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( GUIDED )                             &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                   &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                 &
                !$OMP DEFAULT  ( NONE )                               &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )   &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = nn1,nn2
                do jj       = self % iA(ii),self % iA(ii+1)-1
                   col      = self % jA(jj)
                   raux1    = xx(1,col)
                   raux2    = xx(2,col)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2
                   yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2
                end do
             end do
          end if

       else if( ndofr == 3 .and. ndofc == 3 ) then
          !
          ! NDOF=3
          !
          if( use_openmp ) then
             if( my_schedule == OMP_STATIC ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( STATIC )                              &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2 + self % vA(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2 + self % vA(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + self % vA(1,3,jj) * raux1 + self % vA(2,3,jj) * raux2 + self % vA(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( GUIDED )                              &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2 + self % vA(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2 + self % vA(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + self % vA(1,3,jj) * raux1 + self % vA(2,3,jj) * raux2 + self % vA(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                    &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                  &
                !$OMP DEFAULT  ( NONE )                                &
                !$OMP SHARED   ( self, nn1, nn2, xx, yy )    &
                !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
                do ii = nn1,nn2
                   do jj       = self % iA(ii),self % iA(ii+1)-1
                      col      = self % jA(jj)
                      raux1    = xx(1,col)
                      raux2    = xx(2,col)
                      raux3    = xx(3,col)
                      yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2 + self % vA(3,1,jj) * raux3
                      yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2 + self % vA(3,2,jj) * raux3
                      yy(3,ii) = yy(3,ii) + self % vA(1,3,jj) * raux1 + self % vA(2,3,jj) * raux2 + self % vA(3,3,jj) * raux3
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = nn1,nn2
                do jj       = self % iA(ii),self % iA(ii+1)-1
                   col      = self % jA(jj)
                   raux1    = xx(1,col)
                   raux2    = xx(2,col)
                   raux3    = xx(3,col)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,jj) * raux1 + self % vA(2,1,jj) * raux2 + self % vA(3,1,jj) * raux3
                   yy(2,ii) = yy(2,ii) + self % vA(1,2,jj) * raux1 + self % vA(2,2,jj) * raux2 + self % vA(3,2,jj) * raux3
                   yy(3,ii) = yy(3,ii) + self % vA(1,3,jj) * raux1 + self % vA(2,3,jj) * raux2 + self % vA(3,3,jj) * raux3
                end do
             end do
          end if

       else
          !
          ! NDOF = whatever
          !
          if( use_openmp ) then
             if( my_schedule == OMP_STATIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( STATIC )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( self, nn1, nn2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = nn1,nn2
                   do jj  = self % iA(ii),self % iA(ii+1)-1
                      col = self % jA(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_GUIDED ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( GUIDED )                                         &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( self, nn1, nn2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = nn1,nn2
                   do jj  = self % iA(ii),self % iA(ii+1)-1
                      col = self % jA(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             else if( my_schedule == OMP_DYNAMIC ) then
                !$OMP PARALLEL   DO                                               &
                !$OMP SCHEDULE ( DYNAMIC , my_chunk )                             &
                !$OMP DEFAULT  ( NONE )                                           &
                !$OMP SHARED   ( self, nn1, nn2, ndofr, ndofc, xx, yy ) &
                !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
                do ii = nn1,nn2
                   do jj  = self % iA(ii),self % iA(ii+1)-1
                      col = self % jA(jj)
                      do ll = 1,ndofc
                         raux = xx(ll,col)
                         yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(ll,1:ndofr,jj) * raux
                      end do
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             do ii = nn1,nn2
                do jj  = self % iA(ii),self % iA(ii+1)-1
                   col = self % jA(jj)
                   do ll = 1,ndofc
                      raux = xx(ll,col)
                      yy(1:ndofr,ii) = yy(1:ndofr,ii) + self % vA(ll,1:ndofr,jj) * raux
                   end do
                end do
             end do
          end if
       end if
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
    
    class(mat_csr),                    intent(inout) :: self
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
             !
             ! Eliminate dof of II from other equations (JJ)
             !
             do iz = self % ia(ii),self % ia(ii+1) - 1
                jj = self % ja(iz)
                if( ii /= jj ) then
                   do jz = self % ia(jj),self % ia(jj+1) - 1
                      kk = self % ja(jz)
                      if( kk == ii ) then
                         do jdofn = 1,ndofr
                            self % vA(idofn,jdofn,jz) = 0.0_rp
                         end do
                      end if
                   end do
                end if
             end do
             !
             ! IZDOD: Diagonal
             !
             izdod = self % ia(ii) - 1
             jj    = 0
             do while( jj /= ii )
                izdod = izdod + 1
                jj    = self % ja(izdod)
             end do
             diag(idofn) = self % vA(idofn,idofn,izdod)
             if( abs(diag(idofn)) < epsil ) diag(idofn) = 1.0_rp
             !
             ! Set line to zero and keep diagonal
             !
             do iz = self % ia(ii),self % ia(ii+1) - 1
                do jdofn = 1,ndofc
                   self % vA(jdofn,idofn,iz) = 0.0_rp
                end do
             end do
             self % vA(idofn,idofn,izdod) = diag(idofn)
             if( present(rhs) .and. present(bvess) ) then
                rhs((ii-1)*ndofr+idofn) = bvess((ii-1)*ndofr+idofn) * diag(idofn)
             end if
          end if
       end do
    end do   

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
    
    class(mat_csr),          intent(in)    :: self
    integer(ip),    pointer, intent(inout) :: diagonal(:)
    integer(ip)                            :: ii,jj,iz

    do ii = 1,self % nrows
       iz = self % ia(ii) - 1
       jj = 0
       do while( jj /= ii )
          iz = iz + 1
          jj = self % ja(iz)
       end do
       diagonal(ii) = iz
    end do
       
  end subroutine diagz
     
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Scale
  !> @details Scale a matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine scale(self,dia,SQUARE_ROOT)
    
    class(mat_csr),          intent(inout) :: self
    type(mat_dia),           intent(in)    :: dia
    logical(lg),   optional, intent(in)    :: SQUARE_ROOT
    integer(ip)                            :: ii,jj,iz,idof1

    if( optional_argument(.false.,SQUARE_ROOT) ) then
       do ii = 1,self % nrows
          do iz = self % iA(ii),self % iA(ii+1)-1
             jj = self % jA(iz)
             do idof1 = 1,self % ndof1
                self % va(1:self % ndof2,idof1,iz) = self % va(1:self % ndof2,idof1,iz) &
                     * abs(sqrt(dia % va(idof1,ii)))
             end do
          end do
       end do
    else
       do ii = 1,self % nrows
          do iz = self % iA(ii),self % iA(ii+1)-1
             jj = self % jA(iz)
             do idof1 = 1,self % ndof1
                self % va(1:self % ndof2,idof1,iz) = self % va(1:self % ndof2,idof1,iz) &
                     * dia % va(idof1,ii)
             end do
          end do
       end do
    end if
    
  end subroutine scale
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine diag(self,dia,diagonal,row,val)
    
    class(mat_csr),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,jj,iz,idofn
    integer(ip)                                      :: nrows,ndofn

    nrows = self % get_nrows()
    ndofn = self % ndof2
    
    if( present(diagonal) ) then
       
       do ii = 1,nrows
          iz = self % ia(ii) - 1
          jj = 0
          do while( jj /= ii )
             iz = iz + 1
             jj = self % ja(iz)
          end do
          do idofn = 1,ndofn
             diagonal(idofn,ii) = self % vA(idofn,idofn,iz)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia )
          
          if( .not. associated(dia % vA) ) call dia % alloca((/ndofn,nrows/)) 
          do ii = 1,nrows
             iz = self % ia(ii) - 1
             jj    = 0
             do while( jj /= ii )
                iz = iz + 1
                jj    = self % ja(iz)
             end do
             do idofn = 1,ndofn
                dia % vA(idofn,ii) = self % vA(idofn,idofn,iz)
             end do
          end do
          
       end select
       
    else if( present(row) .and. present(val) ) then
       
       ii = row
       iz = self % ia(ii) - 1
       jj = 0
       do while( jj /= ii )
          iz = iz + 1
          jj = self % ja(iz)
       end do
       do idofn = 1,ndofn
          val(idofn) = self % vA(idofn,idofn,iz)
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

    class(mat_csr), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    real(rp)                      :: dummr
    integer(ip)                   :: ii,jj,iz,ncols
    integer(ip)                   :: idof1,idof2
    real(rp),       pointer       :: aa(:,:)

    anorm = 0.0_rp
    if( self % nz > 0 ) then

       if( wnorm == 'i' .or. wnorm == 'I' ) then

          do ii = 1,self % nrows
             do idof2 = 1,self % ndof2
                dummr = 0.0_rp
                do iz = self % iA(ii),self % iA(ii+1)-1
                   jj =  self % jA(iz)
                   do idof1 = 1,self % ndof1
                      dummr = dummr + abs(self % vA(idof1,idof2,iz))
                   end do
                end do
                anorm = max(anorm,dummr)
             end do
          end do

       else if( wnorm == '1' ) then

          ncols = self % get_ncols()
          allocate(aa(self % ndof1,ncols))
          do ii = 1,ncols
             aa(:,ii) = 0.0_rp
          end do

          do ii = 1,self % nrows
             do idof2 = 1,self % ndof2
                do iz = self % iA(ii),self % iA(ii+1)-1
                   jj = self % jA(iz)
                   do idof1 = 1,self % ndof1
                      aa(idof1,jj) = aa(idof1,jj) + abs(self % vA(idof1,idof2,iz))
                   end do
                end do
             end do
          end do
          do ii = 1,ncols
             anorm = max(anorm,maxval(aa(1:self % ndof1,ii)))
          end do
          deallocate(aa)

       else if( wnorm == 'f' .or. wnorm == 'F' ) then

          do iz = 1,self % nz
             do idof2 = 1,self % ndof2
                do idof1 = 1,self % ndof1
                   anorm = anorm + self % vA(idof1,idof2,iz)*self % vA(idof1,idof2,iz)
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
  !> @date    2021-05-08
  !> @brief   Output
  !> @details Output matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine output(self,FMT,FILENAME,PERM)

    class(mat_csr),                   intent(in) :: self
    character(*),                     intent(in) :: FMT
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: iz,ii,jj,kk,ll,nz
    integer(ip)                                  :: nx,ny,ii_old,jj_old
    integer(ip)                                  :: idof1,idof2,ndof1,ndof2
    integer(ip)                                  :: ie,ke,je
    integer(ip)                                  :: ndofn,nrows,ncols
    integer(4)                                   :: unit4
    character(150)                               :: filename_loc
    logical(lg)                                  :: notfound
    real(rp)                                     :: xx,yy
    !
    ! Dimensions
    !
    ndof1 = self % ndof1
    ndof2 = self % ndof2
    ndofn = self % ndof1 * self % ndof2
    nz    = ndofn * self % nz
    nrows = self % nrows
    ncols = self % get_ncols()
    !
    ! Unit and files
    !
    unit4 = iofile_available_unit(90_ip) ! TODO: possible int8 to int4 conversion

    select case ( upper_case(FMT) )

    case ( 'MATRIX MARKET' , 'MTX' , 'MM' )

       !-----------------------------------------------------------------
       !
       ! Matrix market
       !
       !-----------------------------------------------------------------

       filename_loc = optional_argument('matrix.mtx',FILENAME)
       call add_extension(filename_loc,'mtx')
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       write(unit4,'(a)')         '%%MatrixMarket matrix coordinate real general'
       write(unit4,'(3(1x,i12))') nrows * ndof2, ncols * ndof1,nz       
       do ii = 1,self % nrows
          do idof2 = 1,self % ndof2
             do iz = self % ia(ii),self % ia(ii+1)-1
                jj = self % ja(iz)
                do idof1 = 1,self % ndof1
                   if( present(PERM) ) then
                      write(unit4,'(2(1x,i9),1x,e13.6)') PERM((ii-1)*self % ndof2+idof2),PERM((jj-1)*self % ndof1+idof1),self % vA(idof1,idof2,iz)
                   else
                      write(unit4,'(2(1x,i9),1x,e13.6)') (ii-1)*self % ndof2+idof2,(jj-1)*self % ndof1+idof1,self % vA(idof1,idof2,iz)
                   end if
                end do
             end do
          end do
       end do
       close(unit=unit4)

    case ( 'PAJEK NET' )

       !-----------------------------------------------------------------
       !
       ! Matrix market
       !
       !-----------------------------------------------------------------

       filename_loc = optional_argument('matrix.net',FILENAME)
       call add_extension(filename_loc,'net')
       open(unit=unit4,file=trim(filename_loc),status='unknown')       
       write(unit4,'(a,1x,i7)') '*Vertices ',nrows
       write(unit4,'(a)')       '*arcs'       
       do ii = 1,nrows
          do idof2 = 1,ndof2
             do iz = self % ia(ii),self % ia(ii+1)-1
                jj = self % ja(iz)
                do idof1 = 1,ndof1
                   if( abs(self % vA(idof1,idof2,iz)) > 1.0e-12_rp ) then
                      write(unit4,'(2(1x,i6),1x,e12.6,a)',advance='no') &
                           (ii-1)*ndof2+idof2,(jj-1)*ndof1+idof1,self % vA(idof1,idof2,iz)
                   end if
                end do
             end do
          end do
       end do
       close(unit4)

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
                kk = 0
                ll = self % ia(ii)-1
                do while( ll < self % ia(ii+1)-1 .and. kk == 0 )
                   ll = ll + 1
                   if( self % ja(ll) == jj ) kk = ll
                end do
                do idof1 = 1,ndof1
                   if( kk /= 0 ) then
                      write(unit4,'(1x,e12.6,a)',advance='no') self % vA(idof1,idof2,kk),' ,'
                   else
                      write(unit4,'(1x,e12.6,a)',advance='no') 0.0_rp,' ,'
                   end if
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

       filename_loc = optional_argument('matrix.post.msh',FILENAME)
       call add_extension(filename_loc,'post.msh')
       open(unit=unit4,file=trim(filename_loc),status='unknown')

       nx = nrows * ndof2 + 1
       ny = ncols * ndof1 + 1
       !
       ! Mesh
       !
       write(unit4,'(a)') 'MESH MATRIX dimension 2 Elemtype Quadrilateral Nnode 4'
       kk = 0
       write(unit4,'(a)') 'coordinates'
       do jj = 1,nx
          yy = -real(jj,rp)
          do ii = 1,ny
             xx = real(ii,rp)
             kk = kk + 1
             write(unit4,'(i9,2(1x,e13.6))') kk,xx,yy
          end do
       end do
       write(unit4,'(a)') 'end coordinates'

       ke = 0
       write(unit4,'(a)') 'elements'
       do je = 1,nrows*ndof2
          do ie = 1,ncols*ndof1
             ke = ke + 1
             ii = ny*(je-1) + ie
             write(unit4,'(5(1x,i9))') ke,ii+1,ii,ii+ny,ii+ny+1
          end do
       end do
       write(unit4,'(a)') 'end elements'
       close(unit=unit4)
       !
       ! Matrix coefficients
       !
       filename_loc = optional_argument('matrix.post.res',FILENAME)
       call add_extension(filename_loc,'post.res')
       open(unit=unit4,file=trim(filename_loc),status='unknown')
       write(unit4,*) 'GiD Post Results File 1.0'
       write(unit4,*) ' '
       write(unit4,*) 'GaussPoints GP Elemtype Quadrilateral'
       write(unit4,*) 'Number of Gauss Points: 1'
       write(unit4,*) 'Natural Coordinates: Internal'
       write(unit4,*) 'End GaussPoints'
       write(unit4,*) 'Result MATRIX ALYA  0.00000000E+00 Scalar OnGaussPoints GP'
       write(unit4,*) 'ComponentNames MATRIX'
       write(unit4,*) 'Values'
       kk = 0
       do ii_old = 1,nrows
          do idof2 = 1,ndof2
             if( present(PERM) ) then
                ii = PERM(ii_old)
             else
                ii = ii_old
             end if
             do jj_old = 1,ncols
                if( present(PERM) ) then
                   jj = PERM(jj_old)
                else
                   jj = jj_old
                end if
                iz    = self % ia(ii)-1
                notfound = .true.
                do while( notfound .and. iz < self % ia(ii+1)-1 )
                   iz = iz + 1
                   if( self % ja(iz) == jj ) notfound = .false.
                end do
                do idof1 = 1,ndof1
                   kk = kk + 1
                   if( notfound ) then
                      write(unit4,'(i9,1x,e13.6)') kk,0.0_rp
                   else
                      write(unit4,'(i9,1x,e13.6)') kk,abs(self % vA(idof1,idof2,iz))
                   end if
                end do
             end do
          end do
       end do
       write(unit4,*) 'End Values'
       close(unit=unit4)

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

    class(mat_csr),           intent(in) :: self
    real(rp),       optional, intent(in) :: TOLERANCE
    integer(ip)                          :: ii,iz,lz,jj
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
             do iz = self % ia(ii),self % ia(ii+1)-1
                jj = self % ja(iz)
                loop_lz: do lz = self % ia(jj),self % ia(jj+1)-1
                   if( self % ja(lz) == ii ) then
                      do idof2 = 1,self % ndof2
                         do idof1 = 1,self % ndof1
                            if( idof1 /= idof2 .or. ii /= jj ) then
                               if( abs(self % vA(idof1,idof2,iz)-self % vA(idof2,idof1,lz)) > toler ) then
                                  sym = max(sym,abs(self % vA(idof1,idof2,iz)-self % vA(idof2,idof1,lz)) / vamax)
                               end if
                            end if
                         end do
                      end do
                      exit loop_lz
                   end if
                end do loop_lz
             end do
          end do
       end if

    end if

  end function symmetry

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Scale a vector
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_csr),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: iz,ll,jj

    yy = 0.0_rp
    do iz = self % iA(n),self % iA(n+1)-1
       jj = self % jA(iz)
       do ll = 1,self % ndof2
          yy = yy + self % vA(ll,ndof,iz) * xx((jj-1)*self % ndof2+ll)
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
    
    class(mat_csr), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: iz,id,jd
    integer(ip)                                     :: it,jt
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    do ii = nn1,nn2
       do iz = self % iA(ii),self % iA(ii+1)-1
          jj = self % jA(iz)
          if( jj <= ii ) then
             do id = 1,self % ndof1
                it = (ii-1) * self % ndof1 + id
                do jd = 1,id
                   jt = (jj-1) * self % ndof2+jd
                   yy(it) = yy(it) + self % vA(jd,id,iz) * xx(jt)
                end do
             end do
          end if
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
    
    class(mat_csr), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: iz,id,jd
    integer(ip)                                     :: it,jt
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    do ii = nn1,nn2
       do iz = self % iA(ii),self % iA(ii+1)-1
          jj = self % jA(iz)
          if( jj >= ii ) then
             do id = 1,self % ndof1
                it = (ii-1) * self % ndof1 + id
                do jd = id,self % ndof2
                   jt = (jj-1) * self % ndof2+jd
                   yy(it) = yy(it) + self % vA(jd,id,iz) * xx(jt)
                end do
             end do
          end if
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
    
    class(mat_csr), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    integer(ip)                :: iz

    a = 0.0_rp
    do iz = self % ia(i),self % ia(i+1)-1
       if( self % ja(iz) == j ) then
          a(:,:) = self % vA(:,:,iz)
          return
       end if
    end do
    
  end function get_val

    !----------------------------------------------------------------------
    !>
    !> @author  Guillaume Houzeaux
    !> @date    19/07/2017
    !> @brief   Copy a matrix
    !> @details Copy a matrix Ain to a matrix Aout. If the graph of the
    !>          output matrix is present, then look for corresponding
    !>          positions of Ain
    !>
    !----------------------------------------------------------------------

  subroutine copy(self,mat_in,invp,PARTIAL_MATRIX)

    class(mat_csr),                    intent(inout) :: self            !< Output
    type(mat_csr),                     intent(in)    :: mat_in          !< Input matrix
    integer(ip),    optional, pointer, intent(in)    :: invp(:)         !< Permutation
    logical(lg),    optional,          intent(in)    :: PARTIAL_MATRIX
    integer(ip)                                      :: ii,jj,iz,jz
    integer(ip)                                      :: ii_loc,jj_loc
    integer(ip)                                      :: ndof1,ndof2
    logical(lg)                                      :: full_matrix
    
    ndof1 = min(self % ndof1,mat_in % ndof1)
    ndof2 = min(self % ndof2,mat_in % ndof2)

    full_matrix = .true.
    if( optional_argument(.false.,PARTIAL_MATRIX) ) full_matrix = .false.
    if( present(invp) ) then
       if( associated(invp) ) full_matrix = .false.
    end if

    if( full_matrix ) then
       iz = mat_in %nz
       self % vA(1:ndof2,1:ndof1,1:iz) = mat_in % vA(1:ndof2,1:ndof1,1:iz)
    else      
       if( present(invp) ) then
          do ii = 1,mat_in % nrows
             ii_loc = invp(ii)
             if( ii_loc > 0 ) then
                do iz = mat_in % iA(ii),mat_in % iA(ii+1)-1
                   jj = mat_in % jA(iz)
                   jj_loc = invp(jj)
                   if( jj_loc > 0 ) then
                      jz = self % edge(ii_loc,jj_loc)
                      if( jz /= 0 ) then
                         self % vA(1:ndof2,1:ndof1,jz) = mat_in % vA(1:ndof2,1:ndof1,iz)                         
                      end if
                   end if
                end do
             end if
          end do
       else
          do ii = 1,mat_in % nrows
             do iz = mat_in % iA(ii),mat_in % iA(ii+1)-1
                jj = mat_in % jA(iz)
                jz = self % edge(ii,jj)
                if( jz /= 0 ) &
                     self % vA(1:ndof2,1:ndof1,jz) = mat_in % vA(1:ndof2,1:ndof1,iz)
             end do
          end do
       end if
    end if  
    
  end subroutine copy
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    17/10/2017     
  !> @brief   Find the edge position in a graph
  !> @details Find the edge position IZ in a graph (IA,JA)
  !
  !-----------------------------------------------------------------------

  integer(ip) pure function edge(self,ii,jj) result(iz)

    class(mat_csr),       intent(in)            :: self      !< Matrix
    integer(ip),          intent(in)            :: ii        !< Vertex 1
    integer(ip),          intent(in)            :: jj        !< Vertex 2
    integer(ip)                                 :: kk

    if( ii == 0 .or. jj == 0 ) then
       iz = 0
    else
       iz = 0
       kk = self % ia(ii)
       do while( kk <= self % ia(ii+1)-1 )
          if( self % ja(kk) == jj ) then
             iz = kk
             return
          end if
          kk = kk + 1
       end do
    end if

  end function edge

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/01/1016
  !> @brief   Transform a sum to a linked list
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  subroutine size_to_linked_list(self)
    
    class(mat_csr),       intent(inout) :: self      !< Matrix
    integer(ip)                         :: nn,ii,kk,ll

    nn    = memory_size(self % ia)-1
    if( nn > 0 ) then
       kk    = self % ia(1)
       self % ia(1) = 1 
       do ii = 2,nn+1
          ll     = self % ia(ii)
          self % ia(ii) = self % ia(ii-1) + kk
          kk     = ll
       end do
    end if
    
  end subroutine size_to_linked_list
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/01/1016
  !> @brief   Transform a sum to a linked list
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  function equal(self,csr,TOLERANCE)
    
    class(mat_csr),           intent(inout) :: self      !< Matrix
    class(mat_csr),           intent(inout) :: csr       !< Matrix
    real(rp),       optional, intent(in)    :: TOLERANCE
    integer(ip)                             :: nn,ii,kk,ll
    real(rp)                                :: toler
    logical(lg)                             :: equal
    
    toler = optional_argument(1.0e-16_rp,TOLERANCE)
    equal = .true.
    if( csr % nrows /= self % nrows ) equal = .false.
    if( csr % nz    /= self % nz    ) equal = .false.
    if( csr % nrows > 0 .and. csr % nz > 0 ) then
       if( any(csr % ia /= self % ia) )         equal = .false.
       if( any(csr % ja /= self % ja) )         equal = .false.
       if( any(abs(csr % va-self % va)>toler) ) equal = .false.
    end if
    
  end function equal
  
end module def_mat_csr
!> @}
