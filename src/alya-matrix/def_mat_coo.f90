!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_coo.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   COO Matrix
!> @details COO Matrix class
!>          Required to allocate:
!>          SELF % NZ
!>          SELF % NDOF1
!>          SELF % NDOF2
!>          SELF % XA(SELF % NZ)
!>          SELF % YA(SELF % NZ)
!>          SELF % VA(SELF % NDOF1,SELF % NDOF2,SELF % NZ)
!-----------------------------------------------------------------------

module def_mat_coo

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat,               only : COO_FORMAT
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
  
  type, extends(mat) :: mat_coo
     integer(ip)                      :: nz
     integer(ip),             pointer :: xA(:)
     integer(ip),             pointer :: yA(:)  
     integer(ip),             pointer :: yA_perm(:)
     integer(ip),             pointer :: yA_conc(:)     
     real(rp),    contiguous, pointer :: vA(:,:,:)
   contains
     procedure,               pass    :: init           ! Initialize
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix           
     procedure,               pass    :: diag           ! Diagonal matrix    
     procedure,               pass    :: diagz          ! Diagonal matrix position    
     procedure,               pass    :: norm           ! Matrix norm
     procedure,               pass    :: dirichlet      ! Prescribe a Dirichlet condition
     procedure,               pass    :: get_val        ! Get values i,j
     procedure,               pass    :: output         ! Output 
     procedure,               pass    :: get_nrows      ! Get rows      
     procedure,               pass    :: get_ncols      ! Get columns
     procedure,               pass    :: copy           ! Copy a matrix      
     procedure,               pass    :: move           ! Move a matrix (copy and delete)    
     procedure,               pass    :: compute_perm_and_conc 
     procedure,               pass    :: merge          ! Merge two matrices     
     procedure,               pass    :: mv_row         ! Single row MV product
     procedure,               pass    :: mv_lower       ! (L+D) product
     procedure,               pass    :: mv_upper       ! (U+D) product
     procedure,               pass    :: mv_11
     procedure,               pass    :: mv_12
     procedure,               pass    :: mv_22
     procedure,               pass    :: residual_1
     procedure,               pass    :: residual_2
  end type mat_coo

  character(11), parameter :: vacal = 'def_mat_coo'
  
  public :: mat_coo
  
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

    class(mat_coo), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = COO_FORMAT
    self % nz         = 0
    self % ndof1      = 1
    self % ndof2      = 1
    nullify(self % xA)
    nullify(self % yA)
    nullify(self % yA_perm)
    nullify(self % yA_conc)
    nullify(self % vA)
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Assign an array matrix to mat class vA
  !> @details Assign a rank-1 array to mat class vA
  !> 
  !-----------------------------------------------------------------------
    
  subroutine assign(self,a)
    class(mat_coo),                      intent(inout) :: self
    real(rp),       pointer, contiguous, intent(in)    :: a(:)

    self % vA(1:self % ndof1,1:self % ndof2,1:self % nz) => a
    
  end subroutine assign

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
    
  subroutine alloca(self,param,MEMORY_COUNTER)
    
    class(mat_coo),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: nz_loc,ndof1_loc,ndof2_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nz_loc    = optional_argument(self % nz   ,param,1_ip)
    ndof1_loc = optional_argument(self % ndof1,param,2_ip)
    ndof2_loc = optional_argument(self % ndof2,param,3_ip)
    call memory_alloca(memor_loc,'SELF % XA',vacal,self % xA,nz_loc)
    call memory_alloca(memor_loc,'SELF % YA',vacal,self % yA,nz_loc)
    call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,ndof2_loc,nz_loc)    
    self % ndof1 = ndof1_loc
    self % ndof2 = ndof2_loc       
    self % nz    = nz_loc       
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
    
    class(mat_coo),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % XA',vacal,self % xA)
    call memory_deallo(memor_loc,'SELF % YA',vacal,self % yA)
    call memory_deallo(memor_loc,'SELF % VA',vacal,self % vA)
    if(associated(self%ya_conc)) then
      call memory_deallo(memor_loc,'SELF % YA_CONC',vacal, self%ya_conc)
    end if
    if(associated(self%ya_perm)) then
      call memory_deallo(memor_loc,'SELF % YA_PERM',vacal, self%ya_perm)
    end if
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Get number of rows
  !> @details Get number of rows
  !> 
  !-----------------------------------------------------------------------
    
  pure function get_nrows(self) result(nrows)
    class(mat_coo), intent(in) :: self
    integer(ip)                :: nrows

    if( self % nrows > 0 ) then
       nrows = self % nrows
    else
       if( associated(self % xA) ) then
          nrows = maxval(self % xA)
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
    class(mat_coo), intent(in) :: self
    integer(ip)                :: ncols

    if( self % ncols > 0 ) then
       ncols = self % ncols
    else
       if( associated(self % yA) ) then
          ncols = maxval(self % yA)
       else
          ncols = 0
       end if
    end if
    
  end function get_ncols
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal position
  !> @details Diagonal position in matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diagz(self,diagonal)
    
    class(mat_coo),                    intent(in)    :: self
    integer(ip),    optional, pointer, intent(inout) :: diagonal(:)
    integer(ip)                                      :: ii,iz

    do ii = 1,self % nrows
       iz = 1
       do while( self % xA(iz) /= ii .or. self % yA(iz) /= ii )
          iz = iz + 1
       end do
       diagonal(ii) = iz
    end do
    
  end subroutine diagz
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal
  !> @details Compute the diagonal matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine diag(self,dia,diagonal,row,val)
    
    class(mat_coo),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,iz,idofn

    if( present(diagonal) ) then

       do ii = 1,self % nrows
          iz = 1
          do while( self % xA(iz) /= ii .or. self % yA(iz) /= ii )
             iz = iz + 1
          end do          
          do idofn = 1,self % ndof2
             diagonal(idofn,ii) = self % vA(idofn,idofn,iz)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia)
          
          if( .not. associated(dia % vA) ) then
             dia % nrows = self % nrows
             dia % ndof1 = self % ndof2
             call dia % alloca()
          end if
          
          do ii = 1,self % nrows
             iz = 1
             do while( self % xA(iz) /= ii .or. self % yA(iz) /= ii )
                iz = iz + 1
             end do
             do idofn = 1,self % ndof2
                dia % vA(idofn,ii) = self % vA(idofn,idofn,iz)
             end do
          end do
       end select
       
    else if( present(row) .and. present(val) ) then
       
       ii = row
       iz = 1
       do while( self % xA(iz) /= ii .or. self % yA(iz) /= ii )
          iz = iz + 1
       end do
       do idofn = 1,self % ndof2
          val(idofn) = self % vA(idofn,idofn,iz)          
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
    
    class(mat_coo), intent(in)                      :: self
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
    
    class(mat_coo), intent(in)                      :: self
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
    
    class(mat_coo), intent(in)                      :: self
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
    
    class(mat_coo), intent(in)                      :: self
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
    
    class(mat_coo), intent(in)                      :: self
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

    class(mat_coo), intent(in)                    :: self
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
    integer(ip)                                   :: ii,jj,iz
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

    if( nn1 == 1 .and. nn2 == self % nrows ) then

       !$OMP PARALLEL   DO                                 &
       !$OMP SCHEDULE ( STATIC )                           &
       !$OMP DEFAULT  ( NONE )                             &
       !$OMP SHARED   ( self, nn1, nn2, xx, yy ) &
       !$OMP PRIVATE  ( ii, jj, iz, idofr, idofc )
       do iz = 1,self % nz
          ii = self % xA(iz)
          jj = self % yA(iz)
          do idofc = 1,self % ndof1
             do idofr = 1,self % ndof2
                !$OMP ATOMIC
                yy(idofr,ii) = yy(idofr,ii) + self % vA(idofc,idofr,iz) * xx(idofc,jj)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else      

       !$OMP PARALLEL   DO                                 &
       !$OMP SCHEDULE ( STATIC )                           &
       !$OMP DEFAULT  ( NONE )                             &
       !$OMP SHARED   ( self, nn1, nn2, xx, yy ) &
       !$OMP PRIVATE  ( ii, jj, iz, idofr, idofc )
       do iz = 1,self % nz
          ii = self % xA(iz)
          if( ii >= nn1 .and. ii <= nn2 ) then
             jj = self % yA(iz)
             do idofc = 1,self % ndof1
                do idofr = 1,self % ndof2             
                   !$OMP ATOMIC
                   yy(idofr,ii) = yy(idofr,ii) + self % vA(idofc,idofr,iz) * xx(idofc,jj)
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO
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

    class(mat_coo), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    integer(ip)                   :: ii,jj,iz,ncols
    integer(ip)                   :: idof1,idof2
    real(rp),       pointer       :: aa(:,:)

    anorm = 0.0_rp

    if( self % nz > 0 ) then

       if( wnorm == 'i' .or. wnorm == 'I' ) then

          allocate(aa(self % ndof2,self % nrows))
          do ii = 1,self % nrows
             aa(:,ii) = 0.0_rp
          end do

          do iz = 1,self % nz
             ii = self % xA(iz)
             jj = self % yA(iz)
             do idof2 = 1,self % ndof2
                do idof1 = 1,self % ndof1
                   aa(idof2,ii) = aa(idof2,ii) + abs(self % vA(idof1,idof2,iz))
                end do
             end do
          end do
          do ii = 1,self % nrows
             anorm = max(anorm,maxval(aa(1:self % ndof2,ii)))
          end do
          deallocate(aa)

       else if( wnorm == '1' ) then

          if( self % ncols == 0 ) then
             ncols = maxval(self % yA)
          else
             ncols = self % ncols
          end if
          allocate(aa(self % ndof1,ncols))
          do ii = 1,ncols
             aa(:,ii) = 0.0_rp
          end do

          do iz = 1,self % nz
             ii = self % xA(iz)
             jj = self % yA(iz)
             do idof2 = 1,self % ndof2
                do idof1 = 1,self % ndof1
                   aa(idof1,jj) = aa(idof1,jj) + abs(self % vA(idof1,idof2,iz))
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

    class(mat_coo),                   intent(in) :: self
    character(*),                     intent(in) :: FMT
    character(*),  optional,          intent(in) :: FILENAME
    integer(ip),   optional, pointer, intent(in) :: PERM(:)
    integer(ip)                                  :: iz,ii,jj,kk,nz
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
    unit4 = int(iofile_available_unit(90_ip),4)

    select case ( upper_case(FMT) )

    case ( 'MATRIX MARKET', 'MTX' )

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
       do idof2 = 1,self % ndof2
          do iz = 1,self % nz
             ii = self % xA(iz)
             jj = self % yA(iz)
             do idof1 = 1,self % ndof1
                if( present(PERM) ) then
                   write(unit4,'(2(1x,i9),1x,e13.6)') PERM((ii-1)*self % ndof2+idof2),PERM((jj-1)*self % ndof1+idof1),self % vA(idof1,idof2,iz)
                else
                   write(unit4,'(2(1x,i9),1x,e13.6)') (ii-1)*self % ndof2+idof2,(jj-1)*self % ndof1+idof1,self % vA(idof1,idof2,iz)
                end if
             end do
          end do
       end do
       close(unit=unit4)

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
                iz    = 0
                notfound = .true.
                do while( notfound .and. iz < self % nz )
                   iz = iz + 1
                   if( self % xA(iz) == ii .and. self % yA(iz) == jj ) notfound = .false.
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
  !> @date    2021-01-26
  !> @brief   Copy
  !> @details Copy SELF = SELF2
  !> 
  !-----------------------------------------------------------------------

  subroutine copy(self,self2,MEMORY_COUNTER)

    class(mat_coo),          intent(inout) :: self
    class(mat_coo),          intent(in)    :: self2
    integer(8),   optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                            :: iz
    
    call self % init()
    self % kfl_format = self2 % kfl_format 
    self % nrows      = self2 % nrows      
    self % ncols      = self2 % ncols      
    self % ndof1      = self2 % ndof1      
    self % ndof2      = self2 % ndof2      
    self % nz         = self2 % nz         

    call self % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    if( associated(self % xA) ) then
       do iz = 1,self % nz
          self % xA(iz) = self2 % xA(iz)
       end do
    end if
    
    if( associated(self % yA) ) then
       do iz = 1,self % nz
          self % yA(iz) = self2 % yA(iz)
       end do
    end if
    
    if( associated(self % vA) ) then
       do iz = 1,self % nz
          self % vA(:,:,iz) = self2 % vA(:,:,iz)
       end do
    end if
    
  end subroutine copy
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Move
  !> @details Move SELF <= SELF2 (copy and deallocate)
  !> 
  !-----------------------------------------------------------------------

  subroutine move(self,self2,MEMORY_COUNTER)

    class(mat_coo),          intent(inout)         :: self
    class(mat_coo),          intent(inout)         :: self2
    integer(8),   optional,  intent(inout)         :: MEMORY_COUNTER(2)

    call self  % copy(self2,MEMORY_COUNTER)
    call self2 % deallo(MEMORY_COUNTER)

  end subroutine move

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Merge
  !> @details Merge SELF <= SELF U SELF2
  !>          Matrix should be ordered lexicographically in rows
  !>          self  % xA(:) = 1,1,1,2,2,3,3,3,3
  !>          self2 % xA(:) = 1,1,2
  !>          => self % xA(:) = 1,1,1,2,2,3,3,3,3,4,4,5
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(self,self2,MEMORY_COUNTER,PERMU)

    class(mat_coo),                   intent(inout) :: self
    class(mat_coo),                   intent(in)    :: self2
    integer(8),   optional,           intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),  optional,  pointer, intent(in)    :: permu(:)
    type(mat_coo)                                   :: self1
    integer(ip)                                     :: iz,kz,jz,nn
    integer(ip)                                     :: ii,jj,i1,i2
    integer(ip),             pointer                :: pord1(:)
    integer(ip),             pointer                :: pord2(:)

    if( self % ndof1      /= self2 % ndof1      ) stop
    if( self % ndof2      /= self2 % ndof2      ) stop
    if( self % kfl_format /= self2 % kfl_format ) stop

    
    call self1 % init()
    call self1 % move(self,MEMORY_COUNTER)
    call self  % init()

    self % kfl_format = self1 % kfl_format 
    self % ndof1      = self1 % ndof1
    self % ndof2      = self1 % ndof2
    self % nz         = self1 % nz    + self2 % nz
    self % nrows      = self1 % nrows + self2 % nrows
    self % ncols      = max(self1 % ncols,self2 % ncols)

    call self % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    if( present(permu) ) then
       nn = memory_size(permu)
    else
       nn = 0
    end if
    
    if( nn > 0 .and. present(permu) ) then
       allocate(pord1(self1 % nrows+1))
       allocate(pord2(self2 % nrows+1))
       iz = 1
       do ii = 1,self1 % nrows
          do while( self1 % xA(iz) /= ii )
             iz = iz + 1
          end do
          pord1(ii) = iz
       end do
       pord1(self1 % nrows+1) = self1 % nz + 1
       iz = 1
       do ii = 1,self2 % nrows
          do while( self2 % xA(iz) /= ii )
             iz = iz + 1
          end do
          pord2(ii) = iz
       end do       
       pord2(self2 % nrows+1) = self2 % nz + 1

       kz = 0
       i1 = 0
       i2 = 0
       do ii = 1,self % nrows
          jj = 1
          do while( permu(jj) /= ii )
             jj = jj + 1
          end do
          if( jj <= self1 % nrows ) then
             i1 = i1 + 1
             do iz = pord1(i1),pord1(i1+1)-1
                kz                = kz + 1
                self % xA(kz)     = ii
                self % yA(kz)     = self1 % yA(iz)
                self % vA(:,:,kz) = self1 % vA(:,:,iz)
             end do
          else
             i2 = i2 + 1
             do iz = pord2(i2),pord2(i2+1)-1
                kz                = kz + 1
                self % xA(kz)     = ii
                self % yA(kz)     = self2 % yA(iz)
                self % vA(:,:,kz) = self2 % vA(:,:,iz)
             end do             
          end if
       end do
       deallocate(pord1,pord2)
       
    else
       do iz = 1,self1 % nz
          self % xA(iz)     = self1 % xA(iz) 
          self % yA(iz)     = self1 % yA(iz) 
          self % vA(:,:,iz) = self1 % vA(:,:,iz)
       end do
       iz = self1 % nz
       jz = self1 % nrows
       do kz = 1,self2 % nz
          iz                = iz + 1
          self % xA(iz)     = self2 % xA(kz) + jz
          self % yA(iz)     = self2 % yA(kz)
          self % vA(:,:,iz) = self2 % vA(:,:,kz)
       end do
    end if
    
    call self1 % deallo(MEMORY_COUNTER)

  end subroutine merge

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Scale
  !> @details Scale a matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine scale(self,dia,diagonal1,diagonal)

    class(mat_coo),                    intent(inout) :: self
    class(mat_dia), optional,          intent(in)    :: dia
    real(rp),       optional, pointer, intent(in)    :: diagonal1(:)
    real(rp),       optional, pointer, intent(in)    :: diagonal(:,:)
    integer(ip)                                      :: iz,ii,idofn
    integer(ip)                                      :: ndof1,ndof2

    ndof1 = self % ndof1
    ndof2 = self % ndof2

    if( present(dia) ) then
       !
       ! Diagonal is of type MAT_DIA
       !
       if( dia % ndof1 == self % ndof2 ) then
          do iz = 1,self % nz
             ii = self % xA(iz)
             do idofn = 1,self % ndof2
                self % vA(1:ndof1,idofn,iz) = self % vA(1:ndof1,idofn,iz) / dia % vA(idofn,ii)
             end do
          end do
       else if( dia % ndof1 == 1 ) then
          do iz = 1,self % nz
             ii = self % xA(iz)
             do idofn = 1,self % ndof2
                self % vA(1:ndof1,idofn,iz) = self % vA(1:ndof1,idofn,iz) / dia % vA(1,ii)
             end do
          end do
       end if
       
    else if( present(diagonal1) ) then
       !
       ! Diagonal is of type (:)
       !
       do iz = 1,self % nz
          ii = self % xA(iz)
          do idofn = 1,self % ndof2
             self % vA(1:ndof1,idofn,iz) = self % vA(1:ndof1,idofn,iz) / diagonal1(ii)
          end do
       end do
       
    else if( present(diagonal) ) then
       !
       ! Diagonal is of type (:,:)
       !
       do iz = 1,self % nz
          ii = self % xA(iz)
          do idofn = 1,self % ndof2
             self % vA(1:ndof1,idofn,iz) = self % vA(1:ndof1,idofn,iz) / diagonal(idofn,ii)
          end do
       end do

    end if
    
  end subroutine scale

  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2021-01-26
  !> @brief   Scale
  !> @details Compute ya_perm
  !> 
  !-----------------------------------------------------------------------

  subroutine compute_perm_and_conc(self,MEMORY_COUNTER)


    class(mat_coo), intent(inout) :: self      
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                   :: nz,size_conc,kk,ii,buff,jj
    integer(ip),    pointer       :: ya_buff(:) 
    logical                       :: doub
    integer(8)                    :: memor_loc(2)

    
    

    nullify(self %yA_perm)
    nullify(self %yA_conc)
    nullify(ya_buff)
    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)


    nz = memory_size(self%xA)          
    call memory_alloca(memor_loc,'ya_buff',vacal,ya_buff,nz)
    call memory_alloca(memor_loc,'self%yA_perm',vacal,self%yA_perm,nz)
    ya_buff(:) = -1
    size_conc = 0
    do ii = 1,nz
      buff = self % yA(ii)
      kk = 1
      doub = .false.
      do while (.not.doub .and. kk <= nz )
        if(buff == ya_buff(kk)) then
          doub =.true.
        end if
        kk=kk + 1
      end do
      if(.not.doub) then
        size_conc = size_conc + 1
        ya_buff(size_conc) = buff
        do jj = 1,nz
          if(self %yA(jj) == buff) then
            self %yA_perm(jj) =  size_conc
          end if
        end do
      end if
    end do
    

    if(size_conc>0) then
      call memory_alloca(memor_loc,'ya_conc',vacal,self%yA_conc,size_conc)
      self%yA_conc(1:size_conc) = ya_buff(1:size_conc)
    end if
    
    call memory_deallo(memor_loc,'ya_buff',vacal,ya_buff)    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine compute_perm_and_conc

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume  
  !> @date    2021-04-12
  !> @brief   Mv on a single row
  !> @details Mv on a single row
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function mv_row(self,xx,n,ndof) result(yy)
    
    class(mat_coo),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: iz,ll,jj

    yy = 0.0_rp
    do iz = 1,self % nz
       if( self % xA(iz) == n ) then
          jj = self % yA(iz)
          do ll = 1,self % ndof2
             yy = yy + self % vA(ll,ndof,iz) * xx((jj-1)*self % ndof2+ll)
          end do
       end if
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
    
    class(mat_coo), intent(in)                      :: self
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

    do iz = 1,self % nz
       ii = self % xA(iz)
       if( ii >= nn1 .and. ii <= nn2 ) then
          jj = self % yA(iz)
          if( jj <= ii ) then
             do id = 1,self % ndof1
                it = (ii-1) * self % ndof1 + id
                do jd = 1,id
                   jt = (jj-1) * self % ndof2+jd
                   yy(it) = yy(it) + self % vA(jd,id,iz) * xx(jt)
                end do
             end do
          end if
       end if
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
    
    class(mat_coo), intent(in)                      :: self
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
 
    do iz = 1,self % nz
       ii = self % xA(iz)
       if( ii >= nn1 .and. ii <= nn2 ) then
          jj = self % yA(iz)
          if( jj >= ii ) then
             do id = 1,self % ndof1
                it = (ii-1) * self % ndof1 + id
                do jd = id,self % ndof2
                   jt = (jj-1) * self % ndof2+jd
                   yy(it) = yy(it) + self % vA(jd,id,iz) * xx(jt)
                end do
             end do
          end if
       end if
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
     
    class(mat_coo), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    integer(ip)                :: iz

    a = 0.0_rp
    
   if( i <= self % nrows ) then
       do iz = 1,self % nrows
          if( self % xa(iz) == i ) then
             if( self % ya(iz) == j ) then
                a(:,:) = self % vA(:,:,iz)
                return
             end if
          end if
       end do
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
    
    class(mat_coo),                    intent(inout) :: self
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

end module def_mat_coo
!> @}
