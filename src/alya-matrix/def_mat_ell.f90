!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_mat_ell.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   ELL Matrix
!> @details ELL Matrix class
!>          Required to allocate:
!>          SELF % NDOF1
!>          SELF % NDOF2
!>          SELF % NROWS
!>          SELF % COMAX
!>          SELF % JA(SELF % COMAX,SELF % NROWS)
!>          SELF % VA(SELF % NDOF1,SELF % NDOF2,SELF % COMAX,SELF % NROWS)
!-----------------------------------------------------------------------

module def_mat_ell

  use def_kintyp_basic,      only : ip,rp,lg
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat,               only : ELL_FORMAT
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

  type, extends(mat) :: mat_ell
     integer(ip)                      :: comax
     integer(ip),             pointer :: ja(:,:)        ! ja(nrows,comax)
     real(rp),    contiguous, pointer :: vA(:,:,:,:)
   contains
     procedure,               pass    :: init           ! Initialize the class
     procedure,               pass    :: alloca         ! Allocate   
     procedure,               pass    :: deallo         ! Deallocate           
     procedure,               pass    :: assign         ! Assign a rank-1 matrix 
     procedure,               pass    :: diag           ! Compute diagonal
     procedure,               pass    :: diagz          ! Compute diagonal position
     procedure,               pass    :: set_nrows      ! Set rows      
     procedure,               pass    :: get_nrows      ! Get rows      
     procedure,               pass    :: set_ncols      ! Set columns      
     procedure,               pass    :: get_ncols      ! Get columns
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
  end type mat_ell

  character(11), parameter :: vacal = 'def_mat_ell'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  
  public :: mat_ell
  
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

    class(mat_ell), intent(inout) :: self

    call self % init_mat()

    self % kfl_format = ELL_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    self % nrows      = 0
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
    
    class(mat_ell),            intent(inout) :: self
    integer(ip),    optional,  intent(in)    :: param(:)
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)
    integer(ip)                              :: ndof1_loc,ndof2_loc
    integer(ip)                              :: nrows_loc,comax_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    ndof1_loc = optional_argument(self % ndof1,param,1_ip)
    ndof2_loc = optional_argument(self % ndof2,param,2_ip)
    nrows_loc = optional_argument(self % nrows,param,3_ip)
    comax_loc = optional_argument(self % comax,param,4_ip)
    
    call memory_alloca(memor_loc,'SELF % JA',vacal,self % jA,comax_loc,nrows_loc)
    call memory_alloca(memor_loc,'SELF % VA',vacal,self % vA,ndof1_loc,ndof2_loc,comax_loc,nrows_loc)        
    self % ndof1 = ndof1_loc
    self % ndof2 = ndof2_loc       
    self % nrows = nrows_loc       
    self % comax = comax_loc       
    
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
    
    class(mat_ell),            intent(inout) :: self
    integer(8),     optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                               :: memor_loc(2)

     memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
   
    call memory_deallo(memor_loc,'SELF % JA',vacal,self % jA)
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
    class(mat_ell),                      intent(inout) :: self
    real(rp),      pointer,  contiguous, intent(in)    :: a(:)

    self % vA(1:self % ndof1,1:self % ndof2,1:self % nrows,1:self % comax) => a
    
  end subroutine assign

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-08
  !> @brief   Set number of rows
  !> @details Set number of rows
  !> 
  !-----------------------------------------------------------------------
    
  pure subroutine set_nrows(self)
    class(mat_ell), intent(inout) :: self

    if( associated(self % jA) ) then
       self % nrows = size(self % jA,2)
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
    
  pure subroutine set_ncols(self)
    class(mat_ell), intent(inout) :: self

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
    class(mat_ell), intent(in)    :: self
    integer(ip)                   :: nrows

    if( self % nrows > 0 ) then
       nrows = self % nrows
    else
       if( associated(self % jA) ) then
          nrows = size(self % jA,1)
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
    class(mat_ell), intent(in)    :: self
    integer(ip)                   :: ncols

    if( self % ncols > 0 ) then
       ncols = self % ncols
    else
       if( associated(self % vA) ) then
          ncols = maxval(self % jA)
       else
          ncols = 0
       end if
    end if
    
  end function get_ncols

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-09
  !> @brief   y = b - Ax
  !> @details Compute residual y = b - Ax
  !> 
  !-----------------------------------------------------------------------
  
  subroutine residual_1(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
     
    class(mat_ell), intent(in)                      :: self
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
       ndofc = self % ndof1
       ndofr = self % ndof2
       call mv_go(self,xx,yy,n1,n2,ndofr,ndofc,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
       do ii = 1,self % nrows * ndofr
          yy(ii) = bb(ii) - yy(ii)
       end do
    end if
    
  end subroutine residual_1
  
  subroutine residual_2(self,xx,yy,bb,n1,n2,INITIALIZATION,OPENMP,OPENMP_CHUNK,OPENMP_SCHEDULE)
    
    class(mat_ell), intent(in)                      :: self
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
    
    class(mat_ell), intent(in)                      :: self
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
    
    class(mat_ell), intent(in)                      :: self
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
    
    class(mat_ell), intent(in)                      :: self
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

    class(mat_ell), intent(in)                    :: self
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
    integer(ip)                                   :: col
    integer(ip)                                   :: idof,jdof
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

    if( ndofr == 1 .and. ndofc == 1 ) then
       if( use_openmp ) then
          if( my_schedule == OMP_STATIC ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( STATIC )                        &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( self, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj )
             do ii = n1,n2
                do col = 1,self % comax
                   jj       = self % jA(col,ii)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,col,ii) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == OMP_GUIDED ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( GUIDED )                        &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( self, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj )
             do ii = n1,n2
                do col = 1,self % comax
                   jj       = self % jA(col,ii)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,col,ii) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == OMP_DYNAMIC ) then
             !$OMP PARALLEL   DO                              &
             !$OMP SCHEDULE ( DYNAMIC , my_chunk )            &
             !$OMP DEFAULT  ( NONE )                          &
             !$OMP SHARED   ( self, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj )
             do ii = n1,n2
                do col = 1,self % comax
                   jj       = self % jA(col,ii)
                   yy(1,ii) = yy(1,ii) + self % vA(1,1,col,ii) * xx(1,jj)
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          do ii = n1,n2
             do col = 1,self % comax
                jj       = self % jA(col,ii)
                yy(1,ii) = yy(1,ii) + self % vA(1,1,col,ii) * xx(1,jj)
             end do
          end do
       end if
    else
       if( use_openmp ) then
          if( my_schedule == OMP_STATIC ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( STATIC )                                      &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( self, ndofc, ndofr, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj, idof, jdof )
             do ii = n1,n2
                do col = 1,self % comax
                   jj = self % jA(col,ii)
                   do jdof = 1,ndofc
                      do idof = 1,ndofr
                         yy(idof,ii) = yy(idof,ii) + self % vA(jdof,idof,col,ii) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == OMP_GUIDED ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( GUIDED )                                      &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( self, ndofc, ndofr, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj, idof, jdof )
             do ii = n1,n2
                do col = 1,self % comax
                   jj = self % jA(col,ii)
                   do jdof = 1,ndofc
                      do idof = 1,ndofr
                         yy(idof,ii) = yy(idof,ii) + self % vA(jdof,idof,col,ii) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          else if( my_schedule == OMP_DYNAMIC ) then
             !$OMP PARALLEL   DO                                            &
             !$OMP SCHEDULE ( DYNAMIC , my_chunk )                          &
             !$OMP DEFAULT  ( NONE )                                        &
             !$OMP SHARED   ( self, ndofc, ndofr, n1, n2, xx, yy          ) &
             !$OMP PRIVATE  ( col, ii, jj, idof, jdof )
             do ii = n1,n2
                do col = 1,self % comax
                   jj = self % jA(col,ii)
                   do jdof = 1,ndofc
                      do idof = 1,ndofr
                         yy(idof,ii) = yy(idof,ii) + self % vA(jdof,idof,col,ii) * xx(jdof,jj)
                      end do
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          do ii = n1,n2
             do col = 1,self % comax
                jj = self % jA(col,ii)
                do jdof = 1,ndofc
                   do idof = 1,ndofr
                      yy(idof,ii) = yy(idof,ii) + self % vA(jdof,idof,col,ii) * xx(jdof,jj)
                   end do
                end do
             end do
          end do
       end if
    end if

  end subroutine mv_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Diagonal position
  !> @details Diagonal position in matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diagz(self,diagonal)
    
    class(mat_ell),          intent(in)    :: self
    integer(ip),    pointer, intent(inout) :: diagonal(:)
    integer(ip)                            :: ii,jj

    do ii = 1,self % nrows
       jj = 1          
       do while( self % jA(jj,ii) /= ii )
          jj = jj + 1
       end do
       diagonal(ii) = jj
    end do
    
  end subroutine diagz
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Dirichlet
  !> @details Impose a Dirichlet condition on a matrix
  !> 
  !-----------------------------------------------------------------------

  pure subroutine diag(self,dia,diagonal,row,val)
    class(mat_ell),                    intent(in)    :: self
    class(*),       optional,          intent(inout) :: dia
    real(rp),       optional, pointer, intent(inout) :: diagonal(:,:)
    integer(ip),    optional,          intent(in)    :: row
    real(rp),       optional,          intent(out)   :: val(:)
    integer(ip)                                      :: ii,jj,idofn
    
    if( present(diagonal) ) then
       
       do ii = 1,self % nrows
          jj = 1          
          do while( self % jA(jj,ii) /= ii )
             jj = jj + 1
          end do
          do idofn = 1,self % ndof2
             diagonal(idofn,ii) = self % vA(idofn,idofn,jj,ii)
          end do
       end do
       
    else if( present(dia) ) then
       
       select type ( dia )
       class is ( mat_dia)
          do ii = 1,self % nrows
             jj = 1          
             do while( self % jA(jj,ii) /= ii )
                jj = jj + 1
             end do
             do idofn = 1,self % ndof2
                dia % vA(idofn,ii) = self % vA(idofn,idofn,jj,ii)
             end do
          end do
       end select
       
    else if( present(row) .and. present(val) ) then
       
       ii = row
       jj = 1          
       do while( self % jA(jj,ii) /= ii )
          jj = jj + 1
       end do
       do idofn = 1,self % ndof2
          val(idofn) = self % vA(idofn,idofn,jj,ii)
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

    class(mat_ell), intent(in)    :: self
    character(1),   intent(in)    :: wnorm
    real(rp)                      :: anorm
    integer(ip)                   :: ii,jj,col,ncols
    integer(ip)                   :: idof1,idof2
    real(rp),       pointer       :: aa(:,:)

    anorm = 0.0_rp

    if( self % nrows > 0 ) then

       if( wnorm == 'i' .or. wnorm == 'I' ) then

          allocate(aa(self % ndof2,self % nrows))
          do ii = 1,self % nrows
             aa(:,ii) = 0.0_rp
          end do

          do ii = 1,self % nrows
             do col = 1,self % comax
                jj = self % jA(col,ii)
                do idof2 = 1,self % ndof2
                   do idof1 = 1,self % ndof1
                      aa(idof2,ii) = aa(idof2,ii) + abs(self % vA(idof1,idof2,col,ii))
                   end do
                end do
             end do             
          end do
          do ii = 1,self % nrows
             anorm = max(anorm,maxval(aa(1:self % ndof2,ii)))
          end do
          deallocate(aa)

       else if( wnorm == '1' ) then

          ncols = self % get_ncols()
          allocate(aa(self % ndof1,ncols))
          do ii = 1,ncols
             aa(:,ii) = 0.0_rp
          end do

          do ii = 1,self % nrows
             do col = 1,self % comax
                jj = self % jA(col,ii)
                do idof2 = 1,self % ndof2
                   do idof1 = 1,self % ndof1
                      aa(idof1,jj) = aa(idof1,jj) + abs(self % vA(idof1,idof2,col,ii))
                   end do
                end do
             end do
          end do
          do ii = 1,ncols
             anorm = max(anorm,maxval(aa(1:self % ndof1,ii)))
          end do
          deallocate(aa)

       else if( wnorm == 'f' .or. wnorm == 'F' ) then

          do ii = 1,self % nrows
             do col = 1,self % comax
                jj = self % jA(col,ii)
                do idof2 = 1,self % ndof2
                   do idof1 = 1,self % ndof1
                      anorm = anorm + self % vA(idof1,idof2,col,ii)*self % vA(idof1,idof2,col,ii)
                   end do
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
    
    class(mat_ell),                    intent(in)    :: self
    real(rp),                 pointer, intent(in)    :: xx(:)                   !< Input vector
    integer(ip),                       intent(in)    :: n                       !< Node
    integer(ip),                       intent(in)    :: ndof                    !< dof
    integer(ip)                                      :: ll,jj,col

    yy = 0.0_rp
    do col = 1,self % comax
       jj = self % jA(col,n)
       do ll = 1,self % ndof2
          yy = yy + self % vA(ll,ndof,col,n) * xx((jj-1)*self % ndof2 + ll)
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
    
    class(mat_ell), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,col
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if

    do ii = nn1,nn2
       do col = 1,self % comax
          jj = self % jA(col,ii)
          if( jj <= ii ) then
             do id = 1,self % ndof1
                do jd = 1,id
                   it     = (ii-1) * self % ndof1 + id
                   jt     = (jj-1) * self % ndof2 + jd
                   yy(it) = yy(it) + self % vA(jd,id,col,ii) * xx(jt)
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
    
    class(mat_ell), intent(in)                      :: self
    real(rp),       intent(in),    pointer          :: xx(:)                   !< Input vector
    real(rp),       intent(inout), pointer          :: yy(:)                   !< Output vector
    integer(ip),    intent(in),            optional :: n1                      !< Starting node
    integer(ip),    intent(in),            optional :: n2                      !< Final node
    logical(lg),    intent(in),            optional :: INITIALIZATION          !< If array should be initialized
    logical(lg)                                     :: do_initialize

    integer(ip)                                     :: ii,jj,nn1,nn2
    integer(ip)                                     :: id,jd
    integer(ip)                                     :: it,jt,col
    
    do_initialize = optional_argument(.true. ,INITIALIZATION)
    nn1           = optional_argument(1_ip,n1)
    nn2           = optional_argument(self % nrows,n2)
    
    if( do_initialize ) then
       do ii = (nn1-1) * self % ndof1+1,nn2 * self % ndof2
          yy(ii) = 0.0_rp
       end do
    end if
 
    do ii = nn1,nn2
       do col = 1,self % comax
          jj = self % jA(col,ii)
          if( jj >= ii ) then
             do id = 1,self % ndof1
                do jd = id,self % ndof2
                   it     = (ii-1) * self % ndof1 + id
                   jt     = (jj-1) * self % ndof2 + jd
                   yy(it) = yy(it) + self % vA(jd,id,col,ii) * xx(jt)
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
     
    class(mat_ell), intent(in) :: self
    integer(ip),    intent(in) :: i
    integer(ip),    intent(in) :: j
    real(rp)                   :: a(self % ndof2,self % ndof1)
    integer(ip)                :: col

    a = 0.0_rp
    
    if( i <= self % nrows ) then
       do col = 1,self % comax
          if( self % jA(col,i) == j ) then            
             a(:,:) = self % vA(:,:,col,i)
             return
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
    
    class(mat_ell),                    intent(inout) :: self
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
  
end module def_mat_ell
!> @}
