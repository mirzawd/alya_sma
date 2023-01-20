!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_spmat.f90
!> @author  guillaume
!> @date    2021-01-26
!> @brief   Sparse matrix
!> @details Sparse matrix type
!-----------------------------------------------------------------------

module def_spmat

  use def_kintyp_basic,      only : ip,rp,lg
  use mod_optional_argument, only : optional_argument
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end
  implicit none
  private

  !----------------------------------------------------------------------
  !
  ! Sparse matrix 
  !
  !----------------------------------------------------------------------

  character(9), parameter :: vacal = 'def_spmat'
  integer(ip),  parameter :: COO_FORMAT   = 2

  type :: spmat
     integer(ip)          :: kfl_format
     integer(ip)          :: nrows
     integer(ip)          :: ncols
     integer(ip)          :: ndof1
     integer(ip)          :: ndof2
     integer(ip)          :: nz
     integer(ip), pointer :: iA(:)
     integer(ip), pointer :: jA(:)
     real(rp),    pointer :: vA(:,:,:)
   contains
     procedure,   pass    :: init       ! Initialization           
     procedure,   pass    :: alloca     ! Allocate           
     procedure,   pass    :: deallo     ! Deallocate           
     procedure,   pass    :: copy       ! Copy a matrix      
     procedure,   pass    :: move       ! Move a matrix (copy and delete)     
     procedure,   pass    :: merge      ! Merge two matrices     
     procedure,   pass    :: spmv_10
     procedure,   pass    :: spmv_11
     procedure,   pass    :: spmv_22
     generic              :: spmv =>  & ! Sparse matrix-vector product    
          &                  spmv_10, &
          &                  spmv_11, &
          &                  spmv_22
  end type spmat

  public :: spmat
  public :: COO_FORMAT

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    class(spmat), intent(inout) :: self

    self % kfl_format = COO_FORMAT
    self % ndof1      = 1
    self % ndof2      = 1
    self % nrows      = 0
    self % ncols      = 0
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

  subroutine alloca(self,ndof1,ndof2,nz,MEMORY_COUNTER)

    class(spmat),            intent(inout)         :: self
    integer(ip),  optional,  intent(in)            :: ndof1
    integer(ip),  optional,  intent(in)            :: ndof2
    integer(ip),  optional,  intent(in)            :: nz
    integer(8),   optional,  intent(inout)         :: MEMORY_COUNTER(2)
    integer(8)                                     :: memor_loc(2)
    integer(ip)                                    :: nz_loc,ndof1_loc,ndof2_loc

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    nz_loc    = optional_argument(self % nz   ,nz)
    ndof1_loc = optional_argument(self % ndof1,ndof1)
    ndof2_loc = optional_argument(self % ndof2,ndof2)

    if( nz_loc > 0 ) then
       call memory_alloca(memor_loc,'SPMAT % IA',vacal,self % iA,nz_loc)
       call memory_alloca(memor_loc,'SPMAT % JA',vacal,self % jA,nz_loc)
       call memory_alloca(memor_loc,'SPMAT % VA',vacal,self % vA,ndof1_loc,ndof2_loc,nz_loc)        
       self % ndof1 = ndof1_loc
       self % ndof2 = ndof2_loc       
       self % nz    = nz_loc       
    end if
    !if( nz_loc > 0 ) then
    !   allocate(self % iA(nz_loc))
    !   allocate(self % jA(nz_loc))
    !   allocate(self % vA(ndof1_loc,ndof2_loc,nz_loc))   
    !   self % ndof1 = ndof1_loc
    !   self % ndof2 = ndof2_loc    
    !   self % nz    = nz_loc       
    !end if
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

    class(spmat),            intent(inout)         :: self
    integer(8),   optional,  intent(inout)         :: MEMORY_COUNTER(2)
    integer(8)                                     :: memor_loc(2)

    !if(associated(self % iA) ) deallocate(self % iA)
    !if(associated(self % jA) ) deallocate(self % jA)
    !if(associated(self % vA) ) deallocate(self % vA)
    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    call memory_deallo(memor_loc,'SPMAT % IA',vacal,self % iA)
    call memory_deallo(memor_loc,'SPMAT % JA',vacal,self % jA)
    call memory_deallo(memor_loc,'SPMAT % VA',vacal,self % vA)

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Move
  !> @details Move SELF <= SELF2 (copy and deallocate)
  !> 
  !-----------------------------------------------------------------------

  subroutine move(self,self2,MEMORY_COUNTER)

    class(spmat),            intent(inout)         :: self
    class(spmat),            intent(inout)         :: self2
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
  !>          self  % iA(:) = 1,1,1,2,2,3,3,3,3
  !>          self2 % iA(:) = 1,1,2
  !>          => self % iA(:) = 1,1,1,2,2,3,3,3,3,4,4,5
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(self,self2,MEMORY_COUNTER,PERMU)

    class(spmat),                     intent(inout) :: self
    class(spmat),                     intent(in)    :: self2
    integer(8),   optional,           intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),  optional,  pointer, intent(in)    :: permu(:)
    type(spmat)                                     :: self1
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
          do while( self1 % iA(iz) /= ii )
             iz = iz + 1
          end do
          pord1(ii) = iz
       end do
       pord1(self1 % nrows+1) = self1 % nz + 1
       iz = 1
       do ii = 1,self2 % nrows
          do while( self2 % iA(iz) /= ii )
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
                self % iA(kz)     = ii
                self % jA(kz)     = self1 % jA(iz)
                self % vA(:,:,kz) = self1 % vA(:,:,iz)
             end do
          else
             i2 = i2 + 1
             do iz = pord2(i2),pord2(i2+1)-1
                kz                = kz + 1
                self % iA(kz)     = ii
                self % jA(kz)     = self2 % jA(iz)
                self % vA(:,:,kz) = self2 % vA(:,:,iz)
             end do             
          end if
       end do
       deallocate(pord1,pord2)
       
    else
       do iz = 1,self1 % nz
          self % iA(iz)     = self1 % iA(iz) 
          self % jA(iz)     = self1 % jA(iz) 
          self % vA(:,:,iz) = self1 % vA(:,:,iz)
       end do
       iz = self1 % nz
       jz = self1 % nrows
       do kz = 1,self2 % nz
          iz                = iz + 1
          self % iA(iz)     = self2 % iA(kz) + jz
          self % jA(iz)     = self2 % jA(kz)
          self % vA(:,:,iz) = self2 % vA(:,:,kz)
       end do
    end if
    
    call self1 % deallo(MEMORY_COUNTER)

  end subroutine merge

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-26
  !> @brief   Copy
  !> @details Copy SELF = SELF2
  !> 
  !-----------------------------------------------------------------------

  subroutine copy(self,self2,MEMORY_COUNTER)

    class(spmat),            intent(inout) :: self
    class(spmat),            intent(in)    :: self2
    integer(8),   optional,  intent(inout) :: MEMORY_COUNTER(2)

    call self % init()
    self % kfl_format = self2 % kfl_format 
    self % nrows      = self2 % nrows      
    self % ncols      = self2 % ncols      
    self % ndof1      = self2 % ndof1      
    self % ndof2      = self2 % ndof2      
    self % nz         = self2 % nz         

    call self % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    if( associated(self % iA) ) self % iA = self2 % iA
    if( associated(self % jA) ) self % jA = self2 % jA
    if( associated(self % vA) ) self % vA = self2 % vA

  end subroutine copy

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------

  subroutine spmv_10(self,xx,yy,POINT)

    class(spmat),                          intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx(:)
    real(rp),                               intent(inout) :: yy
    integer(ip),          optional,         intent(in)    :: POINT
    integer(ip)                                           :: ii,iz,ipoin

    if( associated(xx) ) then
       if( present(POINT) ) then
          yy = 0.0_rp
          do iz = 1,self % nz
             ii    = self % iA(iz)
             ipoin = self % jA(iz)
             if( ii == POINT ) then
                yy = yy + self % vA(1,1,iz) * xx(ipoin)
             end if
          end do
       else
          call runend('DEF_SPMAT_METHOD: FOR SINGLE POINT SPMAT, PRESCRIBE POINT')
       end if
    end if

  end subroutine spmv_10
  
  subroutine spmv_11(self,xx,yy,POINT)
    
    class(spmat),                           intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx(:)
    real(rp),                      pointer, intent(inout) :: yy(:)
    integer(ip),          optional,         intent(in)    :: POINT
    integer(ip)                                           :: ii,iz,ipoin,ns

    if( associated(yy) ) then
       if( present(POINT) ) then
          yy(1) = 0.0_rp
          do iz = 1,self % nz
             ii    = self % iA(iz)
             ipoin = self % jA(iz)
             if( ii == POINT ) then
                yy(1) = yy(1) + self % vA(1,1,iz) * xx(ipoin)
             end if
          end do
       else
          ns = int(size(yy,1_ip),ip)
          do ii = 1,ns
             yy(ii) = 0.0_rp
          end do
          do iz = 1,self % nz
             ii    = self % iA(iz)
             ipoin = self % jA(iz)
             if( ii <= ns ) then 
                yy(ii) = yy(ii) + self % vA(1,1,iz) * xx(ipoin)
             end if
          end do
       end if
    end if
    
  end subroutine spmv_11
  
  subroutine spmv_22(self,xx,yy,POINT)

    class(spmat),                           intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx(:,:)
    real(rp),                      pointer, intent(inout) :: yy(:,:)
    integer(ip),          optional,         intent(in)    :: POINT
    integer(ip)                                           :: ii,iz,ipoin,nd,ns

    if( associated(xx) .and. associated(yy) ) then
       nd = min(size(xx,1_ip),size(yy,1_ip))       
    else if( associated(xx) ) then
       nd = size(xx,1_ip)
    else
       nd = 0
    end if
    if( associated(yy) ) then
       ns = size(yy,2_ip)
    else
       ns = 0
    end if

    if( associated(yy) ) then
       do ii = 1,ns
          yy(:,ii) = 0.0_rp
       end do
       if( present(POINT) ) then
          do iz = 1,self % nz
             ii    = self % iA(iz)
             ipoin = self % jA(iz)
             if( ii == POINT ) then
                yy(1:nd,1) = yy(1:nd,1) + self % vA(1,1,iz) * xx(1:nd,ipoin)
             end if
          end do
       else
          do iz = 1,self % nz
             ii    = self % iA(iz)
             ipoin = self % jA(iz)
             if( ii <= ns ) then
                yy(1:nd,ii) = yy(1:nd,ii) + self % vA(1,1,iz) * xx(1:nd,ipoin)
             end if
          end do
       end if
    end if
    
  end subroutine spmv_22

  function spmv_log(ii,ns,POINT) result(res)
    integer(ip),           intent(in) :: ii       
    integer(ip),           intent(in) :: ns      
    integer(ip), optional, intent(in) :: POINT
    integer(ip)                       :: res(2)
    
    if( present(POINT) ) then
       if( ii == POINT ) then
          res(1:2) = 1
       else
          res(1)   = 0
       end if
    else if( ii <= ns ) then
       res(1) = 1
       res(2) = ii
    else
       res(1) = 0
    end if
    
  end function spmv_log
  
end module def_spmat
!> @}
