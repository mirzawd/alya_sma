!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical solvers
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths_arrays.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for arrays
!> @details ToolBox for arrays
!
!-----------------------------------------------------------------------

module mod_maths_arrays

  use def_maths
  implicit none

  interface maths_findloc   
     module procedure &
          maths_findloc_0,&
          maths_findloc_1          
  end interface maths_findloc

  public :: maths_maxloc_nonzero                ! Last non-zero position
  public :: maths_findloc                       ! First position of value in array
  public :: maths_copy_with_perm                ! Copy and permute arrays
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-05-19
  !> @brief   Return position of last non-zero value
  !> @details Examples:
  !>          [ 1 5 3 8 0 0 ] => 4
  !>          [ 0 0 0 ]       => 0
  !>          [ 9 8 3 6 7 ]   => 5
  !>
  !-----------------------------------------------------------------------

  integer(ip) pure function maths_maxloc_nonzero(ll)

    integer(ip), intent(in) :: ll(:)
    integer(ip)             :: ii(1)
    integer(ip)             :: nn

    nn = size(ll,KIND=ip)
    if( ll(nn) /= 0 ) then
       maths_maxloc_nonzero = nn
    else
       ii = minloc(ll,ll==0)    
       if( ii(1) > 0 ) then
          maths_maxloc_nonzero = ii(1)-1
       else if( ll(1) > 0 ) then
          maths_maxloc_nonzero = nn
       else
          maths_maxloc_nonzero = 0
       end if
    end if
    
  end function maths_maxloc_nonzero
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2021-02-15
  !> @brief   Return the position of XX(POS)=VAL
  !> @details Return the position of XX(POS)=VAL
  !>          Return 0 if VAL is not in XX
  !>
  !-----------------------------------------------------------------------

  pure function maths_findloc_0(xx,val) result(pos)

    integer(ip), pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: val
    integer(ip)                      :: kk,pos

    pos = 0
    if( associated(xx) ) then
       kk = 0
       do while( kk < size(xx) )
          kk = kk + 1
          if( xx(kk) == val ) then
             pos = kk
             return
          end if
       end do
    end if
    !
    ! findloc not supported by gcc
    !
    !if( associated(xx) ) then
    !   pos1 = findloc(xx,val)
    !   pos  = pos1(1)
    !else
    !   pos = 0
    !end if
    
  end function maths_findloc_0

  pure function maths_findloc_1(nn,xx,val) result(pos)

    integer(ip),          intent(in) :: nn
    integer(ip),          intent(in) :: xx(:)
    integer(ip),          intent(in) :: val
    integer(ip)                      :: kk,pos

    pos = 0
    do kk = 1,nn
       if( xx(kk) == val ) then
          pos = kk
          return
       end if
    end do
    
  end function maths_findloc_1

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2021-02-15
  !> @brief   Copy an array
  !> @details x1(p1(:)) = x2(p2(:))
  !>
  !-----------------------------------------------------------------------
  
  pure subroutine maths_copy_with_perm(nn,x1,x2,p1,p2)

    integer(ip),                    intent(in)    :: nn
    class(*),              pointer, intent(inout) :: x1(:)
    class(*),              pointer, intent(in)    :: x2(:)
    integer(ip), optional, pointer, intent(in)    :: p1(:)
    integer(ip), optional, pointer, intent(in)    :: p2(:)
    integer(ip)                                   :: ii

    select type ( x1 )
    type is ( real ( kind = rp ) )
       select type ( x2 )
       type is ( real ( kind = rp ) )       
          if( present(p1) ) then
             if( present(p2) ) then
                do ii = 1,nn
                   x1(p1(ii)) = x2(p2(ii))
                end do
             else
                do ii = 1,nn
                   x1(p1(ii)) = x2(ii)
                end do
             end if
          else
             if( present(p2) ) then
                do ii = 1,nn
                   x1(ii) = x2(p2(ii))
                end do
             else
                do ii = 1,nn
                   x1(ii) = x2(ii)
                end do
             end if
          end if
       end select
    end select

  end subroutine maths_copy_with_perm
  
end module mod_maths_arrays
!> @}
