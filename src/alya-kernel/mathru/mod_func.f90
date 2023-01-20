!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_func.f90
!> @author  houzeaux and eduardo perez
!> @date    2020-02-03
!> @brief   Functions
!> @details Generic functions
!-----------------------------------------------------------------------

module mod_func

  use def_kintyp_basic, only : ip,rp
  implicit none

  abstract interface
     real(rp) function func(x,gradx)
       import             :: rp
       real(rp)           :: x(:)
       real(rp), optional :: gradx(:,:)
     end function func
  end interface

  type func_ptr
     procedure(func), pointer, nopass :: f
  end type func_ptr

contains

  real(rp) function func_identity(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_identity = val(1)
  end function func_identity

  real(rp) function func_square(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_square = dot_product(val,val)
  end function func_square

  real(rp) function func_unity(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_unity = 1.0_rp
  end function func_unity

  real(rp) function func_sine(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_sine = sin(val(1))
  end function func_sine

  real(rp) function func_cosine(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_cosine = cos(val(1))
  end function func_cosine

  real(rp) function func_norm(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    func_norm = sqrt(dot_product(val,val))
  end function func_norm

  real(rp) function func_div(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    integer(ip)        :: ii    
    func_div = 0.0_rp
    if( present(grad) ) then
       do ii = 1,min(size(grad,1),size(grad,2))
          func_div = func_div + grad(ii,ii)
       end do
    end if
  end function func_div

  real(rp) function func_dotn(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    integer(ip)        :: ii,nn
    

    func_dotn = 0.0_rp
    nn        = size(grad,2)
    
    if( present(grad) ) then
       do ii = 1,min(size(val,1),size(grad,1))
          func_dotn = func_dotn + val(ii)*grad(ii,nn)
       end do
    end if
  end function func_dotn

  real(rp) function func_rotational_norm(val,grad)
    real(rp)           :: val(:)
    real(rp), optional :: grad(:,:)
    real(rp)           :: rot(3)
    func_rotational_norm = 0.0_rp
    if( present(grad) ) then
       if( size(grad,1) == 2 ) then
          rot(1) = grad(1,2)-grad(2,1)
          func_rotational_norm = abs(rot(1))
       else
          rot(1) = grad(2,3)-grad(3,2)
          rot(2) = grad(3,1)-grad(1,3)
          rot(3) = grad(1,2)-grad(2,1)
          func_rotational_norm = sqrt(dot_product(rot,rot))
       end if
    end if
  end function func_rotational_norm

  subroutine func_initialization(my_func)
    type(func_ptr), intent(inout) :: my_func(:,:)
    integer(ip)                   :: ii,jj
    do jj = 1,size(my_func,2)
       do ii = 1,size(my_func,1)
          my_func(ii,jj) % f => func_unity
       end do
    end do    
  end subroutine func_initialization

end module mod_func
!> @}
