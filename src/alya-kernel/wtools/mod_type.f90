!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_type.f90
!> @author  houzeaux
!> @date    2020-06-17
!> @brief   Return the type of a variable
!> @details Return the type of a variable
!-----------------------------------------------------------------------

module mod_type

  implicit none
  private
  
  interface my_type
     module procedure my_type0,my_type1,my_type2,my_type3,my_type4
  end interface my_type

  public :: my_type
  
contains

  function my_type0(object) result(current_type)

    class(*),         intent(in)  :: object
    character(len=:), allocatable :: current_type

    select type(object)
       
    type is (real(4))
       
       current_type = 'real'
       
    type is (real(8))
       
       current_type = 'real'
       
    type is (integer(4))
       
       current_type = 'integer'
       
    type is (integer(8))
       
       current_type = 'integer'
       
    class default

       current_type = 'unknown'       
       
    end select
    
  end function my_type0

  function my_type1(object) result(current_type)

    class(*),         dimension(:), intent(in) :: object
    character(len=:), allocatable              :: current_type

    select type(object)
       
    type is (real(4))
       
       current_type = 'real'
       
    type is (real(8))
       
       current_type = 'real'
       
    type is (integer(4))
       
       current_type = 'integer'
       
    type is (integer(8))
       
       current_type = 'integer'
       
    end select

  end function my_type1

  function my_type2(object) result(current_type)

    class(*),         dimension(:,:), intent(in) :: object
    character(len=:), allocatable                :: current_type

    select type(object)
       
    type is (real(4))
       
       current_type = 'real'
       
    type is (real(8))
       
       current_type = 'real'
       
    type is (integer(4))
       
       current_type = 'integer'
       
    type is (integer(8))
       
       current_type = 'integer'
       
    class default

       current_type = 'unknown'
       
    end select

  end function my_type2
  
  function my_type3(object) result(current_type)

    class(*),         dimension(:,:,:), intent(in) :: object
    character(len=:), allocatable                  :: current_type

    select type(object)
       
    type is (real(4))
       
       current_type = 'real'
       
    type is (real(8))
       
       current_type = 'real'
       
    type is (integer(4))
       
       current_type = 'integer'
       
    type is (integer(8))
       
       current_type = 'integer'
       
    class default

       current_type = 'unknown'
       
    end select

  end function my_type3

  function my_type4(object) result(current_type)

    class(*),         dimension(:,:,:,:), intent(in) :: object
    character(len=:), allocatable                    :: current_type

    select type(object)
       
    type is (real(4))
       
       current_type = 'real'
       
    type is (real(8))
       
       current_type = 'real'
       
    type is (integer(4))
       
       current_type = 'integer'
       
    type is (integer(8))
       
       current_type = 'integer'
       
    class default

       current_type = 'unknown'
       
    end select

  end function my_type4

end module mod_type
!> @}

