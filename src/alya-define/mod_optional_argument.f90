!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Tools
!> @{
!> @file    mod_optional_argument.f90
!> @author  houzeaux
!> @date    2020-09-18
!> @brief   Optional arguments
!> @details Read some optional arguments or assign a default value
!>         if they are not present
!-----------------------------------------------------------------------

module mod_optional_argument

  use def_kintyp_basic, only : ip,rp,lg
  
  implicit none

  private
  
 interface optional_argument
     module procedure &
          &           optional_argument_ch,    &
          &           optional_argument_lg,    &
          &           optional_argument_lg_1,  &
          &           optional_argument_i4,    &
          &           optional_argument_i4_2,  &
          &           optional_argument_i8,    &
          &           optional_argument_i8_1,  &
          &           optional_argument_i8_2,  &
          &           optional_argument_rp
  end interface optional_argument
  
  public :: optional_argument

contains
  
  pure function optional_argument_ch(DEFAULT,CONDITION) result(my_char)

    character(len=*),            intent(in) :: DEFAULT
    character(len=*), optional,  intent(in) :: CONDITION
    character(len=:), allocatable           :: my_char
    
    if( present(CONDITION) ) then
       my_char = CONDITION
    else       
       my_char = DEFAULT
    end if
    
  end function optional_argument_ch
  
  logical(lg) pure function optional_argument_lg(DEFAULT,CONDITION)

    logical(lg),           intent(in) :: DEFAULT
    logical(lg), optional, intent(in) :: CONDITION

    if( present(CONDITION) ) then
       optional_argument_lg = CONDITION
    else       
       optional_argument_lg = DEFAULT
    end if
    
  end function optional_argument_lg
  
  logical(lg) pure function optional_argument_lg_1(DEFAULT,CONDITION,POSITION)

    logical(lg),                    intent(in) :: DEFAULT
    logical(lg), optional, pointer, intent(in) :: CONDITION(:)
    integer(ip),                    intent(in) :: POSITION

    if( present(CONDITION) ) then
       optional_argument_lg_1 = CONDITION(POSITION)
    else       
       optional_argument_lg_1 = DEFAULT
    end if
    
  end function optional_argument_lg_1
  
  integer(4) pure function optional_argument_i4(DEFAULT,CONDITION)

    integer(4),           intent(in) :: DEFAULT
    integer(4), optional, intent(in) :: CONDITION

    if( present(CONDITION) ) then
       optional_argument_i4 = CONDITION
    else       
       optional_argument_i4 = DEFAULT
    end if
    
  end function optional_argument_i4
  
  integer(4) pure function optional_argument_i4_2(DEFAULT,CONDITION,POSITION)

    integer(4),           intent(in) :: DEFAULT
    integer(4), optional, intent(in) :: CONDITION(:)
    integer(4),           intent(in) :: POSITION

    if( present(CONDITION) ) then
       if( POSITION <= size(CONDITION) ) then
          optional_argument_i4_2 = CONDITION(POSITION)
       else
          optional_argument_i4_2 = DEFAULT          
       end if
    else       
       optional_argument_i4_2 = DEFAULT
    end if
    
  end function optional_argument_i4_2
  
  integer(8) pure function optional_argument_i8_2(DEFAULT,CONDITION,POSITION)

    integer(8),           intent(in) :: DEFAULT
    integer(8), optional, intent(in) :: CONDITION(:)
    integer(8),           intent(in) :: POSITION

    if( present(CONDITION) ) then
       if( POSITION <= int(size(CONDITION),8) ) then
          optional_argument_i8_2 = CONDITION(POSITION)
       else
          optional_argument_i8_2 = DEFAULT          
       end if
    else       
       optional_argument_i8_2 = DEFAULT
    end if
    
  end function optional_argument_i8_2
  
  integer(8) pure function optional_argument_i8(DEFAULT,CONDITION)

    integer(8),           intent(in) :: DEFAULT
    integer(8), optional, intent(in) :: CONDITION

    if( present(CONDITION) ) then
       optional_argument_i8 = CONDITION
    else       
       optional_argument_i8 = DEFAULT
    end if
    
  end function optional_argument_i8
  
  pure function optional_argument_i8_1(DEFAULT,CONDITION) result(out)

    integer(8),           intent(in) :: DEFAULT(:)
    integer(8), optional, intent(in) :: CONDITION(*)
    integer(8)                       :: out(size(DEFAULT,DIM=1))
    integer(ip)                      :: nn
    
    nn = size(DEFAULT)
    
    if( present(CONDITION) ) then
       out(1:nn) = CONDITION(1:nn)
    else       
       out(1:nn) = DEFAULT(1:nn)
    end if
    
  end function optional_argument_i8_1
  
  real(rp) pure function optional_argument_rp(DEFAULT,CONDITION)

    real(rp),           intent(in) :: DEFAULT
    real(rp), optional, intent(in) :: CONDITION

    if( present(CONDITION) ) then
       optional_argument_rp = CONDITION
    else       
       optional_argument_rp = DEFAULT
    end if
    
  end function optional_argument_rp

end module mod_optional_argument
!> @}
