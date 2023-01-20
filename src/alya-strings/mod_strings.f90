!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Tools
!> @{
!> @file    mod_strings.f90
!> @author  houzeaux
!> @date    2020-02-27
!> @brief   Strings
!> @details Operations with strings
!-----------------------------------------------------------------------

module mod_strings

  use def_kintyp_basic, only : ip,rp,lg
  implicit none
  private
  
  character(26),                parameter   :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26),                parameter   :: low = 'abcdefghijklmnopqrstuvwxyz'
  integer(ip),   dimension (:), allocatable :: indexarray
  logical(lg)                               :: CaseSensitive

  interface integer_to_string
     module procedure integer_to_string_4,integer_to_string_8
  end interface integer_to_string

  public :: upper_case
  public :: lower_case
  public :: integer_to_string
  public :: real_to_string
  public :: string_to_integer
  public :: string_sort
  public :: add_extension
  public :: find_substr
  public :: string_continue

contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Sort arrays
  !> @details Subroutine sort uses the quicksort algorithm.
  !>          On input, StringArray is a one-dimensional array of character strings
  !>          to be sorted in ascending lexical order.
  !>          On output, StringArray is the sorted array.
  !>          If the optional argument CaseInsensitive is present and .true.,
  !>          the sort is case-insensitive. If CaseInsensitive is absent or
  !>          if it is .false., the sort is case-sensitive.
  !>          The characters of the elements of the string array are not modified,
  !>          so that if blanks or punctuation characters are to be ignored,
  !>          for instance, this needs to be done before calling sort.
  !> 
  !-----------------------------------------------------------------------

  subroutine string_sort(StringArray, nn, CaseInsensitive_opt)

    character (len = *), dimension (:), intent(inout) :: StringArray
    integer(ip),                        intent(out)   :: nn
    logical(lg),         optional,      intent(in)    :: CaseInsensitive_opt
    integer(ip)                                       :: low, high, k, ii
    ! 
    ! Sort
    !
    if (present(CaseInsensitive_opt)) then
       CaseSensitive = .not. CaseInsensitive_opt
    else
       CaseSensitive = .true.
    end if
    low = 1
    high = int(size(StringArray), ip)
    allocate(indexarray(high))
    indexarray = (/ (k, k = low, high) /)
    call string_quicksort(StringArray, low, high)
    StringArray = StringArray(indexarray)
    deallocate(indexarray)
    !
    ! Merge list
    !
    ii = 0
    nn = 1
    do while( ii < high )
       ii = ii + 1
       if( StringArray(ii) /= StringArray(nn) ) then
          nn = nn + 1
          StringArray(nn) = StringArray(ii)
       end if
    end do

  end subroutine string_sort
  
  recursive subroutine string_quicksort(StringArray, low, high)  
    character (len = *), dimension (:), intent (inout) :: StringArray
    integer(ip),                        intent (in)    :: low, high
    integer(ip)                                        :: pivotlocation
    
    if (low < high) then
      call string_partition(StringArray, low, high, pivotlocation)
      call string_quicksort(StringArray, low, pivotlocation - 1)
      call string_quicksort(StringArray, pivotlocation + 1, high)
   end if
   
  end subroutine string_quicksort
  
  subroutine string_partition(StringArray, low, high, pivotlocation)
    character (len = *), dimension (:), intent (inout) :: StringArray
    integer(ip),                        intent (in)    :: low, high
    integer(ip),                        intent (out)   :: pivotlocation
    integer(ip)                                        :: k, lastsmall
    
    call string_swap(indexarray(low), indexarray((low + high)/2))
    lastsmall = low
    do k = low + 1, high
      if (string_stringComp(StringArray(indexarray(k)), StringArray(indexarray(low)))) then
        lastsmall = lastsmall + 1
        call string_swap(indexarray(lastsmall), indexarray(k))
      end if
    end do
    call string_swap(indexarray(low), indexarray(lastsmall))
    pivotlocation = lastsmall
  end subroutine string_partition
  
  subroutine string_swap(m, n)
    integer(ip), intent (inout) :: m, n
    integer(ip)                 :: temp
    temp = m
    m = n
    n = temp
  end subroutine string_swap
  
  function string_stringComp(p, q) result(lexicalLess)
    character (len = *), intent (in) :: p, q
    logical(lg)                      :: lexicalLess
    integer(ip)                      :: kq, k
    
    if (CaseSensitive) then
      lexicalLess = p < q
    else
      kq = 1
      do k = 1, int(max(len_trim(p), len_trim(q)), ip)
        if (string_UpperCase(p(k:k)) == string_UpperCase(q(k:k)) ) then
          cycle
        else
          kq = k
          exit
        end if
      end do
      lexicalLess = string_UpperCase(p(kq:kq)) < string_UpperCase(q(kq:kq))
    end if
  end function string_stringComp
  
  function string_UpperCase(letter) result(L)
    character (len = *), intent (in) :: letter
    character (len = 1)              :: L
    integer(ip)                      :: k
    
    k = int(index(low, letter), ip)
    if (k > 0) then
      L = cap(k:k)
    else
      L = letter
   end if
   
  end function string_UpperCase
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Integer to string
  !> @details Convert an integer(4) and integer(8) to a string
  !> 
  !-----------------------------------------------------------------------

  function integer_to_string_4(integ) result(intost)

    integer(4),       intent(in)   :: integ
    character(len=:), allocatable  :: intost
    character(20)                  :: intaux

    write(intaux,*) integ  
    intost = trim(adjustl(intaux))

  end function integer_to_string_4

  function integer_to_string_8(integ) result(intost)

    integer(8),       intent(in)   :: integ
    integer(4)                     :: integ4
    character(len=:), allocatable  :: intost
    character(20)                  :: intaux
    
    integ4 = int(integ,4)
    write(intaux,*) integ4    
    intost = trim(adjustl(intaux))

  end function integer_to_string_8

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   String to integer
  !> @details Convert a string to an integer
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function string_to_integer(str,stat) result(integ)

    character(len=*),intent(in)            :: str
    integer(ip),     intent(out), optional :: stat
    integer(4)                             :: stat_loc
    
    read(str,*,iostat=stat_loc) integ
    if( present(stat) ) stat = int(stat_loc,ip)
    
  end function string_to_integer

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Real to string
  !> @details Convert real(rp) to a string
  !> 
  !-----------------------------------------------------------------------

  function real_to_string(realn,REAL_FORMAT) result(retost)

    real(rp)                               :: realn
    character(len=*), intent(in), optional :: REAL_FORMAT
    character(len=:), allocatable          :: retost
    character(20)                          :: reaux
    integer(ip)                            :: ierr

    if( present(REAL_FORMAT) ) then
       write(reaux,REAL_FORMAT,IOSTAT=ierr) realn
       if( ierr /= 0 ) reaux = '0.0'
    else
       write(reaux,'(e19.12)') realn
    end if
    retost = trim(adjustl(reaux))

  end function real_to_string

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to upper case
  !> @details Convert to upper case
  !> 
  !-----------------------------------------------------------------------

  function upper_case (str) result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, int(len_trim(str), ip)
        ic = int(index(low, str(i:i)), ip)
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

  end function upper_case

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to lower case
  !> @details Convert to lower case
  !> 
  !-----------------------------------------------------------------------

  function lower_case (str) result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, int(len_trim(str), ip)
        ic = int(index(cap, str(i:i)), ip)
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

  end function lower_case

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Add extension
  !> @details Add extension
  !> 
  !-----------------------------------------------------------------------

  subroutine add_extension(str,ext)
    
    character(*), intent(inout) :: str
    character(*), intent(in)    :: ext
    
    if( index(str,'.',KIND=ip) == 0 ) str = trim(str)//'.'//trim(ext)

  end subroutine add_extension
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Continue a sentence
  !> @details Continue a sentence
  !> 
  !-----------------------------------------------------------------------
  
  subroutine string_continue(message,link,message_in)

    character(len=:), allocatable, intent(inout) :: message
    character(len=*),              intent(in)    :: link
    character(len=*),              intent(in)    :: message_in

    if( allocated(message) ) then
       if( len_trim(message) == 0 ) then
          message = link // message_in
       else
          message = message // link // message_in
       end if
    else
       message = message_in
    end if

  end subroutine string_continue
  


  integer(ip) function find_substr(str, substr)
    implicit none
    character(len=*), intent(in) :: str, substr
    integer(ip) :: str_pos, substr_pos
    
    find_substr = -1_ip
    substr_pos = 1_ip
    do str_pos = 1_ip, len(str, kind=ip)
      if ( substr(substr_pos:substr_pos)==str(str_pos:str_pos) ) then
          if( find_substr<0_ip ) find_substr = str_pos
          substr_pos = substr_pos + 1_ip
      else
          find_substr = -1_ip
          substr_pos = 1_ip
      end if

      if ( substr_pos>=len(substr, kind=ip) ) exit
    end do

  end function 

  
end module mod_strings
!> @}
