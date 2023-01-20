!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_dictionary.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   dictionary
!> @details dictionary
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

module mod_dictionary
  
  use def_kintyp_basic, only: ip,rp,lg
  
  implicit none

  private
  public    :: dictionary_type

  integer, parameter  :: iphash = ip!16

  type pair_type
     character(len=:), allocatable :: key
     character(len=:), allocatable :: value
  end type pair_type

  type chunk_type
     type(pair_type), allocatable :: pairs(:)
     integer(ip) :: size = 0_ip
     integer(ip) :: count = 0_ip
   contains
     procedure :: find
  end type chunk_type

  type dictionary_type
     type(chunk_type), allocatable :: chunks(:)
     integer(ip) :: size = 0_ip
   contains
     procedure :: init
     procedure :: set
     procedure :: haskey
     procedure :: get
     procedure :: hash => hash_djb2!hash_djb_bis!hash_djb2
     procedure :: deallo
     procedure :: print
  end type dictionary_type
  
  integer(ip), parameter :: dict_size_default  = 1024_ip
  integer(ip), parameter :: chunk_initial_size = 2_ip
  ! with one of the following would be enough
  integer(ip), parameter :: tag_isEmpty    = -2_ip
  integer(ip), parameter :: tag_isNotFound = -3_ip!4_ip

contains

  subroutine init(dict, size)
    class(dictionary_type), intent(out) :: dict
    integer(ip),            intent(in) , optional :: size
    !
    if(present(size)) then
      dict%size = size
      allocate(dict%chunks(size))
    else
      dict%size = dict_size_default
      allocate(dict%chunks(size))
    end if
    !
  end subroutine init
  
  subroutine deallo(dict)
    class(dictionary_type), intent(inout) :: dict
    !
    dict%size = 0
    deallocate(dict%chunks)  ! allocatable things inside should be deallocated automatically
    !
  end subroutine deallo

  subroutine set(dict, k, v)
    class(dictionary_type), intent(inout) :: dict
    character(len=*),       intent(in)    :: k
    character(len=*),       intent(in)    :: v

    type(chunk_type) :: chunk_aux

    integer(ip) :: i, count
    integer(ip) :: h

    h = dict%hash(k)

    count = dict%chunks(h)%find(k)

    if (count == tag_isEmpty) then
       allocate(dict%chunks(h)%pairs(chunk_initial_size))
       dict%chunks(h)%size  = chunk_initial_size
       dict%chunks(h)%count = 1_ip
       
       dict%chunks(h)%pairs(1)%key   = trim(k)
       dict%chunks(h)%pairs(1)%value = trim(v)
       return
    end if

    if (count == tag_isNotFound) then
      if(dict%chunks(h)%count==dict%chunks(h)%size) then !increase size 
        
        chunk_aux%size  = dict%chunks(h)%size*2_ip        
        allocate(chunk_aux%pairs( chunk_aux%size  ))
        do i = 1, dict%chunks(h)%size
           chunk_aux%pairs(i)%key   = dict%chunks(h)%pairs(i)%key
           chunk_aux%pairs(i)%value = dict%chunks(h)%pairs(i)%value
        end do
        
        deallocate(dict%chunks(h)%pairs)
        allocate(dict%chunks(h)%pairs, source=chunk_aux%pairs)
        dict%chunks(h)%size  = chunk_aux%size
        
        deallocate(chunk_aux%pairs)
      end if
      
      dict%chunks(h)%count = dict%chunks(h)%count + 1_ip
      count = dict%chunks(h)%count
    end if

    if (count > 0_ip) then
       dict%chunks(h)%pairs(count)%key   = k
       dict%chunks(h)%pairs(count)%value = v
    end if

  end subroutine set

  function haskey(dict, k) result(is_key_in_dict)
    class(dictionary_type), intent(in) :: dict
    character(len=*), intent(in) :: k
    
    logical(lg) :: is_key_in_dict

    integer(ip) :: h, count

    h = dict%hash(k)
    count = dict%chunks(h)%find(k)

    is_key_in_dict = (count == tag_isEmpty).or.(count == tag_isNotFound)

  end function haskey

  function get(dict, k) result(r)
    class(dictionary_type), intent(in) :: dict
    character(len=*), intent(in) :: k

    character(len=:), allocatable :: r

    integer(ip) :: h, count

    h = dict%hash(k)

    count = dict%chunks(h)%find(k)

    if ( (count == tag_isEmpty) .or. &
         (count == tag_isNotFound) ) then
       r = ''
       return
    end if

    if (count>0_ip) then
       r = dict%chunks(h)%pairs(count)%value
    end if
  end function get

  function hash_djb2(dict, key) result(h_out)  !hash function djb2 (could use sdbm, or lose/lose)
    class(dictionary_type),  intent(in) :: dict
    character(len=*),        intent(in) :: key
    integer(ip) :: h_out
    
    integer(iphash)  :: h
    integer(ip) :: i
    character   :: c
    
    integer(iphash) :: c_i
    integer(iphash), parameter :: h_ini    = 5381_iphash
    integer(iphash), parameter :: h_factor =   33_iphash
    
    h = h_ini
    do i = 1_ip, len(key)
      c = key(i:i)
      c_i = ichar(c)
      
      !h = h*h_factor + c_i
      h = int( modulo(h*h_factor, int(dict%size,iphash) ), ip)
      h = int( modulo( h + c_i  , int(dict%size,iphash) ), ip)
    end do
    
    !h_out = int( modulo(h, int(dict%size,iphash) ) , ip )
    h_out = h
    
    h_out = h_out + 1_ip
    
  end function hash_djb2
  
!   function hash_djb_bis(dict, key) result(h_out)
!     class(dictionary_type),  intent(in) :: dict
!     character(len=*),intent(in) :: key
!     integer(ip) :: h_out
!
!     integer(iphash)  :: h
!
!     integer(ip) :: i
!     integer(iphash) :: c_i
!
!     integer(iphash), parameter :: h_ini    = 5381_iphash
!     integer(iphash), parameter :: h_factor =   33_iphash
!     integer(iphash), parameter :: h_ishft  =    5_iphash
!
!     h = h_ini
!     do i = 1_ip, len(key)
!       c_i = ichar(key(i:i))
!       h = (ishft(h,h_ishft) + h) + c_i
!     end do
!
!     h_out = int( modulo(h, int(dict%size,iphash) ) , ip )
!     h_out = h_out + 1_ip
!
!   end function hash_djb_bis

  subroutine print(dict)
    class(dictionary_type), intent(in) :: dict
    integer(ip) :: i, j, s
    do i = 1, dict%size
       s = dict%chunks(i)%count
       print*, ' '
       print*, 'i: ',i, " from tot: ",dict%size
       print*, 'chunks(i)%size    : ', dict%chunks(i)%size
       print*, 'chunks(i)%count   : ', dict%chunks(i)%count
       if (s > 0) then
          write(*,*) 'chunk   : ', i, ' size ', s
          do j = 1, s
             print*, '-> key      : ', dict%chunks(i)%pairs(j)%key
             print*, '-> value    : ', dict%chunks(i)%pairs(j)%value
          end do
       end if
    end do
  end subroutine print

  function find(chunk, k) result(r)
    class(chunk_type), intent(in) :: chunk
    character(len=*), intent(in) :: k
    integer(ip) :: r
    !
    integer(ip) :: i
    !
    if (chunk%size == 0_ip) then
       r = tag_isEmpty
       return
    end if
    !
    r = tag_isNotFound
    do i = 1_ip, chunk%count
       if (chunk%pairs(i)%key == trim(k)) then
          r = i
          exit
       end if
    end do
    !
  end function find
  
end module mod_dictionary
!> @}








