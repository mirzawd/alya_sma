!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_edict.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   dictionary with edge keys
!> @details dictionary with edge keys
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

module mod_edict
  
  use def_kintyp_basic, only: ip,rp,lg
  
  implicit none

  private
  public    :: edict_type
  !public    :: edge_type

  integer, parameter  :: iphash = ip!16

  type edge_type
    integer(ip) :: n1
    integer(ip) :: n2
  contains
    procedure :: set        => setEdge
    procedure :: isEqualTo  => areEqualEdges
  end type edge_type
  
  type pair_type
     type(edge_type) :: key
     logical(lg) :: value
  end type pair_type

  type chunk_type
     type(pair_type), allocatable :: pairs(:)
     integer(ip) :: size = 0_ip
     integer(ip) :: count = 0_ip
   contains
     procedure :: find
  end type chunk_type

  type edict_type
     type(chunk_type), allocatable :: chunks(:)
     integer(ip) :: size = 0_ip
   contains
     procedure :: init
     procedure :: deallo
     procedure :: isAlloca
     procedure :: setByNodes
     procedure :: setByEdge
     procedure :: haskey
     procedure :: isEdge
     procedure :: areEdge => isEdgeByNodes
     procedure :: get
     procedure :: hash  => hash_djb2!hash_djb_bis!hash_djb2
     procedure :: printDataStructure
     procedure :: printKeys
  end type edict_type
  
  integer(ip), parameter :: chunk_initial_size = 2_ip
  ! with one of the following would be enough
  integer(ip), parameter :: tag_isEmpty    = -1_ip
  integer(ip), parameter :: tag_isNotFound = -2_ip

contains

  
  subroutine setEdge(edge,node1,node2) 
    class(edge_type), intent(out) :: edge
    integer(ip), intent(in) :: node1,node2
    if(node1<node2) then
      edge%n1=node1
      edge%n2=node2
    elseif(node1>node2) then
      edge%n1=node2
      edge%n2=node1
    else
      call runend("edge with two equal nodes..")
    end if
  end subroutine setEdge
  
  function areEqualEdges(edge,edge2) result(isEqual)
    class(edge_type), intent(in) :: edge
    type(edge_type),  intent(in) :: edge2
    logical(lg) :: isEqual
    isEqual = (edge%n1==edge2%n1).and.(edge%n2==edge2%n2)
  end function areEqualEdges

  subroutine init(dict, size)
    class(edict_type), intent(out) :: dict
    integer(ip),            intent(in)  :: size
    !
    dict%size = size
    allocate(dict%chunks(size))
    !
  end subroutine init

  subroutine setByNodes(dict, n1,n2,value)
    class(edict_type), intent(inout) :: dict
    integer(ip),       intent(in)    :: n1,n2
    logical(lg),       intent(in)    :: value
    !
    type(edge_type) :: edge
    !
    call edge%set(n1,n2)
    call dict%setByEdge(edge,value)
    !
  end subroutine setByNodes
  
  subroutine setByEdge(dict, key, value)
    class(edict_type), intent(inout) :: dict
    type(edge_type),   intent(in)    :: key
    logical(lg),       intent(in)    :: value

    type(chunk_type) :: chunk_aux

    integer(ip) :: i, count
    integer(ip) :: h

    h = dict%hash(key)
    count = dict%chunks(h)%find(key)

    if (count == tag_isEmpty) then
       allocate(dict%chunks(h)%pairs(chunk_initial_size))
       dict%chunks(h)%size  = chunk_initial_size
       dict%chunks(h)%count = 1_ip
       
       dict%chunks(h)%pairs(1)%key   = key
       dict%chunks(h)%pairs(1)%value = value
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
       dict%chunks(h)%pairs(count)%key   = key
       dict%chunks(h)%pairs(count)%value = value
    end if

  end subroutine setByEdge

  function haskey(dict, k) result(is_key_in_dict)
    class(edict_type), intent(in) :: dict
    type(edge_type), intent(in) :: k
    
    logical(lg) :: is_key_in_dict

    integer(ip) :: h, count

    h = dict%hash(k)
    count = dict%chunks(h)%find(k)

    is_key_in_dict = (count == tag_isEmpty).or.(count == tag_isNotFound)

  end function haskey

  function isAlloca(dict) result(isAllo)
    class(edict_type), intent(in) :: dict
    logical(lg) :: isAllo
    
    isAllo = dict%size>0_ip
    
  end function isAlloca
  
  function isEdgeByNodes(dict, n1,n2) result(isTrueEdge)
    class(edict_type), intent(in) :: dict
    integer(ip),       intent(in)    :: n1,n2
    logical(lg) :: isTrueEdge
    !
    type(edge_type) :: edge
    !
    call edge%set(n1,n2)
    
    isTrueEdge =  dict%isEdge(edge)
    !
  end function isEdgeByNodes
  
  function isEdge(dict, edge) result(isTrueEdge)
    class(edict_type), intent(in) :: dict
    type(edge_type), intent(in) :: edge
    !
    logical(lg) :: isTrueEdge
    !
    isTrueEdge =  dict%get(edge)
    !
  end function isEdge
  
  function get(dict, k) result(value)
    class(edict_type), intent(in) :: dict
    type(edge_type), intent(in) :: k

    logical(lg) :: value

    integer(ip) :: h, count

    h = dict%hash(k)
    count = dict%chunks(h)%find(k)

    if ( (count == tag_isEmpty) .or. &
         (count == tag_isNotFound) ) then
       value = .false.
       return
    end if

    if (count>0_ip) then
       value = dict%chunks(h)%pairs(count)%value
    end if
  end function get

  function hash_djb2(dict, key) result(h_out)  !hash function djb2 (could use sdbm, or lose/lose)
    class(edict_type),  intent(in) :: dict
    type(edge_type),         intent(in) :: key
    integer(ip) :: h_out
    
    integer(iphash)  :: h
    
    integer(iphash), parameter :: h_ini    = 5381_iphash
    integer(iphash), parameter :: h_factor =   33_iphash
    
    h = h_ini*h_factor + key%n1
    h = h    *h_factor + key%n2
    
    h_out = int( modulo(h, int(dict%size,iphash) ) , ip )
    h_out = h_out + 1_ip
    
  end function hash_djb2
  
  subroutine printDataStructure(dict)
    class(edict_type), intent(in) :: dict
    integer(ip) :: i, j, count_chunk
    do i = 1, dict%size
       count_chunk = dict%chunks(i)%count
       print*, ' '
       print*, 'i: ',i, " from tot: ",dict%size
       print*, 'chunks(i)%size    : ', dict%chunks(i)%size
       print*, 'chunks(i)%count   : ', dict%chunks(i)%count
       if (count_chunk > 0) then
          print*, '  chunk   : ', i, ' size ', count_chunk
          do j = 1, count_chunk
             print*, '  -> edge      : ', dict%chunks(i)%pairs(j)%key%n1, "  ", dict%chunks(i)%pairs(j)%key%n2
             print*, '      -> value    : ',  dict%chunks(i)%pairs(j)%value
          end do
       end if
    end do
  end subroutine printDataStructure
  
  subroutine printKeys(dict)
    class(edict_type), intent(in) :: dict
    integer(ip) :: i, j, count_chunk
    do i = 1, dict%size
       count_chunk = dict%chunks(i)%count
       if (count_chunk > 0) then
          do j = 1, count_chunk
             print*, '  edge      : ', dict%chunks(i)%pairs(j)%key%n1, "  ", dict%chunks(i)%pairs(j)%key%n2
          end do
       end if
    end do
  end subroutine printKeys

  function find(chunk, k) result(r)
    class(chunk_type), intent(in) :: chunk
    type(edge_type), intent(in) :: k
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
      if( k%isEqualTo(chunk%pairs(i)%key) ) then
        r = i
        exit
      end if
    end do
    !
  end function find
  
  subroutine deallo(dict)
    class(edict_type), intent(inout) :: dict
    !
    dict%size = 0
    deallocate(dict%chunks)  ! allocatable things inside should be deallocated automatically
    !
  end subroutine deallo
  
end module mod_edict
!> @}


  
!   function hash_djb_bis(dict, key) result(h_out)
!     class(edict_type),  intent(in) :: dict
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





