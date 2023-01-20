!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_linked_list.f90
!> @author  houzeaux
!> @date    2020-12-06
!> @brief   Linked list
!> @details Linked list
!-----------------------------------------------------------------------

module def_linked_list

  use def_kintyp_basic, only : ip,rp,lg
  implicit none

  private
  integer(ip),              parameter   :: TYPE_INTEGER   = 1
  integer(ip),              parameter   :: TYPE_REAL      = 2
  integer(ip),              parameter   :: TYPE_LOGICAL   = 3
  integer(ip),              parameter   :: TYPE_CHARACTER = 4

  type node_typ
     class(*),              pointer     :: vals 
     class(*),              pointer     :: val1(:)
     class(*),              pointer     :: val2(:,:)
     type(node_typ),        pointer     :: next
     type(node_typ),        pointer     :: prev
     integer(ip)                        :: rank
     integer(ip)                        :: type
     integer(ip)                        :: id
     integer(ip)                        :: assigned
     character(len=:),      allocatable :: name
   contains
     procedure,             pass        :: init => init_node
     procedure,             pass        :: getdim
     procedure,             pass        :: exists
     procedure,             pass        :: sets_ip
     procedure,             pass        :: set1_ip
     procedure,             pass        :: set2_ip
     procedure,             pass        :: sets_rp
     procedure,             pass        :: set1_rp
     procedure,             pass        :: set2_rp
     procedure,             pass        :: sets_lg
     procedure,             pass        :: set1_lg
     procedure,             pass        :: set2_lg
     procedure,             pass        :: sets_ch
     generic                            :: set   => sets_ip,set1_ip,set2_ip,&
          &                                         sets_rp,set1_rp,set2_rp,&
          &                                         sets_lg,set1_lg,set2_lg,&
          &                                         sets_ch
  end type node_typ
  
  type llist_typ
     type(node_typ),        pointer     :: head
     type(node_typ),        pointer     :: tail
     integer(ip)                        :: size
   contains
     procedure,             pass        :: init => init_llist
     procedure,             pass        :: iterate
     procedure,             pass        :: add  
     procedure,             pass        :: deallo  
     procedure,             pass        :: untarget  
     procedure,             pass        :: first 
     procedure,             pass        :: last  
     procedure,             pass        :: nameis
  end type llist_typ

  abstract interface
     subroutine what_to_do(node,ivari)
       import :: ip
       import :: node_typ
       class(node_typ),           intent(inout) :: node
       integer(ip),     optional, intent(in)    :: ivari
     end subroutine what_to_do
  end interface

  public :: node_typ
  public :: llist_typ
  public :: TYPE_INTEGER   
  public :: TYPE_REAL      
  public :: TYPE_LOGICAL   
  public :: TYPE_CHARACTER 
  public :: what_to_do
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Linked list initialization
  !> @details Linked list initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init_llist(self)
    class(llist_typ), intent(inout) :: self
    
    nullify(self % head)
    nullify(self % tail)
    self % size     =  0
    
  end subroutine init_llist
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Node initialization
  !> @details Node initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init_node(self)
    class(node_typ), intent(inout) :: self
    
    nullify(self % next)
    nullify(self % prev)
    nullify(self % vals)
    nullify(self % val1)
    nullify(self % val2)
    self % rank     = -1
    self % type     = -1
    self % assigned =  0

  end subroutine init_node
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Tail and head
  !> @details Point to tail and head
  !> 
  !-----------------------------------------------------------------------

  function first(self)  result(firstnode)
    class(llist_typ), intent(in)         :: self
    type(node_typ),              pointer :: firstnode
    firstnode => self % head
  end function first
  
  function last(self)  result(lastnode)
    class(llist_typ),       intent(in) :: self
    type(node_typ),         pointer    :: lastnode
    lastnode => self % tail
  end function last

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Set value
  !> @details Set value
  !> 
  !-----------------------------------------------------------------------

  subroutine sets_ip(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    integer(ip),     pointer,  intent(in)    :: val
    integer(ip),     optional, intent(out)   :: ierr

    if( associated(val) ) then
       if( present(ierr) ) ierr = 0       
       select type ( v => self % vals )
       type is ( integer(kind=ip) ) ; v = val
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
    
  end subroutine sets_ip
  
  subroutine set1_ip(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    integer(ip),     pointer,  intent(in)    :: val(:)
    integer(ip),     optional, intent(out)   :: ierr
    integer(ip)                              :: n1
    
    if( associated(val) ) then
       if( present(ierr) ) ierr = 0       
       select type ( v => self % val1 )
       type is ( integer(kind=ip) )
          if( size(v) /= size(val) ) then
             n1 = min(size(v),size(val))
             v(1:n1) = val(1:n1)
             if( present(ierr) ) ierr = 1
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
    
  end subroutine set1_ip
  
  subroutine set2_ip(self,val,ierr)
    class(node_typ),            intent(inout) :: self
    integer(ip),      pointer,  intent(in)    :: val(:,:)
    integer(ip),      optional, intent(out)   :: ierr
    integer(ip)                               :: n1,n2
   
    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       
       select type ( v => self % val2 )
       type is ( integer(kind=ip) )
          if( size(v,1) /= size(val,1) .or. size(v,2) /= size(val,2) ) then          
             n1           = min(size(v,1),size(val,1))
             n2           = min(size(v,2),size(val,2))
             v(1:n1,1:n2) = val(1:n1,1:n2)
             if( present(ierr) ) ierr = 1
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
    
  end subroutine set2_ip
  
  subroutine sets_rp(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    real(rp),        pointer,  intent(in)    :: val
    integer(ip),     optional, intent(out)   :: ierr

    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       select type ( v => self % vals )
       type is ( real(kind=rp) ) ; v = val
       end select
    else
       if( present(ierr) ) ierr = 2    
    end if

  end subroutine sets_rp
  
  subroutine set1_rp(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    real(rp),        pointer,  intent(in)    :: val(:)
    integer(ip),     optional, intent(out)   :: ierr
    integer(ip)                              :: n1

    if( associated(val) ) then       
       
       if( present(ierr) ) ierr = 0
       select type ( v => self % val1 )
       type is ( real(kind=rp) )
          if( size(v) /= size(val) ) then
             n1 = min(size(v),size(val))
             v(1:n1) = val(1:n1)
             if( present(ierr) ) ierr = 1
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
    
  end subroutine set1_rp
  
  subroutine set2_rp(self,val,ierr)
    class(node_typ),            intent(inout) :: self
    real(rp),         pointer,  intent(in)    :: val(:,:)
    integer(ip),      optional, intent(out)   :: ierr
    integer(ip)                               :: n1,n2

    if( associated(val) ) then       

       if( present(ierr) ) ierr = 0

       select type ( v => self % val2 )
       type is ( real(kind=rp) )
          if( size(v,1) /= size(val,1) .or. size(v,2) /= size(val,2) ) then          
             n1           = min(size(v,1),size(val,1))
             n2           = min(size(v,2),size(val,2))
             v(1:n1,1:n2) = val(1:n1,1:n2)
             if( present(ierr) ) ierr = 1
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if

  end subroutine set2_rp
  
  subroutine sets_lg(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    logical(lg),     pointer,  intent(in)    :: val
    integer(ip),     optional, intent(out)   :: ierr

    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       select type ( v => self % vals )
       type is ( logical(kind=lg) ) ; v = val
       end select
    else
       if( present(ierr) ) ierr = 2
    end if

  end subroutine sets_lg
  
  subroutine set1_lg(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    logical(lg),     pointer,  intent(in)    :: val(:)
    integer(ip),     optional, intent(out)   :: ierr
    integer(ip)                              :: n1

    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       select type ( v => self % val1 )
       type is ( logical(kind=lg) )
          if( size(v) /= size(val) ) then
             n1 = min(size(v),size(val))
             v(1:n1) = val(1:n1)
             if( present(ierr) ) ierr = 1
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if

  end subroutine set1_lg
  
  subroutine set2_lg(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    logical(lg),     pointer,  intent(in)    :: val(:,:)
    integer(ip),     optional, intent(out)   :: ierr
    integer(ip)                              :: n1,n2
    
    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       select type ( v => self % val2 )
       type is ( logical(kind=lg) )
          if( size(v,1) /= size(val,1) .or. size(v,2) /= size(val,2) ) then          
             n1           = min(size(v,1),size(val,1))
             n2           = min(size(v,2),size(val,2))
             v(1:n1,1:n2) = val(1:n1,1:n2)
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
    
  end subroutine set2_lg
  
  subroutine sets_ch(self,val,ierr)
    class(node_typ),           intent(inout) :: self
    character(len=:), pointer, intent(in)    :: val
    integer(ip),     optional, intent(out)   :: ierr
    integer(ip)                              :: n1,n2

    if( associated(val) ) then
       if( present(ierr) ) ierr = 0
       select type ( v => self % vals )
       type is ( character(len=*) )
          if( len(v) /= len(val) ) then
             n1      = int(len(v),ip)
             n2      = int(len(val),ip)
             n1      = min(n1,n2)
             v(1:n1) = val(1:n1)
          else
             v = val
          end if
       end select
    else
       if( present(ierr) ) ierr = 2
    end if
   
  end subroutine sets_ch
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Get dimension
  !> @details Get dimension
  !> 
  !-----------------------------------------------------------------------

  pure function getdim(self) result(nn)
    class(node_typ), intent(in) :: self
    integer(ip)                 :: nn(2)
    integer(ip)                 :: n1
    
    if( associated(self % vals) ) then
       if( self % TYPE == TYPE_CHARACTER ) then
          select type ( v => self % vals )
          type is ( character(len=*) )
             n1 = int(len(v,KIND=ip),ip)
             nn = (/ n1,0_ip /)
          end select
       else
          nn = (/ 1_ip,0_ip /)
       end if
    else if( associated(self % val1) ) then
       nn = (/ size(self % val1,KIND=IP),0_ip /)
    else if( associated(self % val2) ) then
       nn = (/ size(self % val2,1,KIND=IP) , size(self % val2,2,KIND=IP) /)
    end if
    
  end function getdim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Check if a value exists
  !> @details Check if a value exists
  !> 
  !-----------------------------------------------------------------------
  
  logical(lg) pure function exists(self) 
    
    class(node_typ), intent(in) :: self

    if( self % assigned == 0 ) then
       exists = .false.
    else
       exists = .true.
    end if
    
  end function exists  

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-01
  !> @brief   Deallocate node
  !> @details Deallocate node
  !> 
  !-----------------------------------------------------------------------
  
  subroutine deallo(self)
    class(llist_typ), intent(inout) :: self
    class(node_typ),  pointer       :: node

    node => self % head
    do while( associated(node) )
       if( associated(node % vals) ) deallocate(node % vals)
       if( associated(node % val1) ) deallocate(node % val1)
       if( associated(node % val2) ) deallocate(node % val2) 
       if( allocated (node % name) ) deallocate(node % name)
       node => node % next
    end do
    
  end subroutine deallo
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-01
  !> @brief   Untarget
  !> @details Untarget to values
  !> 
  !-----------------------------------------------------------------------
  
  subroutine untarget(self)
    class(llist_typ), intent(inout) :: self
    class(node_typ),  pointer       :: node

    node => self % head
    do while( associated(node) )
       nullify(node % vals)
       nullify(node % val1)
       nullify(node % val2)
       if( allocated(node % name) ) deallocate(node % name)
       node => node % next
    end do
    
  end subroutine untarget
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-01
  !> @brief   Iterator
  !> @details Iterator
  !> 
  !-----------------------------------------------------------------------

  subroutine iterate(self,what,ivari)

    class(llist_typ),         intent(inout) :: self
    procedure(what_to_do)                   :: what
    integer(ip),     optional, intent(in)   :: ivari
    class(node_typ), pointer                :: current
    class(node_typ), pointer                :: temp
    
    current => self % head
    do while( associated(current) )
       temp => current % next
       call what(current,ivari)
       current => temp
    end do
    
  end subroutine iterate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-01
  !> @brief   Add
  !> @details Add a value
  !> 
  !-----------------------------------------------------------------------

  subroutine add(llist,node) 
    
    class(llist_typ),          intent(inout) :: llist
    class(node_typ),  pointer, intent(inout) :: node
    
    llist % size = llist % size + 1
    
    if( .not. associated(llist % head) ) then
       llist % head        => node
       llist % tail        => node
       nullify(llist % tail % prev)
    else
       llist % tail % prev => llist % tail
       llist % tail % next => node
       llist % tail        => node
    end if
    
    llist % tail % id = llist % size + 1

  end subroutine add

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-01
  !> @brief   Add
  !> @details Add a value
  !> 
  !-----------------------------------------------------------------------

  function nameis(self,name) result(node)

    class(llist_typ), intent(in) :: self
    character(LEN=*), intent(in) :: name
    class(node_typ),  pointer    :: node
    
    node => self % head
    do while( associated(node) )
       if( node % name == name ) return          
       node => node % next
    end do
    nullify(node)

  end function nameis
    
end module def_linked_list
!> @}
