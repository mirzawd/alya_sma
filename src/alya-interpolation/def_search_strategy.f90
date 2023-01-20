!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_search_strategy.f90
!> @author  houzeaux
!> @date    2020-10-02
!> @brief   Class for defining a search strategy
!> @details Class for defining a search strategy
!>
!>                                 Spatial access methods
!>
!>          ------------------------------------------------------------------
!>            Binary space (2 arbitrary divisions)      8 Cartesian divisions
!>                         
!>          ---------------------------------           ----------------------
!>          2 cartesian divisions, successive           
!>          levels are split along different 
!>          dimensions
!>
!>                 Adaptive methods
!>                 ----------------
!>                 Split  divides the set of 
!>                 points into two sets of
!>                 (nearly) equal size
!>
!>          kd-tree            skd-tree                  Octree 
!>          (Non overlapping)  (BIH, overlapping)        (Non overlapping)
!>
!>          https://commons.wikimedia.org/wiki/File:Kd_tree_vs_skd_tree.svg
!>
!>                          SPARSE                          DENSE
!>                   +----------|---------+         +--------------------+  
!>                   |          |         |         |                    |
!>                   | o------o |         |         | o------o           |
!>                   | |      | |         |         | |      |           |
!>          KD-TREE  | o------o |         |         | o------o  ???      |
!>                   |          |         |         |                    |
!>                   |          |o------o |         |      o------o      |
!>                   |          ||      | |         |      |      |      |
!>                   |          |o------o |         |      o------o      |
!>                   +----------|---------+         +--------------------+
!>                   split is unique                objects are duplicated
!>
!>                   +--------|--|---------+         +------|-|-----------+  
!>                   |        |  |         |         |      | |           |
!>                   | o------o  |         |         | o----|-o           |
!>                   | |      |  |         |         | |    | |           |
!>          SKD-TREE | o------o  |         |         | o----|-o           |
!>                   |        |  |         |         |      | |           |
!>                   |        |  o------o  |         |      o-|----o      |
!>                   |        |  |      |  |         |      | |    |      |
!>                   |        |  o------o  |         |      o-|----o      |
!>                   +--------|--|---------+         +------|-|-----------+
!>                   disjoint children               children overlap
!>
!-----------------------------------------------------------------------

module def_search_strategy
  
  use def_kintyp_basic,          only : ip,rp
  use def_search_method,         only : search_method
  use def_search_method,         only : SEARCH_BIN
  use def_search_method,         only : SEARCH_OCTREE 
  use def_search_method,         only : SEARCH_KDTREE 
  use def_search_method,         only : SEARCH_SKDTREE 
  use def_maths_bin,             only : maths_bin
  use def_maths_tree,            only : maths_octree
  use def_maths_tree,            only : maths_kdtree
  use def_maths_tree,            only : maths_skdtree
  use mod_communications_global, only : PAR_SUM
  use mod_strings,               only : integer_to_string
  use mod_exchange,              only : exchange_init
  use mod_exchange,              only : exchange_add
  use mod_exchange,              only : exchange_end
  implicit none
 
  private
  
  type search
     integer(ip)                       :: type      ! Type of search method
     real(rp)                          :: param(10) ! Input parameters
     class(search_method), pointer     :: method    ! Search method class
   contains
     procedure,            pass        :: init
     procedure,            pass        :: set
     procedure,            pass        :: deallo
     procedure,            pass        :: input
     procedure,            pass        :: fill
     procedure,            pass        :: read_data
     procedure,            pass        :: parall
  end type search

  public :: search
  public :: SEARCH_BIN
  public :: SEARCH_OCTREE 
  public :: SEARCH_KDTREE 
  public :: SEARCH_SKDTREE 
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-14
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    class(search) :: self

    self % type = -1
    nullify(self % method)
    
  end subroutine init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-14
  !> @brief   Set the search method
  !> @details Set the search method
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,type_search,param,coord,bobox,name)
    
    class(search),                       intent(inout) :: self
    integer(ip),      optional,          intent(in)    :: type_search
    real(rp),         optional,          intent(in)    :: param(:)
    real(rp),         optional, pointer, intent(in)    :: coord(:,:)   !< Coordinates
    real(rp),         optional, pointer, intent(in)    :: bobox(:,:,:) !< Bounding boxes
    character(len=*), optional,          intent(in)    :: name
    integer(ip)                                        :: nmin
    !
    ! Set search type
    !
    if( present(type_search) ) self % type = type_search
    !
    ! Allocate search method
    !
    select case ( self % type )
    case ( SEARCH_BIN     ) ; allocate( maths_bin     :: self % method )
    case ( SEARCH_OCTREE  ) ; allocate( maths_octree  :: self % method )
    case ( SEARCH_KDTREE  ) ; allocate( maths_skdtree :: self % method )
    case ( SEARCH_SKDTREE ) ; allocate( maths_skdtree :: self % method )
    case default            ; call runend('DEF_SEARCH_STRATEGY: SEARCH METHOD NOT DEFINED: '//integer_to_string(self % type))
    end select
    !
    ! Initialize
    !
    call self % method % init()
    !
    ! Set input
    !
    if( present(param) ) then
       nmin = int(min(size(param),size(self % param)),ip)
       self % param(1:nmin) = param(1:nmin)
    end if
    call self % input(NAME=name)
    !
    ! Fill
    !
    if( present(coord) .or. present(bobox) ) then
       call self % fill(coord,bobox)
    end if
    
  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-14
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    class(search) :: self 

    if( associated(self % method) ) then
       call self % method % deallo()
       deallocate(self % method)
    end if
    
  end subroutine deallo
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-14
  !> @brief   Set input
  !> @details Set input
  !> 
  !-----------------------------------------------------------------------

  subroutine input(self,name)
    class(search)                          :: self    
    character(len=*), optional, intent(in) :: name
   
    select type ( v => self % method )
    class is ( maths_bin     ) ; call v % input(PARAM=self % param,NAME=name)
    class is ( maths_octree  ) ; call v % input(PARAM=self % param,NAME=name)
    class is ( maths_kdtree  ) ; call v % input(PARAM=self % param,NAME=name)
    class is ( maths_skdtree ) ; call v % input(PARAM=self % param,NAME=name)
    class default              ; call runend('DEF_SEARCH_STRATEGY: SEARCH METHOD NOT DEFINED')
    end select
                 
    !select type ( v => self % method )
    !class is ( maths_bin    ) ; print*,'aaaaaaaaaaaaaaaaaaaaaaaa=',v % boxip
    !class is ( maths_octree ) ; print*,'bbbbbbbbbbbbbbbbbbbbbbbbb=',v % limit
    !class default             ; call runend('DEF_SEARCH_STRATEGY: SEARCH METHOD NOT DEFINED')
    !end select
                 
  end subroutine input
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-14
  !> @brief   Fill search method
  !> @details Fill search method
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(self,coord,bobox)
    
    class(search)                              :: self    
    real(rp),    optional, pointer, intent(in) :: coord(:,:)   !< Coordinates
    real(rp),    optional, pointer, intent(in) :: bobox(:,:,:) !< Bounding boxes

    select type ( v => self % method )
    class is ( maths_bin     ) ; call v % fill(coord,bobox)
    class is ( maths_octree  ) ; call v % fill(coord,bobox)
    class is ( maths_kdtree  ) ; call v % fill(coord,bobox)
    class is ( maths_skdtree ) ; call v % fill(coord,bobox)
    class default              ; call runend('DEF_SEARCH_STRATEGY: SEARCH METHOD NOT DEFINED')
    end select
    
  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Read data
  !> @details Read data
  !> 
  !-----------------------------------------------------------------------

  subroutine read_data(self,words,param)

    class(search),        intent(inout) :: self    
    character(len=5),     intent(in)    :: words(:)
    real(rp),             intent(in)    :: param(:)
    
    select case( words(2) )

    case ( 'BIN  ' ) 
       
       self % type = SEARCH_BIN 
       if( words(3) == 'BOXES' ) self % param(1) = param(3)
       self % param(2:3) = param(3)
       
    case ( 'OCTTR' , 'OCT  ' , 'OCTRE' )                

       self % type = SEARCH_OCTREE
       if( words(3) == 'MAXIM' ) self % param(1) = param(3)
       if( words(3) == 'LIMIT' ) self % param(1) = param(3)

    case ( 'SKDTR' , 'KDTRE' )                

       self % type = SEARCH_KDTREE

    end select

  end subroutine read_data

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Broadcast data
  !> @details Broadcast data in parallel 
  !> 
  !-----------------------------------------------------------------------

  subroutine parall(self)

    class(search), intent(inout) :: self    

    call exchange_init()
    call exchange_add(self % type)
    call exchange_add(self % param)
    call exchange_end()
    
  end subroutine parall
  
end module def_search_strategy
!> @}
