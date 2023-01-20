!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to METIS 4 and 5
!> @{
!> @file    mod_alya2metis.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to METIS
!> @details Bridge to usefull functions of METIS 4, 5.0.2 and 5.1.0
!>          WANRING: Integer types of Alya and METIS must correspond.
!>          The following subroutines are interfaced:
!>          METIS_SetDefaultOptions (only METIS5)
!>          METIS_NodeND
!>          METIS_PartGraphRecursive
!>          METIS_PartGraphKway
!>          Manual available at:
!>          http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf
!-----------------------------------------------------------------------

module mod_alya2metis

  use  def_kintyp, only : ip
  use  def_master
  use, intrinsic       :: ISO_C_BINDING
  
  implicit none
    ! 
    ! Size of integers in C same as thos of Alya
    !
#ifdef I8
#define C_SIZE C_LONG
#else
#define C_SIZE C_INT
#endif
  !
  ! METIS parameters
  !
#if   defined V5METIS
    integer(C_SIZE)        :: optio_par(0:40)
    integer(ip), parameter :: METIS_OPTION_NUMBERING=16
#elif defined V51METIS
    integer(C_SIZE)        :: optio_par(0:40)
    integer(ip), parameter :: METIS_OPTION_NUMBERING=17
#elif defined METIS
    integer(C_SIZE)        :: optio_par(1:8)
#else
    integer(ip)            :: optio_par(1:8) 
#endif
    
  private
  !
  ! C interfaces for METIS 5.0.2 and METIS 5.1.0
  !
#if defined V5METIS || defined V51METIS
  interface
     subroutine METIS_SetDefaultOptions(opts) bind(C,name="METIS_SetDefaultOptions")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_SIZE)               :: opts(0:40)
     end subroutine METIS_SetDefaultOptions
  end interface
  interface
     subroutine METIS_NodeND(nvtxs,xadj,adjncy,vwgt,opts,perm,iperm) bind(C,name="METIS_NodeND")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_SIZE)               :: nvtxs
       integer(C_SIZE), dimension(*) :: xadj
       integer(C_SIZE), dimension(*) :: adjncy
       integer(C_SIZE), dimension(*) :: perm
       integer(C_SIZE), dimension(*) :: iperm
       type(C_PTR),     value        :: vwgt
       integer(C_SIZE)               :: opts(0:40)
     end subroutine METIS_NodeND
  end interface
  interface
     subroutine METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,opts,objval,part) bind(C,name="METIS_PartGraphRecursive")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_SIZE)               :: nvtxs
       integer(C_SIZE)               :: ncon
       integer(C_SIZE), dimension(*) :: xadj
       integer(C_SIZE), dimension(*) :: adjncy
       integer(C_SIZE), dimension(*) :: vwgt
       type(C_PTR),     value        :: vsize
       type(C_PTR),     value        :: adjwgt
       integer(C_SIZE)               :: nparts
       type(C_PTR),     value        :: tpwgts
       type(C_PTR),     value        :: ubvec
       integer(C_SIZE)               :: opts(0:40)
       integer(C_SIZE)               :: objval
       integer(C_SIZE), dimension(*) :: part
     end subroutine METIS_PartGraphRecursive
  end interface
  interface
     subroutine METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,opts,objval,part) bind(C,name="METIS_PartGraphKway")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_SIZE)               :: nvtxs
       integer(C_SIZE)               :: ncon
       integer(C_SIZE), dimension(*) :: xadj
       integer(C_SIZE), dimension(*) :: adjncy
       integer(C_SIZE), dimension(*) :: vwgt
       type(C_PTR),     value        :: vsize
       type(C_PTR),     value        :: adjwgt
       integer(C_SIZE)               :: nparts
       type(C_PTR),     value        :: tpwgts
       type(C_PTR),     value        :: ubvec
       integer(C_SIZE)               :: opts(0:40)
       integer(C_SIZE)               :: objval
       integer(C_SIZE), dimension(*) :: part
     end subroutine METIS_PartGraphKway
  end interface
  !
  ! C interfaces for METIS 4.0
  !
#elif defined METIS
  interface
     subroutine METIS_NodeND(nvtxs,xadj,adjncy,npes,opts,perm,iperm) bind(C,name="METIS_NodeND")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_SIZE)               :: nvtxs
       integer(C_SIZE), dimension(*) :: xadj
       integer(C_SIZE), dimension(*) :: adjncy
       integer(C_SIZE)               :: npes
       integer(C_SIZE)               :: opts(1:8)
       integer(C_SIZE), dimension(*) :: perm
       integer(C_SIZE), dimension(*) :: iperm
     end subroutine METIS_NodeND
  end interface
#endif

  public :: alya2metis_initialization  ! Module initialization
  public :: alya2metis_METIS_PartGraph ! Graph partitioning
  public :: alya2metis_METIS_NodeND    ! Node numbering
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-01-07
  !> @brief   METIS Initialization
  !> @details METIS Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2metis_initialization()

#if defined  V5METIS || defined V51METIS
    call METIS_SetDefaultOptions(optio_par)
    optio_par(METIS_OPTION_NUMBERING) = 1_ip
#elif defined METIS
    optio_par = 0_ip
#endif
   
  end subroutine alya2metis_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-01-07
  !> @brief   METIS part graph
  !> @details METIS graph partititioning
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2metis_METIS_PartGraph(nparts,nn,xadj,adj,weight,lpart)

#if defined  V5METIS || defined V51METIS
    integer(C_SIZE), intent(in)             :: nparts
    integer(C_SIZE), intent(in)             :: nn
    integer(C_SIZE), intent(in),    pointer :: xadj(:)
    integer(C_SIZE), intent(in),    pointer :: adj(:)
    integer(C_SIZE), intent(inout), pointer :: weight(:)
    integer(C_SIZE), intent(inout), pointer :: lpart(:)
    integer(C_SIZE)                         :: ncon
    integer(C_SIZE)                         :: edgecut
    type(C_PTR)                             :: vsize
    type(C_PTR)                             :: adjwgt
    type(C_PTR)                             :: tpwgts
    type(C_PTR)                             :: ubvec
#else
    integer(ip),     intent(in)             :: nparts       !< Number of partitions
    integer(ip),     intent(in)             :: nn           !< Size of graph
    integer(ip),     intent(in),    pointer :: xadj(:)  
    integer(ip),     intent(in),    pointer :: adj(:)
    integer(ip),     intent(inout), pointer :: weight(:)
    integer(ip),     intent(inout), pointer :: lpart(:)
    integer(ip)                             :: wedge_par(1)
    integer(ip)                             :: kfl_fortr_par
    integer(ip)                             :: kfl_weigh_par
    integer(ip)                             :: edgecut
#endif
    integer(ip)                             :: ii

    if( nn <  1 ) then
       return
    else if( nn <= 1 .or. nparts <= 1 ) then
       lpart = 1_ip
       return
    else if( nparts > nn ) then
       do ii = 1,nn
          lpart(ii) = ii
       end do
       return
    else

#if defined V5METIS || defined V51METIS    
       !
       ! METIS 5.0.2
       !
       ncon   = 1_ip
       vsize  = C_NULL_PTR
       adjwgt = C_NULL_PTR
       tpwgts = C_NULL_PTR
       ubvec  = C_NULL_PTR
       if( nparts <= 8 ) then
          call METIS_PartGraphRecursive( &
               nn, ncon, xadj, adj, weight,  vsize, &
               adjwgt, nparts, tpwgts, ubvec,       &
               optio_par, edgecut, lpart )
       else
          call METIS_PartGraphKway(&
               nn, ncon, xadj, adj, weight, vsize,  &
               adjwgt, nparts, tpwgts, ubvec,       &
               optio_par, edgecut, lpart )
       end if
       
#elif defined METIS
       !
       ! METIS 4
       !
       kfl_fortr_par = 1_ip
       kfl_weigh_par = 2_ip
       wedge_par     = 0_ip
       optio_par(1)  = 0_ip 
       optio_par(2)  = 3_ip
       optio_par(3)  = 1_ip
       optio_par(4)  = 1_ip
       optio_par(5)  = 0_ip       

       if( nparts <= 8 ) then
          call METIS_PartGraphRecursive( &
               nn, xadj, adj, weight, wedge_par,     &
               kfl_weigh_par, kfl_fortr_par, nparts, &
               optio_par, edgecut, lpart )
       else
          call METIS_PartGraphKway(&
               nn, xadj, adj, weight, wedge_par,     &
               kfl_weigh_par, kfl_fortr_par, nparts, &
               optio_par, edgecut, lpart )
       end if

#else
       !
       ! Nothing
       !
       call runend('TO USE METIS, PLEASE COMPILE WITH PROPER MACRO')

#endif

    end if

  end subroutine alya2metis_METIS_PartGraph

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-01-07
  !> @brief   Bridge to METIS
  !> @details METIS invert the signification of permr and invor with respect to Alya!
  !>          In Alya: NEW = PERMR(OLD) AND OLD = INVPR(NEW). This is why the  permutations
  !>          and inverse permutations are inverted in METIS call
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2metis_METIS_NodeND(nn,xadj,adj,permr,invpr_)

#if defined V5METIS || defined V51METIS    
    integer(C_SIZE), intent(in)                       :: nn
    integer(C_SIZE), intent(in),    pointer           :: xadj(:)
    integer(C_SIZE), intent(in),    pointer           :: adj(:)
    integer(C_SIZE), intent(inout), pointer           :: permr(:)
    integer(C_SIZE), intent(inout), pointer, optional :: invpr_(:)
    type(C_PTR)                                       :: vwgt
    integer(C_SIZE),                pointer           :: invpr(:)
#elif defined  METIS
    integer(C_SIZE), intent(in)                       :: nn
    integer(C_SIZE), intent(in),    pointer           :: xadj(:)
    integer(C_SIZE), intent(in),    pointer           :: adj(:)
    integer(C_SIZE), intent(inout), pointer           :: permr(:)
    integer(C_SIZE), intent(inout), pointer, optional :: invpr_(:)
    integer(C_SIZE)                                   :: kfl_fortr_par
    integer(C_SIZE),                pointer           :: invpr(:) 
#else
    integer(ip),     intent(in)                       :: nn        !< Number of vertices
    integer(ip),     intent(in),    pointer           :: xadj(:)   !< Graph
    integer(ip),     intent(in),    pointer           :: adj(:)
    integer(ip),     intent(inout), pointer           :: permr(:)
    integer(ip),     intent(inout), pointer, optional :: invpr_(:) 
    integer(ip)                                       :: kfl_fortr_par
    integer(ip),                    pointer           :: invpr(:) 
#endif

    if( nn <= 0 ) return

    if( present(invpr_) ) then
       invpr => invpr_
    else
       allocate(invpR(nn))
    end if

    if( nn == 1 ) then
       permr(1) = 1
       invpr(1) = 1
 
    else

#if defined V5METIS || defined V51METIS    
       !
       ! METIS 5.0.2
       !
       vwgt = C_NULL_PTR
       call METIS_NodeND(nn,xadj,adj,vwgt,optio_par,invpr,permr)

#elif defined METIS
       !
       ! METIS 4
       !
       kfl_fortr_par = 1_ip
       optio_par(1)  = 1_ip
       optio_par(2)  = 3_ip
       optio_par(3)  = 2_ip
       optio_par(4)  = 2_ip
       optio_par(5)  = 0_ip
       optio_par(6)  = 1_ip 
       optio_par(7)  = 0_ip
       optio_par(8)  = 1_ip
       call METIS_NodeND(nn,xadj,adj,kfl_fortr_par,optio_par,invpr,permr)

#else    
       !
       ! Nothing
       !
       call runend('TO USE METIS, PLEASE COMPILE WITH PROPER MACRO')

#endif

    end if

    if( .not. present(invpr_) ) then
       deallocate(invpr)
    end if

  end subroutine alya2metis_METIS_NodeND

end module mod_alya2metis
!> @}
