!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling functions
!> @file    mod_couplings.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   Setup coupling
!> @details Setup coupling
!>
!> @{
!------------------------------------------------------------------------

module mod_couplings_setup

  use def_kintyp_basic,      only : ip,rp,lg,r1p,r2p,i1p,i2p
  use def_master,            only : ISEQUEN
  use def_master,            only : INOTMASTER
  use def_master,            only : IMASTER,kfl_paral,lninv_loc
  use mod_strings,           only : integer_to_string
  use def_master,            only : current_code
  use def_master,            only : current_zone
  use def_master,            only : zeror,namda
  use def_master,            only : intost
  use def_domain,            only : lesub,nelem,lnnod
  use def_domain,            only : mnodb,mnode,ndime
  use def_domain,            only : nbono,coord,ltopo,ltype
  use def_domain,            only : lbono,npoin,lnoch
  use def_domain,            only : nboun,lnodb,lnnob
  use def_domain,            only : lelch,meshe,lbset
  use def_domain,            only : nnode,lelbo,lnods
  use def_kermod,            only : ielse,relse,ndivi
  use mod_elmgeo,            only : elmgeo_natural_coordinates
  use mod_elmgeo,            only : elmgeo_natural_coordinates_on_boundaries
  use mod_elsest,            only : elsest_host_element
  use mod_elsest,            only : ELSEST_A_TODA_COSTA
  use mod_parall,            only : PAR_COMM_COLOR_PERM
  use mod_parall,            only : par_part_in_color
  use mod_parall,            only : par_code_zone_subd_to_color
  use mod_parall,            only : color_target
  use mod_parall,            only : color_source
  use mod_parall,            only : par_bin_comin
  use mod_parall,            only : par_bin_comax
  use mod_parall,            only : par_bin_part
  use mod_parall,            only : par_bin_boxes
  use mod_parall,            only : par_bin_size
  use mod_parall,            only : PAR_COMM_CURRENT
  use mod_parall,            only : PAR_COMM_COLOR
  use mod_parall,            only : PAR_COMM_COLOR_ARRAY
  use mod_parall,            only : I_AM_IN_COLOR
  use mod_parall,            only : PAR_MY_CODE_RANK
  use mod_parall,            only : par_part_comin
  use mod_parall,            only : par_part_comax
  use mod_parall,            only : PAR_WORLD_SIZE
  use mod_parall,            only : PAR_MY_WORLD_RANK
  use mod_parall,            only : typ_bin_structure
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_size
  use mod_memory,            only : memory_resize
  use mod_kdtree,            only : typ_kdtree
  use mod_kdtree,            only : kdtree_nearest_boundary
  use mod_htable,            only : hash_t
  use mod_htable,            only : htaini
  use mod_htable,            only : htaadd
  use mod_htable,            only : htades
  use mod_communications,    only : PAR_MIN
  use mod_communications,    only : PAR_SUM
  use mod_communications,    only : PAR_SEND_RECEIVE
  use mod_communications,    only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,    only : PAR_BARRIER
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_START_NON_BLOCKING_COMM
  use mod_communications,    only : PAR_END_NON_BLOCKING_COMM
  use mod_communications,    only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications,    only : PAR_GATHER
  use mod_communications,    only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications,    only : PAR_ALLTOALL
  use def_coupli,            only : typ_color_coupling
  use def_coupli,            only : kdtree_typ
  use def_coupli,            only : mcoup
  use def_coupli,            only : coupling_type
  use def_coupli,            only : ELEMENT_INTERPOLATION
  use def_coupli,            only : BOUNDARY_INTERPOLATION
  use def_coupli,            only : NEAREST_BOUNDARY_NODE
  use def_coupli,            only : NEAREST_ELEMENT_NODE
  use def_coupli,            only : BOUNDARY_VECTOR_PROJECTION
  use def_coupli,            only : GLOBAL_NUMBERING
  use def_coupli,            only : SAME_COORDINATE
  use def_coupli,            only : memor_cou
  use def_coupli,            only : STRESS_PROJECTION
  use def_coupli,            only : PROJECTION
  use def_coupli,            only : BETWEEN_SUBDOMAINS
  use def_coupli,            only : lnodb_cou
  use def_coupli,            only : ltypb_cou
  use def_coupli,            only : lnnob_cou
  use def_coupli,            only : toler_relative_cou
  use def_coupli,            only : toler_absolute_cou
  use mod_iofile,            only : iofile_open_unit
  use mod_iofile,            only : iofile_close_unit
  use mod_iofile,            only : iofile_available_unit
  use mod_iofile,            only : iofile_flush_unit
  use mod_messages,          only : messages_live
  use def_coupli,            only : VALUES_ON_ELEMENTS
  use def_coupli,            only : VALUES_ON_BOUNDARIES
  use def_coupli,            only : VALUES_ON_NODES
  use mod_coupling_toolbox,  only : coupling_toolbox_points_in_partitions
  use mod_par_bin_structure, only : par_bin_structure
  use mod_par_bin_structure, only : par_bin_structure_deallocate
  use mod_par_bin_structure, only : par_bin_structure_initialization
  use mod_par_bin_structure, only : par_bin_structure_partition_bounding_box
  use def_mpi
  use mod_std
#include "def_mpi.inc"
  
  private

  interface COU_INIT_INTERPOLATE_POINTS_VALUES
     module procedure COU_INIT_INTERPOLATE_POINTS_VALUES_OLD,&
          &           COU_INIT_INTERPOLATE_POINTS_VALUES_NEW
  end interface COU_INIT_INTERPOLATE_POINTS_VALUES

  public :: COU_INIT_INTERPOLATE_POINTS_VALUES
  
  integer(ip),   parameter :: my_huge = huge(1_ip)
  character(19), parameter :: vacal='mod_couplings_setup'

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Interpolation communication arrays
  !> @details Interpolation communication arrays
  !>          Input variables:
  !>          COUPLING ............................... Coupling structure
  !>          XX(NDIME,:) ............................ Coordinates of wet points
  !>          COUPLING % GEOME % NUMBER_WET_POINTS ... SIZE(XX,2)
  !>          COUPLING % ITYPE ....................... Vector projection
  !>          COUPLING % KIND ........................ BETWEEN_SUBDOMAINS/BETWEEN_ZONES
  !>          COUPLING % KDTREE ...................... Kdtree of source
  !>
  !>
  !----------------------------------------------------------------------

  subroutine  COU_INIT_INTERPOLATE_POINTS_VALUES_NEW(coupling,COMM,CANDIDATE_SOURCE_NODES, &
       & CANDIDATE_SOURCE_ELEMENTS, CANDIDATE_SOURCE_BOUNDARIES,all_interp_)

    type(typ_color_coupling), intent(inout)          :: coupling
    MY_MPI_COMM   ,           intent(in),   optional :: COMM
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_NODES(:)
    integer(ip), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_ELEMENTS(:)
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_BOUNDARIES(:)
    logical(lg),              intent(in),   optional :: all_interp_
    logical(lg)                                      :: all_interp
    integer(ip)                                      :: icolo_source
    integer(ip)                                      :: icolo_target

    icolo_source = coupling % color_source
    icolo_target = coupling % color_target
    all_interp = .false.
    if(present(all_interp_)) all_interp=all_interp_

    call COU_INIT_INTERPOLATE_POINTS_VALUES_OLD(coupling % wet % coord_wet,icolo_target,icolo_source,coupling,COMM,&
         CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,CANDIDATE_SOURCE_BOUNDARIES,all_interp)   

  end subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_NEW

  subroutine  COU_INIT_INTERPOLATE_POINTS_VALUES_OLD(xx,icolo,jcolo,coupling,COMM,CANDIDATE_SOURCE_NODES, &
       CANDIDATE_SOURCE_ELEMENTS, CANDIDATE_SOURCE_BOUNDARIES,all_interp_)

    integer(ip),              intent(in)             :: icolo
    integer(ip),              intent(in)             :: jcolo
    real(rp),    pointer,     intent(in)             :: xx(:,:)
    type(typ_color_coupling), intent(inout)          :: coupling
    MY_MPI_COMM   ,           intent(in),   optional :: COMM
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_NODES(:)
    integer(ip), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_ELEMENTS(:)
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_BOUNDARIES(:)    
    logical(lg),              intent(in),   optional :: all_interp_
    logical(lg)                                      :: all_interp
    integer(ip)                                      :: ipart,icoup
    integer(ip)                                      :: number_wet_points
    integer(ip)                                      :: cz,cpart
    integer(ip)                                      :: idime,ipoin,itype
    integer(ip)                                      :: PAR_CURRENT_SIZE
    integer(ip)                                      :: PAR_CURRENT_RANK
    MY_MPI_COMM                                      :: PAR_COMM_SAVE
    MY_MPI_COMM                                      :: PAR_COMM_CURRENT4
    type(r2p),   pointer                             :: shapf(:)
    integer(ip), pointer                             :: npoin_send(:)
    integer(ip), pointer                             :: npoin_recv(:)
    integer(ip), pointer                             :: mask_npoin(:)
    type(i1p),   pointer                             :: decision_send(:)
    type(i1p),   pointer                             :: decision_recv(:)
    type(r1p),   pointer                             :: distance_send(:)      ! Distance check
    type(r1p),   pointer                             :: distance_recv(:)      ! Distance check
    type(i1p),   pointer                             :: my_part_to_point(:)
    type(i1p),   pointer                             :: my_point_to_part(:)
    type(r1p),   pointer                             :: coord_send(:)
    type(r1p),   pointer                             :: coord_recv(:)
    type(i1p),   pointer                             :: numer_send(:)
    type(i1p),   pointer                             :: numer_recv(:)
    integer(ip), pointer                             :: lesou(:)
    logical(lg), pointer                             :: lbsou(:)
    logical(lg), pointer                             :: lnsou(:)
    logical(lg), pointer                             :: intersection(:)
    integer(ip), pointer                             :: PAR_WORLD_RANKS(:)
    logical(lg)                                      :: require_distance
    real(rp)                                         :: time0,time1,time2,time3,time4,time5
    integer(4)                                       :: num_intersections4
    type(typ_bin_structure)                          :: bin_structure
    real(rp)                                         :: comin(3),comax(3),delta(3)

    call cputim(time0)
    all_interp = .false.
    if(present(all_interp_)) all_interp=all_interp_
    !
    ! Coupling
    !
    icoup        = coupling % number
    color_target = coupling % color_target
    color_source = coupling % color_source
    itype        = coupling % itype
    if(all_interp) itype = BOUNDARY_INTERPOLATION
    !
    ! Communicator
    !
    PAR_COMM_SAVE = PAR_COMM_CURRENT
    if( present(COMM) ) then
       PAR_COMM_CURRENT = COMM
    else
       PAR_COMM_CURRENT = PAR_COMM_COLOR(icolo,jcolo)
    end if
    PAR_COMM_CURRENT4 = PAR_COMM_CURRENT
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    !
    ! Sizes
    !
    number_wet_points = memory_size(xx,2_ip)
    if( coupling % wet % number_wet_points /= 0 ) then
       if( coupling % wet % number_wet_points /= number_wet_points ) &
            call runend('SOMETHING STRANGE HAPPENS: '//intost(coupling % wet % number_wet_points)//', '//intost(number_wet_points))
    else
       coupling % wet % number_wet_points = number_wet_points
    end if
    !
    ! Nullify and initialize
    !
    nullify(shapf)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(mask_npoin)
    nullify(decision_send)
    nullify(decision_recv)
    nullify(distance_send)
    nullify(distance_recv)
    nullify(my_part_to_point)
    nullify(my_point_to_part)
    nullify(coord_send)
    nullify(coord_recv)
    nullify(numer_send)
    nullify(numer_recv)
    nullify(lesou)
    nullify(lbsou)
    nullify(lnsou)
    nullify(intersection)
    nullify(PAR_WORLD_RANKS)
    comin = 0.0_rp
    comax = 0.0_rp
    cz    = PAR_CURRENT_SIZE
    !
    ! If the distance criterion is required. In the case of global numbering, this
    ! is uncessary, unless meshes are different!
    !
    if(     itype == NEAREST_BOUNDARY_NODE  .or. &
         &  itype == BOUNDARY_INTERPOLATION .or. &
         &  itype == ELEMENT_INTERPOLATION  .or. &
         &  itype == STRESS_PROJECTION      .or. &
         &  itype == PROJECTION             .or. &
         &  itype == NEAREST_ELEMENT_NODE   .or. &
         &  itype == BOUNDARY_VECTOR_PROJECTION ) then
       require_distance = .true.
    else
       require_distance = .false.
    end if

    !---------------------------------------------------------------------------------------------
    !
    ! Source entities
    !
    !---------------------------------------------------------------------------------------------

    call couplings_source_entities(&
       coupling,itype,lesou,lbsou,lnsou,&
       CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,&
       CANDIDATE_SOURCE_BOUNDARIES)

    !---------------------------------------------------------------------------------------------
    !
    ! Allocate memory
    !
    !---------------------------------------------------------------------------------------------

    call memory_alloca(memor_cou,'NPOIN_SEND'      ,vacal,npoin_send,      cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'NPOIN_RECV'      ,vacal,npoin_recv,      cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'DECISION_SEND'   ,vacal,decision_send,   cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'DECISION_RECV'   ,vacal,decision_recv,   cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'MY_PART_TO_POINT',vacal,my_part_to_point,cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'MY_POINT_TO_PART',vacal,my_point_to_part,number_wet_points, 'INITIALIZE')
    call memory_alloca(memor_cou,'SHAPF'           ,vacal,shapf,           cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'INTERSECTION'    ,vacal,intersection,    cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'PAR_WORLD_RANKS' ,vacal,PAR_WORLD_RANKS, cz,                'INITIALIZE',0_ip)

    if( itype == GLOBAL_NUMBERING ) then
       call memory_alloca(memor_cou,'NUMER_SEND',vacal,numer_send,      cz,                'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'NUMER_RECV',vacal,numer_recv,      cz,                'INITIALIZE',0_ip)
    else
       call memory_alloca(memor_cou,'COORD_SEND',vacal,coord_send,      cz,                'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'COORD_RECV',vacal,coord_recv,      cz,                'INITIALIZE',0_ip)
    end if

    if( require_distance ) then
       call memory_alloca(memor_cou,'DISTANCE_SEND',vacal,distance_send,cz,'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'DISTANCE_RECV',vacal,distance_recv,cz,'INITIALIZE',0_ip)
    end if
    !
    ! World ranks. E.g. PAR_WORLD_RANKS(PAR_CURRENT_RANK) = PAR_MY_WORLD_RANK
    !
    if( present(COMM) ) then
       PAR_WORLD_RANKS(PAR_CURRENT_RANK) = PAR_MY_WORLD_RANK
       call PAR_MAX(PAR_CURRENT_SIZE,PAR_WORLD_RANKS(0:PAR_CURRENT_SIZE-1),PAR_COMM_CURRENT4)
    else
       do ipart = 0,PAR_WORLD_SIZE-1
          cpart = PAR_COMM_COLOR_PERM(icolo,jcolo,ipart)
          if( cpart > 0 ) PAR_WORLD_RANKS(cpart) = ipart
       end do
    end if

    !---------------------------------------------------------------------------------------------
    !
    ! Bin structure for partitions
    !
    !---------------------------------------------------------------------------------------------

    call par_bin_structure_initialization(bin_structure)
    call par_bin_structure_partition_bounding_box(&
         comin,comax,delta,&
         RELATIVE_TOLERANCE=toler_relative_cou,&
         ABSOLUTE_TOLERANCE=toler_absolute_cou)
    if( number_wet_points > 0 ) then
       do idime = 1,ndime
          comin(idime) = min(comin(idime),minval(xx(idime,1:number_wet_points)))
          comax(idime) = max(comax(idime),maxval(xx(idime,1:number_wet_points)))
       end do
    end if
    call par_bin_structure(bin_structure,comin,comax,COMM=PAR_COMM_CURRENT,VERBOSE=.false.)
    
    !---------------------------------------------------------------------------------------------
    !
    ! Find list of partitions to which I must send my points
    ! PP:                            : 1 -> NUMBER_WET_POINTS
    ! NPOIN_SEND(CPART)              = number of wet points to send to CPART
    ! MY_POINT_TO_PART(PP) % L(:)    = CPART, list of partitions where PP is in
    ! MY_PART_TO_POINT(CPART) % L(:) = PP,    list of points located in CPART bounding box
    !
    !---------------------------------------------------------------------------------------------

    call coupling_toolbox_points_in_partitions(&
         number_wet_points,ndime,xx,jcolo,bin_structure,coupling % kfl_multi_source,&
         PAR_CURRENT_SIZE,PAR_CURRENT_RANK,PAR_WORLD_RANKS,&
         npoin_send,my_point_to_part,my_part_to_point)

    !---------------------------------------------------------------------------------------------
    !
    ! Intersection arrays
    !
    !---------------------------------------------------------------------------------------------

    call couplings_intersection(&
       itype,PAR_CURRENT_SIZE,PAR_CURRENT_RANK,bin_structure,num_intersections4,intersection)
    
    call cputim(time1)
    coupling % cputim(6) = coupling % cputim(6) + time1-time0

    !---------------------------------------------------------------------------------------------
    !
    ! Coupling coordinates
    !
    !---------------------------------------------------------------------------------------------

    call couplings_coordinates(&
         coupling,itype,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT,&
         require_distance,xx,npoin_recv,numer_recv,coord_recv,decision_recv,distance_recv,&
         npoin_send,numer_send,coord_send,decision_send,distance_send,&
         intersection,my_part_to_point)
    
   call cputim(time2)
    coupling % cputim(3) = coupling % cputim(3) + time2-time1

    !---------------------------------------------------------------------------------------------
    !
    ! Look first for host elements in myslef if I am a source
    !
    !---------------------------------------------------------------------------------------------

    call couplings_in_myself(&
         itype,jcolo,number_wet_points,num_intersections4,&
         PAR_CURRENT_SIZE,PAR_CURRENT_RANK,PAR_COMM_CURRENT,&
         npoin_recv,coord_recv,decision_recv,&
         npoin_send,decision_send,distance_send,&
         lesou,shapf,intersection,my_part_to_point)

    call cputim(time3)
    coupling % cputim(4) = coupling % cputim(4) + time3-time2

    !---------------------------------------------------------------------------------------------
    !
    ! Postprocess
    !
    !---------------------------------------------------------------------------------------------

    call couplings_postprocess(&
       icoup,PAR_CURRENT_SIZE,PAR_CURRENT_RANK,npoin_recv,decision_send)
 
    !---------------------------------------------------------------------------------------------
    !
    ! Search strategy
    !
    !---------------------------------------------------------------------------------------------

    call couplings_search(&
       coupling,itype,number_wet_points,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT,&
       npoin_recv,numer_recv,coord_recv,decision_recv,distance_recv,&
       npoin_send,decision_send,distance_send,&
       lesou,lbsou,lnsou,shapf,&
       intersection,my_part_to_point)

    !---------------------------------------------------------------------------------------------
    !
    ! Decide on the owner
    !
    !---------------------------------------------------------------------------------------------

    call couplings_owner(&
         coupling,itype,number_wet_points,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT,&
         require_distance,npoin_send,decision_send,decision_recv,distance_recv,intersection,&
         my_part_to_point,my_point_to_part)

    call cputim(time4)
    coupling % cputim(4) = coupling % cputim(4) + time4-time3

    !---------------------------------------------------------------------------------------------
    !
    ! Setup communications
    !
    !---------------------------------------------------------------------------------------------

    call couplings_communication(&
         coupling,itype,number_wet_points,PAR_CURRENT_SIZE,&
         npoin_recv,decision_recv,&
         npoin_send,decision_send,&
         shapf,my_part_to_point,my_point_to_part)
    
    !--------------------------------------------------------------------------------------------- 
    !
    ! Allocate coupling % values_frequ for FSI in case of frequency is active for exchange
    !
    !---------------------------------------------------------------------------------------------
    
    if( coupling % frequ_send > 1_ip .or. coupling % frequ_recv > 1_ip ) then
       call memory_alloca(memor_cou,'COUPLING % VALUES_FREQU',vacal,coupling % values_frequ,ndime,coupling % geome % npoin_source,2_ip)
       do ipoin = 1_ip, coupling % geome % npoin_source
          do idime = 1_ip, ndime
             coupling % values_frequ(idime, ipoin, :) = 0_rp
          end do
       end do
    end if

    call cputim(time5)
    coupling % cputim(5) = coupling % cputim(5) + time5-time4
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_cou,'NPOIN_SEND'        ,vacal,npoin_send)
    call memory_deallo(memor_cou,'NPOIN_RECV'        ,vacal,npoin_recv)
    call memory_deallo(memor_cou,'MASK_NPOIN'        ,vacal,mask_npoin)
    call memory_deallo(memor_cou,'DECISION_SEND'     ,vacal,decision_send)
    call memory_deallo(memor_cou,'DECISION_RECV'     ,vacal,decision_recv)
    call memory_deallo(memor_cou,'COORD_SEND'        ,vacal,coord_send)
    call memory_deallo(memor_cou,'COORD_RECV'        ,vacal,coord_recv)
    call memory_deallo(memor_cou,'MY_PART_TO_POINT'  ,vacal,my_part_to_point)
    call memory_deallo(memor_cou,'MY_POINT_TO_PART'  ,vacal,my_point_to_part)
    call memory_deallo(memor_cou,'SHAPF'             ,vacal,shapf)
    call memory_deallo(memor_cou,'DISTANCE_SEND'     ,vacal,distance_send)
    call memory_deallo(memor_cou,'DISTANCE_RECV'     ,vacal,distance_recv)
    call memory_deallo(memor_cou,'INTERSECTION'      ,vacal,intersection)
    call memory_deallo(memor_cou,'PAR_WORLD_RANKS'   ,vacal,PAR_WORLD_RANKS)
    call memory_deallo(memor_cou,'NUMER_SEND'        ,vacal,numer_send)
    call memory_deallo(memor_cou,'NUMER_RECV'        ,vacal,numer_recv)
    if( .not. present(CANDIDATE_SOURCE_ELEMENTS) ) then
       call memory_deallo(memor_cou,'LESOU     '     ,vacal,lesou)
    end if
    if( .not. present(CANDIDATE_SOURCE_NODES) ) then
       call memory_deallo(memor_cou,'LNSOU     '     ,vacal,lnsou)
    end if
    if( .not. present(CANDIDATE_SOURCE_BOUNDARIES) ) then
       call memory_deallo(memor_cou,'LBSOU     '     ,vacal,lbsou)
    end if

    call par_bin_structure_deallocate(bin_structure)
    !
    ! Recover communicator
    !
    PAR_COMM_CURRENT = PAR_COMM_SAVE

  end subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_OLD

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check who owns the same node with same global numbering
  !> @details For each partition whcih have received some points,
  !>          check if I have a node with same global numbering
  !
  !----------------------------------------------------------------------
  subroutine COU_WET_POINTS_GLOBAL_NUMBERING(npoin_test,numer_test,inout_test, ht)

    use mod_htable

    implicit none
    integer(ip),          intent(in)  :: npoin_test                    !< Number of points to test
    integer(ip),          intent(in)  :: numer_test(npoin_test)        !< Global numbering of test points
    integer(ip),          intent(out) :: inout_test(npoin_test)        !< In or outside en element
    type(hash_t),          intent(in)  :: ht

    integer(ip)                       :: kpoin,ipoin_local
    integer(ip)                       :: ipoin_global

    !
    ! obtain lids
    !
    inout_test(:) = 0_ip

    do kpoin = 1,npoin_test
       ipoin_global = abs(numer_test(kpoin))
       ipoin_local  = htalid( ht,ipoin_global)
       if( ipoin_global > 0 ) then
          inout_test(kpoin) = ipoin_local
       endif
    end do

    ! do kpoin = 1,npoin
    !    ipoin_global = lninv_loc(kpoin)
    !    ipoin_local  = htalid( ht,ipoin_global)
    !    if( ipoin_local > 0 ) then
    !       inout_test(ipoin_local) = kpoin
    !    endif
    ! end do

  end subroutine COU_WET_POINTS_GLOBAL_NUMBERING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_HOST_ELEMENTS(ipass,npoin_test,coord_test,dimin_test,shapf_test,inout_test,lesou)

    integer(ip),          intent(in)     :: npoin_test                    !< Number of points to test
    real(rp),             intent(in)     :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),             intent(out)    :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    real(rp),             intent(inout)  :: shapf_test(mnode,npoin_test)  !< Shape function in elements
    integer(ip),          intent(inout)  :: inout_test(npoin_test)        !< In or outside en element
    integer(ip), pointer, intent(in)     :: lesou(:)
    integer(ip),          intent(in)     :: ipass
    integer(ip)                          :: kpoin,ielem
    integer(ip)                          :: iboun,pnodb
    integer(ip)                          :: ielse15,ielse8,ielse14
    real(rp)                             :: coloc(3),deriv(64)
    real(rp)                             :: dista
    integer(ip), pointer                 :: my_lesou(:)

    nullify(my_lesou)
    !
    ! Tell elsest that not all the elements should be considered
    ! Consider element only if LESOU(IELEM) = 1
    !
    ielse14 = ielse(14)
    ielse15 = ielse(15)
    ielse8  = ielse( 8)
    !
    ! Second chance, use A_TODA_COSTA options, but looking at nearest
    ! boundary elements
    !
    if( ipass == 2 ) then
       ielse(8) = ELSEST_A_TODA_COSTA
       call memory_alloca(memor_cou,'MY_LESOU','COU_WET_POINTS_HOST_ELEMENTS',my_lesou,nelem,'INITIALIZE')
       do iboun = 1,nboun
          pnodb = lnnob(iboun)
          ielem = lelbo(iboun)
          my_lesou(ielem) = 1
       end do
       do ielem = 1,nelem
          my_lesou(ielem) = min(my_lesou(ielem),lesou(ielem))
       end do
       ielse(14) = 1
       ielse(15) = 1
    else
       my_lesou => lesou
    end if
    !
    ! Loop over test points
    !
    !---------------------------------------------------------------------------
    !$OMP PARALLEL DO                                                          &
    !$OMP SCHEDULE     ( DYNAMIC , 100 )                                       &
    !$OMP DEFAULT      ( NONE )                                                &
    !$OMP SHARED       ( dimin_test, ielse, relse, meshe, coord_test,          &
    !$OMP                shapf_test, my_lesou, inout_test, ndivi, npoin_test ) &
    !$OMP PRIVATE      ( kpoin, ielem, deriv, coloc, dista )
    !---------------------------------------------------------------------------

    do kpoin = 1,npoin_test

       if( inout_test(kpoin) <= 0 ) then

          ielem = 0
          dimin_test(kpoin) = huge(1.0_rp)

          call elsest_host_element(&
               ielse,relse,1_ip,meshe(ndivi),coord_test(:,kpoin),ielem,&
               shapf_test(:,kpoin),deriv,coloc,dista,my_lesou)

          if( ielem <= 0 ) then
             inout_test(kpoin) = 0
          else
             if( my_lesou(ielem) == 1 ) then
                inout_test(kpoin) = ielem
                dimin_test(kpoin) = dista
             else
                inout_test(kpoin) = 0
             end if
          end if

       else if( inout_test(kpoin) >= huge(1_ip) ) then

          inout_test(kpoin) = my_huge
          dimin_test(kpoin) = huge(1.0_rp)

       end if

    end do
    !$OMP END PARALLEL DO

    ielse( 8) = ielse8
    ielse(14) = ielse14
    ielse(15) = ielse15

    if( ipass == 2 ) then
       ielse(8) = ielse8
       call memory_deallo(memor_cou,'MY_LESOU','COU_WET_POINTS_HOST_ELEMENTS',my_lesou)
    end if

  end subroutine COU_WET_POINTS_HOST_ELEMENTS

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_NEAREST_BOUNDARY_NODE(npoin_test,coord_test,dimin_test,inout_test,lnsou)
    integer(ip), intent(in)          :: npoin_test                    !< Number of points to test
    real(rp),    intent(in)          :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),    intent(out)         :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip), intent(out)         :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), intent(in), pointer :: lnsou(:)
    integer(ip)                      :: ipoin,ibono,ipoin_min
    integer(ip)                      :: idime,kpoin
    real(rp)                         :: dista,dista_min

    do kpoin = 1,npoin_test

       dista_min = huge(1.0_rp)
       ipoin_min = 0

       do ibono = 1,nbono
          ipoin = lbono(ibono)
          if( lnsou(ipoin) ) then
             dista = 0.0_rp
             do idime = 1,ndime
                dista = dista + ( coord_test(idime,kpoin) - coord(idime,ipoin) ) ** 2
             end do
             if( dista < dista_min ) then
                dista_min = dista
                ipoin_min = ipoin
             end if
          end if
       end do

       inout_test(kpoin) = ipoin_min
       dimin_test(kpoin) = dista_min

    end do

  end subroutine COU_WET_POINTS_NEAREST_BOUNDARY_NODE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_NEAREST_ELEMENT_NODE(npoin_test,coord_test,dimin_test,inout_test,lnsou)
    integer(ip), intent(in)          :: npoin_test                    !< Number of points to test
    real(rp),    intent(in)          :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),    intent(out)         :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip), intent(out)         :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), intent(in), pointer :: lnsou(:)
    integer(ip)                      :: ipoin,ipoin_min,kpoin
    real(rp)                         :: dista,dista_min

    do kpoin = 1,npoin_test

       dista_min = huge(1.0_rp)
       ipoin_min = 0

       do ipoin = 1,npoin
          if( lnsou(ipoin) ) then
             dista = dot_product(coord_test(1:ndime,kpoin)-coord(1:ndime,ipoin),coord_test(1:ndime,kpoin)-coord(1:ndime,ipoin))
             if( dista < dista_min ) then
                dista_min = dista
                ipoin_min = ipoin
             end if
          end if
       end do

       inout_test(kpoin) = ipoin_min
       dimin_test(kpoin) = dista_min

    end do

  end subroutine COU_WET_POINTS_NEAREST_ELEMENT_NODE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          Look for a node with the same coordinates
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_SAME_COORDINATE(npoin_test,coord_test,inout_test,lnsou)
    integer(ip),          intent(in)   :: npoin_test                    !< Number of points to test
    real(rp),             intent(in)   :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    integer(ip),          intent(out)  :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), pointer, intent(in)   :: lnsou(:)
    integer(ip)                        :: ipoin,kpoin
    real(rp)                           :: toler

    toler = zeror

    if( ndime == 2 ) then

       if( associated(lnsou) ) then
          do kpoin = 1,npoin_test
             inout_test(kpoin) = 0
             loop_ipoin_2d_1: do ipoin = 1,npoin
                if( lnsou(ipoin) ) then
                   if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < toler ) then
                      if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < toler ) then
                         inout_test(kpoin) = ipoin
                         exit loop_ipoin_2d_1
                      end if
                   end if
                end if
             end do loop_ipoin_2d_1
          end do
       else
          do kpoin = 1,npoin_test
             inout_test(kpoin) = 0
             loop_ipoin_2d_2: do ipoin = 1,npoin
                if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < toler ) then
                   if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < toler ) then
                      inout_test(kpoin) = ipoin
                      exit loop_ipoin_2d_2
                   end if
                end if
             end do loop_ipoin_2d_2
          end do
       end if

    else

       if( associated(lnsou) ) then
          do kpoin = 1,npoin_test
             inout_test(kpoin) = 0
             loop_ipoin_3d_1: do ipoin = 1,npoin
                if( lnsou(ipoin) ) then
                   if( abs(coord_test(1,kpoin)-coord(1,ipoin)) <= toler ) then
                      if( abs(coord_test(2,kpoin)-coord(2,ipoin)) <= toler ) then
                         if( abs(coord_test(3,kpoin)-coord(3,ipoin)) <= toler ) then
                            inout_test(kpoin) = ipoin
                            exit loop_ipoin_3d_1
                         end if
                      end if
                   end if
                end if
             end do loop_ipoin_3d_1
          end do
       else
          do kpoin = 1,npoin_test
             inout_test(kpoin) = 0
             loop_ipoin_3d_2: do ipoin = 1,npoin
                if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < toler ) then
                   if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < toler ) then
                      if( abs(coord_test(3,kpoin)-coord(3,ipoin)) < toler ) then
                         inout_test(kpoin) = ipoin
                         exit loop_ipoin_3d_2
                      end if
                   end if
                end if
             end do loop_ipoin_3d_2
          end do
       end if

    end if

  end subroutine COU_WET_POINTS_SAME_COORDINATE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_HOST_BOUNDARIES(&
       npoin_test,coord_test,dimin_test,inout_test,shapf_test,kdtree_source,lbsou)
    integer(ip),                        intent(in)    :: npoin_test                    !< Number of points to test
    real(rp),                           intent(in)    :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),                           intent(out)   :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip),                        intent(out)   :: inout_test(npoin_test)        !< Node with minimum distance
    real(rp),                           intent(out)   :: shapf_test(mnodb,npoin_test)  !< Shape function in elements
    type(typ_kdtree),                   intent(inout) :: kdtree_source                 !< KD-tree
    logical(lg),      pointer, optional,intent(in)    :: lbsou(:)                      !< List of source boundaries
    integer(ip)                                       :: ipoin,iboun_min
    integer(ip)                                       :: kpoin,inodb
!    integer(ip)                                       :: idime
    integer(ip)                                       :: pnodb,pblty,ifoun,ptopo
    real(rp)                                          :: dista_min,proje(3)
    real(rp)                                          :: deriv(3*mnode)
    real(rp)                                          :: coloc(3),toler
    real(rp)                                          :: bocod(ndime,mnodb)
    logical(lg)                                       :: if_mask
    logical(lg)                                       :: inbox
!    real(rp)                                          :: tolerb

    if_mask = .false.
    if( present(lbsou) ) then
       if( associated(lbsou) ) if_mask = .true.
    end if

    if( kdtree_source % nboun /= 0 ) then

       toler = 0.01_rp

       do kpoin = 1,npoin_test
          !
          ! Check that at least is within the bounding box of boundaries
          !
          inbox=.true.
          !do idime= 1,ndime
          !   tolerb = (kdtree_source % bobox(idime,2)-kdtree_source % bobox(idime,1))*0.01_rp
          !   if(coord_test(idime,kpoin) < kdtree_source % bobox(idime,1)-tolerb) inbox = .false.
          !   if(coord_test(idime,kpoin) > kdtree_source % bobox(idime,2)+tolerb) inbox = .false.
          !enddo
          if( inbox ) then
             dista_min = huge(1.0_rp)
             !
             ! KDTREE returns the boundary number in global numbering because KDTree was constructed 
             ! using a permutation array
             !
             if( if_mask ) then
                call kdtree_nearest_boundary(coord_test(1,kpoin),kdtree_source,iboun_min,dista_min,proje,LMASK=lbsou)
             else
                call kdtree_nearest_boundary(coord_test(1,kpoin),kdtree_source,iboun_min,dista_min,proje)
             end if
             ! if( dista_min < 0 ) print*, "DEBUG: dista min= ",dista_min,coord_test(:,kpoin),kpoin
             dista_min = abs(dista_min)          
             ifoun     = 0
             pblty     = abs(ltypb_cou(iboun_min))
             pnodb     = lnnob_cou(iboun_min)   
             ptopo     = ltopo(pblty)
             do inodb = 1,pnodb
                ipoin = lnodb_cou(inodb,iboun_min)
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
             end do
             !print*,ndime,pblty,ptopo,pnodb

             call elmgeo_natural_coordinates_on_boundaries(&
                  ndime,pblty,pnodb,bocod,         &
                  shapf_test(1,kpoin),deriv,proje, & 
                  coloc,ifoun,toler)             
             !
             ! To avoid problems temporarly
             !

             ifoun = 1_ip

             !
             !if(ifoun/=0) write(*,'(i2,7(1x,e12.6))') kfl_paral,coord_test(1:2,kpoin),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou2,iboun_min))
             !if(ifoun/=0) write(*,'(7(1x,e12.6))') kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),coord(1:ndime,lnodb_cou1,iboun_min))

             if( ifoun == 0 ) then
                !write(6,'(a,3(1x,i5))')             'A=',kfl_paral,kpoin,iboun_min,nboun
                !write(6,'(a,2(1x,i4),3(1x,e12.6))') 'B=',kfl_paral,kpoin,coloc

                !write(6,'(a,i2,7(1x,i5))')          'C=',kfl_paral,kpoin,lnodb_cou1,iboun_min),lnodb_cou2,iboun_min),lnodb_cou3,iboun_min),lnodb_cou4,iboun_min)
                !write(6,'(a,i2,i5,7(1x,e12.6))')    'D=',kfl_paral,kpoin,coord_test(1:3,kpoin)
                !write(6,*) '1: ',coord(1:3,lnodb_cou1,iboun_min))
                !write(6,*) '2: ',coord(1:3,lnodb_cou2,iboun_min))
                !write(6,*) '3: ',coord(1:3,lnodb_cou3,iboun_min))
                !write(6,*) '4: ',coord(1:3,lnodb_cou4,iboun_min))
                !!coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou2,iboun_min))
                inout_test(kpoin) = 0_ip
                dimin_test(kpoin) = dista_min
             else
                inout_test(kpoin) = iboun_min
                dimin_test(kpoin) = dista_min
             end if
          else
             inout_test(kpoin) = 0_ip
             dimin_test(kpoin) = huge(1.0_rp)    
          endif
       end do

    end if

  end subroutine COU_WET_POINTS_HOST_BOUNDARIES

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Search
  !> @details Search elements, boundaries, nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_search(&
       coupling,itype,number_wet_points,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT4,&
       npoin_recv,numer_recv,coord_recv,decision_recv,distance_recv,&
       npoin_send,decision_send,distance_send,&
       lesou,lbsou,lnsou,shapf,&
       intersection,my_part_to_point)

    type(typ_color_coupling),  intent(inout)         :: coupling
    integer(ip),               intent(in)            :: itype
    integer(ip),               intent(in)            :: number_wet_points
    integer(4),                intent(in)            :: num_intersections4 
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    MY_MPI_COMM   ,            intent(in)            :: PAR_COMM_CURRENT4
    type(i1p),   pointer,      intent(inout)         :: numer_recv(:)
    integer(ip), pointer,      intent(inout)         :: npoin_recv(:)
    type(r1p),   pointer,      intent(inout)         :: coord_recv(:)
    type(i1p),   pointer,      intent(inout)         :: decision_recv(:)
    type(r1p),   pointer,      intent(inout)         :: distance_recv(:)      ! Distance check    
    integer(ip), pointer,      intent(inout)         :: npoin_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_send(:)
    type(r1p),   pointer,      intent(inout)         :: distance_send(:)      ! Distance check    
    integer(ip), pointer,      intent(in)            :: lesou(:)
    logical(lg), pointer,      intent(in)            :: lbsou(:)
    logical(lg), pointer,      intent(in)            :: lnsou(:)
    type(r2p),   pointer,      intent(inout)         :: shapf(:)
    logical(lg), pointer,      intent(in)            :: intersection(:)
    type(i1p),   pointer,      intent(in)            :: my_part_to_point(:)

    integer(ip), pointer                             :: check_points(:)
    integer(ip)                                      :: all_points_found
    integer(ip)                                      :: ipass,cpart,pp,kk
    type(hash_t)                                     :: ht
                               
    nullify(check_points)
    all_points_found   = 0
    ipass              = 0
    
    if( coupling % itype == GLOBAL_NUMBERING ) then
       !
       ! Construct hash table
       !
       call htaini( ht, npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
       if(associated(lninv_loc))then
          call htaadd( ht, npoin, lninv_loc)
       else
          if(npoin /= 0_ip) then
             call runend("Cou_init_interpolate_points_values: lninv_loc not associated")
          endif
       endif
       if( ht % nelem /= npoin) call runend("Error: repited elements in lninv_loc")
    endif

    do while( all_points_found == 0 )

       ipass = ipass + 1

       if( IMASTER    ) then

       else

          if( itype == GLOBAL_NUMBERING .or. itype == SAME_COORDINATE ) then
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
          else
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
             call PAR_START_NON_BLOCKING_COMM(2_ip,num_intersections4)
          end if
          !
          ! Check if I own the points
          !       
          if( itype == GLOBAL_NUMBERING ) then
             !
             ! Global numbering
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_GLOBAL_NUMBERING(pp,numer_recv(cpart) % l, decision_send(cpart) % l, ht)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == ELEMENT_INTERPOLATION ) then
             !
             ! Element interpolation
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   if( .not. associated(shapf(cpart)%a) ) &
                        call memory_alloca(memor_cou,'SHAPF % A',vacal,shapf(cpart)%a,mnode,pp)
                   call COU_WET_POINTS_HOST_ELEMENTS(ipass,pp,coord_recv(cpart) % a, distance_send(cpart) % a, shapf(cpart) % a,decision_send(cpart) % l,lesou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == NEAREST_BOUNDARY_NODE ) then
             !
             ! Nearest boundary node
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_NEAREST_BOUNDARY_NODE(pp,coord_recv(cpart) % a,distance_send(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == NEAREST_ELEMENT_NODE ) then
             !
             ! Nearest element node
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_NEAREST_ELEMENT_NODE(pp,coord_recv(cpart) % a,distance_send(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == SAME_COORDINATE ) then
             !
             ! Same coordinate
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_SAME_COORDINATE(pp,coord_recv(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == BOUNDARY_INTERPOLATION .or. &
               &   itype == STRESS_PROJECTION      .or. &
               &   itype == PROJECTION             ) then
             !
             ! Boundary interpolation and Gauss point interpolation
             ! Wet points are nodes in the first case and Gauss points in the second case
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call memory_alloca(memor_cou,'SHAPF % A',vacal,shapf(cpart)%a,mnodb,pp)
                   call COU_WET_POINTS_HOST_BOUNDARIES(pp,coord_recv(cpart) % a,distance_send(cpart) % a,&
                        &                              decision_send(cpart) % l,shapf(cpart) % a,        &
                        &                              coupling % geome % kdtree,lbsou                   )
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do

          else if( itype == BOUNDARY_VECTOR_PROJECTION ) then
             !
             ! Look for the boudnary corssing the projection along a vector
             !
             call runend('MOD_COUPLINGS_SETUP: THIS COUPLING DOES NOT EXIST')

          end if
          if( itype == GLOBAL_NUMBERING .or. itype == SAME_COORDINATE ) then
             call PAR_END_NON_BLOCKING_COMM(1_ip)
          else
             call PAR_END_NON_BLOCKING_COMM(1_ip)
             call PAR_END_NON_BLOCKING_COMM(2_ip)
          end if
       end if
       !
       ! Check if we need a second round
       !
       if( itype == ELEMENT_INTERPOLATION .and. coupling % kfl_toda_costa == 1 ) then

          call memory_alloca(memor_cou,'CHECK_POINTS',vacal,check_points,number_wet_points)
          !
          ! Check if someone has found the point
          !
          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart)
                if( abs(decision_recv(cpart) % l(kk)) > 0 .and. decision_recv(cpart) % l(kk) < my_huge ) then
                   pp = my_part_to_point(cpart) % l(kk)
                   check_points(pp) = 1
                end if
             end do
          end do
          if( number_wet_points > 0 ) then
             all_points_found = minval(check_points)
          else
             all_points_found = 1
          end if
          call PAR_MIN(all_points_found,'IN CURRENT COUPLING')
          !
          ! Tell neighbors which points they should look at
          !
          !if(IMASTER) print*,'popo=',all_points_found
          if( all_points_found == 0 ) then
             !if( ipass == 2 ) call runend('WE ARE IN TROUBLE IN COUPLING')
             if( ipass == 2 ) all_points_found = 1
             call messages_live('WE GO FOR A SECOND ROUND TO FIND HOST ELEMENTS FOR SOME LOST WET NODES','WARNING')
             do cpart = 0,PAR_CURRENT_SIZE-1
                do kk = 1,npoin_send(cpart)
                   pp = my_part_to_point(cpart) % l(kk)
                   if( check_points(pp) == 0 ) then
                      decision_recv(cpart) % l(kk) = -1
                   end if
                end do
             end do
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
             call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( intersection(cpart) ) then
                   call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                end if
             end do
             call PAR_END_NON_BLOCKING_COMM(1_ip)
          end if

          call memory_deallo(memor_cou,'CHECK_POINTS',vacal,check_points)
       else
          all_points_found = 1
       end if

    end do
    if( coupling % itype == GLOBAL_NUMBERING ) then
       !
       ! Dellocate hast table
       !
       call htades( ht )

    endif

    !
    ! Remove -1
    !
    do cpart = 0,PAR_CURRENT_SIZE-1
       do kk = 1,npoin_send(cpart)
          decision_recv(cpart) % l(kk) = max(decision_recv(cpart) % l(kk),0_ip)
          if( decision_recv(cpart) % l(kk) >= my_huge ) decision_recv(cpart) % l(kk) = 0
       end do
       do kk = 1,npoin_recv(cpart)
          decision_send(cpart) % l(kk) = max(decision_send(cpart) % l(kk),0_ip)
          if( decision_send(cpart) % l(kk) >= my_huge ) decision_send(cpart) % l(kk) = 0
       end do
    end do
    
  end subroutine couplings_search

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Search
  !> @details Search elements, boundaries, nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_owner(&
       coupling,itype,number_wet_points,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT4,&
       require_distance,npoin_send,decision_send,decision_recv,distance_recv,intersection,&
       my_part_to_point,my_point_to_part)

    type(typ_color_coupling),  intent(inout)         :: coupling
    integer(ip),               intent(in)            :: itype
    integer(ip),               intent(in)            :: number_wet_points
    integer(4),                intent(in)            :: num_intersections4
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    MY_MPI_COMM   ,            intent(in)            :: PAR_COMM_CURRENT4
    logical(lg),               intent(in)            :: require_distance
    integer(ip), pointer,      intent(in)            :: npoin_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_recv(:)
    type(r1p),   pointer,      intent(in)            :: distance_recv(:)      ! Distance check        
    logical(lg), pointer,      intent(in)            :: intersection(:)
    type(i1p),   pointer,      intent(in)            :: my_part_to_point(:)
    type(i1p),   pointer,      intent(inout)         :: my_point_to_part(:)

    integer(ip), pointer                             :: distance_min_cpart(:) ! Distance check
    real(rp),    pointer                             :: distance_min_value(:) ! Distance check
    integer(ip)                                      :: cpart,kk,pp

    nullify(distance_min_cpart)
    nullify(distance_min_value)

    if( INOTMASTER ) then
       !
       ! Decide the owner
       !
       if( ( itype == GLOBAL_NUMBERING .or. itype == SAME_COORDINATE ) .and. number_wet_points > 0 .and. coupling % kfl_multi_source == 0 ) then
          !
          ! For element interpolation take element which belongs to subdomain
          ! with highest rank
          !
          call memory_alloca(memor_cou,'DISTANCE_MIN_CPART',vacal,distance_min_cpart,number_wet_points)

          distance_min_cpart = 0_ip

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart)
                if( decision_recv(cpart) % l(kk) /= 0 ) then
                   pp = my_part_to_point(cpart) % l(kk)
                   if( cpart > distance_min_cpart(pp) ) then
                      distance_min_cpart(pp) = cpart
                   end if
                end if
             end do
          end do

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart)
                pp = my_part_to_point(cpart) % l(kk)
                if( distance_min_cpart(pp) /= cpart ) decision_recv(cpart) % l(kk) = 0
             end do
          end do

          do pp = 1,number_wet_points
             if( associated(my_point_to_part(pp) % l) ) my_point_to_part(pp) % l(1) = distance_min_cpart(pp)
          end do

          call memory_deallo(memor_cou,'DISTANCE_MIN_CPART',vacal,distance_min_cpart)

       else if( require_distance .and. number_wet_points > 0 ) then
          !
          ! Look for minimum distance
          ! 1. Compute minimum distance for each wet point pp, from 1 to number_wet_points
          ! 2. Put decision_recv(cpart) % l(kk) = 0 when cpart is different from the owner
          ! 3. my_point_to_part(pp) % l(1) = cpart_owner
          !
          call memory_alloca(memor_cou,'DISTANCE_MIN_VALUE',vacal,distance_min_value,number_wet_points)
          call memory_alloca(memor_cou,'DISTANCE_MIN_CPART',vacal,distance_min_cpart,number_wet_points)

          distance_min_value = huge(1.0_rp)

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart)
                if( decision_recv(cpart) % l(kk) /= 0 ) then
                   pp = my_part_to_point(cpart) % l(kk)
                   if( distance_recv(cpart) % a(kk) < distance_min_value(pp) ) then
                      distance_min_value(pp)      = distance_recv(cpart) % a(kk)
                      distance_min_cpart(pp)      = cpart
                   end if
                end if
             end do
          end do

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart)
                pp = my_part_to_point(cpart) % l(kk)
                if( distance_min_cpart(pp) /= cpart ) decision_recv(cpart) % l(kk) = 0
             end do
          end do

          do pp = 1,number_wet_points
             if( associated(my_point_to_part(pp) % l) ) my_point_to_part(pp) % l(1) = distance_min_cpart(pp)
          end do

          call memory_deallo(memor_cou,'DISTANCE_MIN_CPART',vacal,distance_min_cpart)
          call memory_deallo(memor_cou,'DISTANCE_MIN_VALUE',vacal,distance_min_value)

       end if
       !
       ! Send my result back
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
          end if
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)

    end if

  end subroutine couplings_owner

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Look in myself
  !> @details Look for element in myself first
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_in_myself(&
       itype,jcolo,number_wet_points,num_intersections4,&
       PAR_CURRENT_SIZE,PAR_CURRENT_RANK,PAR_COMM_CURRENT4,&
       npoin_recv,coord_recv,decision_recv,&
       npoin_send,decision_send,distance_send,&
       lesou,shapf,intersection,my_part_to_point)

    integer(ip),               intent(in)            :: itype
    integer(ip),               intent(in)            :: jcolo
    integer(ip),               intent(in)            :: number_wet_points
    integer(4),                intent(in)            :: num_intersections4 
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    integer(ip),               intent(in)            :: PAR_CURRENT_RANK
    MY_MPI_COMM   ,            intent(in)            :: PAR_COMM_CURRENT4
    integer(ip), pointer,      intent(in)            :: npoin_recv(:)
    type(r1p),   pointer,      intent(in)            :: coord_recv(:)
    type(i1p),   pointer,      intent(inout)         :: decision_recv(:)
    integer(ip), pointer,      intent(in)            :: npoin_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_send(:)
    type(r1p),   pointer,      intent(inout)         :: distance_send(:)      ! Distance check    
    integer(ip), pointer,      intent(in)            :: lesou(:)
    type(r2p),   pointer,      intent(inout)         :: shapf(:)
    logical(lg), pointer,      intent(in)            :: intersection(:)
    type(i1p),   pointer,      intent(in)            :: my_part_to_point(:)

    logical(lg), pointer                             :: first_myself(:)
    integer(ip)                                      :: ipass,cpart,pp,ii
    integer(ip)                                      :: check_myself_first

    nullify(first_myself)
    
    if( itype == ELEMENT_INTERPOLATION ) then

       check_myself_first = 0
       ipass              = 1

       if( par_part_in_color(PAR_MY_WORLD_RANK,jcolo) .and. I_AM_IN_COLOR(color_source) ) then
          !
          ! Check if I host the wet points
          ! CHECK_MYSLEF_FIRST= number of wet points I have found
          !
          cpart = PAR_CURRENT_RANK
          if( npoin_recv(cpart) /= 0 ) then
             pp = npoin_recv(cpart)
             if( .not. associated(shapf(cpart)%a) ) &
                  call memory_alloca(memor_cou,'SHAPF % A',vacal,shapf(cpart)%a,mnode,pp)
             call COU_WET_POINTS_HOST_ELEMENTS(ipass,pp,coord_recv(cpart) % a, distance_send(cpart) % a, shapf(cpart) % a,decision_send(cpart) % l,lesou)
             check_myself_first = count( decision_send(cpart) % l /= 0 , KIND=ip)
          end if
       end if

       call PAR_MAX(check_myself_first,PAR_COMM_CURRENT4)

       if( check_myself_first > 0 ) then
          !
          ! Mark my wet points and mark lost wet points not to be checked anymore by myself
          ! FIRST_MYSELF(PP) = .TRUE. if I found host element for PP
          !
          call memory_alloca(memor_cou,'FIRST_MYSELF',vacal,first_myself,number_wet_points)
          do ii = 1,npoin_recv(PAR_CURRENT_RANK)
             pp = my_part_to_point(PAR_CURRENT_RANK) % l(ii)
             if( decision_send(PAR_CURRENT_RANK) % l(ii) > 0 ) then
                first_myself(pp) = .true.                            ! I found PP!
             else
                decision_send(PAR_CURRENT_RANK) % l(ii) = my_huge    ! I did not find it, I will not check it anymore
             end if
          end do
          !
          ! Tell people they do not have to look at the wet points I have already found
          !
          do cpart = 0,PAR_CURRENT_SIZE-1
             if( cpart /= PAR_CURRENT_RANK ) then
                do ii = 1,npoin_send(cpart)
                   pp = my_part_to_point(cpart) % l(ii)
                   if( first_myself(pp) ) decision_recv(cpart) % l(ii) = my_huge
                end do
             end if
          end do
          call memory_deallo(memor_cou,'FIRST_MYSELF',vacal,first_myself)
          !
          ! Send my first decision to my neighbor (put result in decision_send
          !
          call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
          call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
          do cpart = 0,PAR_CURRENT_SIZE-1
             if( intersection(cpart) .and. cpart /= PAR_CURRENT_RANK ) then
                call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
             end if
          end do
          call PAR_END_NON_BLOCKING_COMM(1_ip)

       end if

    end if

  end subroutine couplings_in_myself

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Postprocess
  !> @details Postprocess some statistics
  !>          Statistics, to be draw with gnuplot.
  !>          **  is problem name
  !>          *** is the coupling number
  !>          
  !>          set xrange[0.5:]
  !>          set yrange[0.5:]
  !>
  !>          set title  'Number of points rank x receives from rank y'
  !>          set xlabel 'MPI rank x'
  !>          set ylabel 'MPI rank y'
  !>         
  !>          plot '**-histogram-recv-***.cou.res' matrix with image
  !>          set title  'Number of points rank x receives from rank y and checked'
  !>          set xlabel 'MPI rank x'
  !>          set ylabel 'MPI rank y'
  !>          plot '**-histogram-chec-***.cou.res' matrix with image
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_postprocess(&
       icoup,PAR_CURRENT_SIZE,PAR_CURRENT_RANK,npoin_recv,decision_send)
    
    integer(ip),               intent(in)         :: icoup
    integer(ip),               intent(in)         :: PAR_CURRENT_SIZE
    integer(ip),               intent(in)         :: PAR_CURRENT_RANK
    integer(ip), pointer,      intent(in)         :: npoin_recv(:)
    type(i1p),   pointer,      intent(in)         :: decision_send(:)
    
    integer(ip)                                   :: lun_chec,cpart,kpart
    integer(ip)                                   :: lun_recv,ipart
    integer(ip), pointer                          :: num_chec(:)
    integer(ip), pointer                          :: num_chec_gat(:,:)
    integer(ip), pointer                          :: num_recv_gat(:,:)
   
    if( 1 == 2 ) then
       nullify(num_chec)
       nullify(num_chec_gat)
       nullify(num_recv_gat)
       kpart = PAR_CURRENT_SIZE-1
       if( PAR_CURRENT_RANK == 0 ) then
          allocate(num_chec_gat(0:kpart,0:kpart))
          allocate(num_recv_gat(0:kpart,0:kpart))
       end if
       allocate(num_chec(0:kpart))
       num_chec = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( npoin_recv(cpart) /= 0 ) then
             num_chec(cpart) = count(decision_send(cpart) % l == 0, KIND=ip)
          end if
       end do
       call PAR_GATHER(num_chec  ,num_chec_gat,'IN CURRENT COUPLING')
       call PAR_GATHER(npoin_recv,num_recv_gat,'IN CURRENT COUPLING')
       if( PAR_CURRENT_RANK == 0 ) then
          lun_chec = iofile_available_unit()
          call iofile_open_unit(lun_chec,adjustl(trim(namda))//'-histogram-chec-'//trim(integer_to_string(icoup))//'.cou.res')
          lun_recv = iofile_available_unit()
          call iofile_open_unit(lun_recv,adjustl(trim(namda))//'-histogram-recv-'//trim(integer_to_string(icoup))//'.cou.res')
          do cpart = 0,kpart
             write(lun_chec,'(3000(1x,i6))') (num_chec_gat(cpart,ipart),ipart=0,kpart)
             write(lun_recv,'(3000(1x,i6))') (num_recv_gat(cpart,ipart),ipart=0,kpart)
             call iofile_flush_unit(lun_chec)
             call iofile_flush_unit(lun_recv)
          end do
          deallocate(num_chec_gat)
          deallocate(num_recv_gat)
          call iofile_close_unit(lun_chec)
          call iofile_close_unit(lun_recv)
       end if
       deallocate(num_chec)
    end if

  end subroutine couplings_postprocess

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Send receive coordinates
  !> @details Send and receive coordinates
  !>
  !>          1. How many wet points I send to (I) and how many I receive from (J)
  !>         
  !>                                     +-------+
  !>                                 =>  |       |  =>
  !>                                     |       |
  !>                   NPOIN_RECV(I) =>  |       |  => NPOIN_SEND(I)
  !>                                     +-------+
  !>         
  !>           2. Send coordinates and receive coordinates
  !>         
  !>                                     +-------+
  !>                                 =>  |       |  =>
  !>                                     |       |
  !>               COORD_RECV(I) % A =>  |       |  => COORD_SEND(I) % A
  !>                                     +-------+
  !>         
  !>           3. Check if I have what I receive (find host element, nearest point, etc.) and send it back: DECISION_SEND(I) % L
  !>              Receive the decision of the others: DECISION_RECV(I) % L
  !>              In case of nearest point, send the minimal distance I have found as well: DISTANCE_SEND(I) % A
  !>              Receive the minimal distance of others as well: DISTANCE_RECV(I) % A
  !>         
  !>                                     +-------+
  !>                                 <=  |       |  <=
  !>                                     |       |
  !>            DECISION_SEND(I) % L <=  |       |  <= DECISION_RECV(I) % L
  !>                                     +-------+
  !>         
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_coordinates(&
       coupling,itype,num_intersections4,PAR_CURRENT_SIZE,PAR_COMM_CURRENT4,&
       require_distance,xx,npoin_recv,numer_recv,coord_recv,decision_recv,distance_recv,&
       npoin_send,numer_send,coord_send,decision_send,distance_send,&
       intersection,my_part_to_point)

    type(typ_color_coupling),  intent(inout)         :: coupling
    integer(ip),               intent(in)            :: itype
    integer(4),                intent(in)            :: num_intersections4 
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    MY_MPI_COMM   ,            intent(in)            :: PAR_COMM_CURRENT4
    logical(lg),               intent(in)            :: require_distance
    real(rp),    pointer,      intent(in)            :: xx(:,:)
    type(i1p),   pointer,      intent(inout)         :: numer_recv(:)
    integer(ip), pointer,      intent(inout)         :: npoin_recv(:)
    type(r1p),   pointer,      intent(inout)         :: coord_recv(:)
    type(i1p),   pointer,      intent(inout)         :: decision_recv(:)
    type(r1p),   pointer,      intent(inout)         :: distance_recv(:)      ! Distance check    
    integer(ip), pointer,      intent(inout)         :: npoin_send(:)
    type(i1p),   pointer,      intent(inout)         :: numer_send(:)
    type(r1p),   pointer,      intent(inout)         :: coord_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_send(:)
    type(r1p),   pointer,      intent(inout)         :: distance_send(:)  
    logical(lg), pointer,      intent(inout)         :: intersection(:)
    type(i1p),   pointer,      intent(inout)         :: my_part_to_point(:)

    integer(ip)                                      :: ipoin_min_max_send(2)
    integer(ip)                                      :: ipoin_min_max_recv(2)
    integer(ip)                                      :: cpart,kk,pp,ii
    integer(ip)                                      :: idime,ipoin
    integer(ip)                                      :: ipoin_min,ipoin_max
    integer(ip)                                      :: npoin_to_send
    integer(ip)                                      :: npoin_to_recv
    !
    ! Get how many points I should check if I have them
    !
    call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv,PAR_COMM_IN=PAR_COMM_CURRENT4)
    !
    ! In the case of global numbering, reduce the possible number of communications
    ! In fact, the send and receive should be symmetric
    !
    if( coupling % itype == GLOBAL_NUMBERING .and. coupling % kfl_symmetry == 1 ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             if( npoin_send(cpart) == 0 .or. npoin_recv(cpart) == 0 ) then
                npoin_send(cpart)   = 0
                npoin_recv(cpart)   = 0
                intersection(cpart) = .false.
                call memory_deallo(memor_cou,'MY_PART_TO_POINT % L', vacal,my_part_to_point(cpart) % l)
             elseif(1==2) then
                if( npoin_send(cpart) > 0 ) then
                   ipoin_min_max_send(1) = minval(lninv_loc(numer_send(cpart) % l))
                   ipoin_min_max_send(2) = maxval(lninv_loc(numer_send(cpart) % l))
                else
                   ipoin_min_max_send    = 0_rp
                end if
                call PAR_SEND_RECEIVE(2_ip,2_ip,ipoin_min_max_send,ipoin_min_max_recv,'IN CURRENT COUPLING',cpart,'BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
                ipoin_min = max(ipoin_min_max_send(1),ipoin_min_max_recv(1))
                ipoin_max = min(ipoin_min_max_send(2),ipoin_min_max_recv(2))
                if( npoin_send(cpart) > 0_ip ) then
                   kk = 0
                   do ii = 1,npoin_to_send
                      pp    = my_part_to_point(cpart) % l(ii)
                      ipoin = coupling % wet % lpoin_wet(pp)
                      if( lninv_loc(ipoin) >= ipoin_min ) then
                         kk = kk + 1
                         my_part_to_point(cpart) % l(kk) = pp
                      end if
                   end do
                   if( kk /= ii ) then
                      npoin_send(cpart) = kk
                      call memory_resize(memor_cou,'MY_PART_TO_POINT % L', vacal,my_part_to_point(cpart) % l,npoin_send(cpart))
                   end if
                end if
             end if
          end if
       end do
    end if
    !
    ! Save point coordinates
    !
    do cpart = 0,PAR_CURRENT_SIZE-1
       npoin_to_send = npoin_send(cpart)
       npoin_to_recv = npoin_recv(cpart)
       if( itype == GLOBAL_NUMBERING ) then
          call memory_alloca(memor_cou,'NUMER_SEND % L', vacal,numer_send(cpart) % l,npoin_to_send)
          call memory_alloca(memor_cou,'NUMER_RECV % L', vacal,numer_recv(cpart) % l,npoin_to_recv)
       else
          call memory_alloca(memor_cou,'COORD_SEND % A', vacal,coord_send(cpart) % a,npoin_to_send*ndime)
          call memory_alloca(memor_cou,'COORD_RECV % A', vacal,coord_recv(cpart) % a,npoin_to_recv*ndime)
       end if
       call memory_alloca(memor_cou,'DECISION_RECV % L' ,vacal,decision_recv(cpart)  % l,npoin_to_send)
       call memory_alloca(memor_cou,'DECISION_SEND % L' ,vacal,decision_send(cpart)  % l,npoin_to_recv)
       if( require_distance ) then
          call memory_alloca(memor_cou,'DISTANCE_SEND % A',vacal,distance_send(cpart) % a,npoin_to_recv)
          call memory_alloca(memor_cou,'DISTANCE_RECV % A',vacal,distance_recv(cpart) % a,npoin_to_send)
       end if
    end do
    !
    ! Copy points to send buffer
    !
    if( itype == GLOBAL_NUMBERING ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          npoin_to_send = npoin_send(cpart)
          if( npoin_to_send > 0 ) then
             do ii = 1,npoin_to_send
                pp    = my_part_to_point(cpart) % l(ii)
                ipoin = coupling % wet % lpoin_wet(pp)
                numer_send(cpart) % l(ii) = lninv_loc(ipoin)
             end do
          end if
       end do
    else
       do cpart = 0,PAR_CURRENT_SIZE-1
          npoin_to_send = npoin_send(cpart)
          if( npoin_to_send > 0 ) then
             kk = 0
             do ii = 1,npoin_to_send
                pp = my_part_to_point(cpart) % l(ii)
                do idime = 1,ndime
                   kk = kk + 1
                   coord_send(cpart) % a(kk) = xx(idime,pp)
                end do
             end do
          end if
       end do
    end if
    !
    ! Get points coordinates or node numbering
    !
    call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    if( itype == GLOBAL_NUMBERING ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(numer_send(cpart) % l,numer_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
          end if
       end do
    else
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(coord_send(cpart) % a,coord_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN=PAR_COMM_CURRENT4)
          end if
       end do
    end if
    
    call PAR_END_NON_BLOCKING_COMM(1_ip)

  end subroutine couplings_coordinates

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Intersection
  !> @details Compute intersection array
  !>
  !> Bounding box intersection. We will communicate only with these subdomains at first!
  !>
  !> For vector projection, send the wet points to all the partitions as we have no way
  !> to know where the point falls. Idea: consider only the partitions bounding boxes
  !> crossed by the segment formed by the wet node ( x_wet , x_wet + 10000 * projection_vector )
  !>
  !> +-----+------+-------+-----+
  !> |     |      | x_wet |     |
  !> |     |      |    x  |     |
  !> +-----+------+----|--+-----+
  !> |     |      |    |  |     |
  !> |     |      |    |  |     |
  !> +-----+------+----|--+-----+
  !> |     |      |    |  |     |
  !> |     |      |    |  |     |
  !> +-----+------+----|--+-----+
  !> ///////////////// | ////////
  !>                   |
  !>                   | projection_vector
  !>                   |
  !>                   \/
  !>                   x
  !>
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_intersection(&
       itype,PAR_CURRENT_SIZE,PAR_CURRENT_RANK,bin_structure,num_intersections4,intersection)
    
    integer(ip),               intent(in)            :: itype
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    integer(ip),               intent(in)            :: PAR_CURRENT_RANK
    type(typ_bin_structure),   intent(in)            :: bin_structure
    integer(4),                intent(out)           :: num_intersections4 
    logical(lg), pointer,      intent(inout)         :: intersection(:)
    integer(ip)                                      :: cpart,idime
    
    do cpart = 0,PAR_CURRENT_SIZE-1
       intersection(cpart) = .true.
       do idime = 1,ndime
          if(    bin_structure % part_comin(idime,PAR_CURRENT_RANK) > bin_structure % part_comax(idime,cpart) .or. &
               & bin_structure % part_comin(idime,cpart)            > bin_structure % part_comax(idime,PAR_CURRENT_RANK) ) then
             intersection(cpart) = .false.
          end if
       end do
    end do
    if( itype == BOUNDARY_VECTOR_PROJECTION .and. INOTMASTER ) then
       do cpart = 1,PAR_CURRENT_SIZE-1
          intersection(cpart) = .true.
       end do
    end if
    num_intersections4 = int(count(intersection,KIND=ip),KIND=4)

  end subroutine couplings_intersection

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Communication
  !> @details Setup up communiaction arrays for coupling
  !> 
  !-----------------------------------------------------------------------

  
  subroutine couplings_communication(&
       coupling,itype,number_wet_points,PAR_CURRENT_SIZE,&
       npoin_recv,decision_recv,&
       npoin_send,decision_send,&
       shapf,my_part_to_point,my_point_to_part)

    type(typ_color_coupling),  intent(inout)         :: coupling
    integer(ip),               intent(in)            :: itype
    integer(ip),               intent(in)            :: number_wet_points
    integer(ip),               intent(in)            :: PAR_CURRENT_SIZE
    integer(ip), pointer,      intent(in)            :: npoin_recv(:)
    type(i1p),   pointer,      intent(inout)         :: decision_recv(:)
    integer(ip), pointer,      intent(inout)         :: npoin_send(:)
    type(i1p),   pointer,      intent(inout)         :: decision_send(:)
    type(r2p),   pointer,      intent(inout)         :: shapf(:)
    type(i1p),   pointer,      intent(in)            :: my_part_to_point(:)
    type(i1p),   pointer,      intent(in)            :: my_point_to_part(:)

    integer(ip)                                      :: cpart,ii,kk,pp,ierro
    integer(ip)                                      :: irecv,ielem,ineig
    integer(ip)                                      :: isend,kelem,kboun,ksend
    integer(ip)                                      :: krecv,iboun,inodb,ipoin
    integer(ip)                                      :: inode,kpoin,lrecv,nenti
    integer(ip)                                      :: lsend
    integer(ip), pointer                             :: list_neighbors(:,:)
    logical(lg), pointer                             :: list_source_nodes(:)
    
    nullify(list_neighbors)
    nullify(list_source_nodes)
    
    if( INOTMASTER ) then

       !-----------------------------------------------------------------
       !
       ! decision_send(cpart) % l(1:npoin_recv(cpart)) /= 0 => I am in charge of this point OF CPART
       ! decision_recv(cpart) % l(1:npoin_send(cpart)) /= 0 => I need point from CPART
       !
       !-----------------------------------------------------------------
       !
       ! Fill in type
       !
       call memory_alloca(memor_cou,'LIST_NEIGHBORS',vacal,list_neighbors,2_ip,PAR_CURRENT_SIZE,'INITIALIZE',1_ip,0_ip)
       !
       ! Number of subdomains to communicate with
       !
       coupling % commd % nneig = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          do ii = 1,npoin_recv(cpart)
             list_neighbors(1,cpart) = list_neighbors(1,cpart) + min(1_ip,decision_send(cpart) % l(ii))
          end do
          do ii = 1,npoin_send(cpart)
             list_neighbors(2,cpart) = list_neighbors(2,cpart) + min(1_ip,decision_recv(cpart) % l(ii))
          end do
          coupling % commd % nneig = coupling % commd % nneig &
               + min(max(list_neighbors(1,cpart),list_neighbors(2,cpart)),1_ip)
       end do
       !
       ! Send and receives sizes of subdomains
       ! COMMD % LSEND_SIZE(:)
       ! COMMD % LRECV_SIZE(:)
       !
       call memory_alloca(memor_cou,'COUPLING % COMMD % NEIGHTS',   vacal,coupling % commd % neights,   coupling % commd % nneig)
       call memory_alloca(memor_cou,'COUPLING % COMMD % LSEND_SIZE',vacal,coupling % commd % lsend_size,coupling % commd % nneig+1)
       call memory_alloca(memor_cou,'COUPLING % COMMD % LRECV_SIZE',vacal,coupling % commd % lrecv_size,coupling % commd % nneig+1)
       ineig = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          isend = list_neighbors(1,cpart)
          irecv = list_neighbors(2,cpart)
          if( max(isend,irecv) > 0 ) then
             ineig = ineig + 1
             coupling % commd % neights(ineig)    = cpart
             coupling % commd % lsend_size(ineig) = isend
             coupling % commd % lrecv_size(ineig) = irecv
          end if
       end do
       !
       ! Size to list
       !
       ksend                            = coupling % commd % lsend_size(1)
       krecv                            = coupling % commd % lrecv_size(1)
       coupling % commd % lsend_size(1) = 1
       coupling % commd % lrecv_size(1) = 1
       do ineig = 2,coupling % commd % nneig + 1
          lsend                                = coupling % commd % lsend_size(ineig)
          lrecv                                = coupling % commd % lrecv_size(ineig)
          coupling % commd % lsend_size(ineig) = coupling % commd % lsend_size(ineig-1) + ksend
          coupling % commd % lrecv_size(ineig) = coupling % commd % lrecv_size(ineig-1) + krecv
          ksend                                = lsend
          krecv                                = lrecv
       end do
       coupling % commd % lsend_dim = coupling % commd % lsend_size(coupling % commd % nneig + 1) - 1
       coupling % commd % lrecv_dim = coupling % commd % lrecv_size(coupling % commd % nneig + 1) - 1
       !
       ! Order points
       ! KK is the order I receive
       !
       call memory_alloca(memor_cou,'COUPLING % GEOME % STATUS',vacal,coupling % geome % status,coupling % wet % number_wet_points)
       kk = 0

       do cpart = 0,PAR_CURRENT_SIZE-1
          do ii = 1,npoin_send(cpart)
             pp = my_part_to_point(cpart) % l(ii)
             if( associated(my_point_to_part(pp) % l) ) then
                if( my_point_to_part(pp) % l(1) == cpart ) then
                   !
                   ! CPART has wet point pp
                   !
                   kk = kk + 1
                   coupling % geome % status(pp) = kk
                end if
             end if
          end do
       end do
       !
       ! Allocate geometrical information
       ! NENTI= number of entities (elements/boundaries/nodes)
       !
       nenti = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          nenti = nenti + list_neighbors(1,cpart)
       end do

       if( itype == ELEMENT_INTERPOLATION ) then
          !
          ! Element interpolation
          !
          coupling % geome % nelem_source = nenti
          call memory_alloca(memor_cou,'COUPLING % GEOME % SHAPF'       ,vacal,coupling % geome % shapf,mnode, coupling % geome % nelem_source)
          call memory_alloca(memor_cou,'COUPLING % GEOME % LELEM_SOURCE',vacal,coupling % geome % lelem_source,coupling % geome % nelem_source)
          call memory_alloca(memor_cou,'LIST_SOURCE_NODES'              ,vacal,list_source_nodes,npoin)

          kelem = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                ielem = decision_send(cpart) % l(ii)
                if( ielem > 0 ) then
                   kelem = kelem + 1
                   coupling % geome % lelem_source(kelem)  = ielem
                   do inode = 1,mnode
                      coupling % geome % shapf(inode,kelem) = shapf(cpart) % a(inode,ii)
                   end do
                   do inode = 1,lnnod(ielem)
                      ipoin = lnods(inode,ielem)
                      list_source_nodes(ipoin) = .true.
                   end do
                end if
             end do
          end do

       else if( itype == NEAREST_BOUNDARY_NODE      .or. &
            &   itype == NEAREST_ELEMENT_NODE       .or. &
            &   itype == BOUNDARY_VECTOR_PROJECTION .or. &
            &   itype == GLOBAL_NUMBERING           .or. &
            &   itype == SAME_COORDINATE            ) then
          !
          ! Nearest boundary node/global numbering
          !
          coupling % geome % npoin_source = nenti
          call memory_alloca(memor_cou,'COUPLING % GEOME % LPOIN_SOURCE',vacal,coupling % geome % lpoin_source,coupling % geome % npoin_source)

          kpoin = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                ipoin = decision_send(cpart) % l(ii)
                if( ipoin > 0 ) then
                   kpoin = kpoin + 1
                   coupling % geome % lpoin_source(kpoin)  = ipoin
                end if
             end do
          end do

       else if( itype == BOUNDARY_INTERPOLATION ) then       
          !
          ! Boundary interpolation
          !
          coupling % geome % nboun_source = nenti

          call memory_alloca(memor_cou,'COUPLING % GEOME % SHAPF'       ,vacal,coupling % geome % shapf,mnodb, coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'COUPLING % GEOME % LBOUN_SOURCE',vacal,coupling % geome % lboun_source,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'LIST_SOURCE_NODES'              ,vacal,list_source_nodes,npoin)

          kboun = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                iboun = decision_send(cpart) % l(ii)
                if( iboun > 0 ) then
                   kboun = kboun + 1
                   coupling % geome % lboun_source(kboun)  = iboun
                   do inodb = 1,mnodb
                      coupling % geome % shapf(inodb,kboun) = shapf(cpart) % a(inodb,ii)
                   end do
                   do inodb = 1,lnnob_cou(iboun)
                      ipoin = lnodb_cou(inodb,iboun)
                      list_source_nodes(ipoin) = .true.
                   end do
                end if
             end do
          end do

       else if( itype == STRESS_PROJECTION .or. itype == PROJECTION ) then       
          !
          ! Projection
          !
          coupling % geome % nboun_source = nenti
          call memory_alloca(memor_cou,'COUPLING % GEOME % SHAPF'       ,vacal,coupling % geome % shapf,mnodb,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'COUPLING % GEOME % LBOUN_SOURCE',vacal,coupling % geome % lboun_source,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'LIST_SOURCE_NODES'              ,vacal,list_source_nodes,npoin)

          kboun = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                iboun = decision_send(cpart) % l(ii)
                if( iboun > 0 ) then
                   kboun = kboun + 1
                   coupling % geome % lboun_source(kboun)  = iboun
                   do inodb = 1,mnodb
                      coupling % geome % shapf(inodb,kboun) = shapf(cpart) % a(inodb,ii)
                   end do
                   do inodb = 1,lnnob_cou(iboun)
                      ipoin = lnodb_cou(inodb,iboun)
                      list_source_nodes(ipoin) = .true.
                   end do
                   !write(*,'(i2,2(1x,e12.6),i5,i5)') kfl_paral,coord_recv(cpart) % a((ii-1)*2+1),coord_recv(cpart) % a((ii-1)*2+2),lninv_loc(lnodb(1:2,iboun))
                end if
             end do
          end do
       end if
       !
       ! Allocate sched_perm for scheduling permutation
       !
       call memory_alloca(memor_cou,'COUPLING % GEOME % SCHED_PERM',vacal,coupling % geome % sched_perm,coupling % commd % lrecv_dim)
       !
       ! Count number of source nodes if not already done
       !
       if( associated(list_source_nodes) ) then
          coupling % geome % npoin_source = 0
          do ipoin = 1,npoin
             if( list_source_nodes(ipoin) ) coupling % geome % npoin_source = coupling % geome % npoin_source + 1
          end do
          if( associated(coupling % geome % lpoin_source) ) then
             call runend('POINTER ALREADY ASSOCIATED')
          else
             call memory_alloca(memor_cou,'COUPLING % GEOME % LPOIN_SOURCE',vacal,coupling % geome % lpoin_source,coupling % geome % npoin_source)
          end if
          kpoin = 0
          do ipoin = 1,npoin
             if( list_source_nodes(ipoin) ) then
                kpoin = kpoin + 1
                coupling % geome % lpoin_source(kpoin) = ipoin
             end if
          end do
       end if
       !
       ! Allocate memory for previous values
       !
       !call memory_alloca(memor_cou,'VALUES',vacal,coupling % values,coupling % geome % number_wet_points)
       !
       ! Check points that do not have a host partition
       ! If receive:        COUPLI % GEOME % LRECV_DIM values from other subdomains
       ! If need values on: COUPLI % GEOME % NUMBER_WET_POINTS values
       ! COUPLI % GEOME % STATUS(IPOIN) = JPOIN (from 1 to COUPLI % GEOME % LRECV_DIM)
       !
       ierro = 0
       ii    = 0
       do pp = 1,number_wet_points
          if( associated(my_point_to_part(pp) % l) ) then
             if( my_point_to_part(pp) % l(1) == 0 ) then
                ierro = 1
             else
                ii = ii + 1
             end if
          else
             ierro = 1
          end if
       end do

    end if
    
    call memory_deallo(memor_cou,'LIST_NEIGHBORS'    ,vacal,list_neighbors)
    call memory_deallo(memor_cou,'LIST_SOURCE_NODES' ,vacal,list_source_nodes)

  end subroutine couplings_communication

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-22
  !> @brief   Source entities
  !> @details Determin source entities
  !>          LESOU, LBSOU, LNSOU: Source elements/boundaries/nodes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_source_entities(&
       coupling,itype,lesou,lbsou,lnsou,&
       CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,&
       CANDIDATE_SOURCE_BOUNDARIES)
    
    type(typ_color_coupling),          intent(in)             :: coupling
    integer(ip),                       intent(in)             :: itype
    integer(ip),              pointer, intent(inout)          :: lesou(:)
    logical(lg),              pointer, intent(inout)          :: lbsou(:)
    logical(lg),              pointer, intent(inout)          :: lnsou(:)
    logical(lg),              pointer, intent(in),   optional :: CANDIDATE_SOURCE_NODES(:)
    integer(ip),              pointer, intent(in),   optional :: CANDIDATE_SOURCE_ELEMENTS(:)
    logical(lg),              pointer, intent(in),   optional :: CANDIDATE_SOURCE_BOUNDARIES(:)
    integer(ip)                                               :: ielem,inode,ipoin

    if( INOTMASTER ) then

       if( present(CANDIDATE_SOURCE_ELEMENTS) ) then
          lesou => CANDIDATE_SOURCE_ELEMENTS
       else
          if(     itype == NEAREST_ELEMENT_NODE  .or. &
               &  itype == NEAREST_BOUNDARY_NODE .or. &
               &  itype == ELEMENT_INTERPOLATION ) then
             call memory_alloca(memor_cou,'LESOU',vacal,lesou,nelem)
             if( coupling % kind == BETWEEN_SUBDOMAINS ) then
                do ielem = 1,nelem
                   if( lesub(ielem) == coupling % subdomain_source ) then
                      lesou(ielem) = 1
                   end if
                end do
             else
                if( nelem > 0 ) lesou = 1
             end if
          end if
       end if

       if( present(CANDIDATE_SOURCE_BOUNDARIES) ) then
          lbsou => CANDIDATE_SOURCE_BOUNDARIES
       end if

       if( present(CANDIDATE_SOURCE_NODES) ) then
          lnsou => CANDIDATE_SOURCE_NODES
       else if( associated(lesou) ) then
          if(     itype == NEAREST_ELEMENT_NODE .or. &
               &  itype == NEAREST_BOUNDARY_NODE ) then
             call memory_alloca(memor_cou,'LNSOU',vacal,lnsou,npoin)
             do ielem = 1,nelem
                if( lesou(ielem) == 1 ) then
                   do inode = 1,lnnod(ielem)
                      ipoin = lnods(inode,ielem)
                      lnsou(ipoin) = .true.
                   end do
                end if
             end do
          end if
       end if

    end if

  end subroutine couplings_source_entities
  
end module mod_couplings_setup
!> @}
