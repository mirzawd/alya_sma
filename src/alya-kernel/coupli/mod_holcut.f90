!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Coupling
!> @{
!> @name    ToolBox for elemental and general geometrical operations
!> @file    mod_elmgeo.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for graphs and renumbering.
!> @details ToolBox for graphs and renumbering. Uses METIS_NodeND,
!>          (Node dissection) for renumbering
!
!-----------------------------------------------------------------------

module mod_holcut
  use def_master,         only : INOTMASTER,zeror
  use def_master,         only : IMASTER
  use def_master,         only : lninv_loc,kfl_paral,leinv_loc
  use def_master,         only : current_code
  use def_kintyp,         only : ip,rp,lg,r1p,r2p,i1p 
  use def_elmtyp,         only : NOHOL,NOFRI,NOFEM,ELHOL,BOFRI
  use def_elmgeo,         only : element_type
  use def_domain,         only : lnods,lnodb,lboel,lelbo
  use def_domain,         only : lesub,nelem,lnnod
  use def_domain,         only : mnodb,mnode,nboun,ndime,ndimb
  use def_domain,         only : nbono,coord,ltopo,ltype
  use def_domain,         only : lbono,npoin,lboch,lnnob
  use def_domain,         only : ltypb,nnode
  use def_domain,         only : lnoch,lelch
  use def_coupli,         only : kdtree_typ
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_type
  use def_coupli,         only : ON_CHIMERA_MESH
  use def_coupli,         only : ON_EMBEDDED_MESH
  use def_coupli,         only : memor_cou
  use def_coupli,         only : nboun_cou
  use def_coupli,         only : lnodb_cou
  use def_coupli,         only : ltypb_cou
  use def_coupli,         only : lboch_cou
  use def_coupli,         only : lnnob_cou
  use def_coupli,         only : lboel_cou
  use def_coupli,         only : lelbo_cou
  use def_coupli,         only : number_of_holes
  use def_coupli,         only : typ_coupling_wet
  use mod_maths,          only : maths_mapping_coord_to_3d
  use mod_maths,          only : maths_in_box
  use mod_parall,         only : par_part_in_color
  use mod_parall,         only : par_code_zone_subd_to_color
  use mod_parall,         only : color_target
  use mod_parall,         only : color_source
  use mod_parall,         only : par_bin_comin
  use mod_parall,         only : par_bin_comax
  use mod_parall,         only : par_bin_part
  use mod_parall,         only : par_bin_boxes
  use mod_parall,         only : par_bin_size
  use mod_parall,         only : PAR_COMM_CURRENT
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_parall,         only : PAR_MY_WORLD_RANK
  use mod_parall,         only : par_part_comin
  use mod_parall,         only : par_part_comax
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_resize
  use mod_memory,         only : memory_copy
  use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_kdtree,         only : typ_kdtree
  use mod_kdtree,         only : kdtree_construct
  use mod_kdtree,         only : kdtree_initialize
  use mod_kdtree,         only : kdtree_nearest_boundary
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_COMM_SPLIT
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_couplings,      only : I_HAVE_A_FRINGE_ELEMENT
  use def_kintyp_mesh_basic
  use def_mpi
#include "def_mpi.inc"

  implicit none
  private

  type(mesh_type_basic)        :: mesh_free_bou
  type(typ_kdtree)             :: kdtree_chimi
  real(rp)                     :: hole_box_min(3), hole_box_max(3)
  logical(lg),         pointer :: lmask(:) 
  integer(ip),         pointer :: lwetno(:) ! list of wet nodes
  logical(lg)                  :: i_am_embed
  logical(lg)                  :: i_am_backg
  logical(lg)                  :: i_am_patch
  logical(lg)                  :: i_am_source
  logical(lg)                  :: i_am_target

  public :: cou_holcut
  public :: cou_holcut_multicode

contains 

  !-----------------------------------------------------------------------
  !
  !> @author  Beatriz Eguzkitza
  !> @date    03/10/2014
  !> @brief   ToolBox for hole-cutting process
  !> @details ToolBox for hole-cutting process in Chimera-type problems
  !
  !-----------------------------------------------------------------------

  subroutine cou_holcut()

    integer(ip)            :: jcoup,ipoin,idime,ielem,inode,pnode,icoup
    integer(ip)            :: iboun,jboun,iboun_wet,iboun_min,inodb,pnodb,pblty
    integer(ip)            :: ielty,iface,nboun_ori,nboun_hole
    integer(ip)            :: cont,kk,cont_free,cont_hole
    integer(ip)            :: cpart,cz,nb_send,nb_recv
    integer(ip)            :: code_target,subdomain_target
    integer(ip)            :: subdomain_source
    real(rp)               :: time1,time2
    real(rp)               :: dista_min,proje(3),coord_test(3)
    MY_MPI_COMM            :: PAR_COMM_SAVE
    integer(ip)            :: PAR_CURRENT_RANK
    integer(ip)            :: PAR_CURRENT_SIZE
    logical(lg)            :: i_am_in_source
    logical(lg)            :: i_am_in_source_bound
    logical(lg)            :: i_am_in_target
    logical(lg)            :: i_am_a_fringe_element
    integer(ip), pointer   :: PAR_WORLD_RANKS(:)
    integer(ip), pointer   :: what_send_patch_bound(:)
    integer(ip), pointer   :: what_recv_patch_bound(:)
    integer(ip), pointer   :: lnfri(:)
    integer(ip), pointer   :: what_send_npoin(:)
    integer(ip), pointer   :: what_recv_npoin(:)
    integer(ip), pointer   :: lnodb_aux(:,:)
    integer(ip), pointer   :: ltypb_aux(:)
    integer(ip), pointer   :: lboch_aux(:)
    integer(ip), pointer   :: lnnob_aux(:)
    integer(ip), pointer   :: lboel_aux(:,:)
    integer(ip), pointer   :: lelbo_aux(:)
    logical(lg), pointer   :: lnsou(:)
    logical(lg), pointer   :: intersection(:)
    logical(lg), pointer   :: whos_in_target(:)
    logical(lg), pointer   :: whos_in_source_bound(:)
    type(r2p),   pointer   :: send_coord(:) 
    type(r2p),   pointer   :: recv_coord(:) 
    type(i1p),   pointer   :: permu(:)
    type(r1p),   pointer   :: send_nhole_distance(:)
    type(r1p),   pointer   :: recv_nhole_distance(:)
    type(r2p),   pointer   :: send_nhole_projection(:)
    type(r2p),   pointer   :: recv_nhole_projection(:)
    logical(lg)            :: make_a_hole
    logical(lg)            :: there_is_a_hole
    logical(lg)            :: i_am_in_chimera
    !
    ! Nullify
    ! 
    nullify(PAR_WORLD_RANKS)
    nullify(what_send_patch_bound)
    nullify(what_recv_patch_bound)
    nullify(lnfri)
    nullify(what_send_npoin)
    nullify(what_recv_npoin)
    nullify(lnodb_aux)
    nullify(ltypb_aux)
    nullify(lboch_aux)
    nullify(lnnob_aux)
    nullify(lboel_aux)
    nullify(lelbo_aux)
    nullify(lnsou)
    nullify(intersection)
    nullify(whos_in_target)
    nullify(whos_in_source_bound)
    nullify(send_coord) 
    nullify(recv_coord) 
    nullify(permu)
    nullify(send_nhole_distance)
    nullify(recv_nhole_distance)
    nullify(send_nhole_projection)
    nullify(recv_nhole_projection)

    !--------------------------------------------------------------------
    !
    ! Define what to do
    !
    !--------------------------------------------------------------------

    i_am_in_chimera = .false.
    there_is_a_hole = .false.
    do icoup = 1,mcoup
       color_target = coupling_type(icoup) % color_target
       color_source = coupling_type(icoup) % color_source
       if( coupling_type(icoup) % where_type == ON_CHIMERA_MESH ) then
          if( I_AM_IN_COLOR(color_target) ) i_am_in_chimera = .true.
          there_is_a_hole = .true.
       end if
    end do
    if( IMASTER ) i_am_in_chimera = .false.
    if( .not. there_is_a_hole ) return

    !--------------------------------------------------------------------
    !
    ! Identify hole nodes LNOCH
    !
    !--------------------------------------------------------------------

    do icoup = 1,mcoup

       color_target = coupling_type(icoup) % color_target
       color_source = coupling_type(icoup) % color_source

       if(    coupling_type(icoup) % where_type == ON_CHIMERA_MESH .and. &
            ( I_AM_IN_COLOR(color_target) .or. I_AM_IN_COLOR(color_source) ) ) then
          make_a_hole = .true.
       else
          make_a_hole = .false.
       end if

       if( make_a_hole ) then

          call cputim(time1)
          code_target          = coupling_type(icoup) % code_target
          subdomain_source     = coupling_type(icoup) % subdomain_source
          subdomain_target     = coupling_type(icoup) % subdomain_target
          jcoup                = coupling_type(icoup) % mirror_coupling !Source:background / Target:patch
          i_am_in_source_bound = .false.
          !
          ! Communicator
          !
          PAR_COMM_SAVE    = PAR_COMM_CURRENT
          PAR_COMM_CURRENT = PAR_COMM_COLOR(color_target,color_source) 
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)

          if( I_AM_IN_COLOR(color_target) .or. I_AM_IN_COLOR(color_source) ) then
             !
             ! Allocate memory
             !
             cz = PAR_CURRENT_SIZE
             call memory_alloca(memor_cou,'INTERSECTION'         ,'cou_holcut',intersection   ,       cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'PAR_WORLD_RANKS'      ,'cou_holcut',PAR_WORLD_RANKS,       cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHAT_SEND_NPOIN'      ,'cou_holcut',what_send_npoin,       cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHAT_RECV_NPOIN'      ,'cou_holcut',what_recv_npoin,       cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHOS_IN_SOURCE_BOUND' ,'cou_holcut',whos_in_source_bound , cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHOS_IN_TARGET'       ,'cou_holcut',whos_in_target ,       cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHAT_SEND_PATCH_BOUND','cou_holcut',what_send_patch_bound, cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'WHAT_RECV_PATCH_BOUND','cou_holcut',what_recv_patch_bound, cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'SEND_COORD'           ,'cou_holcut',send_coord,            cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'RECV_COORD'           ,'cou_holcut',recv_coord,            cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'PERMU'                ,'cou_holcut',permu,                 cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'SEND_NHOLE_DISTANCE'  ,'cou_holcut',send_nhole_distance,   cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'RECV_NHOLE_DISTANCE'  ,'cou_holcut',recv_nhole_distance,   cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'SEND_NHOLE_PROJECTION','cou_holcut',send_nhole_projection, cz,  'INITIALIZE',0_ip)
             call memory_alloca(memor_cou,'RECV_NHOLE_PROJECTION','cou_holcut',recv_nhole_projection, cz,  'INITIALIZE',0_ip)
             !
             !if (coupling_type(icoup) % where_type == ON_CHIMERA_MESH  )then !Source:patch/ Target:background  !.and. INOTMASTER
             !
             ! Compute the bounding box of the patch, two strategies (should be optimized)
             ! Stategy 1
             ! 1. Patch CPUs compute their bounding box
             ! 2. All reduce between the patch (maybe between the ones who have the interface only)
             ! 3. One patch CPU broadcasts the result to the background CPUs
             ! Strategy 2 (the one implemented)
             ! 1. Patch CPUs compute their bounding box
             ! 2. All reduce in the coupling communicator
             !
             hole_box_min(1:3) =  huge(1.0_rp)
             hole_box_max(1:3) = -huge(1.0_rp)

             if( jcoup == 0 ) then

                call runend('COU_HOLCUT: NO MIRROR COUPLING FOUND')

             else if( jcoup /= 0 .and. I_AM_IN_COLOR(color_source) .and. INOTMASTER ) then

                if( coupling_type(jcoup) % wet % nboun_wet > 0 ) i_am_in_source_bound = .true.

                do iboun_wet = 1,coupling_type(jcoup) % wet % nboun_wet
                   iboun = coupling_type(jcoup) % wet % lboun_wet(iboun_wet) !!!COMPROBAR
                   pblty = abs(ltypb(iboun))
                   pnodb = nnode(pblty)
                   do inodb = 1,pnodb
                      ipoin = lnodb(inodb,iboun)
                      hole_box_min(1:ndime) = min( hole_box_min(1:ndime) , coord(1:ndime,ipoin) ) 
                      hole_box_max(1:ndime) = max( hole_box_max(1:ndime) , coord(1:ndime,ipoin) ) 
                   end do
                end do
                call kdtree_initialize(kdtree_chimi)
                call kdtree_construct(&
                     nboun,npoin,lnodb,ltypb,coord,kdtree_chimi,&
                     coupling_type(jcoup) % wet % lboun_wet)  
             end if
             call PAR_MIN(ndime,hole_box_min,'IN CURRENT COUPLING')
             call PAR_MAX(ndime,hole_box_max,'IN CURRENT COUPLING')
             !
             ! shoudl be optimized
             !    (tengo un array que me da si las particiones son del source o del target)
             !
             if( IMASTER ) then
                i_am_in_source = .false.
                i_am_in_target = .false.
             else
                i_am_in_source = I_AM_IN_COLOR(color_source)
                i_am_in_target = I_AM_IN_COLOR(color_target)
             end if
             if( i_am_in_target ) then
                do idime = 1,ndime
                   if(    hole_box_min(idime)                     >= par_part_comax(idime,PAR_MY_WORLD_RANK) .or. &
                        & par_part_comin(idime,PAR_MY_WORLD_RANK) >= hole_box_max(idime) ) then
                      i_am_in_target = .false.
                   end if
                end do
             end if
             call PAR_ALLGATHER(i_am_in_source_bound,whos_in_source_bound,1_4,'IN CURRENT COUPLING') 
             call PAR_ALLGATHER(i_am_in_target,      whos_in_target,      1_4,'IN CURRENT COUPLING')
             !
             ! INTERSECTION(CPART)
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if(    ( i_am_in_target .and. whos_in_source_bound(cpart) ) .or. &
                     & ( i_am_in_source_bound .and. whos_in_target(cpart) ) ) then
                   intersection(cpart) = .true.
                else
                   intersection(cpart) = .false.
                end if
             end do
             !
             ! Check background nodes located above the patch mesh
             !
             if( i_am_in_target .and. INOTMASTER ) then
                call memory_alloca(memor_cou,'LNSOU','cou_holcut',lnsou,npoin)
                do ielem = 1,nelem
                   if( lesub(ielem) == subdomain_target) then
                      do inode = 1,lnnod(ielem)
                         ipoin = lnods(inode,ielem)
                         lnsou(ipoin) = .true.
                      end do
                   end if
                end do
             end if

             do cpart = 1,PAR_CURRENT_SIZE-1
                if( intersection(cpart) ) then
                   if( i_am_in_target .and. whos_in_source_bound(cpart) ) then !!!Guillaume,he sustituido por if i_am_patch_bound(cpart)
                      do ipoin = 1,npoin
                         if( lnsou(ipoin) ) then
                            cont = 0
                            do idime = 1,ndime
                               if ( (coord(idime,ipoin)-hole_box_min(idime) > -zeror) .and.(coord(idime,ipoin)-hole_box_max(idime) < zeror) ) cont = cont + 1
                            end do
                            if(cont==ndime) then
                               what_send_npoin(cpart) = what_send_npoin(cpart) + 1
                            end if
                         end if
                      end do
                   end if
                   nb_send = 0
                   nb_recv = 0
                   if( i_am_in_target .and. whos_in_source_bound(cpart) ) nb_send = 1
                   if( i_am_in_source .and. whos_in_target(cpart) )       nb_recv = 1
                   call PAR_SEND_RECEIVE(nb_send,nb_recv,what_send_npoin(cpart:cpart),what_recv_npoin(cpart:cpart),'IN CURRENT COUPLING',cpart)
                   !
                   ! Compute and send/recv coords
                   !                   
                   if( nb_send == 1 ) then
                      call memory_alloca(memor_cou,'SEND_COORD % A',           'cou_holcut',send_coord(cpart)%a,ndime,what_send_npoin(cpart))   
                      call memory_alloca(memor_cou,'PERMU % L',                'cou_holcut',permu(cpart)%l,what_send_npoin(cpart))   
                      call memory_alloca(memor_cou,'RECV_NHOLE_DISTANCE % A',  'cou_holcut',recv_nhole_distance(cpart)   % a,what_send_npoin(cpart))   
                      call memory_alloca(memor_cou,'RECV_NHOLE_PROJECTION % A','cou_holcut',recv_nhole_projection(cpart) % a,ndime,what_send_npoin(cpart))   
                   end if
                   if( nb_recv == 1 ) then
                      call memory_alloca(memor_cou,'RECV_COORD % A',           'cou_holcut',recv_coord(cpart)%a,ndime,what_recv_npoin(cpart))     
                      call memory_alloca(memor_cou,'SEND_NHOLE_DISTANCE % A',  'cou_holcut',send_nhole_distance(cpart)   % a,what_recv_npoin(cpart)) 
                      call memory_alloca(memor_cou,'SEND_NHOLE_PROJECTION % A','cou_holcut',send_nhole_projection(cpart) % a,ndime,what_recv_npoin(cpart)) 
                   end if
                   if( i_am_in_target .and. whos_in_source_bound(cpart) ) then !!!Guillaume,he sustituido por if i_am_patch_bound(cpart)
                      kk = 0
                      do ipoin =1,npoin
                         if( lnsou(ipoin) ) then !soy un punto del background
                            cont = 0
                            do idime = 1,ndime  !optimizar este bucle
                               if ( (coord(idime,ipoin)-hole_box_min(idime) > -zeror) .and.(coord(idime,ipoin)-hole_box_max(idime) < zeror) ) cont = cont + 1
                            end do
                            if(cont==ndime) then
                               kk = kk + 1
                               send_coord(cpart) % a(1:ndime,kk) = coord(1:ndime,ipoin)
                               permu(cpart) % l(kk)              = ipoin
                            end if
                         end if
                      end do
                   end if
                   if( nb_send == 1 .or. nb_recv == 1 ) then
                      call PAR_SEND_RECEIVE(send_coord(cpart) % a,recv_coord(cpart) % a,'IN CURRENT COUPLING',cpart)
                   end if
                end if
             end do

             if( i_am_in_source_bound ) then
                do cpart = 0,PAR_CURRENT_SIZE-1
                   do ipoin = 1,what_recv_npoin(cpart)
                      do idime=1,ndime
                         coord_test(idime) = recv_coord(cpart)%a(idime,ipoin)
                      end do
                      dista_min = huge(1.0_rp)
                      call kdtree_nearest_boundary(coord_test(1),kdtree_chimi,iboun_min,dista_min,proje)  
                      send_nhole_distance(cpart)   % a(ipoin)         = dista_min 
                      send_nhole_projection(cpart) % a(1:ndime,ipoin) = proje(1:ndime) 
                      !
                      ! PARA CONSIDERAR SOLAPE:Llamar con ielem= lboel(inode+1,iboun_min) la subru elmgeo_where_is sin ulitmo arg:llist
                      !
                   end do
                end do
             end if

             do cpart = 0,PAR_CURRENT_SIZE-1
                if( intersection(cpart) ) then
                   call PAR_SEND_RECEIVE(send_nhole_distance(cpart)   % a,recv_nhole_distance(cpart)   % a,'IN CURRENT COUPLING ',cpart)
                   call PAR_SEND_RECEIVE(send_nhole_projection(cpart) % a,recv_nhole_projection(cpart) % a,'IN CURRENT COUPLING ',cpart)
                end if
             end do

             if( INOTMASTER .and. i_am_in_target ) then

                do ipoin = 1,npoin
                   dista_min = huge(1.0_rp)
                   if( lnsou(ipoin) ) then
                      do cpart = 0, PAR_CURRENT_SIZE-1
                         kk = 0
                         do while ( kk < what_send_npoin(cpart) )
                            kk = kk + 1
                            if( ipoin == permu(cpart)%l(kk) ) then
                               !
                               !choose abs(smalest distance) and decide if the point is inside 
                               !
                               if( abs(recv_nhole_distance(cpart) % a(kk)) < abs(dista_min) ) then
                                  dista_min = recv_nhole_distance(cpart) % a(kk)
                                  !bvess_defor_ker(1:ndime,ipoin)     = recv_nhole_projection(cpart) % a(1:ndime,kk)
                                  !kfl_fixno_defor_ker(1:ndime,ipoin) = 1
                               end if
                            end if
                         end do
                      end do
                      if( dista_min <= 0.0_rp ) lnoch(ipoin) = NOHOL
                   end if
                end do
             end if
             !
             ! Deallocate memory
             !
             call memory_deallo(memor_cou,'RECV_NHOLE_PROJECTION','cou_holcut',  recv_nhole_projection)
             call memory_deallo(memor_cou,'SEND_NHOLE_PROJECTION','cou_holcut',  send_nhole_projection)
             call memory_deallo(memor_cou,'RECV_NHOLE_DISTANCE'  ,'cou_holcut',  recv_nhole_distance)
             call memory_deallo(memor_cou,'SEND_NHOLE_DISTANCE'  ,'cou_holcut',  send_nhole_distance)
             call memory_deallo(memor_cou,'PERMU'                ,'cou_holcut',                permu)
             call memory_deallo(memor_cou,'RECV_COORD'           ,'cou_holcut',           recv_coord)
             call memory_deallo(memor_cou,'SEND_COORD'           ,'cou_holcut',           send_coord)
             call memory_deallo(memor_cou,'WHAT_RECV_PATCH_BOUND','cou_holcut',what_recv_patch_bound)
             call memory_deallo(memor_cou,'WHAT_SEND_PATCH_BOUND','cou_holcut',what_send_patch_bound)
             call memory_deallo(memor_cou,'WHOS_IN_TARGET'       ,'cou_holcut',       whos_in_target)
             call memory_deallo(memor_cou,'WHOS_IN_SOURCE_BOUND' ,'cou_holcut', whos_in_source_bound)
             call memory_deallo(memor_cou,'WHAT_RECV_NPOIN'      ,'cou_holcut',      what_recv_npoin)
             call memory_deallo(memor_cou,'WHAT_SEND_NPOIN'      ,'cou_holcut',      what_send_npoin)
             call memory_deallo(memor_cou,'LNSOU'                ,'cou_holcut',                lnsou)
             call memory_deallo(memor_cou,'PAR_WORLD_RANKS'      ,'cou_holcut',      par_world_ranks)
             call memory_deallo(memor_cou,'INTERSECTION'         ,'cou_holcut',         intersection)

          end if
          call cputim(time2)
          coupling_type(icoup) % cputim(10) = coupling_type(icoup) % cputim(10) + time2 - time1
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Construct hole
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then

       if( i_am_in_chimera ) then
          !
          ! Mark the hole elements: with all nodes nohol
          !
          do ielem =1,nelem
             pnode = lnnod(ielem)
             cont  = 0_ip
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                if( lnoch(ipoin) == NOHOL ) then
                   cont = cont + 1
                end if
             end do
             if( cont == pnode ) then
                lelch(ielem) = ELHOL
             end if
          end do
          !
          ! saco pelos
          ! 
          do ielem = 1,nelem
             if( lelch(ielem) == ELHOL )then
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   lnoch(ipoin) = 10
                end do
             end if
          end do

       end if

       call PAR_INTERFACE_NODE_EXCHANGE(lnoch,'MAX','IN THE WORLD')
       call memory_alloca(memor_cou,'LNFRI','cou_holcut',lnfri,npoin)   

       if( i_am_in_chimera ) then
          !
          ! LNOCH
          !
          do ipoin = 1,npoin
             if( lnoch(ipoin) == NOHOL ) then
                lnoch(ipoin) = NOFEM
             else if( lnoch(ipoin) == 10 ) then
                lnoch(ipoin) = NOHOL
             end if
          end do
          !
          ! Mark fringe nodes
          !
          do ielem = 1,nelem
             pnode     = lnnod(ielem) 
             cont_hole = 0_ip
             cont_free = 0_ip
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                if( lnoch(ipoin) == NOHOL ) cont_hole = cont_hole + 1
                if( lnoch(ipoin) == NOFEM ) cont_free = cont_free + 1
             end do
             if( cont_hole /= 0 .and. cont_free /= 0 ) then
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   if( lnoch(ipoin) == NOHOL ) lnfri(ipoin) = 1  
                end do
             end if
          end do
       end if

       call PAR_INTERFACE_NODE_EXCHANGE(lnfri,'MAX','IN THE WORLD')

       if( i_am_in_chimera ) then
          do ielem = 1,nelem
             do inode = 1,lnnod(ielem)
                ipoin = lnods(inode,ielem)
                if( lnfri(ipoin) == 1 ) lnoch(ipoin) = NOFRI
             end do
          end do
          !
          ! Put LTYPE < 0 for holes
          !
          do ielem = 1,nelem
             if( lelch(ielem) == ELHOL ) then
                ltype(ielem) = -abs(ltype(ielem))
             end if
          end do
       end if
       !
       ! Merge LNOD into LNODB_COU with the new interface
       !
       if( I_HAVE_A_FRINGE_ELEMENT() ) then

          call memory_alloca(memor_cou,'LNODB_AUX','cou_holcut',lnodb_aux,mnodb,nelem) 
          call memory_alloca(memor_cou,'LTYPB_AUX','cou_holcut',ltypb_aux,nelem)       
          call memory_alloca(memor_cou,'LBOCH_AUX','cou_holcut',lboch_aux,nelem)       
          call memory_alloca(memor_cou,'LNNOB_AUX','cou_holcut',lnnob_aux,nelem)       
          call memory_alloca(memor_cou,'LBOEL_AUX','cou_holcut',lboel_aux,mnodb,nelem)       
          call memory_alloca(memor_cou,'LELBO_AUX','cou_holcut',lelbo_aux,nelem)       

          nboun_hole = 0
          do ielem = 1,nelem
             i_am_a_fringe_element = .false.
             if( lelch(ielem) == ELHOL ) then
                ielty = abs(ltype(ielem))
                pnode = lnnod(ielem)
                inode = 0
                do while( inode < pnode )
                   inode = inode + 1
                   ipoin = lnods(inode,ielem)
                   if( lnoch(ipoin) == NOFRI ) then
                      inode = pnode
                      i_am_a_fringe_element = .true.
                   end if
                end do
                if( i_am_a_fringe_element ) then
                   do iface = 1,element_type(ielty) % number_faces
                      cont = 0
                      do inodb = 1,element_type(ielty) % node_faces(iface)
                         if( lnoch(lnods(element_type(ielty) % list_faces(inodb,iface),ielem)) == NOFRI ) then
                            !if(lnoch(faces(ielem)%l(inodb,iface)) == NOFRI )then
                            cont = cont + 1
                         end if
                      end do
                      if ( cont == element_type(ielty) % node_faces(iface) ) then
                         nboun_hole = nboun_hole + 1
                         ltypb_aux(nboun_hole) = element_type(ielty) % type_faces(iface)
                         lboch_aux(nboun_hole) = BOFRI
                         lnnob_aux(nboun_hole) = element_type(ielty) % node_faces(iface)
                         do inodb = 1,element_type(ielty) % node_faces(iface)
                            inode = element_type(ielty) % list_faces(inodb,iface) 
                            lboel_aux(inodb,nboun_hole) = inode
                            lnodb_aux(inodb,nboun_hole) = lnods(inode,ielem)  !!sobredimensiono  lnodb_aux?
                         end do
                         lelbo_aux(nboun_hole) = ielem 
                      end if
                   end do
                end if
             end if
          end do
          !
          ! Reallocate memory for boundary mesh
          !
          nboun_ori = nboun_cou
          nboun_cou = nboun_ori + nboun_hole

          if( nboun_hole /= 0 ) then

             if( number_of_holes == 0 ) then
                nullify(lnodb_cou)
                nullify(ltypb_cou)
                nullify(lboch_cou)
                nullify(lnnob_cou)
                nullify(lboel_cou)
                nullify(lelbo_cou)
                call memory_copy(memor_cou,'LNODB_COU','cou_holcut',lnodb,lnodb_cou,'DO_NOT_DEALLOCATE')
                call memory_copy(memor_cou,'LTYPB_COU','cou_holcut',ltypb,ltypb_cou,'DO_NOT_DEALLOCATE')
                call memory_copy(memor_cou,'LBOCH_COU','cou_holcut',lboch,lboch_cou,'DO_NOT_DEALLOCATE')
                call memory_copy(memor_cou,'LNNOB_COU','cou_holcut',lnnob,lnnob_cou,'DO_NOT_DEALLOCATE')                
                call memory_copy(memor_cou,'LBOEL_COU','cou_holcut',lboel,lboel_cou,'DO_NOT_DEALLOCATE')                
                call memory_copy(memor_cou,'LELBO_COU','cou_holcut',lelbo,lelbo_cou,'DO_NOT_DEALLOCATE')                
             end if

             number_of_holes = number_of_holes + 1

             call memory_resize(memor_cou,'LNODB_COU','cou_holcut',lnodb_cou,mnodb,nboun_cou)
             call memory_resize(memor_cou,'LTYPB_COU','cou_holcut',ltypb_cou,nboun_cou)
             call memory_resize(memor_cou,'LBOCH_COU','cou_holcut',lboch_cou,nboun_cou)
             call memory_resize(memor_cou,'LNNOB_COU','cou_holcut',lnnob_cou,nboun_cou)
             call memory_resize(memor_cou,'LBOEL_COU','cou_holcut',lboel_cou,mnodb,nboun_cou)
             call memory_resize(memor_cou,'LELBO_COU','cou_holcut',lelbo_cou,nboun_cou)
             !
             ! Merge new boundary mesh
             !         
             iboun = nboun_ori
             do jboun = 1,nboun_hole
                iboun = iboun + 1
                ltypb_cou(iboun)         = ltypb_aux(jboun)
                lboch_cou(iboun)         = lboch_aux(jboun)
                lnnob_cou(iboun)         = lnnob_aux(jboun)
                lnodb_cou(1:mnodb,iboun) = lnodb_aux(1:mnodb,jboun)
                lboel_cou(1:mnodb,iboun) = lboel_aux(1:mnodb,jboun)
                lelbo_cou(iboun)         = lelbo_aux(jboun)
             end do

          end if
          !
          ! Deallocate memory
          !
          call memory_deallo(memor_cou,'LNODB_AUX'            ,'cou_holcut',            lnodb_aux) 
          call memory_deallo(memor_cou,'LTYPB_AUX'            ,'cou_holcut',            ltypb_aux)       
          call memory_deallo(memor_cou,'LBOCH_AUX'            ,'cou_holcut',            lboch_aux)       
          call memory_deallo(memor_cou,'LNNOB_AUX'            ,'cou_holcut',            lnnob_aux)       
          call memory_deallo(memor_cou,'LBOEL_AUX'            ,'cou_holcut',            lboel_aux)       
          call memory_deallo(memor_cou,'LELBO_AUX'            ,'cou_holcut',            lelbo_aux)       

       end if

       call memory_deallo(memor_cou,'LNFRI'                ,'cou_holcut',                lnfri)

    end if

  end subroutine cou_holcut


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    21/03/2022
  !> @brief   Cut hole for multi-code coupled problems 
  !> @details Cut hole for multi-code coupled problems, ie: Embeded Finite
  !>          Element Method Coupling Technique (EFECT).
  !
  !-----------------------------------------------------------------------
  subroutine cou_holcut_multicode()

    implicit none

    logical(lg)         :: make_a_hole
    integer(ip)         :: icoup

    make_a_hole = embedded_exists()

    ! Make hole
    if( make_a_hole ) then

       do icoup = 1,mcoup
          call get_flags( icoup )

          if( i_am_embed ) then
             call setup_communicators()
             call allocate_variables()
             call build_kd_tree( icoup )

             if( i_am_backg ) then
                call make_hole( icoup )
                call build_hole_boundary()
                call merge_hole_boundary()
             end if

             call deallocate_variables()

          end if
       end do
    end if

  end subroutine cou_holcut_multicode


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Get flags for specific coupling in EFECT method
  !> @details Get flags for specific coupling in EFECT method
  !
  !-----------------------------------------------------------------------
  subroutine get_flags(icoup)

    implicit none

    integer(ip), intent(in)  :: icoup

    i_am_embed  = .false.
    i_am_source = .false.
    i_am_target = .false.

    if( current_code == coupling_type(icoup) % code_source) i_am_source = .true.
    if( current_code == coupling_type(icoup) % code_target) i_am_target = .true.

    if( coupling_type(icoup) % where_type == ON_EMBEDDED_MESH ) then
        i_am_embed = .true.
        if( i_am_source ) i_am_patch = .true. ! source is solid patch
        if( i_am_target ) i_am_backg = .true. ! target is fluid background
    else
        if( i_am_source ) i_am_backg = .true.
        if( i_am_target ) i_am_patch = .true.
    end if

  end subroutine get_flags


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Check if embedded coupling exists in order to make a hole
  !> @details Check if embedded coupling exists in order to make a hole
  !
  !-----------------------------------------------------------------------
  function embedded_exists()

    implicit none

    integer(ip) :: icoup
    logical(lg) :: embedded_exists

    embedded_exists = .false. 
    do icoup = 1,mcoup
       if( coupling_type(icoup) % where_type == ON_EMBEDDED_MESH) embedded_exists = .true.
    end do

  end function embedded_exists


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    26/03/2022
  !> @brief   Setup MPI communicators for multi-code hole-cutting
  !> @details Setup MPI communicators for multi-code hole-cutting
  !
  !-----------------------------------------------------------------------
  subroutine setup_communicators()

     implicit none

  end subroutine setup_communicators


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Build kd-tree
  !> @details Build kd-tree for EFECT
  !
  !-----------------------------------------------------------------------
  subroutine build_kd_tree(icoup)

     implicit none

     integer(ip),            intent(in) :: icoup
     type(typ_coupling_wet), pointer    :: mirror_wet 
     integer(ip)                        :: jcoup_mirror
     !
     ! Initialize variables 
     !
     jcoup_mirror =  coupling_type(icoup) % mirror_coupling
     mirror_wet   => coupling_type(jcoup_mirror) % wet

  end subroutine build_kd_tree 


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @brief   Create holes on wet elements 
  !> @details Set background wet elements to be "hole" elements. 
  !>          TODO: Add hair-removal.
  !
  !-----------------------------------------------------------------------
  subroutine make_hole(icoup)

     use mod_par_element_loop,   only : par_element_loop 

     implicit none

     integer(ip), intent(in)         :: icoup
     type(typ_coupling_wet), pointer :: mirror_wet 
     integer(ip)                     :: kpoin, ipoin, cont
     integer(ip)                     :: pnode, ielem, inode
     integer(ip)                     :: jcoup_mirror

     ! Nullifications
     nullify(mirror_wet)

     ! Get wet mirror coupling 
     jcoup_mirror =  coupling_type(icoup) % mirror_coupling
     mirror_wet   => coupling_type(jcoup_mirror) % wet
 
     if ( INOTMASTER ) then
  
        ! Save wet nodes in LIST_WET_NODES and mark interior wet nodes as NOHOL in LNOCH 
        do kpoin = 1,mirror_wet % npoin_wet
           ipoin = mirror_wet % lpoin_wet(kpoin)
           ! Save wet node in list
           lwetno(ipoin) = 1
           ! Mark interior wet node as NOHOL
           if (mirror_wet % kfl_fringe_wetnodes(ipoin) < 0) then
              lnoch(ipoin) = NOHOL
           end if
        end do

        ! Parallel exchange of LNOCH
        call PAR_INTERFACE_NODE_EXCHANGE(lnoch, 'MAX','IN THE WORLD')
        call PAR_INTERFACE_NODE_EXCHANGE(lwetno,'MAX','IN THE WORLD')

        ! Mark the hole elements with all nodes NOHOL as ELHOL
        do ielem = 1,nelem
           pnode = lnnod(ielem)
           cont  = 0_ip
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              if( lwetno(ipoin) == 1 ) then
                 cont = cont + 1
              end if
           end do
           if( cont == pnode ) then
              lelch(ielem) =  ELHOL
              ltype(ielem) = -abs(ltype(ielem))
           end if
        end do
     end if

     ! Redefine arrays for hybrid parallelization 
     call par_element_loop()

  end subroutine make_hole


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    25/04/2022
  !> @brief   Build mesh structure for hole boundary 
  !> @details Build mesh structure for boundary of hole in background mesh
  !
  !-----------------------------------------------------------------------
  subroutine build_hole_boundary()

     use def_kermod, only      : ndivi
     use def_domain, only      : meshe

     implicit none

     type(mesh_type_basic)    :: mesh_hole
     type(mesh_type_basic)    :: mesh_free

     call mesh_free     % init('MY_MESH')
     call mesh_free_bou % init('BOUNDARY')

     call create_mask( lmask )
     call mesh_free     % extract( meshe(ndivi) , lmask )
     call mesh_free_bou % boundary_mesh( mesh_free )
     call mesh_free_bou % grandpa_permutation( mesh_free )

     call mesh_free     % deallo()
     call mesh_hole     % deallo()

  end subroutine build_hole_boundary


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    25/04/2022
  !> @brief   Merge hole boundary to existing mesh boundary
  !> @details Merge hole boundary to existing mesh boundary
  !
  !-----------------------------------------------------------------------
  subroutine merge_hole_boundary()

    implicit none

    integer(ip)             :: ipoin,ielem,inode,pnode
    integer(ip)             :: iboun,jboun,inodb,ielty
    integer(ip)             :: iface,nboun_ori,nboun_hole
    integer(ip)             :: cont
    logical(lg)             :: i_am_a_fringe_element
    integer(ip), pointer    :: lnodb_aux(:,:)
    integer(ip), pointer    :: ltypb_aux(:)
    integer(ip), pointer    :: lboch_aux(:)
    integer(ip), pointer    :: lnnob_aux(:)
    integer(ip), pointer    :: lboel_aux(:,:)
    integer(ip), pointer    :: lelbo_aux(:)

    ! Nullifications
    nullify(lnodb_aux)
    nullify(ltypb_aux)
    nullify(lboch_aux)
    nullify(lnnob_aux)
    nullify(lboel_aux)
    nullify(lelbo_aux)

    ! Allocations
    call memory_alloca(memor_cou,'LNODB_AUX','merge_hole_boundary',lnodb_aux,mnodb,nelem) 
    call memory_alloca(memor_cou,'LTYPB_AUX','merge_hole_boundary',ltypb_aux,nelem)       
    call memory_alloca(memor_cou,'LBOCH_AUX','merge_hole_boundary',lboch_aux,nelem)       
    call memory_alloca(memor_cou,'LNNOB_AUX','merge_hole_boundary',lnnob_aux,nelem)       
    call memory_alloca(memor_cou,'LBOEL_AUX','merge_hole_boundary',lboel_aux,mnodb,nelem)       
    call memory_alloca(memor_cou,'LELBO_AUX','merge_hole_boundary',lelbo_aux,nelem)       

    nboun_hole = 0
    do ielem = 1,nelem
       i_am_a_fringe_element = .false.
       if( lelch(ielem) == ELHOL ) then
          ielty = abs(ltype(ielem))
          pnode = lnnod(ielem)
          inode = 0
          do while( inode < pnode )
             inode = inode + 1
             ipoin = lnods(inode,ielem)
             if( lnoch(ipoin) == NOFRI ) then
                inode = pnode
                i_am_a_fringe_element = .true.
             end if
          end do
          if( i_am_a_fringe_element ) then
             do iface = 1,element_type(ielty) % number_faces
                cont = 0
                do inodb = 1, element_type(ielty) % node_faces(iface)
                   inode = element_type(ielty) % list_faces(inodb,iface)
                   if( lnoch(lnods(inode,ielem)) == NOFRI ) then
                      cont = cont + 1
                   end if
                end do
                if ( cont == element_type(ielty) % node_faces(iface) ) then
                   nboun_hole = nboun_hole + 1
                   ltypb_aux(nboun_hole) = element_type(ielty) % type_faces(iface)
                   lboch_aux(nboun_hole) = BOFRI
                   lnnob_aux(nboun_hole) = element_type(ielty) % node_faces(iface) 
                   do inodb = 1, element_type(ielty) % node_faces(iface)
                      inode                       = element_type(ielty) % list_faces(inodb,iface) 
                      lboel_aux(inodb,nboun_hole) = inode
                      lnodb_aux(inodb,nboun_hole) = lnods(inode,ielem)
                   end do
                   lelbo_aux(nboun_hole) = ielem 
                end if
             end do
          end if
       end if
    end do
    !
    ! Reallocate memory for boundary mesh
    !
    nboun_ori = nboun_cou
    nboun_cou = nboun_ori + nboun_hole

    if( nboun_hole /= 0 ) then

       if( number_of_holes == 0 ) then
          nullify(lnodb_cou)
          nullify(ltypb_cou)
          nullify(lboch_cou)
          nullify(lnnob_cou)
          nullify(lboel_cou)
          nullify(lelbo_cou)
          call memory_copy(memor_cou,'LNODB_COU','merge_hole_boundary',lnodb,lnodb_cou,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_cou,'LTYPB_COU','merge_hole_boundary',ltypb,ltypb_cou,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_cou,'LBOCH_COU','merge_hole_boundary',lboch,lboch_cou,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_cou,'LNNOB_COU','merge_hole_boundary',lnnob,lnnob_cou,'DO_NOT_DEALLOCATE')                
          call memory_copy(memor_cou,'LBOEL_COU','merge_hole_boundary',lboel,lboel_cou,'DO_NOT_DEALLOCATE')                
          call memory_copy(memor_cou,'LELBO_COU','merge_hole_boundary',lelbo,lelbo_cou,'DO_NOT_DEALLOCATE')                
       end if

       number_of_holes = number_of_holes + 1

       call memory_resize(memor_cou,'LNODB_COU','merge_hole_boundary',lnodb_cou,mnodb,nboun_cou)
       call memory_resize(memor_cou,'LTYPB_COU','merge_hole_boundary',ltypb_cou,nboun_cou)
       call memory_resize(memor_cou,'LBOCH_COU','merge_hole_boundary',lboch_cou,nboun_cou)
       call memory_resize(memor_cou,'LNNOB_COU','merge_hole_boundary',lnnob_cou,nboun_cou)
       call memory_resize(memor_cou,'LBOEL_COU','merge_hole_boundary',lboel_cou,mnodb,nboun_cou)
       call memory_resize(memor_cou,'LELBO_COU','merge_hole_boundary',lelbo_cou,nboun_cou)
       !
       ! Merge new boundary mesh
       !         
       iboun = nboun_ori
       do jboun = 1,nboun_hole
          iboun = iboun + 1
          ltypb_cou(iboun)         = ltypb_aux(jboun)
          lboch_cou(iboun)         = lboch_aux(jboun)
          lnnob_cou(iboun)         = lnnob_aux(jboun)
          lnodb_cou(1:mnodb,iboun) = lnodb_aux(1:mnodb,jboun)
          lboel_cou(1:mnodb,iboun) = lboel_aux(1:mnodb,jboun)
          lelbo_cou(iboun)         = lelbo_aux(jboun)
       end do

    end if
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_cou,'LNODB_AUX','merge_hole_boundary',lnodb_aux) 
    call memory_deallo(memor_cou,'LTYPB_AUX','merge_hole_boundary',ltypb_aux)       
    call memory_deallo(memor_cou,'LBOCH_AUX','merge_hole_boundary',lboch_aux)       
    call memory_deallo(memor_cou,'LNNOB_AUX','merge_hole_boundary',lnnob_aux)       
    call memory_deallo(memor_cou,'LBOEL_AUX','merge_hole_boundary',lboel_aux)       
    call memory_deallo(memor_cou,'LELBO_AUX','merge_hole_boundary',lelbo_aux)       

  end subroutine merge_hole_boundary


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    25/04/2022
  !> @brief   Create mask to identify hole elements
  !> @details Create mask to identify hole elements
  !
  !-----------------------------------------------------------------------
  subroutine create_mask(lmask)

    implicit none

    logical(lg), intent(inout) :: lmask(:)
    integer(ip)                :: ielem

    do ielem = 1,nelem
       if( lelch(ielem) /= ELHOL ) then
          lmask(ielem) = .true.
           lmask(ielem) = .true. 
          lmask(ielem) = .true.
       end if
    end do

  end subroutine create_mask


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    25/04/2022
  !> @brief   Allocate variables for multi-code hole-cutting
  !> @details Allocate variables for multi-code hole-cutting, that is the
  !>          "cou_holcut_multicode" subroutine.
  !
  !-----------------------------------------------------------------------
  subroutine allocate_variables()

    implicit none

    if( i_am_backg ) then
       ! Nullifications
       nullify( lmask )
       nullify( lwetno )

       ! Allocations
       call memory_alloca(memor_cou,'LWETNO','allocate_variables',lwetno,npoin)
       call memory_alloca(memor_cou,'LMASK', 'allocate_variables',lmask, nelem)

       ! Initializations
       if( npoin > 0 ) lwetno = 0_ip
       if( nelem > 0 ) lmask  = .false.
    end if

  end subroutine allocate_variables


  !-----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    25/04/2022
  !> @brief   Deallocate variables for multi-code hole-cutting
  !> @details Deallocate variables for multi-code hole-cutting, that is
  !>          the "cou_holcut_multicode" subroutine.
  !
  !-----------------------------------------------------------------------
  subroutine deallocate_variables()

    implicit none

    if( i_am_backg ) then
       call memory_deallo(memor_cou,'LWETNO','deallocate_variables',lwetno)
       call memory_deallo(memor_cou,'LMASK', 'deallocate_variables',lmask)
    end if

  end subroutine deallocate_variables
  

end module mod_holcut
!> @}
