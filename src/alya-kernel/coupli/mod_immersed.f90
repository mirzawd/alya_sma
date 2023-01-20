!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Immersed coupling functions
!> @file    mod_immersed.f90
!> @author  David Oks
!> @date    13/04/2022
!> @brief   ToolBox for immersed couplings
!> @details ToolBox for immersed couplings
!>
!> @{
!------------------------------------------------------------------------

module mod_immersed

  use def_kintyp_basic,   only : ip, rp, lg
  use def_master,         only : ISEQUEN
  use def_master,         only : current_code, current_zone
  use def_master,         only : INOTMASTER
  use def_domain,         only : npoin, nelem, ndime, coord, ltype, lnoch, lelch
  use def_kermod,         only : ielse
  use mod_parall,         only : color_target, color_source
  use mod_memory,         only : memory_alloca, memory_deallo
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_type
  use def_coupli,         only : ON_IMMERSED_MESH
  use def_coupli,         only : ON_EMBEDDED_MESH
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : typ_coupling_wet
  use def_coupli,         only : memor_cou
  use def_coupli,         only : kfl_efect
  use mod_messages,       only : messages_live

  implicit none

  private
  
  character(13), parameter :: vacal='mod_immersed'

  public :: cou_update_couplings
  public :: cou_find_nearest_neighbor
  public :: cou_check_if_fringe_nodes_are_necessary
  public :: cou_get_nodes_in_surf
  public :: cou_get_wet_elements
  public :: cou_get_fringe_nodes
  
contains

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2020-12-21
   !> @brief   Update couplings for multi-code immersed-type problems
   !> @details Update couplings for multi-code immersed-type problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_update_couplings()

     use mod_coupling_memory
     use mod_communications, only : PAR_MAX
     use mod_holcut,         only : cou_holcut_multicode

     logical(lg) :: kfl_immersed = .false.
     logical(lg) :: i_am_backg   = .false.
     logical(lg) :: i_am_patch   = .false.
     integer(ip) :: icoup
     integer(ip) :: icoup_immers, icoup_mirror
     integer(ip) :: code_backg, code_patch
     integer(ip) :: istat = 1_ip
     
     !
     ! Check if coupling is of immersed type 
     !
     do icoup = 1,mcoup
        if( coupling_type(icoup) % where_type == ON_IMMERSED_MESH .or. &
          & coupling_type(icoup) % where_type == ON_EMBEDDED_MESH ) then
           kfl_immersed = .true.
           icoup_immers = icoup
           icoup_mirror = coupling_type(icoup) % mirror_coupling
           code_backg   = coupling_type(icoup) % code_source
           code_patch   = coupling_type(icoup) % code_target
        end if
     end do
     if( current_code == code_backg ) i_am_backg = .true.
     if( current_code == code_patch ) i_am_patch = .true.
     !
     ! Immersed-type coupling: 
     !
     if ( kfl_immersed ) then

        call messages_live('UPDATE COUPLINGS','START SECTION')
        call cou_reset_elem_types(i_am_backg)                             ! Reset element & node types
        call cou_update_coords('UPDATE', i_am_patch)                      ! Update patch coords to deformed configuration
        call cou_reset_par_bin_structure()                                ! Reset parallel bin structures
        call cou_save_previous_wetnodes()                                 ! Save last timestep wetnodes
        call COU_REINITIALIZATION_ALL_COUPLINGS(WET=.true.,INPUT=.false.) ! Reinitialize coupling data-structures
        call cou_define_wet_geometry()                                    ! Redefine wet geometry
        current_zone = 0_ip
        call cou_initialize_coupling(istat)                               ! Reinitialize coouplings
        call cou_update_coords('RESET', i_am_patch)                       ! Reset patch coords to reference configuration
        call cou_holcut_multicode()                                       ! Create hole in background mesh
        call cou_get_fresh_nodes()                                        ! Get freshly uncovered nodes in background
        if( kfl_efect ) call domarr(2_ip)                                 ! Recalculate all arrays dependent on mesh
        call messages_live('UPDATE COUPLINGS','END SECTION')

     end if

   end subroutine cou_update_couplings

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Update coords of patch for multi-code immersed-type problems
   !> @details Update coords of patch for multi-code immersed-type problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_update_coords(task, i_am_patch)

      use def_master,                  only : ID_KERMOD
      use def_master,                  only : ITER_K 
      use def_master,                  only : mem_modul
      use def_master,                  only : displ
      use mod_memory,                  only : memory_copy
      use def_domain,                  only : coord
      use def_domain,                  only : cooin
      use def_domain,                  only : ndime

      implicit none

      character(*),             intent(in) :: task
      logical(lg),              intent(in) :: i_am_patch
      integer(ip)                          :: ipoin, idime

      if( INOTMASTER .and. i_am_patch ) then

         if( task == 'UPDATE' ) then
            !
            ! Allocate COOIN: Initial coordinates
            !
            if( .not. associated(cooin) ) then
               call memory_copy(mem_modul(1:2,ID_KERMOD),'COOIN',vacal,coord,cooin,'DO_NOT_DEALLOCATE')
            end if
            !
            ! Update COORD to update configuration
            !
            if( associated(displ) ) then
               do ipoin = 1,npoin
                  do idime = 1,ndime
                     coord(idime,ipoin) = cooin(idime,ipoin) + displ(idime,ipoin,ITER_K)
                  end do
               end do
            end if

         else if( task == 'RESET' ) then
            !
            ! Reset COORD to reference configuration
            !
            do ipoin = 1,npoin
               do idime = 1,ndime
                  coord(idime,ipoin) = cooin(idime,ipoin)
               end do
            end do
         end if

      end if

   end subroutine cou_update_coords

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Reset parallel bin structures after coupling update
   !> @details Reset parallel bin structures after coupling update
   !>
   !-----------------------------------------------------------------------

   subroutine cou_reset_par_bin_structure()

      use mod_elsest,                  only : elsest_deallocate
      use mod_elsest,                  only : elsest_initialization
      use mod_par_bin_structure,       only : par_bin_structure
      use mod_par_bin_structure,       only : par_bin_structure_deallocate
      use mod_communications,          only : PAR_BARRIER

      call PAR_BARRIER()

      if ( INOTMASTER ) call elsest_deallocate(ielse)
      call elsini()
      call par_bin_structure_deallocate()
      call par_bin_structure()

   end subroutine cou_reset_par_bin_structure


   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Save previous wet nodes 
   !> @details Save last time-step wet nodes for multi-code immersed-type
   !>          problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_save_previous_wetnodes()
      !
      ! Save previous wet-nodes
      !
      type(typ_coupling_wet), pointer :: wet 
      integer(ip)                     :: icoup, kpoin

      do icoup = 1,mcoup

         nullify(wet)
         wet => coupling_type(icoup) % wet

         nullify(wet % lpoin_wet_prev)
         if (INOTMASTER) then
            if (wet % npoin_wet > 0_ip) then
               if( associated(wet % lpoin_wet_prev) )then
                  nullify( wet % lpoin_wet_prev )
                  call memory_deallo(memor_cou,'COUPLING % WET % LPOIN_WET_PREV','cou_save_previous_wetnodes',wet % lpoin_wet_prev)
               end if
               call memory_alloca(memor_cou,'COUPLING % WET % LPOIN_WET_PREV','cou_save_previous_wetnodes', &
                                  wet % lpoin_wet_prev,wet % npoin_wet)
               wet % npoin_wet_prev = wet % npoin_wet
               do kpoin = 1,wet % npoin_wet_prev
                  wet % lpoin_wet_prev(kpoin) = wet % lpoin_wet(kpoin)
               end do
            end if
         end if

      end do

   end subroutine cou_save_previous_wetnodes

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Get nodes newly transformed from wet to normal nodes 
   !> @details Get nodes newly transformed from wet to normal nodes for
   !>          multi-code immersed-type problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_get_fresh_nodes()
      !
      ! Get "fresh" nodes, that is, nodes that in last timestep where wet nodes but now are regular nodes 
      !
      type(typ_coupling_wet), pointer :: wet 
      logical(lg),            pointer :: lpoin_fresh(:)
      integer(ip)                     :: icoup, npoin_fresh, ipoin, jpoin, kpoin, mpoin, ipoin_fresh

      do icoup = 1,mcoup

         nullify(wet)
         nullify(lpoin_fresh)

         wet => coupling_type(icoup) % wet

         if (INOTMASTER) then
            if (wet % npoin_wet_prev > 0_ip) then

               call memory_alloca(memor_cou,'LPOIN_FRESH','cou_get_fresh_nodes',lpoin_fresh,wet % npoin_wet_prev)

               ! Build boolean array of length npoin_wet_prev which is True for fresh nodes
               npoin_fresh = 0_ip
               do kpoin = 1,wet % npoin_wet_prev
                  ipoin = wet % lpoin_wet_prev(kpoin)
                  lpoin_fresh(kpoin) = .true.
                  mpoin = 0_ip
                  jpoin = 0_ip
                  do while ( (mpoin .ne. wet % npoin_wet) .and. lpoin_fresh(kpoin) )
                     mpoin = mpoin + 1
                     jpoin = wet % lpoin_wet(mpoin)
                     lpoin_fresh(kpoin) = (ipoin .ne. jpoin) 
                  end do
                  if (lpoin_fresh(kpoin)) npoin_fresh = npoin_fresh + 1_ip
               end do
               wet % npoin_fresh = npoin_fresh
              
               ! Allocate list in which fresh nodes will be saved 
               if ( associated(wet % lpoin_fresh) ) then
                  nullify( wet % lpoin_fresh )
                  call memory_deallo(memor_cou,'COUPLING % WET % LPOIN_FRESH','cou_get_fresh_nodes',wet % lpoin_fresh)
               end if
               call memory_alloca(memor_cou,'LPOIN_FRESH','cou_get_fresh_nodes',wet % lpoin_fresh,npoin_fresh)

               ! Build array which saves nodal indices of fresh nodes 
               ipoin_fresh = 0
               do kpoin = 1,wet % npoin_wet_prev
                  if (lpoin_fresh(kpoin)) then
                     ipoin_fresh = ipoin_fresh + 1
                     wet % lpoin_fresh(ipoin_fresh) = wet % lpoin_wet_prev(kpoin)
                  end if
               end do

               call memory_deallo(memor_cou,'LPOIN_FRESH','cou_get_fresh_nodes',lpoin_fresh)

            end if
         end if

      end do

   end subroutine cou_get_fresh_nodes

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @date    2022-03-23
   !> @brief   Find nearest non-wet neighbor node to a given wetnode
   !> @details Find nearest non-wet neighbor node to a given wetnode for
   !>          multi-code immersed-type problems
   !>
   !-----------------------------------------------------------------------

   subroutine cou_find_nearest_neighbor(jpoin_neighbor, ipoin_mine, icoup, iswet)

     use def_domain, only : c_dom, r_dom

     implicit none

     integer(ip), intent(out) :: jpoin_neighbor
     integer(ip), intent(in)  :: ipoin_mine 
     integer(ip), intent(in)  :: icoup 
     logical(lg), intent(in)  :: iswet
     logical(lg)              :: jpoin_was_wet
     integer(ip)              :: jpoin, kpoin, izdom
     real(rp)                 :: disti, distj

     disti = 1e16_rp

     ! Loop over neighbors
     do izdom = r_dom(ipoin_mine), r_dom(ipoin_mine+1)-1
        jpoin = c_dom(izdom)

        ! Check if node was wet in previous timestep
        jpoin_was_wet = .false.
        neighbors: do kpoin = 1,coupling_type(icoup) % wet % npoin_wet_prev
            if (jpoin == coupling_type(icoup) % wet % lpoin_wet_prev(kpoin)) then
               jpoin_was_wet = .true.
               exit neighbors
            end if
        end do neighbors

        ! Compute squared distance to neighbor
        distj = sum((coord(:,jpoin) - coord(:,ipoin_mine))**2)

        ! If is smaller and satisfies wet/dry node condition => keep
        if ( (distj < disti) .and. (jpoin_was_wet .eqv. iswet) ) then
           jpoin_neighbor = jpoin
           disti          = distj
        end if

     end do

   end subroutine cou_find_nearest_neighbor

   !-----------------------------------------------------------------------
   !>
   !> @author  David Oks
   !> @brief   Check if I have to calculate fringe nodes of wet nodes
   !> @details Check if I have to calculate fringe nodes of wet nodes for
   !>          Embedded Finite Element Coupling Technique
   !>
   !-----------------------------------------------------------------------

   subroutine cou_check_if_fringe_nodes_are_necessary(icoup)

     implicit none

     integer(ip),      intent(in)   :: icoup 

     ! Activate fringe wetnodes flag if they are necessary     
     if( coupling_type(icoup) % where_type == ON_EMBEDDED_MESH ) then
        coupling_type(icoup) % wet % kfl_get_fringe = .true.
        kfl_efect                                   = .true.
     end if

   end subroutine cou_check_if_fringe_nodes_are_necessary

   !-----------------------------------------------------------------------
   !
   !> @author  David Oks
   !> @brief   Sets fixities on fringe nodes 
   !> @details Sets fixities on fringe nodes
   !
   !-----------------------------------------------------------------------

   subroutine cou_set_fixities_on_fringe_nodes(icoup)
   
     use def_elmtyp,         only : NOHOL
     use def_master,         only : solve 

     implicit none
   
     integer(ip), intent(in)    :: icoup 
     integer(ip)                :: ipoin

     if ( INOTMASTER ) then
        !
        ! Set fringe node fixities to Dirichlet type in all dimensions
        !
        do ipoin = 1,npoin
           if (coupling_type(icoup) % wet % kfl_fringe_wetnodes(ipoin) == 1_ip) then
              solve(1) % block_array(1) % kfl_fixno(1:ndime,ipoin) = 99_ip
           end if
        end do

     end if
   
   end subroutine cou_set_fixities_on_fringe_nodes

   !-----------------------------------------------------------------------
   !
   !> @author  David Oks
   !> @brief   Sets fixities on hole nodes 
   !> @details Sets fixities on hole nodes
   !
   !-----------------------------------------------------------------------

   subroutine cou_set_fixities_on_hole_nodes()
   
     use def_elmtyp,         only : NOHOL
     use def_master,         only : solve 

     implicit none
   
     integer(ip)                :: ipoin

     if ( INOTMASTER ) then
        !
        ! Set hole node fixities to Dirichlet type in all dimensions
        !
        do ipoin = 1,npoin
           if( lnoch(ipoin) == NOHOL ) then
              solve(1) % kfl_fixno(1:ndime,ipoin) = 99_ip
           end if
        end do

     end if
   
   end subroutine cou_set_fixities_on_hole_nodes
  
   !-----------------------------------------------------------------------
   !
   !> @author  David Oks
   !> @brief   Reset element types 
   !> @details Reset element types if using Embedded Finite Element Coupling
   !>          Technique
   !
   !-----------------------------------------------------------------------

   subroutine cou_reset_elem_types(i_am_backg)
   
     use def_elmtyp,         only : NOHOL, NOFEM, ELHOL, ELFEM

     implicit none
   
     logical(lg), intent(in) :: i_am_backg
     integer(ip)             :: ielem, ipoin

     ! Recover original element types 
     if( kfl_efect .and. i_am_backg ) then
        if ( INOTMASTER ) then
           do ielem = 1,nelem
              lelch(ielem) = ELFEM              ! TODO: Do as Tobal & save initial elem types
              ltype(ielem) = abs(ltype(ielem))
           end do
           do ipoin=1,npoin
              lnoch(ipoin) = NOFEM
           end do
        end if
        call PAR_INTERFACE_NODE_EXCHANGE(lnoch,'MIN','IN THE WORLD')
     end if
   
   end subroutine cou_reset_elem_types
  
   subroutine cou_get_nodes_in_surf(mesh, dista, lnins)

     use def_kintyp_mesh_basic

     implicit none
     class(mesh_type_basic),           intent(in)    :: mesh
     real(rp),                pointer, intent(in)    :: dista(:)
     integer(ip),             pointer, intent(inout) :: lnins(:)
     integer(ip)                                     :: ipoin

     do ipoin = 1, mesh % npoin
        if( dista(ipoin) <= 0.0_rp ) lnins(ipoin) = 1_ip
     end do

     call PAR_INTERFACE_NODE_EXCHANGE(lnins, 'MAX', 'IN THE WORLD')

   end subroutine cou_get_nodes_in_surf

   subroutine cou_get_wet_elements(mesh, lnins, lewet)

     use def_kintyp_mesh_basic

     implicit none
     class(mesh_type_basic),           intent(in)    :: mesh
     integer(ip),             pointer, intent(in)    :: lnins(:)
     integer(ip),             pointer, intent(inout) :: lewet(:)
     integer(ip)                                     :: ielem, ipoin, inode

     do ielem = 1, mesh % nelem
        do inode = 1, mesh % mnode
           ipoin = mesh % lnods(inode,ielem)
           if( lnins(ipoin) == 1_ip ) then
              lewet(ielem) = 1_ip
           end if
        end do
     end do

   end subroutine cou_get_wet_elements

   subroutine cou_get_fringe_nodes(mesh, lewet, lnins, lnfri)

     use def_kintyp_mesh_basic

     implicit none
     class(mesh_type_basic),           intent(in)    :: mesh
     integer(ip),             pointer, intent(in)    :: lewet(:)
     integer(ip),             pointer, intent(in)    :: lnins(:)
     integer(ip),             pointer, intent(inout) :: lnfri(:)
     integer(ip)                                     :: ielem, ipoin, inode

     do ielem = 1, mesh % nelem
        if( lewet(ielem) == 1_ip ) then
           do inode = 1, mesh % mnode
              ipoin = mesh % lnods(inode,ielem)
              if( lnins(ipoin) /= 1_ip ) lnfri(ipoin) = 1_ip
           end do
        end if
     end do

   end subroutine cou_get_fringe_nodes

end module mod_immersed
