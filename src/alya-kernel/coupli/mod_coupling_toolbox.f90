!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    mod_coupling_toolbox.f90
!> @author  houzeaux
!> @date    2018-04-24
!> @brief   Toolbox ffor the coupling
!> @details Toolbox for setting a up coupling
!>
!-----------------------------------------------------------------------

module mod_coupling_toolbox

  use def_kintyp,            only : ip,rp,i1p,lg
  use def_master,            only : INOTMASTER
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_size
  use mod_parall,            only : PAR_COMM_COLOR_PERM
  use mod_parall,            only : par_bin_comin
  use mod_parall,            only : par_bin_comax
  use mod_parall,            only : par_bin_part
  use mod_parall,            only : par_bin_boxes
  use mod_parall,            only : par_bin_size
  use mod_parall,            only : par_part_comin
  use mod_parall,            only : par_part_comax
  use mod_parall,            only : par_part_in_color
  use def_coupli,            only : memor_cou
  use mod_parall,            only : typ_bin_structure
  implicit none
  private

  public :: coupling_toolbox_points_in_partitions

contains

   !-----------------------------------------------------------------------
   !>
   !> @author  houzeaux -  borrell
   !> @date    2018-11-26
   !> @brief   Find list of partitions to which I must send my points
   !> @details Find list of partitions to which I must send my points
   !>          MY_POINT_TO_PART(PP) % L(:)    = list of partitions where PP is in
   !>          MY_PART_TO_POINT(CPART) % L(:) = list of points located in CPART bounding box
   !>
   !-----------------------------------------------------------------------

   subroutine coupling_toolbox_points_in_partitions(&
         nn,ndime,xx,jcolo,bin_structure,kfl_multiple_source,&
         PAR_CURRENT_SIZE,PAR_CURRENT_RANK,PAR_WORLD_RANKS,&
         npoin_send,my_point_to_part,my_part_to_point)

      use mod_parall, only : PAR_CPU_TO_COLOR
      use mod_maths,  only : maths_mapping_coord_to_3d
      use mod_maths,  only : maths_in_box

      integer(ip),             intent(in)    :: nn                  !< Number of points
      integer(ip),             intent(in)    :: ndime               !< Space dimension
      real(rp),    pointer,    intent(in)    :: xx(:,:)             !< Point coordinates
      integer(ip),             intent(in)    :: jcolo               !< Source coupling color
      type(typ_bin_structure), intent(inout) :: bin_structure       !< Bin structure
      integer(ip),             intent(in)    :: kfl_multiple_source !< Communicator size
      integer(ip),             intent(in)    :: PAR_CURRENT_SIZE    !< Communicator size
      integer(ip),             intent(in)    :: PAR_CURRENT_RANK    !< My current rank
      integer(ip), pointer,    intent(in)    :: PAR_WORLD_RANKS(:)  !< World ranks
      integer(ip), pointer,    intent(inout) :: npoin_send(:)       !< Points to send
      type(i1p),   pointer,    intent(inout) :: my_point_to_part(:) !< List of partitions the point belongs to
      type(i1p),   pointer,    intent(inout) :: my_part_to_point(:) !< List of partitions the point belongs to
      integer(ip)                            :: pp,ii,jj,kk
      integer(ip)                            :: ksize,kpart
      integer(ip)                            :: ipart,cpart
      integer(ip), pointer                   :: aux_myp2p(:)
      logical(lg), pointer                   :: aux_ptinc(:)

      character(100), PARAMETER :: vacal = "coupling_toolbox_points_in_partitions"

      !
      ! Allocate if necessary
      !
      if( nn > 0 ) then
         if( .not. associated(npoin_send) ) &
            call memory_alloca(memor_cou,'NPOIN_SEND',vacal,npoin_send,PAR_CURRENT_SIZE,'INITIALIZE',0_ip)
         if( .not. associated(my_point_to_part) )&
            call memory_alloca(memor_cou,'MY_POINT_TO_PART',vacal,my_point_to_part,nn,'INITIALIZE')
         if( .not. associated(my_part_to_point) )&
            call memory_alloca(memor_cou,'MY_PART_TO_POINT',vacal,my_part_to_point,PAR_CURRENT_SIZE,'INITIALIZE',0_ip)
         nullify(aux_myp2p,aux_ptinc)
         call memory_alloca(memor_cou,'AUX_MYP2P',vacal,aux_myp2p,PAR_CURRENT_SIZE)
         call memory_alloca(memor_cou,'AUX_PTINC',vacal,aux_ptinc,PAR_CURRENT_SIZE,'INITIALIZE',0_ip)
         if( associated(PAR_CPU_TO_COLOR) .and. jcolo /= 0 ) then
            do kpart = 0, PAR_CURRENT_SIZE-1
               ipart = PAR_WORLD_RANKS(kpart)
               aux_ptinc(kpart)=par_part_in_color(ipart,jcolo)
            enddo
         endif
      end if
      !
      ! MY_POINT_TO_PART(PP) % L(:) = list of partitions where PP is in
      !
      do pp = 1,nn
         call maths_mapping_coord_to_3d(ndime,bin_structure % boxes,bin_structure % comin,bin_structure % comax,xx(1:ndime,pp),ii,jj,kk)
         if( ii /= 0 ) then
            ksize = 0
            do kpart = 1,bin_structure % size(ii,jj,kk)
               cpart = bin_structure % part(ii,jj,kk) % l(kpart)
               if( kfl_multiple_source /= -1 .or. cpart /= PAR_CURRENT_RANK ) then
                  if( maths_in_box(ndime,xx(1:ndime,pp),bin_structure % part_comin(1:ndime,cpart),bin_structure % part_comax(1:ndime,cpart)) ) then
                     !
                     ! Take those partitions only in color JCOLO
                     !
                     if( associated(PAR_CPU_TO_COLOR) .and. jcolo /= 0 ) then
                        if( aux_ptinc(cpart) ) then
                           ksize = ksize + 1
                           npoin_send(cpart) = npoin_send(cpart) + 1
                           aux_myp2p(ksize) = cpart
                        end if
                     else
                        ksize = ksize + 1
                        npoin_send(cpart) = npoin_send(cpart) + 1
                        aux_myp2p(ksize) = cpart
                     end if
                  end if
               end if
            end do
            if( ksize > 0 ) then
               call memory_alloca(memor_cou,'MY_POINT_TO_PART % L','cou_init_interpolate_points_values',my_point_to_part(pp) % l,ksize)
               my_point_to_part(pp) % l(1:ksize) = aux_myp2p(1:ksize)
            end if
         end if
      end do
      !
      ! MY_PART_TO_POINT(CPART) % L(:) = list of points located in CPART bounding box
      !
      do cpart = 0,PAR_CURRENT_SIZE-1
         if( npoin_send(cpart) > 0 ) then
            call memory_alloca(memor_cou,'MY_PART_TO_POINT % L','cou_init_interpolate_points_values',my_part_to_point(cpart) % l,npoin_send(cpart))
         end if
         npoin_send(cpart) = 0
      end do
      do pp = 1,nn
         do kk = 1,memory_size(my_point_to_part(pp) % l)
            cpart = my_point_to_part(pp) % l(kk)
            npoin_send(cpart) = npoin_send(cpart) + 1
            my_part_to_point(cpart) % l(npoin_send(cpart)) = pp
         end do
      end do
      if( nn > 0 ) then
         call memory_deallo(memor_cou,'AUX_MYP2P',vacal,aux_myp2p)
         call memory_deallo(memor_cou,'AUX_PTINC',vacal,aux_ptinc)
      endif

   end subroutine coupling_toolbox_points_in_partitions

end module mod_coupling_toolbox
!> @}
