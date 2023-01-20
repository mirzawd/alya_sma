!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_outset()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_outset
  ! NAME 
  !    exm_outset
  ! DESCRIPTION
  !    Compute and write results on sets:
  !    - Element, boundary and node sets:
  !      1. INTRA
  !      2. EXTRA
  !      3. RECOV
  ! USES
  ! USED BY
  !    exm_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_iofile
  use mod_output_postprocess, only : output_postprocess_node_sets_parall
  use mod_output_postprocess, only : output_postprocess_element_sets_parall
  implicit none
  integer(ip) :: inset
  integer(ip) :: ieset
  integer(ip) :: ivar
  integer(ip), dimension(4), parameter :: vars_woese_to_normalize = [1,2,3,4]

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then

         if( INOTMASTER ) then
            do inset = 1,nnset
               call exm_nodset(lnsec(inset),inset)
            end do
         end if
         call output_postprocess_node_sets_parall()
      end if
  end if


   !----------------------------------------------------------------------
   !
   ! Element sets
   !
   !----------------------------------------------------------------------

   if( maxval(postp(1) % npp_setse) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            do ieset = 1,neset
               call exm_elmset(lesec(ieset),ieset)
            end do
         end if
         !
         ! Parall
         !
         call output_postprocess_element_sets_parall()

         !
         ! Averaged/Norm of some variables
         !
         if( INOTSLAVE ) then
            do ieset = 1,neset
               do ivar = 1,size(vars_woese_to_normalize)
                  if( postp(1) % npp_setse( vars_woese_to_normalize(ivar) ) /= 0 .and. veset(postp(1) % nvaes+1,ieset) > 0.0_rp ) then
                     veset( vars_woese_to_normalize(ivar),ieset ) = veset( vars_woese_to_normalize(ivar),ieset )/veset(postp(1) % nvaes+1,ieset)
                  end if
               end do
            end do
         end if         
      end if
   end if



end subroutine exm_outset
