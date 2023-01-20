!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outset()
  !-----------------------------------------------------------------------
  !****f* partis/chm_outset
  ! NAME
  !    chm_outset
  ! DESCRIPTION
  !    Compute and write results on sets
  ! USES
  ! USED BY
  !    chm_output
  !***
  !-----------------------------------------------------------------------
!  use def_parame
  use def_master,        only : INOTMASTER, INOTSLAVE, postp, vbset, cutim, ittim, mitim, timef, conce
  use def_domain,        only : nbset, neset, nnset, lesec, lbsec, lnsec
  use def_kintyp,        only : ip, rp
  use def_chemic,        only : nclas_chm
  use mod_output_postprocess, only : output_postprocess_node_sets_parall
  use mod_output_postprocess, only : output_postprocess_boundary_sets_parall
  use mod_output_postprocess, only : output_postprocess_element_sets_parall
  implicit none
  integer(ip) :: ieset,ibset,inset,iclas

  external    :: chm_elmset
  external    :: chm_bouset

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse)>0 ) then
      if(( mod(ittim, postp(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            do ieset = 1,neset
               call chm_elmset(lesec(ieset),ieset)
            end do
         end if
         call output_postprocess_element_sets_parall()
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepboset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then

         do ibset = 1,nbset
            call chm_bouset(lbsec(ibset),ibset)
         end do
         call output_postprocess_boundary_sets_parall()
         !
         ! Normalize results
         !
         if( INOTSLAVE ) then
            do ibset = 1,nbset
               do iclas=1,min(8_ip,nclas_chm)
                  !
                  ! <Yk> = int_S rho*u* Y_k dS / int_S rho*u dS, for k = 1,...,8
                  !
                  if(postp(1) % npp_setsb(2)/=0.and.vbset(1,ibset) /= 0.0_rp) then
                     vbset( iclas+1_ip,ibset) = vbset( iclas+1_ip,ibset) / vbset( 1,ibset)
                  end if

               end do
            end do
         end if
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn)>0 ) then
      if(( mod(ittim, postp(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            do inset=1,nnset
               if(lnsec(inset)/=0) then
                  iclas=1
                  if(postp(1) % npp_setsn(1)/=0) postp(1) % vnset(1,inset)=conce(lnsec(inset),iclas,1)
               end if
            end do
         end if
         call output_postprocess_node_sets_parall()
      end if
  end if

end subroutine chm_outset
